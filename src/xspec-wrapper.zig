const std = @import("std");
const io = @import("io.zig");
const util = @import("utils.zig");
const lineprofile = @import("line-profile.zig");

const NPARAMS = 2;
const REFINEMENT = 10;
var allocator = std.heap.c_allocator;

// singleton
var profile: ?lineprofile.LineProfileTable(NPARAMS, f32) = null;
var g_grid: []f32 = undefined;
var r_grid: []f32 = undefined;
var flux_cache: []f32 = undefined;

fn setup() !void {
    profile = try io.readFitsFile("rel_table.fits", allocator);
    errdefer profile.?.deinit();

    // build r grid
    var ritt = util.RangeIterator(f32).init(3.0, 50.0, 3000);
    r_grid = try ritt.drain(allocator);
    errdefer allocator.free(r_grid);

    // build some g grid, it will be refined anyway
    var gitt = util.RangeIterator(f32).init(0.1, 1.5, 200);
    g_grid = try gitt.drain(allocator);
    errdefer allocator.free(g_grid);

    // allocate output
    flux_cache = try allocator.alloc(f32, g_grid.len - 1);
}

fn refine_grid(comptime T: type, grid: []const T) void {
    if (grid.len * REFINEMENT != g_grid.len) {
        if (!allocator.resize(g_grid, grid.len * REFINEMENT)) {
            @panic("Failed to refine energy grid.");
        }
        if (!allocator.resize(flux_cache, grid.len * REFINEMENT - 1)) {
            @panic("Failed to refine flux cache.");
        }
    }
    util.refine_grid(T, grid, g_grid, REFINEMENT);
}

export fn kerrlineprofile(
    energy_ptr: *const f32,
    n_flux: c_int,
    parameters_ptr: *const f32,
    spectrum: c_int,
    flux_ptr: *f32,
    flux_variance_ptr: *f32,
    init_ptr: *const u8,
) callconv(.C) void {
    // unused
    _ = init_ptr;
    _ = flux_variance_ptr;
    _ = spectrum;

    const N = @intCast(usize, n_flux);

    // stack allocate
    var parameters: [NPARAMS]f32 = undefined;
    std.mem.copy(
        f32,
        &parameters,
        @ptrCast([*]const f32, parameters_ptr)[0..NPARAMS],
    );
    // convert inclination to mu
    parameters[1] = std.math.cos(parameters[1]);

    // convert to slices
    const energy = @ptrCast([*]const f32, energy_ptr)[0 .. N + 1];
    var flux = @ptrCast([*]f32, flux_ptr)[0..N];

    // refine our grid
    refine_grid(f32, energy);
    // convert from E to g
    for (g_grid) |*g| g.* = g.* / 6.4;

    // do we need to do first time setup?
    var lp = profile orelse blk: {
        setup() catch @panic("COULD NOT INITIALIZE MODEL.");
        break :blk profile.?;
    };

    // zero the flux cache
    for (flux_cache) |*f| f.* = 0;

    // do the integration
    var itf = lp.interpolate_parameters(parameters);
    itf.integrate(r_grid, g_grid, flux_cache);

    // rebin for output
    for (0..flux.len) |i| {
        // ensure it is zeroed
        flux[i] = 0;
        // sum up the grid values
        for (0..REFINEMENT) |k| {
            flux[i] += flux_cache[i + k];
        }
        // normalize to counts per bin
        const e_midpoint = 0.5 * (energy[i + 1] + energy[i]);
        flux[i] = flux[i] / e_midpoint;
    }
    // normalize output by area
    util.normalize(f32, flux);
}
