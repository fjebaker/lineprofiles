const std = @import("std");
const io = @import("io.zig");
const util = @import("utils.zig");
const lineprofile = @import("line-profile.zig");

const NPARAMS = 2;
const NMODELPARAMS = NPARAMS + 1;
const REFINEMENT = 5;
const MODELPATH = "rel_table.fits";
var allocator = std.heap.c_allocator;

// singleton
var profile: ?lineprofile.LineProfileTable(NPARAMS, f32) = null;
var g_grid: []f32 = undefined;
var r_grid: []f32 = undefined;
var flux_cache: []f32 = undefined;

fn setup() !void {
    profile = try io.readFitsFile(NPARAMS, f32, MODELPATH, allocator);
    errdefer profile.?.deinit();

    std.debug.print("kerrlineprofile: Read in {d} transfer function tables.\n", .{profile.?.transfer_functions.len});

    // build r grid
    var ritt = util.RangeIterator(f32).init(3.0, 50.0, 2000);
    r_grid = try ritt.drain(allocator);
    errdefer allocator.free(r_grid);

    // build some g grid, it will be refined anyway
    var gitt = util.RangeIterator(f32).init(0.1, 1.5, 200);
    g_grid = try gitt.drain(allocator);
    errdefer allocator.free(g_grid);

    // allocate output
    flux_cache = try allocator.alloc(f32, g_grid.len - 1);

    std.debug.print("kerrlineprofile: Finished setup.\n", .{});
}

fn refine_grid(comptime T: type, grid: anytype, norm: T) void {
    const N = (grid.len - 1) * REFINEMENT;
    if (N != g_grid.len) {
        allocator.free(g_grid);
        g_grid = allocator.alloc(T, N) catch {
            @panic("Failed to refine energy grid.");
        };
        allocator.free(flux_cache);
        flux_cache = allocator.alloc(T, N - 1) catch {
            @panic("Failed to refine flux cache.");
        };
    }
    util.refine_grid(T, grid, g_grid, REFINEMENT, norm);
}

export fn kerrlineprofile(
    // all inputs are double precision
    energy_ptr: *const f64,
    n_flux: c_int,
    parameters_ptr: *const f64,
    spectrum: c_int,
    flux_ptr: *f64,
    flux_variance_ptr: *f64,
    init_ptr: *const u8,
) callconv(.C) void {
    // unused
    _ = init_ptr;
    _ = flux_variance_ptr;
    _ = spectrum;

    // do we need to do first time setup?
    var lp = profile orelse blk: {
        setup() catch @panic("COULD NOT INITIALIZE MODEL.");
        break :blk profile.?;
    };

    const N = @intCast(usize, n_flux);
    // convert to slices
    const energy = @ptrCast([*]const f64, energy_ptr)[0 .. N + 1];
    var flux = @ptrCast([*]f64, flux_ptr)[0..N];
    var pslice = @ptrCast([*]const f64, parameters_ptr)[0..NMODELPARAMS];

    // stack allocate
    var parameters: [NMODELPARAMS]f32 = undefined;

    // a, incl, Eline
    for (0..NMODELPARAMS) |i| parameters[i] = @floatCast(f32, pslice[i]);

    // convert inclination to mu
    parameters[1] = @cos(std.math.degreesToRadians(f32, parameters[1]));

    // refine our grid
    refine_grid(f32, energy, parameters[2]);

    // zero the flux cache
    for (flux_cache) |*f| f.* = 0;

    // do the integration
    var itf = lp.interpolate_parameters(parameters[0..NPARAMS].*);
    itf.integrate(r_grid, g_grid, flux_cache);

    // rebin for output
    var j: usize = 0;
    for (0..flux.len) |i| {
        // ensure it is zeroed
        flux[i] = 0;
        // sum up the grid values
        for (0..REFINEMENT) |_| {
            flux[i] += @floatCast(f64, flux_cache[j]);
            j += 1;
            if (j == flux_cache.len) {
                break;
            }
        }
        // normalize to counts per bin
        const e_midpoint = 0.5 * (energy[i + 1] + energy[i]);
        flux[i] = flux[i] / e_midpoint;
    }
    // normalize output by area
    util.normalize(f64, flux);
}
