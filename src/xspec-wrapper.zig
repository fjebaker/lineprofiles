const std = @import("std");
const io = @import("io.zig");
const util = @import("utils.zig");
const lineprofile = @import("line-profile.zig");

const NPARAMS = 2;
var allocator = std.heap.c_allocator;

// singleton
var profile: ?lineprofile.LineProfileTable(NPARAMS, f32) = null;
var g_grid: []f32 = undefined;
var r_grid: []f32 = undefined;
var flux_cache: []f32 = undefined;

fn setup() !*lineprofile.LineProfileTable(NPARAMS, f32) {
    profile = try io.readFitsFile("rel_table.fits", allocator);
    errdefer profile.deinit();

    // build r grid
    var ritt = util.RangeIterator(f32).init(3.0, 50.0, 3000);
    r_grid = try ritt.drain(allocator);
    errdefer allocator.free(r_grid);

    // build g grid
    var gitt = util.RangeIterator(f32).init(0.1, 1.5, 200);
    g_grid = try gitt.drain(allocator);
    errdefer allocator.free(g_grid);

    // allocate output
    flux_cache = try allocator.alloc(f32, g_grid.len - 1);
}

export fn modelfunc(
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
    std.mem.copy(f32, parameters, parameters_ptr[0..NPARAMS]);

    // convert to slices
    const energy = energy_ptr[0 .. N + 1];
    var flux = flux_ptr[0..N];
    _ = energy;
    _ = flux;

    // do we need to do first time setup?
    var lp = profile orelse setup() catch @panic("COULD NOT INITIALIZE MODEL.");

    // zero the flux cache
    for (flux_cache) |*f| f.* = 0;

    // do the integration
    var itf = lp.interpolate_parameters(parameters);
    itf.integrate(r_grid, g_grid, flux_cache);

    // rebin for output
}
