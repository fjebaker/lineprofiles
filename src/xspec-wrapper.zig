const std = @import("std");
const io = @import("io.zig");
const util = @import("utils.zig");
const lineprofile = @import("line-profile.zig");

// number of parameters in the table
const NPARAMS = 2;
const MODELPATH = "kerr-transfer-functions.fits";

// refinement for the energy grid
const REFINEMENT = 5;

var allocator = std.heap.c_allocator;

// singleton
var profile: ?lineprofile.LineProfileTable(NPARAMS, f32) = null;
var g_grid: []f32 = undefined;
var r_grid: []f32 = undefined;
var flux_cache: []f32 = undefined;

fn Parameters(comptime T: type) type {
    return struct {
        const Self = @This();
        a: T,
        inclination: T,
        eline: T,
        rmin: T,
        rmax: T,
        pub fn from_ptr(ptr: *const f64) Parameters(T) {
            const N = @typeInfo(Self).Struct.fields.len;
            var slice = @ptrCast([*]const f64, ptr)[0..N];
            return .{
                .a = @floatCast(T, slice[0]),
                // convert to cos(i)
                .inclination = @cos(std.math.degreesToRadians(
                    T,
                    @floatCast(T, slice[1]),
                )),
                .eline = @floatCast(T, slice[2]),
                .rmin = @floatCast(T, slice[3]),
                .rmax = @floatCast(T, slice[4]),
            };
        }
        pub fn table_parameters(self: Self) [NPARAMS]T {
            return [_]T{ self.a, self.inclination };
        }
    };
}

fn set_radial_grid(comptime T: type, rmin: T, rmax: T) void {
    util.inverse_grid_inplace(T, r_grid, rmin, rmax);
}

fn setup() !void {
    profile = try io.readFitsFile(NPARAMS, f32, MODELPATH, allocator);
    errdefer profile.?.deinit();

    std.debug.print(
        "kerrlineprofile: Read in {d} transfer function tables.\n",
        .{profile.?.transfer_functions.len},
    );

    // build r grid
    r_grid = try util.inverse_grid(f32, allocator, 3.0, 50.0, 2000);
    errdefer allocator.free(r_grid);

    // build some g grid, it will be refined anyway
    var gitt = util.RangeIterator(f32).init(0.1, 1.5, 200);
    g_grid = try gitt.drain(allocator);
    errdefer allocator.free(g_grid);

    // allocate output
    flux_cache = try allocator.alloc(f32, g_grid.len - 1);

    std.debug.print("kerrlineprofile: Finished one-time setup.\n", .{});
}

fn refine_grid(comptime T: type, grid: anytype, norm: T) void {
    const N = (grid.len - 1) * REFINEMENT;
    if (N != g_grid.len) {
        allocator.free(g_grid);
        g_grid = allocator.alloc(T, N) catch {
            @panic("kerrlineprofile: Failed to refine energy grid.");
        };
        allocator.free(flux_cache);
        flux_cache = allocator.alloc(T, N - 1) catch {
            @panic("kerrlineprofile: Failed to refine flux cache.");
        };
    }
    util.refine_grid(T, grid, g_grid, REFINEMENT, norm);
}

fn integrate_lineprofile(
    comptime T: type,
    energy: []const f64,
    flux: []f64,
    params: Parameters(T),
) void {
    // do we need to do first time setup?
    var lp = profile orelse blk: {
        setup() catch |e| {
            std.debug.print("kerrlineprofile: error: {!}\n", .{e});
            @panic("kerrlineprofile: fatal: COULD NOT INITIALIZE MODEL.");
        };
        break :blk profile.?;
    };

    // refine our grids
    refine_grid(T, energy, params.eline);
    set_radial_grid(T, params.rmin, params.rmax);

    // zero the flux cache
    for (flux_cache) |*f| f.* = 0;

    // do the integration
    var itf = lp.interpolate_parameters(params.table_parameters());
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

    const N = @intCast(usize, n_flux);
    // convert to slices
    const energy = @ptrCast([*]const f64, energy_ptr)[0 .. N + 1];
    var flux = @ptrCast([*]f64, flux_ptr)[0..N];

    var parameters = Parameters(f32).from_ptr(parameters_ptr);
    integrate_lineprofile(f32, energy, flux, parameters);
}
