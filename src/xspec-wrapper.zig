const std = @import("std");
const io = @import("io.zig");
const util = @import("utils.zig");
const lineprofile = @import("line-profile.zig");
const emissivity = @import("emissivity.zig");
const convolution = @import("convolution.zig");

// number of parameters in the table
const NPARAMS = 2;
const MODELPATH = "./kerr-transfer-functions.fits";

pub fn getModelPath() ![]const u8 {
    return MODELPATH;
}

// refinement for the energy grid
const REFINEMENT = 5;

// convolution has a fixed bin width for the model
// that is to say, the model is evaluated on a grid
// of g values in [0,2], where the bin width is:
pub const CONVOLUTION_RESOLUTION = 1e-3;

var allocator = std.heap.c_allocator;

// singleton
var profile: ?lineprofile.LineProfileTable(NPARAMS, f32) = null;
var g_grid: []f32 = undefined;
var r_grid: []f32 = undefined;
var flux_cache: []f32 = undefined;
// convolution specifics, need to be double precision to mimic
// the XSPEC behaviour
var conv_g_grid: ?[]f64 = null;
var conv_flux: []f64 = undefined;

fn set_radial_grid(comptime T: type, rmin: T, rmax: T) void {
    util.inverse_grid_inplace(T, r_grid, rmin, rmax);
}

fn setup() !void {
    profile = try io.readFitsFile(NPARAMS, f32, 30, try getModelPath(), allocator);
    errdefer profile.?.deinit();

    std.debug.print(
        "kerrlineprofile: Read in {d} transfer function tables.\n",
        .{profile.?.transfer_functions.len},
    );

    // build r grid
    r_grid = try util.inverse_grid(f32, allocator, 1.0, 50.0, 2000);
    errdefer allocator.free(r_grid);

    // build some g grid, it will be refined anyway
    var gitt = util.RangeIterator(f32).init(0.1, 1.5, 200);
    g_grid = try gitt.drain(allocator);
    errdefer allocator.free(g_grid);

    // allocate output
    flux_cache = try allocator.alloc(f32, g_grid.len - 1);

    std.debug.print("kerrlineprofile: Finished one-time setup.\n", .{});
}

fn convolve_setup() !void {
    // build the convolution energy grid
    var itt = util.RangeIterator(f64).init(
        CONVOLUTION_RESOLUTION,
        2.0,
        @trunc((2.0 - CONVOLUTION_RESOLUTION) / CONVOLUTION_RESOLUTION),
    );
    var conv_energy = try itt.drain(allocator);
    errdefer allocator.free(conv_energy);

    // temporary flux array
    var model_flux = try allocator.alloc(f64, conv_energy.len - 1);
    errdefer allocator.free(model_flux);

    // assign to singletons
    conv_g_grid = conv_energy;
    conv_flux = model_flux;
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
    params: anytype,
    emis: anytype,
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
    itf.integrate(r_grid, g_grid, flux_cache, emis);

    // rebin for output
    var j: usize = 0;
    var total_flux: f64 = 0;
    for (0..flux.len) |i| {
        // ensure it is zeroed
        flux[i] = 0;
        // sum up the grid values
        for (0..REFINEMENT) |_| {
            flux[i] += @as(f64, @floatCast(flux_cache[j]));
            j += 1;
            if (j == flux_cache.len) {
                break;
            }
        }
        // normalize to counts per bin
        const e_midpoint = 0.5 * (energy[i + 1] + energy[i]);
        flux[i] = flux[i] / e_midpoint;
        total_flux += flux[i];
    }
    // normalize output by area
    for (flux) |*f| f.* /= total_flux;
}

fn kerr_convolve(
    energy: []const f64,
    flux: []f64,
    parameters: anytype,
    emis: anytype,
) void {
    const conv_energy = conv_g_grid orelse blk: {
        convolve_setup() catch |e| {
            std.debug.print("kerrlineprofile: error: {!}\n", .{e});
            @panic("kerrlineprofile: fatal: COULD NOT ALLOCATE CONVOLUTION MEMORY.");
        };
        break :blk conv_g_grid.?;
    };

    // invoke the model
    integrate_lineprofile(f32, conv_energy, conv_flux, parameters, emis);

    // make a copy of the flux to use in calculating Toeplitz
    var flux_copy = allocator.dupe(f64, flux) catch |e| {
        std.debug.print("kerrlineprofile: error: {!}\n", .{e});
        @panic("kerrlineprofile: fatal: COULD NOT DUPE FLUX ARRAY.");
    };
    defer allocator.free(flux_copy);

    // convolve and write result directly into output array
    convolution.convolve(f64, flux, conv_energy, conv_flux, energy, flux_copy);
}

// model paramter definitions

fn Parameters(comptime T: type) type {
    return struct {
        const Self = @This();
        a: T,
        inclination: T,
        eline: T,
        alpha: T,
        rmin: T,
        rmax: T,
        pub fn from_ptr(ptr: *const f64) Parameters(T) {
            const N = 6;
            var slice = @as([*]const f64, @ptrCast(ptr))[0..N];
            return .{
                .a = @as(T, @floatCast(slice[0])),
                // convert to cos(i)
                .inclination = @cos(std.math.degreesToRadians(
                    T,
                    @as(T, @floatCast(slice[1])),
                )),
                .eline = @as(T, @floatCast(slice[2])),
                .alpha = @as(T, @floatCast(slice[3])),
                .rmin = @as(T, @floatCast(slice[4])),
                .rmax = @as(T, @floatCast(slice[5])),
            };
        }
        pub fn from_ptr_conv(ptr: *const f64) Parameters(T) {
            const N = 5;
            var slice = @as([*]const f64, @ptrCast(ptr))[0..N];
            return .{
                .a = @as(T, @floatCast(slice[0])),
                // convert to cos(i)
                .inclination = @cos(std.math.degreesToRadians(
                    T,
                    @as(T, @floatCast(slice[1])),
                )),
                // fixed for convolution models
                .eline = 1,
                .alpha = @as(T, @floatCast(slice[2])),
                .rmin = @as(T, @floatCast(slice[3])),
                .rmax = @as(T, @floatCast(slice[4])),
            };
        }
        pub fn table_parameters(self: Self) [NPARAMS]T {
            return [_]T{ self.a, self.inclination };
        }
    };
}

fn LinEmisParameters(comptime T: type, comptime Nbins: comptime_int) type {
    return struct {
        const Self = @This();
        a: T,
        inclination: T,
        eline: T,
        rmin: T,
        rmax: T,
        rcutoff: T,
        alpha: T,
        weights: [Nbins]T,

        fn read_weights(slice: []const f64, first_index: usize) [Nbins]T {
            // read in the emissivity weights
            var weights: [Nbins]T = undefined;
            for (slice[first_index..], 0..) |w, i| {
                weights[i] = @as(T, @floatCast(w));
            }
            return weights;
        }

        pub fn from_ptr(ptr: *const f64) Self {
            const N = 7 + Nbins;
            var slice = @as([*]const f64, @ptrCast(ptr))[0..N];
            return .{
                .a = @as(T, @floatCast(slice[0])),
                // convert to cos(i)
                .inclination = @cos(std.math.degreesToRadians(
                    T,
                    @as(T, @floatCast(slice[1])),
                )),
                .eline = @as(T, @floatCast(slice[2])),
                .rmin = @as(T, @floatCast(slice[3])),
                .rmax = @as(T, @floatCast(slice[4])),
                .rcutoff = @as(T, @floatCast(slice[5])),
                .alpha = @as(T, @floatCast(slice[6])),
                .weights = read_weights(slice, 7),
            };
        }

        pub fn from_ptr_conv(ptr: *const f64) Self {
            const N = 6 + Nbins;
            var slice = @as([*]const f64, @ptrCast(ptr))[0..N];
            return .{
                .a = @as(T, @floatCast(slice[0])),
                // convert to cos(i)
                .inclination = @cos(std.math.degreesToRadians(
                    T,
                    @as(T, @floatCast(slice[1])),
                )),
                .eline = 1,
                .rmin = @as(T, @floatCast(slice[2])),
                .rmax = @as(T, @floatCast(slice[3])),
                .rcutoff = @as(T, @floatCast(slice[4])),
                .alpha = @as(T, @floatCast(slice[5])),
                .weights = read_weights(slice, 6),
            };
        }

        pub fn table_parameters(self: Self) [NPARAMS]T {
            return [_]T{ self.a, self.inclination };
        }
    };
}

// model dispatchers

inline fn kerr_lin_emisN(
    comptime Nemis: comptime_int,
    energy_ptr: *const f64,
    n_flux: c_int,
    parameters_ptr: *const f64,
    spectrum: c_int,
    flux_ptr: *f64,
    flux_variance_ptr: *f64,
    init_ptr: *const u8,
) void {
    // unused
    _ = spectrum;
    _ = flux_variance_ptr;
    _ = init_ptr;

    const N = @as(usize, @intCast(n_flux));
    // convert to slices
    const energy = @as([*]const f64, @ptrCast(energy_ptr))[0 .. N + 1];
    var flux = @as([*]f64, @ptrCast(flux_ptr))[0..N];

    const parameters = LinEmisParameters(f32, Nemis).from_ptr(parameters_ptr);
    const lin_emis = emissivity.LinInterpEmissivity(f32, Nemis).init(
        parameters.weights,
        parameters.rmin,
        parameters.rmax,
        parameters.alpha,
    );
    integrate_lineprofile(f32, energy, flux, parameters, lin_emis);
}

inline fn kerr_conv_emisN(
    comptime Nemis: comptime_int,
    energy_ptr: *const f64,
    n_flux: c_int,
    parameters_ptr: *const f64,
    spectrum: c_int,
    flux_ptr: *f64,
    flux_variance_ptr: *f64,
    init_ptr: *const u8,
) void {
    // unused
    _ = spectrum;
    _ = flux_variance_ptr;
    _ = init_ptr;

    const N = @as(usize, @intCast(n_flux));
    // convert to slices
    const energy = @as([*]const f64, @ptrCast(energy_ptr))[0 .. N + 1];
    var flux = @as([*]f64, @ptrCast(flux_ptr))[0..N];

    const parameters = LinEmisParameters(f32, Nemis).from_ptr_conv(parameters_ptr);
    const lin_emis = emissivity.LinInterpEmissivity(f32, Nemis).init(
        parameters.weights,
        parameters.rmin,
        parameters.rmax,
        parameters.alpha,
    );

    kerr_convolve(energy, flux, parameters, lin_emis);
}

// XSPEC API

pub export fn kline(
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
    _ = spectrum;
    _ = flux_variance_ptr;
    _ = init_ptr;

    const N = @as(usize, @intCast(n_flux));
    // convert to slices
    const energy = @as([*]const f64, @ptrCast(energy_ptr))[0 .. N + 1];
    var flux = @as([*]f64, @ptrCast(flux_ptr))[0..N];

    const parameters = Parameters(f32).from_ptr(parameters_ptr);
    const fixed_emis = emissivity.PowerLawEmissivity(f32).init(-parameters.alpha);

    integrate_lineprofile(f32, energy, flux, parameters, fixed_emis);
}

pub export fn kconv(
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
    _ = spectrum;
    _ = flux_variance_ptr;
    _ = init_ptr;

    const N = @as(usize, @intCast(n_flux));
    // convert to slices
    const energy = @as([*]const f64, @ptrCast(energy_ptr))[0 .. N + 1];
    var flux = @as([*]f64, @ptrCast(flux_ptr))[0..N];

    const parameters = Parameters(f32).from_ptr_conv(parameters_ptr);
    const fixed_emis = emissivity.PowerLawEmissivity(f32).init(-parameters.alpha);

    kerr_convolve(energy, flux, parameters, fixed_emis);
}

pub export fn kline5(
    energy_ptr: *const f64,
    n_flux: c_int,
    parameters_ptr: *const f64,
    spectrum: c_int,
    flux_ptr: *f64,
    flux_variance_ptr: *f64,
    init_ptr: *const u8,
) callconv(.C) void {
    return kerr_lin_emisN(5, energy_ptr, n_flux, parameters_ptr, spectrum, flux_ptr, flux_variance_ptr, init_ptr);
}

pub export fn kconv5(
    energy_ptr: *const f64,
    n_flux: c_int,
    parameters_ptr: *const f64,
    spectrum: c_int,
    flux_ptr: *f64,
    flux_variance_ptr: *f64,
    init_ptr: *const u8,
) callconv(.C) void {
    return kerr_conv_emisN(5, energy_ptr, n_flux, parameters_ptr, spectrum, flux_ptr, flux_variance_ptr, init_ptr);
}
