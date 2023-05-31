const std = @import("std");
const io = @import("io.zig");
const util = @import("utils.zig");
const lineprofile = @import("line-profile.zig");
const emissivity = @import("emissivity.zig");

// number of parameters in the table
const NPARAMS = 2;
const MODELPATH = "./kerr-transfer-functions.fits";

// refinement for the energy grid
const REFINEMENT = 5;

var allocator = std.heap.c_allocator;

// singleton
var profile: ?lineprofile.LineProfileTable(NPARAMS, f32) = null;
var g_grid: []f32 = undefined;
var r_grid: []f32 = undefined;
var flux_cache: []f32 = undefined;

fn set_radial_grid(comptime T: type, rmin: T, rmax: T) void {
    util.inverse_grid_inplace(T, r_grid, rmin, rmax);
}

fn setup() !void {
    profile = try io.readFitsFile(NPARAMS, f32, 30, MODELPATH, allocator);
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
            flux[i] += @floatCast(f64, flux_cache[j]);
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

const CONVOLUTION_RESOLUTION = 1e-3;
fn _kerr_convolve(
    comptime T: type,
    energy: []const T,
    flux: []T,
    parameters: anytype,
    emis: anytype,
) !void {
    // build the convolution energy grid
    var itt = util.RangeIterator(f64).init(
        CONVOLUTION_RESOLUTION,
        2.0,
        @trunc((2.0 - CONVOLUTION_RESOLUTION) / CONVOLUTION_RESOLUTION),
    );
    var conv_energy = try itt.drain(allocator);
    defer allocator.free(conv_energy);

    // temporary flux array
    var model_flux = try allocator.alloc(f64, conv_energy.len - 1);
    defer allocator.free(model_flux);

    // invoke the model
    integrate_lineprofile(f32, conv_energy, model_flux, parameters, emis);

    // shift energy into g
    var gs = try allocator.alloc(f64, energy.len);
    defer allocator.free(gs);

    var flux_copy = try allocator.dupe(f64, flux);
    defer allocator.free(flux_copy);

    _convolve(f64, flux, conv_energy, model_flux, energy, flux_copy);
}

fn kerr_convolve(
    comptime T: type,
    energy: []const T,
    flux: []T,
    parameters: anytype,
    emis: anytype,
) void {
    _kerr_convolve(T, energy, flux, parameters, emis) catch |e| {
        std.debug.print("kerrlineprofile: ERR: {any}\n", .{e});
        @panic("kerrlineprofile: fatal error");
    };
}

/// convolves a and b on domains g and x respectively by rebinning
/// if necessary, assumes output is on domain x
pub fn _convolve(
    comptime T: type,
    output: []T,
    g: []const T,
    a: []const T,
    x: []const T,
    b: []const T,
) void {
    std.debug.assert(output.len == b.len);

    // find window of non-zero
    const window_start = blk: {
        for (0..a.len) |i| if (a[i] > 0)
            break :blk if (i > 0) i - 1 else i;
        break :blk 0;
    };
    const window_end = blk: {
        for (0..a.len) |i| {
            const j = a.len - i - 1;
            if (a[j] > 0)
                break :blk if (j < a.len - 1) j + 1 else j;
        }
        break :blk 0;
    };

    const g_min = g[window_start];
    const g_max = g[window_end];

    for (0..output.len) |i| {
        output[i] = 0;
        const avg = 2 / (x[i + 1] + x[i]);
        for (0..b.len) |j| {
            // output bin extremes
            const low = x[j] * avg;
            const high = x[j + 1] * avg;

            // skip if outside of the window
            if (high < g_min or low > g_max) continue;

            // integrate the window
            for (window_start..window_end) |w| {
                const bin_low = g[w];
                const bin_high = g[w + 1];

                // check different cases
                // 1. check if bin is not in output bin
                //    bin_low |---| bin_high     low |---| high
                if (bin_low > high or bin_high < low) {
                    continue;
                }
                // 2. check if bin is striding
                //    bin_low |---| low |.....| bin_high |---| high
                // or
                //    low |---| bin_low |.....| high |---| bin_high
                else if ((bin_low < low and bin_high < high) or
                    (bin_low > low and bin_high > high))
                {
                    // integrate stride
                    const overlap = @max(low - bin_low, bin_high - high);
                    output[i] += a[w] * b[j] * overlap;
                }
                // 3. check if bin is bigger than output bin
                //    bin_low |...| low |.....| high |...| bin_high
                // 4. bin must be contained
                //    low |---| bin_low |.....| bin_high |---| high
                else {
                    output[i] += a[w] * b[j] * CONVOLUTION_RESOLUTION;
                }
            }
        }
    }
}

// model definitions

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

pub export fn kerr_line_profile(
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

    const N = @intCast(usize, n_flux);
    // convert to slices
    const energy = @ptrCast([*]const f64, energy_ptr)[0 .. N + 1];
    var flux = @ptrCast([*]f64, flux_ptr)[0..N];

    const parameters = Parameters(f32).from_ptr(parameters_ptr);
    const fixed_emis = emissivity.PowerLawEmissivity(f32).init(-3);

    integrate_lineprofile(f32, energy, flux, parameters, fixed_emis);
}

pub export fn kerr_conv_profile(
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

    const N = @intCast(usize, n_flux);
    // convert to slices
    const energy = @ptrCast([*]const f64, energy_ptr)[0 .. N + 1];
    var flux = @ptrCast([*]f64, flux_ptr)[0..N];

    const parameters = Parameters(f32).from_ptr(parameters_ptr);
    const fixed_emis = emissivity.PowerLawEmissivity(f32).init(-3);

    kerr_convolve(f64, energy, flux, parameters, fixed_emis);
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

        pub fn from_ptr(ptr: *const f64) Self {
            const N = 7 + Nbins;
            var slice = @ptrCast([*]const f64, ptr)[0..N];

            // read in the emissivity weights
            var weights: [Nbins]T = undefined;
            for (slice[7..], 0..) |w, i| {
                weights[i] = @floatCast(T, w);
            }

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
                .rcutoff = @floatCast(T, slice[5]),
                .alpha = @floatCast(T, slice[6]),
                .weights = weights,
            };
        }
        pub fn table_parameters(self: Self) [NPARAMS]T {
            return [_]T{ self.a, self.inclination };
        }
    };
}

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

    const N = @intCast(usize, n_flux);
    // convert to slices
    const energy = @ptrCast([*]const f64, energy_ptr)[0 .. N + 1];
    var flux = @ptrCast([*]f64, flux_ptr)[0..N];

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

    const N = @intCast(usize, n_flux);
    // convert to slices
    const energy = @ptrCast([*]const f64, energy_ptr)[0 .. N + 1];
    var flux = @ptrCast([*]f64, flux_ptr)[0..N];

    const parameters = LinEmisParameters(f32, Nemis).from_ptr(parameters_ptr);
    const lin_emis = emissivity.LinInterpEmissivity(f32, Nemis).init(parameters.weights, parameters.rmin, parameters.rmax, parameters.alpha);

    kerr_convolve(f64, energy, flux, parameters, lin_emis);
}

pub export fn kerr_lin_emis5(
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

pub export fn kerr_conv_emis5(
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
