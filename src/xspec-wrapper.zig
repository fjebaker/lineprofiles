const std = @import("std");
const io = @import("io.zig");
const util = @import("utils.zig");
const lineprofile = @import("line-profile.zig");
const emissivity = @import("emissivity.zig");
const convolution = @import("convolution.zig");

// number of parameters in the table
const NPARAMS = 2;
const DATA_DIR_ENV_VAR = "KLINE_PROF_DATA_DIR";
const MODELPATH = "kerr-transfer-functions.fits";
const MSG_PREFIX = "kerrlineprofile: ";
const REFINEMENT = 4;

const FINE_G_BINS = 2000;
const RADIAL_BINS = 1500;

fn debugPrint(comptime fmt: []const u8, args: anytype) void {
    if (!verbose) return;

    var out_writer = std.fs.File.stdout().writer(&.{});
    const out = &out_writer.interface;

    out.print(fmt, args) catch |err| {
        // try to debug log, else no worries
        std.debug.print(MSG_PREFIX ++ "Error in debug printing: {any}\n", .{err});
    };
}

pub fn getModelPath() ![]const u8 {
    return data_file orelse {
        const root = std.process.getEnvVarOwned(
            allocator,
            DATA_DIR_ENV_VAR,
        ) catch |err|
            if (err == error.EnvironmentVariableNotFound)
                try allocator.dupe(u8, ".")
            else
                return err;

        defer allocator.free(root);
        data_file = try std.fs.path.join(allocator, &.{ root, MODELPATH });

        debugPrint("Using table path '{s}'\n", .{data_file.?});

        return data_file.?;
    };
}

var verbose = true;

// convolution has a fixed bin width for the model
// that is to say, the model is evaluated on a grid
// of g values in [0,2], where the bin width is:
pub const CONVOLUTION_RESOLUTION = 1e-3;

var data_file: ?[]const u8 = null;
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

    debugPrint(
        "Read in {d} transfer function tables.\n",
        .{profile.?.transfer_functions.len},
    );

    // build r grid
    r_grid = try util.inverse_grid(f32, allocator, 1.0, 50.0, RADIAL_BINS);
    errdefer allocator.free(r_grid);

    // build some g grid, it will be refined anyway
    var gitt = util.RangeIterator(f32).init(0.1, 1.5, 200);
    g_grid = try gitt.drain(allocator);
    errdefer allocator.free(g_grid);

    // allocate output
    flux_cache = try allocator.alloc(f32, g_grid.len - 1);

    debugPrint("Finished one-time setup.\n", .{});
}

fn convolve_setup() !void {
    // build the convolution energy grid
    var itt = util.RangeIterator(f64).init(
        CONVOLUTION_RESOLUTION,
        2.0,
        @trunc((2.0 - CONVOLUTION_RESOLUTION) / CONVOLUTION_RESOLUTION),
    );
    const conv_energy = try itt.drain(allocator);
    errdefer allocator.free(conv_energy);

    // temporary flux array
    const model_flux = try allocator.alloc(f64, conv_energy.len - 1);
    errdefer allocator.free(model_flux);

    // assign to singletons
    conv_g_grid = conv_energy;
    conv_flux = model_flux;
}

fn refine_grid(comptime T: type, grid: anytype, norm: T) !void {
    // don't do more than 2000
    if (FINE_G_BINS != g_grid.len) {
        allocator.free(g_grid);
        g_grid = try allocator.alloc(T, FINE_G_BINS);
        allocator.free(flux_cache);
        flux_cache = try allocator.alloc(T, FINE_G_BINS - 1);
    }
    var itt = util.RangeIterator(f32).init(@floatCast(grid[0]), @floatCast(grid[grid.len - 1]), FINE_G_BINS);
    var i: usize = 0;
    const inorm = 1 / norm;
    while (itt.next()) |g| {
        g_grid[i] = g * inorm;
        i += 1;
    }
}

fn integrate_lineprofile(
    comptime T: type,
    energy: []const f64,
    flux: []f64,
    params: anytype,
    emis: anytype,
) !void {
    // do we need to do first time setup?
    var lp = profile orelse blk: {
        try setup();
        break :blk profile.?;
    };

    // refine our grids
    try refine_grid(T, energy, params.eline);
    set_radial_grid(
        T,
        @max(r_grid[0], params.rmin),
        @min(params.rmax, r_grid[r_grid.len - 1]),
    );

    // zero the flux cache
    for (flux_cache) |*f| f.* = 0;

    // do the integration
    var itf = lp.interpolate_parameters(params.table_parameters());
    itf.integrate(r_grid, g_grid, flux_cache, emis);

    const inorm = 1 / params.eline;

    // rebin for output
    var j: usize = 0;
    var total_flux: f64 = 0;
    for (0..flux.len) |i| {
        // ensure it is zeroed
        flux[i] = 0;
        // sum up the grid values
        var summed_flux: f32 = 0;
        while (g_grid[j] < energy[i] * inorm) : (j += 1) {
            summed_flux += flux_cache[j];
        }
        flux[i] += @floatCast(summed_flux);
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
) !void {
    const conv_energy = conv_g_grid orelse blk: {
        try convolve_setup();
        break :blk conv_g_grid.?;
    };

    // invoke the model
    try integrate_lineprofile(f32, conv_energy, conv_flux, parameters, emis);

    // make a copy of the flux to use in calculating Toeplitz
    const flux_copy = try allocator.dupe(f64, flux);
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
            const slice = @as([*]const f64, @ptrCast(ptr))[0..N];
            return .{
                .a = @as(T, @floatCast(slice[0])),
                // convert to cos(i)
                .inclination = @cos(std.math.degreesToRadians(
                    @as(T, @floatCast(slice[1])),
                )),
                .eline = @as(T, @floatCast(slice[2])),
                .rmin = @as(T, @floatCast(slice[3])),
                .rmax = @as(T, @floatCast(slice[4])),
                .alpha = @as(T, @floatCast(slice[5])),
            };
        }
        pub fn from_ptr_conv(ptr: *const f64) Parameters(T) {
            const N = 5;
            const slice = @as([*]const f64, @ptrCast(ptr))[0..N];
            return .{
                .a = @as(T, @floatCast(slice[0])),
                // convert to cos(i)
                .inclination = @cos(std.math.degreesToRadians(
                    @as(T, @floatCast(slice[1])),
                )),
                // fixed for convolution models
                .eline = 1,
                .rmin = @as(T, @floatCast(slice[2])),
                .rmax = @as(T, @floatCast(slice[3])),
                .alpha = @as(T, @floatCast(slice[4])),
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
        rcut: T,
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
            const slice = @as([*]const f64, @ptrCast(ptr))[0..N];
            return .{
                .a = @as(T, @floatCast(slice[0])),
                // convert to cos(i)
                .inclination = @cos(std.math.degreesToRadians(
                    @as(T, @floatCast(slice[1])),
                )),
                .eline = @as(T, @floatCast(slice[2])),
                .rmin = @as(T, @floatCast(slice[3])),
                .rmax = @as(T, @floatCast(slice[4])),
                .rcut = @as(T, @floatCast(slice[5])),
                .alpha = @as(T, @floatCast(slice[6])),
                .weights = read_weights(slice, 7),
            };
        }

        pub fn from_ptr_conv(ptr: *const f64) Self {
            const N = 6 + Nbins;
            const slice = @as([*]const f64, @ptrCast(ptr))[0..N];
            return .{
                .a = @as(T, @floatCast(slice[0])),
                // convert to cos(i)
                .inclination = @cos(std.math.degreesToRadians(
                    @as(T, @floatCast(slice[1])),
                )),
                .eline = 1,
                .rmin = @as(T, @floatCast(slice[2])),
                .rmax = @as(T, @floatCast(slice[3])),
                .rcut = @as(T, @floatCast(slice[4])),
                .alpha = @as(T, @floatCast(slice[5])),
                .weights = read_weights(slice, 6),
            };
        }

        pub fn table_parameters(self: Self) [NPARAMS]T {
            return [_]T{ self.a, self.inclination };
        }

        // pub fn checkParameters(self: Self) void {
        //     self.a = std.math.clamp(self.a, -0.998, 0.998);
        //     self.inclination = std.math.clamp(self.a, 2, 88);
        //     self.rmin = std.math.clamp(self.rmin, 2, 88);
        // }
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
    const flux = @as([*]f64, @ptrCast(flux_ptr))[0..N];

    const parameters = if (Nemis == 1)
        Parameters(f32).from_ptr(parameters_ptr)
    else
        LinEmisParameters(f32, Nemis).from_ptr(parameters_ptr);

    const emis = if (Nemis == 1)
        emissivity.PowerLawEmissivity(f32).init(-parameters.alpha)
    else
        emissivity.LinInterpEmissivity(f32, Nemis).init(
            parameters.weights,
            parameters.rmin,
            parameters.rcut,
            parameters.alpha,
        );

    integrate_lineprofile(f32, energy, flux, parameters, emis) catch |e| {
        if (@errorReturnTrace()) |t| std.debug.dumpStackTrace(t.*);
        switch (e) {
            error.FileNotFound => std.log.err(
                "{s}Failed to find table model. Set the {s} environment variable to where {s} is located",
                .{ MSG_PREFIX, DATA_DIR_ENV_VAR, MODELPATH },
            ),
            else => std.log.err("{s}{any}", .{ MSG_PREFIX, e }),
        }
        std.log.err("{s}MODEL OUTPUT SET TO ZERO", .{MSG_PREFIX});
        @memset(flux, 0);
    };
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
    const flux = @as([*]f64, @ptrCast(flux_ptr))[0..N];

    const parameters = if (Nemis == 1)
        Parameters(f32).from_ptr_conv(parameters_ptr)
    else
        LinEmisParameters(f32, Nemis).from_ptr_conv(parameters_ptr);

    const emis = if (Nemis == 1)
        emissivity.PowerLawEmissivity(f32).init(-parameters.alpha)
    else
        emissivity.LinInterpEmissivity(f32, Nemis).init(
            parameters.weights,
            parameters.rmin,
            parameters.rcut,
            parameters.alpha,
        );

    kerr_convolve(energy, flux, parameters, emis) catch |e| {
        switch (e) {
            error.FileNotFound => std.log.err(
                "{s}Failed to find table model. Set the {s} environment variable to where {s} is located",
                .{ MSG_PREFIX, DATA_DIR_ENV_VAR, MODELPATH },
            ),
            else => std.log.err("{s} {any}", .{ MSG_PREFIX, e }),
        }
        std.log.err("{s}MODEL OUTPUT SET TO ZERO", .{MSG_PREFIX});
        @memset(flux, 0);
    };
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
) callconv(.c) void {
    return kerr_lin_emisN(
        1,
        energy_ptr,
        n_flux,
        parameters_ptr,
        spectrum,
        flux_ptr,
        flux_variance_ptr,
        init_ptr,
    );
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
) callconv(.c) void {
    return kerr_conv_emisN(
        1,
        energy_ptr,
        n_flux,
        parameters_ptr,
        spectrum,
        flux_ptr,
        flux_variance_ptr,
        init_ptr,
    );
}

pub export fn kline5(
    energy_ptr: *const f64,
    n_flux: c_int,
    parameters_ptr: *const f64,
    spectrum: c_int,
    flux_ptr: *f64,
    flux_variance_ptr: *f64,
    init_ptr: *const u8,
) callconv(.c) void {
    return kerr_lin_emisN(
        5,
        energy_ptr,
        n_flux,
        parameters_ptr,
        spectrum,
        flux_ptr,
        flux_variance_ptr,
        init_ptr,
    );
}

pub export fn kconv5(
    energy_ptr: *const f64,
    n_flux: c_int,
    parameters_ptr: *const f64,
    spectrum: c_int,
    flux_ptr: *f64,
    flux_variance_ptr: *f64,
    init_ptr: *const u8,
) callconv(.c) void {
    return kerr_conv_emisN(
        5,
        energy_ptr,
        n_flux,
        parameters_ptr,
        spectrum,
        flux_ptr,
        flux_variance_ptr,
        init_ptr,
    );
}

pub export fn kconv10(
    energy_ptr: *const f64,
    n_flux: c_int,
    parameters_ptr: *const f64,
    spectrum: c_int,
    flux_ptr: *f64,
    flux_variance_ptr: *f64,
    init_ptr: *const u8,
) callconv(.c) void {
    return kerr_conv_emisN(
        10,
        energy_ptr,
        n_flux,
        parameters_ptr,
        spectrum,
        flux_ptr,
        flux_variance_ptr,
        init_ptr,
    );
}

fn smokeTestModel(domain: []const f64, comptime f: anytype, params: []const f64) !void {
    const output = try std.testing.allocator.alloc(f64, domain.len - 1);
    defer std.testing.allocator.free(output);
    var f_var: f64 = 0;
    const init_ptr = "";
    f(
        @ptrCast(domain.ptr),
        @intCast(output.len),
        @ptrCast(params.ptr),
        0,
        @ptrCast(output.ptr),
        &f_var,
        @ptrCast(init_ptr),
    );
}

fn exampleDomain(alloc: std.mem.Allocator, size: usize) ![]f64 {
    var itt = util.RangeIterator(f64).init(0.1, 10, size);
    return try itt.drain(alloc);
}

test "smoke test all" {
    verbose = false;
    const domain = try exampleDomain(std.testing.allocator, 100);
    defer std.testing.allocator.free(domain);

    try smokeTestModel(domain, kline, &[_]f64{ 0.998, 60, 6.4, 3.0, 1.0, 400.0 });
    try smokeTestModel(domain, kline, &[_]f64{ 0.998, 0.1, 6.4, 3.0, 1.0, 400.0 });
    try smokeTestModel(domain, kline, &[_]f64{ 0.998, 89.0, 6.4, 3.0, 1.0, 400.0 });
    try smokeTestModel(domain, kline, &[_]f64{ 0.998, 60.0, 6.4, 0.0, 1.0, 400.0 });
    try smokeTestModel(domain, kline, &[_]f64{ 0.998, 60.0, 6.4, 0.0, 0.0, 400.0 });

    try smokeTestModel(domain, kconv5, &[_]f64{ 0.998, 60.0, 1.0, 400.0, 100.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0 });
    // with zero emissivities
    try smokeTestModel(domain, kconv5, &[_]f64{ 0.998, 60.0, 1.0, 400.0, 100.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0 });
    try smokeTestModel(domain, kconv10, &[_]f64{
        0.998,
        60.0,
        1.0,
        400.0,
        50.0,
        3.0,
        3.0,
        3.0,
        3.0,
        3.0,
        3.0,
        3.0,
        3.0,
        3.0,
        3.0,
        3.0,
    });
}
