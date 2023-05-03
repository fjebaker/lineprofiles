const std = @import("std");
const io = @import("io.zig");

const xfunc = @import("transfer-functions.zig");
const util = @import("utils.zig");
const emissivity = @import("emissivity.zig");

const TABLE_FILE = "kerr-transfer-functions.fits";

const CHOICE = 52;
const NUM_FMT = "{d: <8.6}\t";

pub fn printPlot(
    comptime T: type,
    tfs: []*const xfunc.TransferFunction(T),
    itfs: []*const xfunc.InterpolatingTransferFunction(T),
) !void {
    const stream = std.io.getStdOut().writer();
    const n = if (tfs.len > 0)
        tfs[0].lower_branch[0].len
    else
        itfs[0].tf.lower_branch[0].len;

    var grid = util.relline_gstar_grid(f32);
    for (0..n) |i| {
        const x = grid[i];
        try stream.print(NUM_FMT, .{x});
        for (itfs) |tf| {
            const y1 = tf.cache_lower[i];
            const y2 = tf.cache_upper[i];
            try stream.print(NUM_FMT, .{y1});
            try stream.print(NUM_FMT, .{y2});
        }
        for (tfs) |tf| {
            const y1 = tf.lower_branch[CHOICE][i];
            const y2 = tf.upper_branch[CHOICE][i];
            try stream.print(NUM_FMT, .{y1});
            try stream.print(NUM_FMT, .{y2});
        }
        try stream.print("\n", .{});
    }
}

pub fn interpolation_test(allocator: std.mem.Allocator) !void {
    // std.debug.print("Starting read...\n", .{});
    var data = try io.readFitsFile(2, f32, TABLE_FILE, allocator);
    defer data.deinit();
    // std.debug.print("Done.\n", .{});

    var params = [2]f32{ 0.99, 0.421 };
    var itf = data.interpolate_parameters(params);
    // data.interpolated_cache[0].assignFrom(&data.transfer_functions[428]);
    // var itf = data.interpolated_transfer_function;

    itf.stage_radius(4.616);

    try printPlot(
        f32,
        &[_]*const xfunc.TransferFunction(f32){
            &data.transfer_functions[428],
        },
        &[_]*const xfunc.InterpolatingTransferFunction(f32){
            &itf,
        },
    );
}

pub fn integrate(allocator: std.mem.Allocator, filepath: [:0]const u8) !void {
    const stream = std.io.getStdOut().writer();

    std.debug.print("Starting read...\n", .{});
    var data = try io.readFitsFile(2, f32, 30, filepath, allocator);
    defer data.deinit();
    var data2 = try io.readFitsFile(2, f32, 20, "relline.fits", allocator);
    defer data2.deinit();
    std.debug.print("Done.\n", .{});

    // build r grid
    var r_grid = try util.inverse_grid(f32, allocator, 1.0, 50.0, 2000);
    defer allocator.free(r_grid);

    // build g grid
    var gitt = util.RangeIterator(f32).init(0.0, 1.8, 200);
    var g_grid = try gitt.drain(allocator);
    defer allocator.free(g_grid);

    // build refined grid
    const refinement = 5;
    var fine_grid = try allocator.alloc(f32, refinement * (g_grid.len - 1));
    defer allocator.free(fine_grid);
    util.refine_grid(f32, g_grid, fine_grid, refinement, 1.0);

    // allocate output
    var flux = try allocator.alloc(f32, g_grid.len - 1);
    defer allocator.free(flux);
    var fflux1 = try allocator.alloc(f32, fine_grid.len - 1);
    defer allocator.free(fflux1);
    // and zero it
    for (flux) |*f| f.* = 0;
    for (fflux1) |*f| f.* = 0;

    var flux2 = try allocator.alloc(f32, g_grid.len - 1);
    defer allocator.free(flux2);
    var fflux2 = try allocator.alloc(f32, fine_grid.len - 1);
    defer allocator.free(fflux2);
    // and zero it
    for (flux2) |*f| f.* = 0;
    for (fflux2) |*f| f.* = 0;

    var params = [2]f32{ 0.998, 0.90 };
    var itf = data.interpolate_parameters(params);
    var itf2 = data2.interpolate_parameters(params);

    var fixed_emis = emissivity.StepFunctionEmissivity(f32, 5).init([_]f32{ 1, 0, 0, 0, 0 }, 1.0, 50.0);
    std.debug.print("bins: {any}\n", .{fixed_emis.cutoffs});
    const emis = fixed_emis.emissivity();

    itf.integrate(r_grid, fine_grid, fflux1, emis);
    itf2.integrate(r_grid, fine_grid, fflux2, emis);

    var j: usize = 0;
    for (0..flux.len) |i| {
        for (0..refinement) |_| {
            flux[i] += fflux1[j];
            flux2[i] += fflux2[j];
            j += 1;
            if (j == fflux1.len) break;
        }
        // normalize to counts per bin
        const e_midpoint = 0.5 * (g_grid[i + 1] + g_grid[i]);
        flux[i] = flux[i] / e_midpoint;
        flux2[i] = flux2[i] / e_midpoint;
    }

    util.normalize(f32, flux);
    util.normalize(f32, flux2);

    for (0..flux.len) |i| {
        try stream.print("{d}\t{d}\t{d}\n", .{ g_grid[i], flux[i], flux2[i] });
    }
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var alloc = gpa.allocator();

    var args = try std.process.argsAlloc(alloc);
    defer std.process.argsFree(alloc, args);

    var filepath = if (args.len < 2)
        "kerr-transfer-functions.fits"
    else
        args[1];

    std.debug.print("filepath: {s}\n", .{filepath});

    // try interpolation_test(alloc);
    try integrate(alloc, filepath);
}
