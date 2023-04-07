const std = @import("std");
const io = @import("io.zig");

const xfunc = @import("transfer-functions.zig");
const util = @import("utils.zig");

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
    std.debug.print("Starting read...\n", .{});
    var data = try io.readFitsFile("../relline/rel_table.fits", allocator);
    defer data.deinit();
    std.debug.print("Done.\n", .{});

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

pub fn integrate(allocator: std.mem.Allocator) !void {
    const stream = std.io.getStdOut().writer();

    std.debug.print("Starting read...\n", .{});
    var data = try io.readFitsFile("../relline/rel_table.fits", allocator);
    defer data.deinit();
    std.debug.print("Done.\n", .{});

    // build r grid
    var ritt = util.RangeIterator(f32).init(3.0, 50.0, 1000);
    var r_grid = try ritt.drain(allocator);
    defer allocator.free(r_grid);

    // build g grid
    var gitt = util.RangeIterator(f32).init(0.1, 1.5, 100);
    var g_grid = try gitt.drain(allocator);
    defer allocator.free(g_grid);

    // allocate output
    var flux = try allocator.alloc(f32, g_grid.len - 1);
    defer allocator.free(flux);
    // and zero it
    for (flux) |*f| f.* = 0;

    var params = [2]f32{ 0.99, 0.5 };
    var itf = data.interpolate_parameters(params);

    itf.integrate(r_grid, g_grid, flux);

    for (0..flux.len) |i| {
        try stream.print("{d}\t{d}\n", .{ g_grid[i], flux[i] });
    }
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var alloc = gpa.allocator();

    // try interpolation_test(alloc);
    try integrate(alloc);
}
