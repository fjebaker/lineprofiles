const std = @import("std");
const io = @import("io.zig");

const xfunc = @import("transfer-functions.zig");
const util = @import("utils.zig");

pub fn printPlot(comptime T: type, tfs: []*const xfunc.TransferFunction(T)) !void {
    const CHOICE = 80;
    const stream = std.io.getStdOut().writer();
    const n = tfs[0].lower_branch[0].len;

    const NUM_FMT = "{d: <8.6}\t";

    var grid = util.relline_gstar_grid(f32);
    for (0..n) |i| {
        const x = grid[i];
        try stream.print(NUM_FMT, .{x});
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
    var data = try io.readFitsFile("../relline/rel_table.fits", allocator);
    defer data.deinit();

    var params = [2]f32{ 0.998, 0.276 };
    var itf = data.interpolate_parameters(params);
    itf.stage_radius(1.509);

    try printPlot(
        f32,
        &[_]*const xfunc.TransferFunction(f32){
            itf.tf,
            &data.transfer_functions[588],
        },
    );
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    var alloc = gpa.allocator();

    try interpolation_test(alloc);
}
