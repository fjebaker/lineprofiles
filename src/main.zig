const std = @import("std");
const io = @import("io.zig");

const util = @import("utils.zig");

pub fn print_for_plotting(allocator: std.mem.Allocator) !void {
    // var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    // defer _ = gpa.deinit();
    // var alloc = gpa.allocator();

    var data = try io.readFitsFile("../relline/rel_table.fits", allocator);

    defer data.deinit();
    var itf = data.interpolated_transfer_function;
    itf.stage_radius(10.0);

    for (itf.cache_lower, 0..) |l, i| {
        std.debug.print("{d}, {d}\n", .{ i, l });
    }
}

test {
    _ = util;
    var allocator = std.testing.allocator;
    try print_for_plotting(allocator);
    // var data = io.readFitsFile("../relline/rel_table.fits", allocator) catch |err| {
    //     std.debug.print("{any}\n", .{err});
    //     return err;
    // };
    // defer data.deinit();

    // var itf = data.interpolated_transfer_function;
    // std.debug.print("{any}\n", .{itf.radii});
    // std.debug.print("{any}\n", .{itf.cache_lower});
    // std.debug.print("INTERPOLATING\n", .{});
    // itf.stage_radius(10.0);
    // std.debug.print("{any}\n", .{itf.cache_lower});
    // std.debug.print("N cache {d}, N branch {d}\n", .{ itf.cache_lower.len, itf.tf.lower_branch[0].len });
    // std.debug.print("{any}\n", .{itf.tf.lower_branch[20]});
}
