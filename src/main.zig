const std = @import("std");
const io = @import("io.zig");
const util = @import("utils.zig");
const tf = @import("transfer-functions.zig");
const lp = @import("line-profile.zig");

test "main" {
    _ = tf;
    _ = lp;
}

test "index-lookups" {
    var lineprof = io.readFitsFile(2, f32, 30, "./kerr-transfer-functions.fits", std.testing.allocator) catch |e| {
        std.debug.print("{any}\n", .{e});
        return;
    };
    defer lineprof.deinit();

    try std.testing.expectEqual(lineprof.parameter_indices_to_table_index([2]usize{ 27, 12 }), 552);
    try std.testing.expectEqual(lineprof.parameter_indices_to_table_index([2]usize{ 0, 12 }), 12);

    const t1 = lineprof.find_tables([2]usize{ 27, 12 });
    try std.testing.expectEqualSlices(usize, &[4]usize{ 532, 552, 531, 551 }, &t1);

    const t2 = lineprof.find_tables([2]usize{ 0, 1 });
    try std.testing.expectEqualSlices(usize, &[4]usize{ 1, 1, 0, 0 }, &t2);

    const parameters = [2]f32{ 0.92, 0.4 };
    const p1 = lineprof.find_parameter_indices(parameters);
    try std.testing.expectEqualSlices(usize, &[2]usize{ 24, 8 }, &p1);

    const w1 = lineprof.determine_parameter_weight(&parameters, 0, 24);
    const w2 = lineprof.determine_parameter_weight(&parameters, 1, 8);
    std.debug.print("w1: {any}, w2: {any}\n", .{ w1, w2 });
}
