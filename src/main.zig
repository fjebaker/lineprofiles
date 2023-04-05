const std = @import("std");
const io = @import("io.zig");

test {
    var allocator = std.testing.allocator;
    io.readFitsFile("../relline/rel_table.fits", allocator) catch |err| {
        std.debug.print("{any}\n", .{err});
    };
}
