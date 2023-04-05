const std = @import("std");
const zfits = @import("zfits");

const lineprof = @import("line-profile.zig");
const TransferFunction = lineprof.TransferFunction(f32);

pub fn readFitsFile(path: []const u8, allocator: std.mem.Allocator) !void {
    var f = try zfits.FITS.initFromFile(path);
    defer f.deinit();

    const alpha_hdu = try f.getHDU(2);
    var alphas = try alpha_hdu.BinaryTable.getColumnTyped(f32, 1, allocator, .{});
    defer allocator.free(alphas);

    const mu_hdu = try f.getHDU(3);
    var mus = try mu_hdu.BinaryTable.getColumnTyped(f32, 1, allocator, .{});
    defer allocator.free(mus);

    var list = try std.ArrayList(TransferFunction).initCapacity(allocator, f.num_hdus - 4);
    defer list.deinit();
    defer for (list.items) |*tf| {
        tf.free(allocator);
    };

    for (4..f.num_hdus) |i| {
        std.debug.print("{d}\n", .{i});
        const table = (try f.getHDU(i)).BinaryTable;
        _ = table;
    }
}
