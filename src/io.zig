const std = @import("std");
const zfits = @import("zfits");

const lineprof = @import("line-profile.zig");
const xfunc = @import("transfer-functions.zig");

const TransferFunction = xfunc.TransferFunction(f32);
const LineProfileTable = lineprof.LineProfileTable(2, f32);

pub fn readTransferFunctions(
    comptime T: type,
    allocator: std.mem.Allocator,
    f: *zfits.FITS,
    offset: usize,
) ![]xfunc.TransferFunction(T) {
    var list = try std.ArrayList(TransferFunction).initCapacity(allocator, f.num_hdus - offset);
    errdefer list.deinit();
    errdefer for (list.items) |*tf| tf.free(allocator);

    // read in all of the transfer function data
    for (offset..f.num_hdus) |i| {
        const table = (try f.getHDU(i)).BinaryTable;

        var upper_f = try table.getColumnVectorTyped(f32, 4, allocator);
        errdefer allocator.free(upper_f);
        errdefer for (upper_f) |line| allocator.free(line);

        var lower_f = try table.getColumnVectorTyped(f32, 5, allocator);
        errdefer allocator.free(lower_f);
        errdefer for (lower_f) |line| allocator.free(line);

        var radii = try table.getColumnTyped(f32, 1, allocator, .{});
        errdefer allocator.free(radii);

        var gmins = try table.getColumnTyped(f32, 2, allocator, .{});
        errdefer allocator.free(gmins);

        var gmaxs = try table.getColumnTyped(f32, 3, allocator, .{});
        errdefer allocator.free(gmaxs);

        list.appendAssumeCapacity(.{
            .upper_branch = upper_f,
            .lower_branch = lower_f,
            .gmin = gmins,
            .gmax = gmaxs,
            .radii = radii,
        });
    }
    return list.toOwnedSlice();
}

pub fn readFitsFile(path: []const u8, allocator: std.mem.Allocator) !LineProfileTable {
    var f = try zfits.FITS.initFromFile(path);
    errdefer f.deinit();

    // todo: n generic parameters
    const alpha_hdu = try f.getHDU(2);
    var alphas = try alpha_hdu.BinaryTable.getColumnTyped(f32, 1, allocator, .{});
    errdefer allocator.free(alphas);

    const mu_hdu = try f.getHDU(3);
    var mus = try mu_hdu.BinaryTable.getColumnTyped(f32, 1, allocator, .{});
    errdefer allocator.free(mus);

    // todo: get the actual gstars
    var gstars = try allocator.dupe(f32, &[_]f32{ 0.0, 1.0 });
    errdefer allocator.free(gstars);

    // unpack the list
    var tf = try readTransferFunctions(f32, allocator, &f, 4);
    errdefer allocator.free(tf);
    errdefer for (tf) |*t| t.free(allocator);

    return LineProfileTable.init(
        allocator,
        gstars,
        tf,
        [2][]f32{ alphas, mus },
    );
}
