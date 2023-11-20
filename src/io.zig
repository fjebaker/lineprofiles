const std = @import("std");
const zigfitsio = @import("zfitsio");

const lineprof = @import("line-profile.zig");
const xfunc = @import("transfer-functions.zig");
const util = @import("utils.zig");

pub fn readTransferFunctions(
    comptime T: type,
    comptime DataType: type,
    allocator: std.mem.Allocator,
    f: *zigfitsio.FITS,
    offset: usize,
) ![]xfunc.TransferFunction(T) {
    if (T != DataType) @compileError("Type mistmatch.");

    var list = try std.ArrayList(xfunc.TransferFunction(T)).initCapacity(allocator, f.num_hdus - offset + 1);
    errdefer list.deinit();
    errdefer for (list.items) |*tf| tf.free(allocator);

    // read in all of the transfer function data
    for (offset..f.num_hdus + 1) |i| {
        const table = (try f.getHDU(i)).BinaryTable;

        var upper_f = try table.getColumnVectorTyped(DataType, 4, allocator);
        errdefer allocator.free(upper_f);
        errdefer for (upper_f) |line| allocator.free(line);

        var lower_f = try table.getColumnVectorTyped(DataType, 5, allocator);
        errdefer allocator.free(lower_f);
        errdefer for (lower_f) |line| allocator.free(line);

        var radii = try table.getColumnTyped(DataType, 1, allocator, .{});
        errdefer allocator.free(radii);

        var gmins = try table.getColumnTyped(DataType, 2, allocator, .{});
        errdefer allocator.free(gmins);

        var gmaxs = try table.getColumnTyped(DataType, 3, allocator, .{});
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

pub fn readFitsFile(
    comptime NParams: comptime_int,
    comptime T: type,
    comptime Ngstar: comptime_int,
    path: []const u8,
    allocator: std.mem.Allocator,
) !lineprof.LineProfileTable(NParams, T) {
    var f = try zigfitsio.FITS.initFromFile(path);
    errdefer f.deinit();

    // todo: n generic parameters
    const alpha_hdu = try f.getHDU(2);
    var alphas = try alpha_hdu.BinaryTable.getColumnTyped(T, 1, allocator, .{});
    errdefer allocator.free(alphas);

    const incl = try f.getHDU(3);
    var angle = try incl.BinaryTable.getColumnTyped(T, 1, allocator, .{});
    errdefer allocator.free(angle);

    var gstars = try allocator.dupe(T, &util.gstar_grid(T, Ngstar));
    errdefer allocator.free(gstars);

    // unpack the list
    var tf = try readTransferFunctions(T, f32, allocator, &f, 4);
    errdefer allocator.free(tf);
    errdefer for (tf) |*t| t.free(allocator);

    return lineprof.LineProfileTable(NParams, T).init(
        allocator,
        gstars,
        tf,
        [2][]T{ alphas, angle },
    );
}
