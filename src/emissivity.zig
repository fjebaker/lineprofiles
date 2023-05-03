const std = @import("std");
const util = @import("utils.zig");

pub fn PowerLawEmissivity(comptime T: type) type {
    return struct {
        const Self = @This();
        index: T,
        pub fn init(index: T) Self {
            return .{ .index = index };
        }
        pub fn emissivity(self: *const Self, r: T) T {
            return std.math.pow(T, r, self.index);
        }
    };
}

pub fn StepFunctionEmissivity(comptime T: type, comptime Nbins: comptime_int) type {
    return struct {
        const Self = @This();
        weights: [Nbins]T,
        cutoffs: [Nbins]T,
        pub fn init(weights: [Nbins]T, rmin: T, rmax: T) Self {
            var cutoffs: [Nbins]T = undefined;
            var itt = util.RangeIterator(T).init(std.math.log10(rmin), std.math.log10(rmax), Nbins + 1);
            // throw away minimum
            _ = itt.next();
            for (&cutoffs) |*c| {
                c.* = std.math.pow(T, 10, itt.next().?);
            }
            return .{ .weights = weights, .cutoffs = cutoffs };
        }

        pub fn emissivity(self: *const Self, r: T) T {
            const i = if (util.find_first_geq(T, &self.cutoffs, r, 0)) |vi| vi.index else Nbins - 1;
            return self.weights[i];
        }
    };
}
