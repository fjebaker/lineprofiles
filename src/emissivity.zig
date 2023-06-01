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

pub fn LinInterpEmissivity(comptime T: type, comptime Nbins: comptime_int) type {
    return struct {
        const Self = @This();
        knots: [Nbins]T,
        cutoffs: [Nbins]T,
        alpha: T,

        /// Knots are defines as the edges between which linear interpolations are reconstructed in log space.
        ///
        ///   k1   k2   k3   kn
        ///   |    |    |    |
        ///   |    x -- x    |
        ///   |  / |   linear interpolation between knots
        ///   | /  |    |  \ |
        ///   x    |    |   \
        ///   |    |    |    \   regular power law
        ///
        ///  log r
        ///
        /// This function will not extrapolate beyond rmin and rmax. Beyond kn, the emissivity model assumes
        /// a standard power law with fixed emissivity alpha, renormalized to intersect with the knot at kn.
        pub fn init(knots: [Nbins]T, rmin: T, rmax: T, alpha: T) Self {
            var cutoffs: [Nbins]T = undefined;
            var itt = util.RangeIterator(T).init(std.math.log10(rmin), std.math.log10(rmax), Nbins);
            for (&cutoffs) |*c| {
                c.* = std.math.pow(T, 10, itt.next().?);
            }
            return .{ .knots = knots, .cutoffs = cutoffs, .alpha = -alpha };
        }
        pub fn emissivity(self: *const Self, r: T) T {
            const i = if (util.find_first_geq(T, &self.cutoffs, r, 0)) |vi| vi.index else {
                // constant power law
                return self.knots[Nbins - 1] * std.math.pow(T, r, self.alpha) / std.math.pow(T, self.cutoffs[Nbins - 1], self.alpha);
            };
            if (i == 0) {
                return self.knots[0];
            }

            const r1 = std.math.log10(self.cutoffs[i]);
            const r0 = std.math.log10(self.cutoffs[i - 1]);
            const y1 = std.math.log10(self.knots[i]);
            const y0 = std.math.log10(self.knots[i - 1]);
            const weight = (std.math.log10(r) - r0) / (r1 - r0);

            return std.math.pow(T, 10, (y1 - y0) * weight + y0);
        }
    };
}
