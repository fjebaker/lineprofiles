const std = @import("std");
const util = @import("utils.zig");

pub fn TransferFunction(comptime T: type) type {
    return struct {
        const Self = @This();
        // we will store [r][g*] -> f, since we will be accessing
        // over the latter more often, thus needs to be contiguous
        upper_branch: [][]T,
        lower_branch: [][]T,
        gmin: []T,
        gmax: []T,
        radii: []T,

        pub fn copy(self: *const Self, allocator: std.mem.Allocator) !Self {
            var ub = try util.dupe2d(T, allocator, self.upper_branch);
            errdefer util.free2d(T, allocator, ub);

            var lb = try util.dupe2d(T, allocator, self.lower_branch);
            errdefer util.free2d(T, allocator, lb);

            var gmin = try allocator.dupe(T, self.gmin);
            errdefer allocator.free(gmin);

            var gmax = try allocator.dupe(T, self.gmax);
            errdefer allocator.free(gmax);

            var radii = try allocator.dupe(T, self.radii);
            errdefer allocator.free(radii);

            return .{
                .upper_branch = ub,
                .lower_branch = lb,
                .gmin = gmin,
                .gmax = gmax,
                .radii = radii,
            };
        }

        pub fn free(self: *Self, allocator: std.mem.Allocator) void {
            util.free2d(T, allocator, self.upper_branch);
            util.free2d(T, allocator, self.lower_branch);
            allocator.free(self.gmin);
            allocator.free(self.gmax);
            allocator.free(self.radii);
        }

        pub fn assignFrom(self: *Self, other: *const Self) void {
            for (self.upper_branch, 0..) |ub, j| {
                std.mem.copy(T, ub, other.upper_branch[j]);
            }
            for (self.lower_branch, 0..) |lb, j| {
                std.mem.copy(T, lb, other.lower_branch[j]);
            }
            std.mem.copy(T, self.gmin, other.gmin);
            std.mem.copy(T, self.gmax, other.gmax);
            std.mem.copy(T, self.radii, other.radii);
        }

        pub fn interpolateBetween(
            self: *const Self,
            other: *const Self,
            out: *Self,
            weight: T,
        ) void {
            for (self.upper_branch, 0..) |ub, j| {
                util.linear_interpolate_arrays(
                    T,
                    weight,
                    ub,
                    other.upper_branch[j],
                    out.upper_branch[j],
                );
            }
            for (self.lower_branch, 0..) |lb, j| {
                util.linear_interpolate_arrays(
                    T,
                    weight,
                    lb,
                    other.lower_branch[j],
                    out.lower_branch[j],
                );
            }
            util.linear_interpolate_arrays(
                T,
                weight,
                self.gmin,
                other.gmin,
                out.gmin,
            );
            util.linear_interpolate_arrays(
                T,
                weight,
                self.gmax,
                other.gmax,
                out.gmax,
            );
            util.linear_interpolate_arrays(
                T,
                weight,
                self.radii,
                other.radii,
                out.radii,
            );
        }
    };
}

pub fn InterpolatingTransferFunction(comptime T: type) type {
    return struct {
        const Self = @This();
        // points to a cached branch
        tf: *TransferFunction(T),
        // cache of the current row of the data
        // into which we write the interpolations
        cache_upper: []T,
        cache_lower: []T,
        cache_gmin: T = 0,
        cache_gmax: T = 0,
        // we don't own these, they are referenced from
        // the parent
        gstars: []const T,
        // store the last high indices for cached searches
        last_r_index: usize = 0,
        last_g_index: usize = 0,

        fn stage_index(self: *Self, index: usize) void {
            for (self.cache_lower, 0..) |*f, i| f.* = self.tf.lower_branch[index][i];
            for (self.cache_upper, 0..) |*f, i| f.* = self.tf.upper_branch[index][i];
        }

        pub fn stage_radius(self: *Self, r: T) void {
            // radii are reversed
            const start = if (r > self.tf.radii[self.last_r_index])
                0
            else
                self.last_r_index;
            const loc = util.find_first_leq(T, self.tf.radii, r, start) orelse {
                // stage the minimal possible radius
                self.stage_index(self.tf.radii.len - 1);
                self.last_r_index = self.tf.radii.len - 1;
                return;
            };
            const ilow = if (loc.index > 1) loc.index - 1 else {
                // stage the minimal possible radius
                self.stage_index(0);
                self.last_r_index = 0;
                return;
            };

            // clamp to avoid extrapolation
            var r_clamped = std.math.clamp(
                r,
                self.tf.radii[self.tf.radii.len - 1],
                self.tf.radii[0],
            );

            // interpolant factor: since everything is for the same radius
            // we build it as (r - r0) / (r1 - r0) to then be multiplied by (f1 - f0)
            const r0 = self.tf.radii[ilow];
            const factor = (r_clamped - r0) / (loc.value - r0);
            self.stage_interpolated(loc.index, factor);
            self.last_r_index = loc.index;
        }

        fn stage_interpolated(self: *Self, index: usize, factor: T) void {
            for (self.cache_lower, 0..) |*f, i| {
                const f0 = self.tf.lower_branch[index - 1][i];
                f.* = (self.tf.lower_branch[index][i] - f0) * factor + f0;
            }

            for (self.cache_upper, 0..) |*f, i| {
                const f0 = self.tf.upper_branch[index - 1][i];
                f.* = (self.tf.upper_branch[index][i] - f0) * factor + f0;
            }

            const gmax0 = self.tf.gmax[index - 1];
            const gmax1 = self.tf.gmax[index];
            self.cache_gmax = (gmax1 - gmax0) * factor + gmax0;

            const gmin0 = self.tf.gmin[index - 1];
            const gmin1 = self.tf.gmin[index];
            self.cache_gmin = (gmin1 - gmin0) * factor + gmin0;
        }

        inline fn get_gstar_interpolating_factor(
            self: *const Self,
            gstar: T,
        ) ?util.ValueIndex(T) {
            const start = if (gstar >= self.gstars[self.last_g_index])
                self.last_g_index
            else
                0;
            // should not be possible for gstar to not be in [0,1]
            const loc = util.find_first_geq(T, self.gstars, gstar, start) catch unreachable;
            const ilow = if (loc.index > 1) loc.index - 1 else return null;
            const g0 = self.gstars[ilow];
            // interpolant factor
            const factor = (gstar - g0) / (self.gstars[loc.index] - g0);
            self.last_g_index = loc.index;
            return .{ .index = loc.index, .value = factor };
        }

        pub fn branches_at_gstar(
            self: *const Self,
            gstar: T,
        ) struct { lower: T, upper: T } {
            const loc = self.get_gstar_interpolating_factor(gstar) orelse {
                return .{
                    .lower = self.cache_lower[0],
                    .upper = self.cache_upper[0],
                };
            };
            // interpolate lower branch
            const f0lower = self.cache_lower[loc.index - 1];
            const lower = loc.factor * (self.cache_lower[loc.index] - f0lower) + f0lower;
            // interpolate upper branch
            const f0upper = self.cache_upper[loc.index - 1];
            const upper = loc.factor * (self.cache_upper[loc.index] - f0upper) + f0upper;
            return .{ .lower = lower, .upper = upper };
        }

        pub fn init(
            gstars: []const T,
            allocator: std.mem.Allocator,
            tf: *TransferFunction(T),
        ) !Self {
            const n = tf.lower_branch[0].len;
            var cache = try allocator.alloc(T, 2 * n);
            return .{
                .tf = tf,
                .cache_upper = cache[0..n],
                .cache_lower = cache[n..],
                .gstars = gstars,
            };
        }

        pub fn free(self: *Self, allocator: std.mem.Allocator) void {
            // only the upper cache has the pointer
            allocator.free(self.cache_upper);
        }
    };
}
