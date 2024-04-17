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

pub fn LinInterpEmissivity(comptime T: type, comptime N: comptime_int) type {
    return struct {
        const Self = @This();

        powers: [N]T,
        radii: [N]T,
        coeffs: [N]T,
        alpha: T,

        fn calculateLogCoefficients(coeffs: []T, radii: []const T, powers: []const T, alpha: T) void {
            // do the edge case
            coeffs[coeffs.len - 1] = radii[radii.len - 1] * (powers[powers.len - 1] - alpha);

            for (0..radii.len - 2) |j| {
                // run in reverse, i goes form radii.len - 1 to 1
                const i = radii.len - j - 2;
                const delta_p = powers[i] - powers[i + 1];
                coeffs[i] = coeffs[i + 1] + (radii[i] * delta_p);
            }
        }

        /// Given some power indices, an inner and outer radius, and the
        /// regular behaviour, returns a `LinInterpEmissivity` that
        /// interpolates the emissivity function appropriately.
        pub fn init(powers: [N]T, rmin: T, rmax: T, alpha: T) Self {
            var radii: [N]T = undefined;
            var coeffs: [N]T = undefined;
            var itt = util.RangeIterator(T).init(std.math.log10(rmin), std.math.log10(rmax) / 2.5, N + 1);
            _ = itt.next();

            for (&radii) |*r| r.* = itt.next().?;
            calculateLogCoefficients(&coeffs, &radii, &powers, alpha);

            // exponentiate everything back to regular values
            for (&radii) |*r| r.* = std.math.pow(T, 10, r.*);
            for (&coeffs) |*c| c.* = std.math.pow(T, 10, c.*);

            const em: Self = .{
                .powers = powers,
                .radii = radii,
                .alpha = alpha,
                .coeffs = coeffs,
            };
            return em;
        }

        pub fn emissivity(self: *const Self, r: T) T {
            const i = if (util.find_first_geq(T, &self.radii, r, 0)) |vi| vi.index else {
                // constant power law
                return std.math.pow(T, r, -self.alpha);
            };
            if (i == 0) {
                return std.math.pow(T, r, -self.powers[0]);
            }

            const pow = self.powers[i];
            const coeff = self.coeffs[i];
            return coeff * std.math.pow(T, r, -pow);
        }
    };
}

test "emissivity interpolations" {
    const em = LinInterpEmissivity(f32, 4).init(
        [_]f32{ 5, 1, 3, 2 },
        1.0,
        50.0,
        3,
    );
    _ = em;
}
