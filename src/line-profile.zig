const std = @import("std");

pub fn TransferFunction(comptime T: type) type {
    return struct {
        const Self = @This();
        // we will store [r][g*] -> f, since we will be accessing
        // over the latter more often, thus needs to be contiguous
        upper_branch: [][]T,
        lower_branch: [][]T,
        gmin: T,
        gmax: T,
        pub fn free(self: *Self, alloctor: std.mem.Allocator) void {
            for (self.upper_branch) |ptr| {
                alloctor.free(ptr);
            }
            alloctor.free(self.upper_branch);
            for (self.lower_branch) |ptr| {
                alloctor.free(ptr);
            }
            alloctor.free(self.lower_branch);
        }
    };
}

fn ValueIndex(comptime T: type) type {
    return struct { index: usize, value: T };
}

fn find_first_geq(comptime T: type, haystack: []T, item: T, start: usize) ?ValueIndex(T) {
    for (haystack[start..], start..) |t, i| {
        if (t >= item) return i;
    }
    return null;
}

const InterpolationError = error{IdenticalDomain};

fn linear_interpolation(comptime T: type, x: T, x0: T, x1: T, y0: T, y1: T) InterpolationError!T {
    if (x0 == x1) return .IdenticalDomain;
    const factor = (y1 - y0) / (x1 - x0);
    return (x - x0) * factor + y0;
}

fn InterpolatingTransferFunction(comptime T: type) type {
    return struct {
        const Self = @This();
        tf: TransferFunction(T),
        // cache of the current row of the data
        // into which we write the interpolations
        cache_upper: []T,
        cache_lower: []T,
        // we don't own these, they are referenced from
        // the parent
        radii: []T,
        gstars: []T,
        // store the last high indices for cached searches
        last_r_index: usize,
        last_g_index: usize,

        fn stage_index(self: *Self, index: usize) void {
            for (self.cache_lower, 0..) |*g, i| g.* = self.tf.lower_branch[index][i];
            for (self.cache_upper, 0..) |*g, i| g.* = self.tf.upper_branch[index][i];
        }

        pub fn stage_radius(self: *Self, r: T) void {
            const start = if (r >= self.radii[self.last_r_index]) self.last_r_index else 0;
            const loc = find_first_geq(T, self.radii, r, start) orelse {
                // stage the greatest possible radius
                self.stage_index(self.radii.len - 1);
                self.last_r_index = self.radii.len - 1;
                return;
            };
            const ilow = if (loc.index > 1) loc.index - 1 else {
                // stage the minimal possible radius
                self.stage_index(0);
                self.last_r_index = 0;
                return;
            };

            // interpolant factor: since everything is for the same radius
            // we build it as (r - r0) / (r1 - r0) to then be multiplied by (f1 - f0)
            const r0 = self.radii[ilow];
            const factor = (r - r0) / (loc.value - r0);
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
        }

        inline fn get_gstar_interpolating_factor(self: *const Self, gstar: T) ?ValueIndex(T) {
            const start = if (gstar >= self.gstars[self.last_g_index]) self.last_g_index else 0;
            // should not be possible for gstar to not be in [0,1]
            const loc = find_first_geq(T, self.gstars, gstar, start) catch unreachable;
            const ilow = if (loc.index > 1) loc.index - 1 else return null;
            const g0 = self.gstars[ilow];
            // interpolant factor
            const factor = (gstar - g0) / (self.gstars[loc.index] - g0);
            self.last_g_index = loc.index;
            return .{ .index = loc.index, .value = factor };
        }

        pub fn branches_at_gstar(self: *const Self, gstar: T) struct { lower: T, upper: T } {
            const loc = self.get_gstar_interpolating_factor(gstar) orelse {
                return .{ .lower = self.cache_lower[0], .upper = self.cache_upper[0] };
            };
            // interpolate lower branch
            const f0lower = self.cache_lower[loc.index - 1];
            const lower = loc.factor * (self.cache_lower[loc.index] - f0lower) + f0lower;
            // interpolate upper branch
            const f0upper = self.cache_upper[loc.index - 1];
            const upper = loc.factor * (self.cache_upper[loc.index] - f0upper) + f0upper;
            return .{ .lower = lower, .upper = upper };
        }

        pub fn free(self: *Self, allocator: std.mem.Allocator) void {
            self.tf.free(allocator);
            allocator.free(self.cache_lower);
            allocator.free(self.cache_upper);
        }
    };
}

/// keeps all of the read in data in a structure
fn LineProfileTable(comptime NParams: comptime_int, comptime T: type) type {
    return struct {
        const Self = @This();
        const TFunc = TransferFunction(T);

        allocator: std.mem.Allocator,
        // the transfer function data maps g* -> f for different radii
        // these must be the same for all tables of data
        // if data for a specific radius does not exist, it is zeroed
        gstars: []T,
        radii: []T,
        // there is a matrix of values for each parameter combination
        // there will be Np1 * Np2 * ... * Npn many frames, scaling with the
        // number of parameters, and how its been interpolated
        transfer_functions: []TFunc,
        // we store one array for each parameter with all of its values
        parameters: [NParams][]T,

        // since the transfer function could be laid out in memory in an
        // inconvenient way, we do all of the interpolations first,
        // and then return a `TranferFunction`, storing indexes to the
        // frames that were used to generate the table, since it is likely
        // we'll just have to reinterpoalte, and can save the parameter lookup.
        last_upper_frame: [2 * NParams]usize = .{0 ** (2 * NParams)},
        last_lower_frame: [2 * NParams]usize = .{0 ** (2 * NParams)},
        // that does also mean we store the parameters which strided the last used
        // frames, such that the needed parameters P are
        //       last_upper_parameters >= P >= last_lower_parameters
        // which means we need 2 frames pointers for each parameter
        last_lower_parameters: [NParams]usize = .{0 ** NParams},
        // last_upper_parameters is implied as last_lower_parameters[i] + 1

        // then, when we receive some new parameters, we can optimize the search for
        // the new index, and if the same index is used, we already have the frames over
        // which to redo the interpolation

        // dont want to have to allocate new memory each time, so we keep one frame special
        // into which we write the new interpolation and return a copy which does not need
        // to be freed (since we own the memory)
        interpolated_cache: InterpolatingTransferFunction(T),

        pub fn interpolate_parameters(
            self: *Self,
            parameters: [NParams]T,
        ) TFunc {
            _ = parameters;

            return self.interpolated_cache;
        }

        pub fn init(
            alloc: std.mem.Allocator,
            gstars: []T,
            radii: []T,
            transfer_functions: []TFunc,
            parameters: [NParams][]T,
        ) Self {
            return .{
                .alloctor = alloc,
                .gstars = gstars,
                .radii = radii,
                .transfer_functions = transfer_functions,
                .parameters = parameters,
            };
        }

        // standard free everything approach
        pub fn deinit(self: *Self) void {
            self.allocator.free(self.gstars);
            self.allocator.free(self.radii);
            // free the parameter arrays
            for (self.parameters) |*p| {
                self.allocator.free(p);
            }
            // free the transfer functions
            for (self.transfer_functions) |*tf| {
                tf.free(self.allocator);
            }
            self.allocator.free(self.transfer_functions);
            // free the interpolant cache
            self.interpolated_cache.free(self.alloctor);
        }
    };
}

fn Integrator(comptime T: type) type {
    return struct {
        const Self = @This();
        // unlikely that the intervals in x are evenly spaced
        // so if we can mandate the spacing, use a gaussian
        // quadrature rule
        pub fn integrate_gauss(self: *Self, a: T, b: T, y: []const T) T {
            _ = self;
            _ = a;
            _ = b;
            _ = y;
        }
        // if they are equally spaced
        pub fn integrate_trapezoid(self: *Self, x: []const T, y: []const T) T {
            _ = self;
            _ = x;
            _ = y;
        }
        pub fn integrate_simpsons(self: *Self, x: []const T, y: []const T) T {
            _ = self;
            _ = x;
            _ = y;
        }
        // todo: romberg
    };
}

/// facilitates the integration over arbitrary parameter, energy
/// and flux arrays
fn MultiEmissivityLineProfile(
    comptime NParams: comptime_int,
    comptime T: type,
) type {
    return struct {
        table: LineProfileTable(NParams, T),
        integrator: Integrator(T),
    };
}
