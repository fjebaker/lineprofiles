const std = @import("std");
const util = @import("utils.zig");

pub fn Emissivity(comptime T: type) type {
    return struct {
        const Self = @This();
        const EmissProto = fn (*anyopaque, r: T) T;

        state: *anyopaque,
        emissivity_function: *const EmissProto,

        pub fn init(state: anytype, comptime emissivity_function: anytype) Self {
            const State = @TypeOf(state);
            const alignment = @typeInfo(State).Pointer.alignment;

            const gen = struct {
                fn _emissivity_function(ptr: *anyopaque, r: T) T {
                    const _state = @ptrCast(State, @alignCast(alignment, ptr));
                    return @call(.always_inline, emissivity_function, .{ _state, r });
                }
            };

            return .{
                .state = state,
                .emissivity_function = gen._emissivity_function,
            };
        }

        pub fn call(self: *const Self, r: T) T {
            return self.emissivity_function(self.state, r);
        }
    };
}

pub fn PowerLawEmissivity(comptime T: type) type {
    return struct {
        const Self = @This();

        index: T,

        pub fn init(index: T) Self {
            return .{ .index = index };
        }

        pub fn emissivity(self: *Self) Emissivity(T) {
            return Emissivity(T).init(self, Self.emissivity_function);
        }

        fn emissivity_function(self: *Self, r: T) T {
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

        pub fn emissivity(self: *Self) Emissivity(T) {
            return Emissivity(T).init(self, Self.emissivity_function);
        }

        fn emissivity_function(self: *Self, r: T) T {
            const i = if (util.find_first_geq(T, &self.cutoffs, r, 0)) |vi| vi.index else Nbins - 1;
            return self.weights[i];
        }
    };
}
