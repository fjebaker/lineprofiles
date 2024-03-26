const std = @import("std");

pub fn product(
    comptime T: type,
    array: []const T,
) T {
    var prod: T = 1;
    for (array) |v| {
        prod *= v;
    }
    return prod;
}

pub fn convert(
    comptime T: type,
    comptime C: type,
    allocator: std.mem.Allocator,
    source: []const C,
) ![]T {
    var out = try allocator.alloc(T, source.len);
    for (0..source.len) |i| {
        out[i] = @as(T, @floatCast(source[i]));
    }
    return out;
}

pub fn g_to_gstar(
    g: anytype,
    gmin: @TypeOf(g),
    gmax: @TypeOf(g),
) @TypeOf(g) {
    return (g - gmin) / (gmax - gmin);
}

pub fn gstar_to_g(
    gstar: anytype,
    gmin: @TypeOf(gstar),
    gmax: @TypeOf(gstar),
) @TypeOf(gstar) {
    return (gmax - gmin) * gstar + gmin;
}

pub fn linear_interpolation_weighted(
    comptime T: type,
    weight: T,
    y0: T,
    y1: T,
) T {
    // weight is (x - x0) / (x1 - x0)
    return (y1 - y0) * weight + y0;
}

pub fn linear_interpolate_arrays(
    comptime T: type,
    weight: T,
    a1: []const T,
    a2: []const T,
    out: []T,
) void {
    std.debug.assert(a1.len == a2.len);
    std.debug.assert(a2.len == out.len);
    for (0..a1.len) |i| {
        out[i] = linear_interpolation_weighted(T, weight, a1[i], a2[i]);
    }
}

pub fn dupe2d(
    comptime T: type,
    allocator: std.mem.Allocator,
    source: []const []T,
) ![][]T {
    var list = try std.ArrayList([]T).initCapacity(allocator, source.len);
    errdefer list.deinit();
    errdefer for (list.items) |item| allocator.free(item);
    for (source) |item| {
        list.appendAssumeCapacity(try allocator.dupe(T, item));
    }
    return list.toOwnedSlice();
}

pub fn free2d(comptime T: type, allocator: std.mem.Allocator, arr: [][]T) void {
    for (arr) |slice| {
        allocator.free(slice);
    }
    allocator.free(arr);
}

pub fn ValueIndex(comptime T: type) type {
    return struct { index: usize, value: T };
}

pub fn find_first_geq(
    comptime T: type,
    haystack: []const T,
    item: T,
    start: usize,
) ?ValueIndex(T) {
    for (haystack[start..], start..) |t, i| {
        if (t >= item) return .{
            .index = i,
            .value = t,
        };
    }
    return null;
}
pub fn find_first_leq(
    comptime T: type,
    haystack: []const T,
    item: T,
    start: usize,
) ?ValueIndex(T) {
    for (haystack[start..], start..) |t, i| {
        if (t <= item) return .{
            .index = i,
            .value = t,
        };
    }
    return null;
}

pub fn RangeIterator(comptime T: type) type {
    return struct {
        const Self = @This();
        delta: T,
        current: T,
        remaining: usize,

        pub fn next(self: *Self) ?T {
            if (self.remaining > 0) {
                const v = self.current;
                self.current += self.delta;
                self.remaining -= 1;
                return v;
            } else return null;
        }

        pub fn init(min: T, max: T, N: usize) Self {
            const delta = (max - min) / @as(T, @floatFromInt(N - 1));
            return .{
                .delta = delta,
                .remaining = N,
                .current = min,
            };
        }

        pub fn drain(self: *Self, allocator: std.mem.Allocator) ![]T {
            var out = try allocator.alloc(T, self.remaining);
            var i: usize = 0;
            while (self.next()) |v| {
                out[i] = v;
                i += 1;
            }
            return out;
        }
    };
}

test "range-iterator" {
    // 32/8 = 4
    var itt = RangeIterator(f32).init(8, 32, 4);
    try std.testing.expectEqual(itt.delta, 8);
    var last: f32 = itt.next().?;
    try std.testing.expectEqual(last, 8);
    var count: usize = 1;
    while (itt.next()) |v| {
        last = v;
        count += 1;
    }
    try std.testing.expectEqual(count, 4);
    try std.testing.expectEqual(last, 32);
}

pub fn inverse_grid_inplace(
    comptime T: type,
    grid: []T,
    min: T,
    max: T,
) void {
    var itt = RangeIterator(T).init(1 / max, 1 / min, grid.len);
    for (grid) |*g| {
        g.* = 1 / itt.next().?;
    }
    std.mem.reverse(T, grid);
}

pub fn inverse_grid(
    comptime T: type,
    alloc: std.mem.Allocator,
    min: T,
    max: T,
    N: usize,
) ![]T {
    const grid = try alloc.alloc(T, N);
    inverse_grid_inplace(T, grid, min, max);
    return grid;
}

pub fn normalize(comptime T: type, arr: []T) void {
    var sum: T = 0;
    for (arr) |v| sum += v;
    for (arr) |*v| v.* = v.* / sum;
}

pub fn refine_grid(comptime T: type, grid: anytype, fine_grid: []T, N: usize, norm: T) void {
    std.debug.assert((grid.len - 1) * N == fine_grid.len);
    const inorm = 1 / norm;

    var j: usize = 0;
    for (1..grid.len) |i| {
        const v0 = @as(T, @floatCast(grid[i - 1]));
        const v1 = @as(T, @floatCast(grid[i]));

        // interpolate the grid
        for (0..N) |k| {
            const factor = @as(T, @floatFromInt(k)) / @as(T, @floatFromInt(N));
            const v = factor * v1 + (1 - factor) * v0;
            fine_grid[j] = v * inorm;
            j += 1;
        }
    }
    // TODO: note that this is technically missing the very last
    // edge in grid
}

pub const H = 5e-4;
pub fn gstar_grid(comptime T: type, comptime N: comptime_int) [N]T {
    const g1 = H;
    const g2 = 1.0 - H;
    var grid: [N]T = undefined;

    var itt = RangeIterator(T).init(g1, g2, N);
    for (0..grid.len) |i| {
        const v = itt.next().?;
        grid[i] = v;
    }

    return grid;
}
