const std = @import("std");
const xspec = @import("./xspec-wrapper.zig");
const utils = @import("./utils.zig");

const TestSetup = struct {
    allocator: std.mem.Allocator,
    energy: []f64,
    n_flux: i64,
    parameters: []f64,
    flux: []f64,

    pub fn init(alloc: std.mem.Allocator, emin: f64, emax: f64, N: usize) !TestSetup {
        var erange = utils.RangeIterator(f64).init(emin, emax, N);
        const energy = try erange.drain(alloc);
        errdefer alloc.free(energy);

        const flux = try alloc.alloc(f64, energy.len - 1);
        errdefer alloc.free(flux);
        for (flux) |*f| f.* = 0;

        const parameters = try alloc.dupe(f64, &[_]f64{ 0.9, 80, 6.4, 3.0, 1.0, 50.0 });
        errdefer alloc.free(parameters);

        return .{
            .allocator = alloc,
            .energy = energy,
            .n_flux = @as(i64, @intCast(N - 1)),
            .parameters = parameters,
            .flux = flux,
        };
    }

    pub fn setParams(self: *@This(), params: []const f64) !void {
        self.allocator.free(self.parameters);
        self.parameters = try self.allocator.dupe(f64, params);
    }

    pub fn deinit(self: *@This()) void {
        self.allocator.free(self.energy);
        self.allocator.free(self.flux);
        self.allocator.free(self.parameters);
    }

    pub fn call(self: *@This(), comptime f: anytype) void {
        const start = std.time.microTimestamp();
        f(
            @as(*const f64, @ptrCast(self.energy.ptr)),
            @as(c_int, @intCast(self.n_flux)),
            @as(*const f64, @ptrCast(self.parameters.ptr)),
            @as(c_int, 0),
            @as(*f64, @ptrCast(self.flux.ptr)),
            undefined,
            undefined,
        );
        const delta = std.time.microTimestamp() - start;
        std.debug.print(
            "   elapsed time: {s}\n",
            .{std.fmt.fmtDuration(@intCast(delta * 1000))},
        );
    }
};

pub fn main() !void {
    const allocator = std.heap.c_allocator;

    { // additive model
        std.debug.print(" 1. Additive\n", .{});
        var setup = try TestSetup.init(allocator, 0.1, 12.0, 300);
        defer setup.deinit();
        try setup.setParams(&[_]f64{ 0.998, 40, 6.4, 3.0, 1.0, 50.0 });

        setup.call(xspec.kline);
    }

    { // convolutional model
        std.debug.print(" 2. Convolutional\n", .{});
        var setup = try TestSetup.init(allocator, 0.1, 12.0, 300);
        defer setup.deinit();
        try setup.setParams(&[_]f64{ 0.998, 40, 3.0, 1.0, 50.0 });

        // init delta flux
        setup.flux[150] = 1;
        setup.flux[200] = 1;
        setup.call(xspec.kconv);
        utils.normalize(f64, setup.flux);
    }

    { // convolutional model
        std.debug.print(" 3. Convolutional Five\n", .{});
        var setup = try TestSetup.init(allocator, 0.1, 12.0, 300);
        defer setup.deinit();

        var prng = std.Random.DefaultPrng.init(0);

        // init delta flux
        setup.flux[150] = 1;
        setup.flux[200] = 1;

        while (true) {
            try setup.setParams(&[_]f64{
                prng.random().float(f64) * 0.998,
                prng.random().float(f64) * 90,
                prng.random().float(f64) + 1.0,
                prng.random().float(f64) * 500,
                prng.random().float(f64) * 5.0,
                prng.random().float(f64) * 5.0,
                prng.random().float(f64) * 5.0,
                prng.random().float(f64) * 5.0,
                prng.random().float(f64) * 5.0,
                prng.random().float(f64) * 5.0,
            });
            std.debug.print(">> Params: {any}", .{setup.parameters});
            setup.call(xspec.kconv5);
        }

        utils.normalize(f64, setup.flux);
    }
}
