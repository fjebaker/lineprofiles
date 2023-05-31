const std = @import("std");
const xspec = @import("./xspec-wrapper.zig");
const utils = @import("./utils.zig");

fn plotxy(
    alloc: std.mem.Allocator,
    comptime T: type,
    x: []const T,
    y: []const T,
    filename: []const u8,
) !void {
    var buf = std.ArrayList(u8).init(alloc);
    defer buf.deinit();

    var writer = buf.writer();
    for (x, y) |i, j| {
        try writer.print("{d} {d}\n", .{ i, j });
    }

    var cmd = std.ChildProcess.init(
        &[_][]const u8{ "graph", "-T", "png" },
        alloc,
    );
    cmd.stdin_behavior = .Pipe;
    cmd.stdout_behavior = .Pipe;
    cmd.stderr_behavior = .Pipe;

    var output = std.ArrayList(u8).init(alloc);
    defer output.deinit();

    var err = std.ArrayList(u8).init(alloc);
    defer err.deinit();

    try cmd.spawn();

    try cmd.stdin.?.writer().writeAll(buf.items);
    cmd.stdin.?.close();
    cmd.stdin = null;

    try cmd.collectOutput(&output, &err, 1 << 16);

    _ = try cmd.wait();

    var outfile = try std.fs.cwd().createFile(filename, .{});
    try outfile.writeAll(output.items);
}

const TestSetup = struct {
    allocator: std.mem.Allocator,
    energy: []f64,
    n_flux: i64,
    parameters: []f64,
    flux: []f64,

    pub fn init(alloc: std.mem.Allocator, emin: f64, emax: f64, N: usize) !TestSetup {
        var erange = utils.RangeIterator(f64).init(emin, emax, N);
        var energy = try erange.drain(alloc);
        errdefer alloc.free(energy);

        var flux = try alloc.alloc(f64, energy.len - 1);
        errdefer alloc.free(flux);
        for (flux) |*f| f.* = 0;

        var parameters = try alloc.dupe(f64, &[_]f64{ 0.9, 80, 6.4, 1.0, 50.0 });
        errdefer alloc.free(parameters);

        return .{
            .allocator = alloc,
            .energy = energy,
            .n_flux = @intCast(i64, N - 1),
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
        f(
            @ptrCast(*const f64, self.energy.ptr),
            @intCast(c_int, self.n_flux),
            @ptrCast(*const f64, self.parameters.ptr),
            @as(c_int, 0),
            @ptrCast(*f64, self.flux.ptr),
            undefined,
            undefined,
        );
    }
};

pub fn main() !void {
    var allocator = std.heap.c_allocator;

    try plotxy(
        allocator,
        f64,
        &[_]f64{ 1, 2, 3 },
        &[_]f64{ 2, 4, 8 },
        "hello.png",
    );

    { // additive model
        var setup = try TestSetup.init(allocator, 0.1, 12.0, 300);
        defer setup.deinit();
        try setup.setParams(&[_]f64{ 0.9, 80, 6.4, 1.0, 50.0 });

        setup.call(xspec.kerr_line_profile);
        try plotxy(
            allocator,
            f64,
            setup.energy[0 .. setup.energy.len - 1],
            setup.flux,
            "additive-example.png",
        );
    }

    { // additive model
        var setup = try TestSetup.init(allocator, 0.1, 12.0, 300);
        defer setup.deinit();
        try setup.setParams(&[_]f64{
            0.9,
            80,
            6.4,
            1.0,
            50.0,
            20.0,
            3.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
        });

        // init delta flux
        setup.flux[100] = 1;
        setup.call(xspec.kerr_conv_profile);
        try plotxy(
            allocator,
            f64,
            setup.energy[0 .. setup.energy.len - 1],
            setup.flux,
            "convolve-example.png",
        );
    }
}
