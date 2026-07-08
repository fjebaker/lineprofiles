const std = @import("std");
const clippy = @import("clippy");

const io = @import("io.zig");
const util = @import("utils.zig");
const tf = @import("transfer-functions.zig");
const lp = @import("line-profile.zig");
const xspec = @import("xspec-wrapper.zig");
const emissivity = @import("emissivity.zig");

test "main" {
    _ = tf;
    _ = lp;
}

test "index-lookups" {
    var lineprof = io.readFitsFile(2, f32, 30, "./kerr-transfer-functions.fits", std.testing.allocator) catch |e| {
        std.debug.print("{any}\n", .{e});
        return;
    };
    defer lineprof.deinit();

    try std.testing.expectEqual(lineprof.parameter_indices_to_table_index([2]usize{ 27, 12 }), 552);
    try std.testing.expectEqual(lineprof.parameter_indices_to_table_index([2]usize{ 0, 12 }), 12);

    const t1 = lineprof.find_tables([2]usize{ 27, 12 });
    try std.testing.expectEqualSlices(usize, &[4]usize{ 532, 552, 531, 551 }, &t1);

    const t2 = lineprof.find_tables([2]usize{ 0, 1 });
    try std.testing.expectEqualSlices(usize, &[4]usize{ 1, 1, 0, 0 }, &t2);

    const parameters = [2]f32{ 0.92, 0.4 };
    const p1 = lineprof.find_parameter_indices(parameters);
    try std.testing.expectEqualSlices(usize, &[2]usize{ 24, 8 }, &p1);

    const w1 = lineprof.determine_parameter_weight(&parameters, 0, 24);
    const w2 = lineprof.determine_parameter_weight(&parameters, 1, 8);
    std.debug.print("w1: {any}, w2: {any}\n", .{ w1, w2 });
}

const Kconv5Arguments = clippy.Arguments(
    &.{},
);

const EmissivityArguments = clippy.Arguments(
    &.{
        .{
            .arg = "--indices indexes",
            .help = "Comma seperated emissivity index used in the fixed mode.",
            .required = true,
        },
        .{
            .arg = "--alpha alpha",
            .help = "Powerlaw index.",
            .argtype = f64,
            .default = "3.0",
        },
        .{
            .arg = "--rcut rcut",
            .argtype = f64,
            .default = "100.0",
            .help = "Where the cutoff radius.",
        },
    },
);

const HelpArguments = clippy.Arguments(&.{
    .{ .arg = "--help", .help = "Print help message" },
});

pub const Commands = clippy.Commands(union(enum) {
    kconv5: Kconv5Arguments,
    emissivity: EmissivityArguments,
});

pub fn main() !void {
    var gpa: std.heap.DebugAllocator(.{}) = .init;
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const raw_args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, raw_args);

    var itt = clippy.ArgumentIterator.init(raw_args[1..]);
    var help_itt = itt;

    var stdout = std.fs.File.stdout();
    var stdout_writer = stdout.writer(&.{});
    const writer = &stdout_writer.interface;

    // first we parse to see if help was given
    const help_parsed = try HelpArguments.initParseAll(&help_itt, .{ .forgiving = true });

    if (help_parsed.help) {
        try Commands.writeHelp(writer, .{});
        return std.process.cleanExit();
    }

    const cmd = try Commands.initParseAll(&itt, .{});

    switch (cmd) {
        .kconv5 => |args| {
            _ = args;
        },
        .emissivity => |args| {
            var fixed_indices: [4]f64 = .{0} ** 4;
            var index: usize = 0;

            var token_itt = std.mem.tokenizeScalar(u8, args.indices, ',');
            while (token_itt.next()) |token| {
                fixed_indices[index] = try std.fmt.parseFloat(f64, token);
                index += 1;
            }

            if (index != 4) {
                return try clippy.throwError(
                    error.TooFewIndices,
                    "At least 4 indices are required to `--indices`. Only {d} were passed",
                    .{index},
                );
            }

            const emiss = emissivity.fixedEmissivity(
                f64,
                fixed_indices,
                args.rcut,
                args.alpha,
            );

            var radii = util.RangeIterator(f64).init(std.math.log10(1.0), std.math.log10(1e4), 1000);

            while (radii.next()) |logr| {
                const r = std.math.pow(f64, 10, logr);
                const em = emiss.emissivity(r);
                // _ = em;
                try writer.print("{d} {d}\n", .{ r, em });
            }
        },
    }
    try writer.flush();
}
