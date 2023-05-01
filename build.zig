const std = @import("std");

const zigFITSTIO = @import("./vendor/zigFITSIO/zigFITSIO.zig");

const HEADER_PATH = "include";
// Although this function looks imperative, note that its job is to
// declaratively construct a build graph that will be executed by an external
// runner.
pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const zigfitsio = zigFITSTIO.create(b, target);
    const build_xspec = b.option(bool, "xspec", "Build XSPEC shared library.") orelse false;
    if (build_xspec) {
        const xspec = b.addStaticLibrary(.{
            .name = "xslineprofiles",
            .root_source_file = .{ .path = "src/xspec-wrapper.zig" },
            .target = target,
            .optimize = optimize,
        });
        zigfitsio.link(xspec);
        b.installArtifact(xspec);
    } else {
        const xspec = b.addSharedLibrary(.{
            .name = "xslineprofiles",
            .root_source_file = .{ .path = "src/xspec-wrapper.zig" },
            .target = target,
            .optimize = optimize,
            .version = .{ .major = 0, .minor = 2 },
        });
        zigfitsio.link(xspec);
        b.installArtifact(xspec);
    }

    // create a temporary executable
    const debug_exe = b.option(bool, "exe", "Build a CLI executable (intended for debug).") orelse false;
    if (debug_exe) {
        const exe = b.addExecutable(.{
            .name = "test-exe",
            .root_source_file = .{ .path = "src/exe.zig" },
            .target = target,
            .optimize = optimize,
        });
        zigfitsio.link(exe);
        b.installArtifact(exe);
    }

    // Creates a step for unit testing.
    var main_tests = b.addTest(.{
        .root_source_file = .{ .path = "src/main.zig" },
        .target = target,
        .optimize = optimize,
    });
    zigfitsio.link(main_tests);

    const test_step = b.step("test", "Run library tests");
    const run_step = b.addRunArtifact(main_tests);
    test_step.dependOn(&run_step.step);
}
