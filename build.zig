const std = @import("std");

const HEADER_PATH = "include";
// Although this function looks imperative, note that its job is to
// declaratively construct a build graph that will be executed by an external
// runner.
pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const zfits = b.dependency("zfits", .{ .target = target, .optimize = optimize });

    const lib = b.addStaticLibrary(.{
        .name = "lineprofile-integrator",
        .root_source_file = .{ .path = "src/main.zig" },
        .target = target,
        .optimize = optimize,
    });
    lib.addModule("zfits", zfits.module("zfits"));
    lib.addIncludePath(HEADER_PATH);
    lib.linkLibrary(zfits.artifact("cfitsio"));

    lib.install();

    // create a temporary executable
    const exe = b.addExecutable(.{
        .name = "profl",
        .root_source_file = .{ .path = "src/exe.zig" },
        .target = target,
        .optimize = optimize,
    });
    exe.addModule("zfits", zfits.module("zfits"));
    exe.addIncludePath(HEADER_PATH);
    exe.linkLibrary(zfits.artifact("cfitsio"));
    exe.install();

    // Creates a step for unit testing.
    const main_tests = b.addTest(.{
        .root_source_file = .{ .path = "src/main.zig" },
        .target = target,
        .optimize = optimize,
    });
    main_tests.addModule("zfits", zfits.module("zfits"));
    main_tests.linkLibrary(zfits.artifact("cfitsio"));
    main_tests.addIncludePath(HEADER_PATH);

    const test_step = b.step("test", "Run library tests");
    test_step.dependOn(&main_tests.run().step);
}
