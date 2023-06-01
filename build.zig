const std = @import("std");
const zigFITSTIO = @import("./vendor/zigFITSIO/zigFITSIO.zig");

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
            .version = .{ .major = 0, .minor = 1 },
        });
        zigfitsio.link(xspec);
        b.installArtifact(xspec);
    }

    var main_tests = b.addExecutable(.{
        .name = "plots",
        .root_source_file = .{ .path = "src/test.zig" },
        .target = target,
        .optimize = optimize,
    });
    zigfitsio.link(main_tests);

    const test_step = b.step("plots", "Create debugging plots");
    const run_step = b.addRunArtifact(main_tests);
    const ins_step = b.addInstallArtifact(main_tests);
    test_step.dependOn(&ins_step.step);
    test_step.dependOn(&run_step.step);
}
