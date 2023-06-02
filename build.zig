const std = @import("std");
const zigFITSTIO = @import("./vendor/zigFITSIO/zigFITSIO.zig");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const zigfitsio = zigFITSTIO.create(b, target);

    const xspec = b.addStaticLibrary(.{
        .name = "xsklineprofiles",
        .root_source_file = .{ .path = "src/xspec-wrapper.zig" },
        .target = target,
        .optimize = optimize,
    });
    zigfitsio.link(xspec);

    const xspec_step = b.step("xspec", "Build the XSPEC static library");
    const install_xspec = b.addInstallArtifact(xspec);
    xspec_step.dependOn(&install_xspec.step);

    var plots = b.addExecutable(.{
        .name = "plots",
        .root_source_file = .{ .path = "src/test.zig" },
        .target = target,
        .optimize = optimize,
    });
    zigfitsio.link(plots);

    const plots_step = b.step("plots", "Create debugging plots");
    const install_plots = b.addInstallArtifact(plots);
    const run_plots = b.addRunArtifact(plots);
    plots_step.dependOn(&install_plots.step);
    plots_step.dependOn(&run_plots.step);
}
