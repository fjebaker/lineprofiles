const std = @import("std");
pub fn addTracy(b: *std.Build, step: *std.Build.Step.Compile) !void {
    step.want_lto = false;
    step.addCSourceFile(
        .{
            .file = b.path("../tracy/public/TracyClient.cpp"),
            .flags = &.{
                "-DTRACY_ENABLE=1",
                "-fno-sanitize=undefined",
            },
        },
    );
    step.linkLibCpp();
}

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const zfits = b.dependency("zfits", .{ .target = target, .optimize = optimize });
    const zfitsio = zfits.module("zfitsio");

    const tracy = b.option(bool, "tracy", "Compile against tracy client") orelse false;
    var opts = b.addOptions();
    opts.addOption(bool, "tracy_enabled", tracy);

    const xspec = b.addLibrary(.{
        .name = "xsklineprofiles",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/xspec-wrapper.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zfitsio", .module = zfitsio },
                .{ .name = "options", .module = opts.createModule() },
            },
        }),
    });

    const xspec_step = b.step("xspec", "Build the XSPEC static library");
    const install_xspec = b.addInstallArtifact(xspec, .{});
    xspec_step.dependOn(&install_xspec.step);

    const plots = b.addExecutable(.{
        .name = "plots",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/plots.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zfitsio", .module = zfitsio },
                .{ .name = "options", .module = opts.createModule() },
            },
        }),
    });

    if (tracy) {
        try addTracy(b, plots);
    }

    const plots_step = b.step("plots", "Create debugging plots");
    const install_plots = b.addInstallArtifact(plots, .{});
    // const run_plots = b.addRunArtifact(plots);
    plots_step.dependOn(&install_plots.step);
    // plots_step.dependOn(&run_plots.step);

    const fuzz = b.addExecutable(.{
        .name = "fuzz",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/fuzz.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "zfitsio", .module = zfitsio }},
        }),
    });

    const fuzz_step = b.step("fuzz", "Create debugging fuzz");
    const fuzz_install = b.addInstallArtifact(fuzz, .{});
    const fuzz_run = b.addRunArtifact(fuzz);
    fuzz_step.dependOn(&fuzz_install.step);
    fuzz_step.dependOn(&fuzz_run.step);

    const unit_tests = b.addTest(.{
        .root_module = xspec.root_module,
    });

    const run_unit_tests = b.addRunArtifact(unit_tests);
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_unit_tests.step);
}
