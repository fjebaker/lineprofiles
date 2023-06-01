const std = @import("std");

/// convolves a and b on domains g and x respectively by rebinning
/// if necessary, assumes output is on domain x
/// approximately follows a method of constructing a Toeplitz matrix on
/// an (irregular) grid, but without physically storing the matrix
pub fn convolve(
    comptime T: type,
    output: []T,
    g: []const T,
    a: []const T,
    x: []const T,
    b: []const T,
) void {
    std.debug.assert(output.len == b.len);

    // find window of non-zero
    const window_start = blk: {
        for (0..a.len) |i| if (a[i] > 0)
            break :blk if (i > 0) i - 1 else i;
        break :blk 0;
    };
    const window_end = blk: {
        for (0..a.len) |i| {
            const j = a.len - i - 1;
            if (a[j] > 0)
                break :blk if (j < a.len - 1) j + 1 else j;
        }
        break :blk 0;
    };

    // find the extremal
    const g_min = g[window_start];
    const g_max = g[window_end];

    // zero the output
    for (output) |*o| o.* = 0;

    for (0..output.len) |i| {
        const avg = 2 / (x[i + 1] + x[i]);
        for (0..output.len) |j| {
            // output bin extremes shifted
            const low = x[j] * avg;
            const high = x[j + 1] * avg;

            // skip if outside of the window
            if (high < g_min or low > g_max) continue;

            // integrate the window
            const weight = convolution_weight(T, window_start, window_end, g, a, low, high);

            // accumulate the dot product
            output[j] += weight * b[i];
        }
    }
}

fn convolution_weight(
    comptime T: type,
    start: usize,
    end: usize,
    g: []const T,
    a: []const T,
    low: T,
    high: T,
) T {
    var weight: T = 0;
    for (start..end) |w| {
        const bin_low = g[w];
        const bin_high = g[w + 1];
        const bin_width = bin_high - bin_low;
        // check different cases to calculate overlap amount
        // 1. check if bin is not in output bin
        //    bin_low |---| bin_high     low |---| high
        if (bin_high < low) {
            continue;
        } else if (bin_low > high) {
            // break the loop early
            break;
        }
        var overlap: T = 1;
        // 2. check if bin is bigger than output bin
        //    bin_low |...| low |.....| high |...| bin_high
        if (bin_low < low and bin_high > high) {
            overlap = (high - low) / bin_width;
        }
        // 3. check if bin is striding
        //    bin_low |---| low |.....| bin_high |---| high
        else if (bin_low < low and bin_high < high) {
            overlap = (bin_high - low) / bin_width;
        }
        // or
        //    low |---| bin_low |.....| high |---| bin_high
        else if (bin_low > low and bin_high > high) {
            overlap = (high - bin_low) / bin_width;
        }
        // 4. bin must be contained, and overlap is 1
        //    low |---| bin_low |.....| bin_high |---| high
        weight += a[w] * bin_width * overlap;
    }
    return weight;
}
