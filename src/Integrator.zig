const std = @import("std");
const xfunc = @import("transfer-functions.zig");

pub fn trapezoid_integration_weight(
    comptime T: type,
    domain: []const T,
    index: usize,
) T {
    if (index == 0) {
        return domain[1] - domain[0];
    }
    if (index == domain.len - 1) {
        return domain[index] - domain[index - 1];
    }
    // TODO: check this
    return domain[index + 1] - domain[index - 1];
}

// unlikely that the intervals in x are evenly spaced
// so if we can mandate the spacing, use a gaussian
// quadrature rule
pub fn integrate_gauss(comptime T: type, a: T, b: T, y: []const T) T {
    _ = a;
    _ = b;
    _ = y;
}
// if they are equally spaced
pub fn integrate_trapezoid(comptime T: type, x: []const T, y: []const T) T {
    _ = x;
    _ = y;
}
pub fn integrate_simpsons(comptime T: type, x: []const T, y: []const T) T {
    _ = x;
    _ = y;
}
// todo: romberg

