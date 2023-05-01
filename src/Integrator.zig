const std = @import("std");
const xfunc = @import("transfer-functions.zig");

// gauss 7 weights
const GAUSS_X = [_]comptime_float{
    9.4910791234275852452618968404809e-01,
    7.415311855993944398638647732811e-01,
    4.0584515137739716690660641207707e-01,
    0.0,
};
const GAUSS_W = [_]comptime_float{
    1.2948496616886969327061143267787e-01,
    2.797053914892766679014677714229e-01,
    3.8183005050511894495036977548818e-01,
    4.1795918367346938775510204081658e-01,
};
const KRONROD_X = [_]comptime_float{
    -9.9145537112081263920685469752598e-01,
    -9.4910791234275852452618968404809e-01,
    -8.6486442335976907278971278864098e-01,
    -7.415311855993944398638647732811e-01,
    -5.8608723546769113029414483825842e-01,
    -4.0584515137739716690660641207707e-01,
    -2.0778495500789846760068940377309e-01,
    0.0,
};
const KRONROD_W = [_]comptime_float{
    2.2935322010529224963732008059913e-02,
    6.3092092629978553290700663189093e-02,
    1.0479001032225018383987632254189e-01,
    1.4065325971552591874518959051021e-01,
    1.6900472663926790282658342659795e-01,
    1.9035057806478540991325640242055e-01,
    2.0443294007529889241416199923466e-01,
    2.0948214108472782801299917489173e-01,
};

pub fn integrate_gauss(
    comptime T: type,
    pointer: anytype,
    comptime f: fn (@TypeOf(pointer), T) T,
    a: T,
    b: T,
) T {
    const scale = (b - a) / 2;
    var sum: T = 0;

    // alternating signs
    inline for (GAUSS_X[0..3], GAUSS_W[0..3]) |xi, wi| {
        const xp = (xi + 1) * scale + a;
        const xm = (-xi + 1) * scale + a;
        const w = wi * scale;

        sum += (f(pointer, xp) + f(pointer, xm)) * w;
    }

    // central point
    const x0 = scale + a;
    const w0 = GAUSS_W[3];

    sum += f(pointer, x0) * w0;
    return sum;
}

pub fn integrate_kronrod(
    comptime T: type,
    pointer: anytype,
    comptime f: fn (@TypeOf(pointer), T) T,
    a: T,
    b: T,
) T {
    const scale = (b - a) / 2;
    var sum: T = 0;

    // alternating signs
    inline for (KRONROD_X[0..7], KRONROD_W[0..7]) |xi, wi| {
        const xp = (xi + 1) * scale + a;
        const xm = (-xi + 1) * scale + a;
        const w = wi * scale;

        sum += (f(pointer, xp) + f(pointer, xm)) * w;
    }

    // central point
    const x0 = scale + a;
    const w0 = KRONROD_W[7];

    sum += f(pointer, x0) * w0;
    return sum;
}

pub fn trapezoid_integration_weight(
    comptime T: type,
    domain: []const T,
    index: usize,
) T {
    // first index
    if (index == 0) {
        return domain[1] - domain[0];
    }
    // last index
    if (index == domain.len - 1) {
        return domain[index] - domain[index - 1];
    }
    // TODO: check this
    return domain[index + 1] - domain[index - 1];
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

