const std = @import("std");
const xfunc = @import("transfer-functions.zig");
const util = @import("utils.zig");

const InterpolatingTransferFunction = xfunc.InterpolatingTransferFunction;
const Integrator = @import("Integrator.zig");

/// keeps all of the read in data in a structure
pub fn LineProfileTable(comptime NParams: comptime_int, comptime T: type) type {
    return struct {
        const Self = @This();
        const TFunc = xfunc.TransferFunction(T);

        allocator: std.mem.Allocator,

        // the transfer function data maps g* -> f for different radii
        // these must be the same for all tables of data
        // if data for a specific radius does not exist, it is zeroed
        gstars: []T,

        // there is a matrix of values for each parameter combination
        // there will be Np1 * Np2 * ... * Npn many frames, scaling with the
        // number of parameters, and how its been interpolated
        transfer_functions: []TFunc,

        // we store one array for each parameter with all of its values
        parameters: [NParams][]T,
        n_parameters: [NParams]usize,

        // since the transfer function could be laid out in memory in an
        // inconvenient way, we do all of the interpolations first,
        // and then return a `TranferFunction`, storing indexes to the
        // frames that were used to generate the table, since it is likely
        // we'll just have to reinterpoalte, and can save the parameter lookup.
        last_upper_frame: [2 * NParams]usize = .{0} ** (2 * NParams),
        last_lower_frame: [2 * NParams]usize = .{0} ** (2 * NParams),

        // that does also mean we store the parameters which strided the last used
        // frames, such that the needed parameters P are
        //       last_upper_parameters >= P >= last_lower_parameters
        // which means we need 2 frames pointers for each parameter
        last_lower_parameters: [NParams]usize = .{0} ** NParams,

        // last_upper_parameters is implied as last_lower_parameters[i] + 1

        // then, when we receive some new parameters, we can optimize the search for
        // the new index, and if the same index is used, we already have the frames over
        // which to redo the interpolation

        // dont want to have to allocate new memory each time, so we keep one frame special
        // into which we write the new interpolation and return a copy which does not need
        // to be freed (since we own the memory)
        interpolated_cache: []xfunc.TransferFunction(T),
        interpolated_transfer_function: InterpolatingTransferFunction(T),

        /// Convert the indices of the parameters to an index that may be used
        /// to fetch the corresponding transfer function from `self.transfer_functions`
        /// For parameters p,q, it assumes the table is laid out as
        ///   (p1, q1), (p1, q2), ... (p1, qn), (p2, q1), (p2, q2), ...
        pub fn parameter_indices_to_table_index(self: *const Self, param_indices: [NParams]usize) usize {
            var index: usize = 0;
            for (0..NParams) |i| {
                const stride: usize = if (i < NParams - 1)
                    util.product(usize, self.n_parameters[i + 1 ..])
                else
                    1;
                index += param_indices[i] * stride;
            }
            return index;
        }

        pub fn find_parameter_indices(self: *Self, parameters: [NParams]T) [NParams]usize {
            var base_indices: [NParams]usize = undefined;
            for (0..NParams) |i| {
                const last_i = self.last_lower_parameters[i];

                // todo: parameters are probably going to be pretty close to eachother
                // so might be better to switch the direction of search given the last
                // one, instead of looping back to 0 if it's less than
                const search_up = parameters[i] >= self.parameters[i][last_i];

                const start = if (search_up) last_i else 0;
                const new_i = util.find_first_geq(
                    T,
                    self.parameters[i],
                    parameters[i],
                    start,
                );

                base_indices[i] = if (new_i) |vi| vi.index else 0;
            }
            return base_indices;
        }

        /// Multi-linear interpolation of the parameters. Interpolates first in
        /// q and then in p (in 2d):
        ///
        /// (p1,q0)      (p1,q1)
        ///    \      B  /
        ///     x ----+ x
        ///     ^       |
        ///     |     c |  <- new chosen parameters
        ///     |       |
        ///     x ----+ x
        ///    /      A  \
        /// (p0,q0)    (p0,q1)
        ///  base
        ///
        /// The interpolation over some value y is then
        ///   k = (q1 - q0)
        ///   A = (y(p0, q1) - y(p0, q0)) / k
        ///   B = (y(p1, q1) - y(p1, q0)) / k
        ///   y(c) = A + (B - A) / (p1 - p0)
        ///
        /// The extension to higher dimensions is hopefully apparent.
        pub fn interpolate_parameters(
            self: *Self,
            parameters: [NParams]T,
        ) InterpolatingTransferFunction(T) {
            const base_indices = self.find_parameter_indices(parameters);

            // find which tables we are interpolating between
            const table_indices = self.find_tables(base_indices);

            // now that we know which tables we want to interpolate over
            // we just need to know the weights
            var weights: [NParams]?T = undefined;
            for (0..NParams) |i| {
                weights[i] = self.determine_parameter_weight(
                    &parameters,
                    i,
                    base_indices[i],
                );
            }

            // then interpolate
            // do a few base cases by hand
            switch (NParams) {
                // simple linear case
                1 => {
                    self.interpolate(
                        &self.transfer_functions[table_indices[0]],
                        &self.transfer_functions[table_indices[1]],
                        weights[0],
                        0,
                    );
                },
                // slightly more complex bi-linear case
                2 => {
                    // first (p0,q0) to (p0,q1)
                    self.interpolate(
                        &self.transfer_functions[table_indices[0]],
                        &self.transfer_functions[table_indices[1]],
                        weights[0],
                        0,
                    );
                    // then (p1,q0) to (p1,q1)
                    self.interpolate(
                        &self.transfer_functions[table_indices[2]],
                        &self.transfer_functions[table_indices[3]],
                        weights[0],
                        1,
                    );
                    // and then between those two
                    self.interpolate(
                        &self.interpolated_cache[1],
                        &self.interpolated_cache[0],
                        weights[1],
                        0,
                    );
                },
                // too complex for now tri-linear and above case
                else => @panic("Not implemented yet."),
            }
            return self.interpolated_transfer_function;
        }

        fn interpolate(
            self: *Self,
            from: *const xfunc.TransferFunction(T),
            to: *const xfunc.TransferFunction(T),
            weight: ?T,
            cache_index: usize,
        ) void {
            if (weight) |w| {
                from.interpolateBetween(
                    to,
                    &self.interpolated_cache[cache_index],
                    w,
                );
            } else {
                self.interpolated_cache[cache_index].assignFrom(from);
            }
        }

        pub fn determine_parameter_weight(
            self: *const Self,
            parameters: []const T,
            parameter_index: usize,
            index: usize,
        ) ?T {
            // no interpolation needs to happen here
            if (index == 0) return null;
            const p0 = self.parameters[parameter_index][index - 1];
            const p1 = self.parameters[parameter_index][index];
            return (parameters[parameter_index] - p0) / (p1 - p0);
        }

        pub fn find_tables(
            self: *const Self,
            parameter_indices: [NParams]usize,
        ) [NParams * 2]usize {
            // mutable copy of the parameter_indices
            var pindices: [NParams]usize = undefined;
            std.mem.copyForwards(usize, &pindices, &parameter_indices);

            // store the table striding the chosen parameter interpolation
            var table_indices: [NParams * 2]usize = undefined;

            // TODO: this only works for NParams == 2
            table_indices[1] = self.parameter_indices_to_table_index(pindices);
            pindices[0] -|= 1;
            table_indices[0] = self.parameter_indices_to_table_index(pindices);
            pindices[1] -|= 1;
            table_indices[2] = self.parameter_indices_to_table_index(pindices);
            if (pindices[0] < parameter_indices[0]) {
                pindices[0] += 1;
            }
            table_indices[3] = self.parameter_indices_to_table_index(pindices);
            return table_indices;
        }

        pub fn init(
            alloc: std.mem.Allocator,
            gstars: []T,
            transfer_functions: []TFunc,
            parameters: [NParams][]T,
        ) !Self {
            // create the cache
            // needs to be on the heap so we can point to it
            // from the interpolation structs
            var cache = try alloc.alloc(xfunc.TransferFunction(T), NParams);
            errdefer alloc.free(cache);
            errdefer for (cache) |*tf| tf.free(alloc);

            for (0..NParams) |i| {
                cache[i] = try transfer_functions[0].copy(alloc);
            }
            // make a single interpolated transfer function regardless of
            // NParams
            const itf = try InterpolatingTransferFunction(T).init(
                gstars,
                alloc,
                // point to the cached function
                &cache[0],
            );

            // count how many of each parameter we have
            var n_params: [NParams]usize = undefined;
            for (0..NParams) |i| {
                n_params[i] = parameters[i].len;
            }

            return .{
                .allocator = alloc,
                .gstars = gstars,
                .transfer_functions = transfer_functions,
                .parameters = parameters,
                .n_parameters = n_params,
                .interpolated_cache = cache,
                .interpolated_transfer_function = itf,
            };
        }

        // standard free everything approach
        pub fn deinit(self: *Self) void {
            self.allocator.free(self.gstars);
            // free the parameter arrays
            for (self.parameters) |p| {
                self.allocator.free(p);
            }
            // free the transfer functions
            for (self.transfer_functions) |*tf| {
                tf.free(self.allocator);
            }
            self.allocator.free(self.transfer_functions);
            // free the interpolant cache
            self.interpolated_transfer_function.free(self.allocator);
            for (self.interpolated_cache) |*tf| tf.free(self.allocator);
            self.allocator.free(self.interpolated_cache);
        }
    };
}

/// facilitates the integration over arbitrary parameter, energy
/// and flux arrays
pub fn MultiEmissivityLineProfile(
    comptime NParams: comptime_int,
    comptime T: type,
) type {
    return struct {
        table: LineProfileTable(NParams, T),
    };
}
