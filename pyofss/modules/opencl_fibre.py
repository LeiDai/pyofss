
"""
    Copyright (C) 2013 David Bolt

    This file is part of pyofss.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array

from string import Template

OPENCL_OPERATIONS = Template("""
#ifdef cl_khr_fp64 // Khronos extension
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64) // AMD extension
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
    #error "Double precision floating point is not supported"
#endif

#define PYOPENCL_DEFINE_CDOUBLE
#include <pyopencl-complex.h>

__kernel void cl_cache(__global c${dorf}_t* factor,
                       const ${dorf} h) {
    int gid = get_global_id(0);

    factor[gid] = c${dorf}_exp(h * factor[gid]);
}

__kernel void cl_linear_cached(__global c${dorf}_t* A,
                               __global c${dorf}_t* factor) {
    int gid = get_global_id(0);

    A[gid] = c${dorf}_mul(A[gid], factor[gid]);
}

__kernel void cl_linear(__global c${dorf}_t* A,
                        __global c${dorf}_t* factor,
                        const ${dorf} h ) {
    int gid = get_global_id(0);

    // IMPORTANT: Cannot just multiply two complex numbers!
    // Must use the appropriate function (e.g. cfloat_mul).
    A[gid] = c${dorf}_mul(A[gid], c${dorf}_exp(h * factor[gid]));
}

c${dorf}_t cl_square_abs(const c${dorf}_t a) {
    return c${dorf}_mul(a, c${dorf}_conj(a));
}

__kernel void cl_nonlinear(__global c${dorf}_t* A,
                           const ${dorf} gamma,
                           const ${dorf} h) {
    int gid = get_global_id(0);

    const c${dorf}_t im_gamma = (c${dorf}_t)(0.0, h * gamma);

    A[gid] = c${dorf}_mul(A[gid],
                          c${dorf}_mul(im_gamma, cl_square_abs(A[gid])));
}

__kernel void cl_sum(__global c${dorf}_t* A,
                     const ${dorf} A_factor,
                     __global c${dorf}_t* B,
                     const ${dorf} B_factor) {
    int gid = get_global_id(0);

    A[gid] *= A_factor;
    A[gid] += B_factor * B[gid];
}

__kernel void cl_copy(__global c${dorf}_t* A,
                      __global c${dorf}_t* B) {
    int gid = get_global_id(0);

    A[gid] = B[gid];
}

""")


class OpenclFibre(object):
    """
    This optical module is similar to Fibre, but uses PyOpenCl (Python
    bindings around OpenCL) to generate parallelised code.
    """
    def __init__(self, name="clfibre", dorf="float"):
        self.queue = None
        self.np_float = None
        self.np_complex = None
        self.prg = None
        self.cl_initialise(dorf)

        self.name = name

        self.buf_field = None
        self.buf_temp = None
        self.buf_interaction = None
        self.buf_factor = None

    def __call__(self, domain, field):
        pass

    def cl_initialise(self, dorf="float"):
        """ Initialise opencl related parameters. """
        float_conversions = {"float": np.float32, "double": np.float64}
        complex_conversions = {"float": np.complex64, "double": np.complex128}

        self.np_float = float_conversions[dorf]
        self.np_complex = complex_conversions[dorf]

        for platform in cl.get_platforms():
            if platform.name == "NVIDIA CUDA":
                print("Using compiler optimisations suitable for Nvidia GPUs")
                compiler_options = "-cl-mad-enable -cl-fast-relaxed-math"
            else:
                compiler_options = ""

        ctx = cl.create_some_context(interactive=False)
        self.queue = cl.CommandQueue(ctx)

        substitutions = {"dorf": dorf}
        code = OPENCL_OPERATIONS.substitute(substitutions)
        self.prg = cl.Program(ctx, code).build(options=compiler_options)

    @staticmethod
    def print_device_info():
        """ Output information on each OpenCL platform and device. """
        for platform in cl.get_platforms():
            print "=" * 60
            print "Platform information:"
            print "Name: ", platform.name
            print "Profile: ", platform.profile
            print "Vender: ", platform.vendor
            print "Version: ", platform.version

            for device in platform.get_devices():
                print "-" * 60
                print "Device information:"
                print "Name: ", device.name
                print "Type: ", cl.device_type.to_string(device.type)
                print "Memory: ", device.global_mem_size // (1024 ** 2), "MB"
                print "Max clock speed: ", device.max_clock_frequency, "MHz"
                print "Compute units: ", device.max_compute_units

            print "=" * 60

    def send_arrays_to_device(self, field, field_temp, field_interaction,
                              factor):
        """ Move numpy arrays onto compute device. """
        self.buf_field = cl_array.to_device(
            self.queue, field.astype(self.np_complex))
        self.buf_temp = cl_array.to_device(
            self.queue, field_temp.astype(self.np_complex))
        self.buf_interaction = cl_array.to_device(
            self.queue, field_interaction.astype(self.np_complex))

        self.buf_factor = cl_array.to_device(
            self.queue, factor.astype(self.np_complex))

if __name__ == "__main__":
    # Compare simulations using Fibre and OpenclFibre modules.
    OCL_FIBRE = OpenclFibre(dorf="float")
    OCL_FIBRE.print_device_info()