
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

from pyfft.cl import Plan
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
                           const ${dorf} stepsize) {
        int gid = get_global_id(0);

        factor[gid] = c${dorf}_exp(stepsize * factor[gid]);
    }

    __kernel void cl_linear_cached(__global c${dorf}_t* field,
                                   __global c${dorf}_t* factor) {
        int gid = get_global_id(0);

        field[gid] = c${dorf}_mul(field[gid], factor[gid]);
    }

    __kernel void cl_linear(__global c${dorf}_t* field,
                            __global c${dorf}_t* factor,
                            const ${dorf} stepsize) {
        int gid = get_global_id(0);

        // IMPORTANT: Cannot just multiply two complex numbers!
        // Must use the appropriate function (e.g. cfloat_mul).
        field[gid] = c${dorf}_mul(field[gid],
                                  c${dorf}_exp(stepsize * factor[gid]));
    }

    c${dorf}_t cl_square_abs(const c${dorf}_t element) {
        return c${dorf}_mul(element, c${dorf}_conj(element));
    }

    __kernel void cl_nonlinear(__global c${dorf}_t* field,
                               const ${dorf} gamma,
                               const ${dorf} stepsize) {
        int gid = get_global_id(0);

        const c${dorf}_t im_gamma = (c${dorf}_t)(0.0, stepsize * gamma);

        field[gid] = c${dorf}_mul(
            field[gid], c${dorf}_mul(im_gamma, cl_square_abs(field[gid])));
    }

    __kernel void cl_sum(__global c${dorf}_t* first_field,
                         const ${dorf} first_factor,
                         __global c${dorf}_t* second_field,
                         const ${dorf} second_factor) {
        int gid = get_global_id(0);

        first_field[gid] *= first_factor;
        first_field[gid] += (second_factor * second_field[gid]);
    }

    __kernel void cl_copy(__global c${dorf}_t* first_field,
                          __global c${dorf}_t* second_field) {
        int gid = get_global_id(0);

        first_field[gid] = second_field[gid];
    }
""")


class OpenclFibre(object):
    """
    This optical module is similar to Fibre, but uses PyOpenCl (Python
    bindings around OpenCL) to generate parallelised code.
    """
    def __init__(self, total_samples, dorf, length=None, total_steps=None,
                 name="ocl_fibre"):
        self.name = name

        self.queue = None
        self.np_float = None
        self.np_complex = None
        self.prg = None
        self.cl_initialise(dorf)

        self.plan = Plan(total_samples, queue=self.queue)

        self.buf_field = None
        self.buf_temp = None
        self.buf_interaction = None
        self.buf_factor = None

        self.shape = None
        self.plan = None

        self.cached_factor = False
        # Force usage of cached version of function:
        self.cl_linear = self.cl_linear_cached

        self.length = length
        self.total_steps = total_steps

    def __call__(self, domain, field):
        # Setup plan for calculating fast Fourier transforms:
        self.plan = Plan(domain.total_samples, queue=self.queue)

        field_temp = np.empty_like(field)
        field_interaction = np.empty_like(field)

        from pyofss.modules.linearity import Linearity
        dispersion = Linearity(beta=[0.0, 0.0, 0.0, 1.0], sim_type="default")
        factor = dispersion(domain)

        self.send_arrays_to_device(
            field, field_temp, field_interaction, factor)

        stepsize = self.length / self.total_steps
        zs = np.linspace(0.0, self.length, self.total_steps + 1)

        #start = time.clock()
        for z in zs[:-1]:
            self.cl_rk4ip(self.buf_field, self.buf_temp,
                          self.buf_interaction, self.buf_factor, stepsize)
        #stop = time.clock()

        #cl_result = self.buf_field.get()
        #print("cl_result: %e" % ((stop - start) / 1000.0))

        return self.buf_field.get()

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

    def send_arrays_to_device(self, field, field_temp,
                              field_interaction, factor):
        """ Move numpy arrays onto compute device. """
        self.shape = field.shape

        self.buf_field = cl_array.to_device(
            self.queue, field.astype(self.np_complex))
        self.buf_temp = cl_array.to_device(
            self.queue, field_temp.astype(self.np_complex))
        self.buf_interaction = cl_array.to_device(
            self.queue, field_interaction.astype(self.np_complex))

        self.buf_factor = cl_array.to_device(
            self.queue, factor.astype(self.np_complex))

    def cl_copy(self, first_buffer, second_buffer):
        """ Copy contents of one buffer into another. """
        self.prg.cl_copy(self.queue, self.shape, None,
                         first_buffer.data, second_buffer.data)

    def cl_linear(self, field_buffer, stepsize, factor_buffer):
        """ Linear part of step. """
        self.plan.execute(field_buffer.data, inverse=True)
        self.prg.cl_linear(self.queue, self.shape, None, field_buffer.data,
                           factor_buffer.data, self.np_float(stepsize))
        self.plan.execute(field_buffer.data)

    def cl_linear_cached(self, field_buffer, stepsize, factor_buffer):
        """ Linear part of step (cached version). """
        if (self.cached_factor is False):
            print "Caching factor"
            self.prg.cl_cache(self.queue, self.shape, None,
                              factor_buffer.data, self.np_float(stepsize))
            self.cached_factor = True

        self.plan.execute(field_buffer.data, inverse=True)
        self.prg.cl_linear_cached(self.queue, self.shape, None,
                                  field_buffer.data, factor_buffer.data)
        self.plan.execute(field_buffer.data)

    def cl_nonlinear(self, field_buffer, stepsize, gamma=100.0):
        """ Nonlinear part of step. """
        self.prg.cl_nonlinear(self.queue, self.shape, None, field_buffer.data,
                              self.np_float(gamma), self.np_float(stepsize))

    def cl_sum(self, first_buffer, first_factor, second_buffer, second_factor):
        """ Calculate weighted summation. """
        self.prg.cl_sum(self.queue, self.shape, None,
                        first_buffer.data, self.np_float(first_factor),
                        second_buffer.data, self.np_float(second_factor))

    def cl_rk4ip(self, field, field_temp, field_interaction, factor, stepsize):
        """ Runge-Kutta in the interaction picture method using OpenCL. """
        inv_six = 1.0 / 6.0
        inv_three = 1.0 / 3.0
        half_step = 0.5 * stepsize

        self.cl_copy(field_temp, field)
        self.cl_linear(field, half_step, factor)

        self.cl_copy(field_interaction, field)
        self.cl_nonlinear(field_temp, stepsize)
        self.cl_linear(field_temp, half_step, factor)

        self.cl_sum(field, 1.0, field_temp, inv_six)
        self.cl_sum(field_temp, 0.5, field_interaction, 1.0)
        self.cl_nonlinear(field_temp, stepsize)

        self.cl_sum(field, 1.0, field_temp, inv_three)
        self.cl_sum(field_temp, 0.5, field_interaction, 1.0)
        self.cl_nonlinear(field_temp, stepsize)

        self.cl_sum(field, 1.0, field_temp, inv_three)
        self.cl_sum(field_temp, 1.0, field_interaction, 1.0)
        self.cl_linear(field_interaction, half_step, factor)

        self.cl_linear(field, half_step, factor)
        self.cl_nonlinear(field_temp, stepsize)

        self.cl_sum(field, 1.0, field_temp, inv_six)

if __name__ == "__main__":
    # Compare simulations using Fibre and OpenclFibre modules.
    from pyofss import Domain, System, Gaussian, Fibre
    from pyofss import temporal_power, double_plot, labels

    import time

    TS = 4096
    GAMMA = 100.0
    STEPS = 800
    LENGTH = 0.1

    DOMAIN = Domain(bit_width=30.0, samples_per_bit=TS)

    SYS = System(DOMAIN)
    SYS.add(Gaussian("gaussian", peak_power=1.0, width=1.0))
    SYS.add(Fibre("fibre", beta=[0.0, 0.0, 0.0, 1.0], gamma=GAMMA,
                  length=LENGTH, total_steps=STEPS, method="RK4IP"))

    start = time.clock()
    SYS.run()
    stop = time.clock()
    NO_OCL_DURATION = (stop - start) / 1000.0
    NO_OCL_OUT = SYS.fields["fibre"]

    sys = System(DOMAIN)
    sys.add(Gaussian("gaussian", peak_power=1.0, width=1.0))
    sys.add(OpenclFibre(TS, dorf="float", length=LENGTH, total_steps=STEPS))

    start = time.clock()
    sys.run()
    stop = time.clock()
    OCL_DURATION = (stop - start) / 1000.0
    OCL_OUT = sys.fields["ocl_fibre"]

    NO_OCL_POWER = temporal_power(NO_OCL_OUT)
    OCL_POWER = temporal_power(OCL_OUT)
    DELTA_POWER = NO_OCL_POWER - OCL_POWER

    MEAN_RELATIVE_ERROR = np.mean(np.abs(DELTA_POWER))
    MEAN_RELATIVE_ERROR /= np.amax(temporal_power(NO_OCL_OUT))

    print("Run time without OpenCL: %e" % NO_OCL_DURATION)
    print("Run time with OpenCL: %e" % OCL_DURATION)
    print("Mean relative error: %e" % MEAN_RELATIVE_ERROR)

    # Expect both plots to appear identical:
    double_plot(SYS.domain.t, NO_OCL_POWER, SYS.domain.t, OCL_POWER,
                x_label=labels["t"], y_label=labels["P_t"],
                X_label=labels["t"], Y_label=labels["P_t"])
