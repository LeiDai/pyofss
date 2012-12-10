
"""
    Copyright (C) 2011, 2012  David Bolt

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
from scipy import linalg

from storage import Storage
from solver import Solver


class Stepper(object):
    """
    :param Uint traces: Number of ouput trace to use
    :param double local_error: Relative local error required in adaptive method
    :param string method: ODE solver method to use
    :param object f: Derivative function to be solved
    :param double length: Length to integrate over
    :param Uint total_steps: Number of steps to use for ODE integration

    method:
      * EULER -- Euler method;
      * MIDPOINT -- Midpoint method;
      * RK4 -- Fourth order Runge-Kutta method;
      * BS -- Bogacki-Shampine method;
      * RKF -- Runge-Kutta-Fehlberg method;
      * CK -- Cash-Karp method;
      * DP -- Dormand-Prince method;
      * SS_SIMPLE -- Simple split-step method;
      * SS_SYMMETRIC -- Symmetric split-step method;
      * SS_REDUCED -- Reduced split-step method;
      * SS_AGRAWAL -- Agrawal (iterative) split-step method;
      * SS_SYM_MIDPOINT -- Symmetric split-step method (MIDPOINT for nonlinear)
      * SS_SYM_RK4 -- Symmetric split-step method (RK4 for nonlinear);
      * RK4IP -- Runge-Kutta in the interaction picture method.

      Each method may use an adaptive stepper by prepending an 'A' to the name.

    traces:
      * 0 -- Store A for each succesful step;
      * 1 -- Store A at final value (length) only;
      * >1 -- Store A for each succesful step then use interpolation to get A
         values for equally spaced z-values, calculated using traces.
    """
    def __init__(self, traces=1, local_error=1.0e-6, method="RK4",
                 f=None, length=1.0, total_steps=100):
        self.traces = traces
        self.local_error = local_error

        # Check if adaptive stepsize is required:
        if method.upper().startswith('A'):
            self.adaptive = True
            self.method = method[1:]
        else:
            self.adaptive = False
            self.method = method

        #~print "Using {0} method".format( self.method )

        # Delegate method and function to solver
        self.solver = Solver(self.method, f)
        self.step = self.solver

        self.length = length
        self.total_steps = total_steps

        # Use a list of tuples ( z, A(z) ) for dense output if required:
        self.storage = Storage()

        # Store constants for adaptive method:
        self.total_attempts = 100
        self.steps_max = 50000

        self.safety = 0.9
        self.max_factor = 10.0
        self.min_factor = 0.2

        # Store local error of method:
        self.eta = self.solver.errors[self.method.lower()]

        self.A_out = None

    def __call__(self, A):
        """ Delegate to appropriate function, adaptive- or standard-stepper """

        self.storage.reset_fft_counter()

        if self.adaptive:
            return self.adaptive_stepper(A)
        else:
            return self.standard_stepper(A)

    def standard_stepper(self, A):
        """ Take a fixed number of steps, each of equal length """
        #~print( "Starting ODE integration with fixed step-size... " ),

        # Initialise:
        self.A_out = A

        # Require an initial step-size:
        h = self.length / self.total_steps

        # Construct mesh points for z:
        zs = np.linspace(0.0, self.length, self.total_steps + 1)

        # Construct mesh points for traces:
        if self.traces != self.total_steps:
            trace_zs = np.linspace(0.0, self.length, self.traces + 1)

        # Make sure to store the initial A if more than one trace is required:
        if self.traces != 1:
            self.storage.append(zs[0], self.A_out)

        # Start at z = 0.0 and repeat until z = length - h (inclusive),
        # i.e. z[-1]
        for z in zs[:-1]:
            # Currently at L = z
            if self.solver.embedded:
                self.A_out, A_other = self.step(self.A_out, z, h)
            else:
                self.A_out = self.step(self.A_out, z, h)
            # Now at L = z + h

            # If multiple traces required, store A_out at each relavant z
            # value:
            if self.traces != 1:
                self.storage.append(z + h, self.A_out)

        # Store total number of fft and ifft operations that were used:
        self.storage.store_current_fft_count()

        # Need to interpolate dense output to grid points set by traces:
        if self.traces > 1 and (self.traces != self.total_steps):
            self.storage.interpolate_As_for_z_values(trace_zs)

        return self.A_out

    @staticmethod
    def relative_local_error(A_fine, A_coarse):
        """ Calculate an estimate of the relative local error """

        norm_fine = linalg.norm(A_fine)

        # Avoid possible divide by zero:
        if norm_fine != 0.0:
            return linalg.norm(A_fine - A_coarse) / norm_fine
        else:
            return linalg.norm(A_fine - A_coarse)

    def adaptive_stepper(self, A):
        """ Take multiple steps, with variable length, until target reached """

        #~print( "Starting ODE integration with adaptive step-size... " ),

        # Initialise:
        self.A_out = A
        z = 0.0

        # Require an initial step-size which will be adapted by the routine:
        if self.traces > self.total_steps:
            h = self.length / self.traces
        else:
            h = self.length / self.total_steps

        # Constants used for approximation of solution using local
        # extrapolation:
        f_eta = np.power(2, self.eta - 1.0)
        f_alpha = f_eta / (f_eta - 1.0)
        f_beta = 1.0 / (f_eta - 1.0)

        # Calculate z-values at which to save traces.
        if self.traces > 1:
            # zs contains z-values for each trace, as well as the initial
            # trace:
            zs = np.linspace(0.0, self.length, self.traces + 1)

        # Store initial trace:
        if self.traces != 1:
            self.storage.append(z, self.A_out)

        # Limit the number of steps in case of slowly converging runs:
        for s in range(1, self.steps_max):
            # If step-size takes z our of range [0.0, length], then correct it:
            if (z + h) > self.length:
                h = self.length - z

            # Take an adaptive step:
            for ta in range(0, self.total_attempts):
                h_half = 0.5 * h
                z_half = z + h_half

                # Calculate A_fine and A_coarse internally if using an
                # embedded method. Otherwise use method of step-doubling:
                if self.solver.embedded:
                    A_fine, A_coarse = self.step(self.A_out, z, h)
                else:
                    # Calculate fine solution using two steps of size h_half:
                    A_half = self.step(self.A_out, z, h_half)
                    A_fine = self.step(A_half, z_half, h_half)
                    # Calculate coarse solution using one step of size h:
                    A_coarse = self.step(self.A_out, z, h)

                # Calculate an estimate of relative local error:
                delta = self.relative_local_error(A_fine, A_coarse)

                # Store current stepsize:
                h_temp = h

                # Adjust stepsize for next step:
                if delta > 0.0:
                    error_ratio = (self.local_error / delta)
                    factor = \
                        self.safety * np.power(error_ratio, 1.0 / self.eta)
                    h = h_temp * min(self.max_factor,
                                     max(self.min_factor, factor))
                else:
                    # Error approximately zero, so use largest stepsize
                    # increase:
                    h = h_temp * self.max_factor

                if delta < 2.0 * self.local_error:
                    # Successful step, so increment z h_temp (which is the
                    # stepsize that was used for this step):
                    z += h_temp

                    if self.solver.embedded:
                        # Accept the higher order method:
                        self.A_out = A_fine
                    else:
                        # Use local extrapolation to form a higher order
                        # solution:
                        self.A_out = f_alpha * A_fine - f_beta * A_coarse

                    # Store data on current z and stepsize used for each
                    # succesful step:
                    self.storage.step_sizes.append((z, h_temp))

                    # Most dense storage (stores a trace for each successful
                    # step):
                    if self.traces != 1:
                        self.storage.append(z, self.A_out)

                    break  # Successful attempt at step, move on to next step.

                # Otherwise error was too large, continue with next attempt.

            else:
                raise Exception("Failed to set suitable step-size")

            # If the desired z has been reached, then finish:
            if z >= self.length:
                # Store total number of fft and ifft operations that were used:
                self.storage.store_current_fft_count()

                # Interpolate dense output to uniformly-spaced z values:
                if self.traces > 1:
                    self.storage.interpolate_As_for_z_values(zs)

                return self.A_out

        raise Exception("Failed to complete with maximum steps allocated")

if __name__ == "__main__":
    """
    Exact solution: A(z) = 0.5 * ( 5.0 * exp(-2.0 * z) - 3.0 * exp(-4.0 * z) )
    A(0) = 1.0
    A(0.5) = 0.71669567807368684
    Numerical solution (RK4, total_steps = 5):      0.71668876283331295
    Numerical solution (RK4, total_steps = 50):     0.71669567757603803
    Numerical solution (RK4, total_steps = 500):    0.71669567807363854
    Numerical solution (RKF, total_steps = 5):      0.71669606109336026
    Numerical solution (RKF, total_steps = 50):     0.71669567807672185
    Numerical solution (RKF, total_steps = 500):    0.71669567807368773
    """
    import matplotlib.pyplot as plt

    def simple(A, z):
        """ Just a simple function. """
        return 3.0 * np.exp(-4.0 * z) - 2.0 * A

    stepper = Stepper(f=simple, length=0.5, total_steps=50,
                      method="RKF", traces=50)
    A_in = 1.0
    A_out = stepper(A_in)
    print "A_out = %.17f" % (A_out)

    x = stepper.storage.z
    y = stepper.storage.As

    title = r'$\frac{dA}{dz} + 2A = 3 e^{-4z}$'
    plt.title(r'Numerical integration of ODE:' + title)
    plt.xlabel('z')
    plt.ylabel('A(z)')
    plt.plot(x, y, label='RKF: 50 steps')
    plt.legend()
    plt.show()
