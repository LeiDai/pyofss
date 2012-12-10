
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


# Define exceptions
class SolverError(Exception):
    pass


class NoDerivativeFunctionError(SolverError):
    pass


class Solver(object):
    """
    :param string method: Name of solver method to be used
    :param object f: Derivative function

    Solver consists of both explicit and embedded routines for ODE integration.

    .. note::
      For embedded methods:
      Even though the calculated error applies to A_coarse, it has become
      common to use a weighted combination of A_coarse and A_fine to produce
      a higher order result (local extrapolation).
      It is for this reason that the embedded methods only return A_fine
      rather than A_coarse, or even both (A_coarse, A_fine). To be strict,
      return A_coarse in the methods.
    """
    # The following ssfm_solvers have been replaced with ss_solvers:
    ssfm_solvers = ["ssfm", "ssfm_reduced", "ssfm_sym",
                    "ssfm_sym_midpoint", "ssfm_sym_rk4", "ssfm_sym_rkf"]
    explicit_solvers = ["euler", "midpoint", "rk4"]
    embedded_solvers = ["bs", "rkf", "ck", "dp"]
    ss_solvers = ["ss_simple", "ss_symmetric", "ss_reduced",
                  "ss_agrawal", "ss_sym_midpoint", "ss_sym_rk4"]
    other_solvers = ["rk4ip", ]

    # For embedded methods, use local error of the lower order method:
    errors = {"euler": 2, "midpoint": 3, "rk4": 5,
              "bs": 2, "rkf": 4, "ck": 4, "dp": 4, "ss_simple": 2,
              "ss_symmetric": 3, "ss_reduced": 2, "ss_agrawal": 3,
              "ss_sym_midpoint": 3, "ss_sym_rk4": 5, "rk4ip": 5}

    def __init__(self, method="rk4", f=None):
        if method.lower() in self.embedded_solvers:
            self.embedded = True
        else:
            self.embedded = False

        # Do not use a default for getattr. Better to raise an exception.
        self.method = getattr(self, method.lower())

        # Make sure there is a function/functor attached to f:
        if f is None:
            raise NoDerivativeFunctionError(
                "Require a function for ODE integration")
        else:
            self.f = f

    def __call__(self, A, z, h):
        """ Return A_fine, calculated by method. """
        return self.method(A, z, h, self.f)

    @staticmethod
    def euler(A, z, h, f):
        """ Euler method """
        k0 = h * f(A, z)

        return A + k0

    @staticmethod
    def midpoint(A, z, h, f):
        """ Midpoint method """
        k0 = h * f(A, z)
        k1 = h * f(A + 0.5 * k0, z + 0.5 * h)

        return A + k1

    @staticmethod
    def rk4(A, z, h, f):
        """ Runge-Kutta fourth-order method """
        k0 = h * f(A, z)
        k1 = h * f(A + 0.5 * k0, z + 0.5 * h)
        k2 = h * f(A + 0.5 * k1, z + 0.5 * h)
        k3 = h * f(A + k2, z + h)

        return A + (k0 + 2.0 * (k1 + k2) + k3) / 6.0

    @staticmethod
    def bs(A, z, h, f):
        """ Bogacki-Shampine method (local orders: three and two) """
        p = [1.0 / 2.0, 3.0 / 4.0, 2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0]
        q = [1.0 / 2.0, 3.0 / 4.0, 1.0]
        r = [2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0]
        s = [7.0 / 24.0, 1.0 / 4.0, 1.0 / 3.0, 1.0 / 8.0]

        k0 = h * f(A, z)
        k1 = h * f(A + p[0] * k0, z + q[0] * h)
        k2 = h * f(A + p[1] * k1, z + q[1] * h)
        k3 = h * f(A + p[2] * k0 + p[3] * k1 + p[4] * k2, z + q[2] * h)

        A_fine = A + r[0] * k0 + r[1] * k1 + r[2] * k2
        A_coarse = A + s[0] * k0 + s[1] * k1 + s[2] * k2 + s[3] * k3

        return A_fine, A_coarse

    @staticmethod
    def rkf(A, z, h, f):
        """ Runge-Kutta-Fehlberg method (local orders: four and five)"""
        # List all constants used in this method (Butcher tableau):
        p = [1.0 / 4.0, 3.0 / 32.0, 9.0 / 32.0, 1932.0 / 2197.0,
             -7200.0 / 2197.0, 7296.0 / 2197.0, 439.0 / 216.0, -8.0,
             3680.0 / 513.0, -845.0 / 4104.0, -8.0 / 27.0, 2.0,
             -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0]
        q = [1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5]
        r = [25.0 / 216.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0]
        s = [16.0 / 135.0, 6656.0 / 12825.0, 28561.0 / 56430.0,
             -9.0 / 50.0, 2.0 / 55.0]

        k0 = h * f(A, z)
        k1 = h * f(A + p[0] * k0, z + q[0] * h)
        k2 = h * f(A + p[1] * k0 + p[2] * k1, z + q[1] * h)
        k3 = h * f(A + p[3] * k0 + p[4] * k1 + p[5] * k2, z + q[2] * h)
        k4 = h * f(A + p[6] * k0 + p[7] * k1 + p[8] * k2 +
                   p[9] * k3, z + q[3] * h)
        k5 = h * f(A + p[10] * k0 + p[11] * k1 + p[12] * k2 +
                   p[13] * k3 + p[14] * k4, z + q[4] * h)

        # 4th order (A_coarse) and 5th order (A_fine) approximations:
        A_coarse = A + r[0] * k0 + r[1] * k2 + r[2] * k3 + r[3] * k4
        A_fine = A + s[0] * k0 + s[1] * k2 + s[2] * k3 + s[3] * k4 + s[4] * k5

        return A_fine, A_coarse

    @staticmethod
    def ck(A, z, h, f):
        """ Cash-Karp method (local order: four and five) """
        p = [1.0 / 5.0, 3.0 / 40.0, 9.0 / 40.0, 3.0 / 10.0, -9.0 / 10.0,
             6.0 / 5.0, -11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0,
             1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0,
             44275.0 / 110592.0, 253.0 / 4096.0]
        q = [1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0]
        r = [37.0 / 378.0, 250.0 / 621.0, 125.0 / 594.0, 512.0 / 1771.0]
        s = [2825.0 / 27648.0, 18575.0 / 48384.0, 13525.0 / 55296.0,
             277.0 / 14336.0, 1.0 / 4.0]

        k0 = h * f(A, z)
        k1 = h * f(A + p[0] * k0, z + q[0] * h)
        k2 = h * f(A + p[1] * k0 + p[2] * k1, z + q[1] * h)
        k3 = h * f(A + p[3] * k0 + p[4] * k1 + p[5] * k2, z + q[2] * h)
        k4 = h * f(A + p[6] * k0 + p[7] * k1 + p[8] * k2 +
                   p[9] * k3, z + q[3] * h)
        k5 = h * f(A + p[10] * k0 + p[11] * k1 + p[12] * k2 +
                   p[13] * k3 + p[14] * k4, z + q[4] * h)

        A_coarse = A + r[0] * k0 + r[1] * k2 + r[2] * k3 + r[3] * k5
        A_fine = A + s[0] * k0 + s[1] * k2 + s[2] * k3 + s[3] * k4 + s[4] * k5

        return A_fine, A_coarse

    @staticmethod
    def dp(A, z, h, f):
        """ Dormand-Prince method (local orders: four and five) """
        p = [1.0 / 5.0, 3.0 / 40.0, 9.0 / 40.0, 44.0 / 45.0, -56.0 / 15.0,
             32.0 / 9.0, 19372.0 / 6561.0, -25360.0 / 2187.0,
             64448.0 / 6561.0, -212.0 / 729.0, 9017.0 / 3168.0, -355.0 / 33.0,
             46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 35.0 / 384.0,
             500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0]

        q = [1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0]

        r = [5179.0 / 57600.0, 7571.0 / 16695.0, 393.0 / 640.0,
             -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0]

        s = [35.0 / 384.0, 500.0 / 1113.0, 125.0 / 192.0,
             -2187.0 / 6784.0, 11.0 / 84.0]

        k0 = h * f(A, z)
        k1 = h * f(A + p[0] * k0, z + q[0] * h)
        k2 = h * f(A + p[1] * k0 + p[2] * k1, z + q[1] * h)
        k3 = h * f(A + p[3] * k0 + p[4] * k1 + p[5] * k2, z + q[2] * h)
        k4 = h * f(A + p[6] * k0 + p[7] * k1 + p[8] * k2 + p[9] * k3,
                   z + q[3] * h)
        k5 = h * f(A + p[10] * k0 + p[11] * k1 + p[12] * k2 +
                   p[13] * k3 + p[14] * k4, z + q[4] * h)
        k6 = h * f(A + p[15] * k0 + p[16] * k2 + p[17] * k3 +
                   p[18] * k5 + p[19] * k5, z + q[5] * h)

        A_coarse = A + r[0] * k0 + r[1] * k2 + r[2] * k3 + \
            r[3] * k4 + r[4] * k5 + r[5] * k6
        A_fine = A + s[0] * k0 + s[1] * k2 + s[2] * k3 + s[3] * k4 + s[4] * k5

        return A_fine, A_coarse

    @staticmethod
    def ss_simple(A, z, h, f):
        """ Simple split-step method """
        # Alternative:
        # A_L = f.linear(A, h)
        # return f.nonlinear(A, h, A_L)

        A_N = f.nonlinear(A, h, A)

        return f.linear(A_N, h)

    @staticmethod
    def ss_symmetric(A, z, h, f):
        """ Symmetric split-step method """
        # Alternative:
        # A_N = f.nonlinear(A, 0.5 * h, A)
        # A_L = f.linear(A_N, h)
        # return f.nonlinear(A, 0.5 * h, A_L)

        A_L = f.linear(A, 0.5 * h)
        A_N = f.nonlinear(A, h, A_L)

        return f.linear(A_N, 0.5 * h)

    @staticmethod
    def ss_reduced(A, z, h, f):
        """ Reduced split-step method """
        # Alternative:
        # A_N = f.nonlinear(A, 0.5 * h, A)
        # A_L = f.linear(A_N, h)
        # return f.nonlinear(A_L, 0.5 * h, A_L)

        A_L = f.linear(A, 0.5 * h)
        A_N = f.nonlinear(A_L, h, A_L)

        return f.linear(A_N, 0.5 * h)

    @staticmethod
    def ss_agrawal(A, z, h, f):
        """ Agrawal (iterative) split-step method """
        Y_0 = A
        A_L = f.linear(A, 0.5 * h)

        A_N1 = np.exp(0.5 * h *
                      (np.where(abs(Y_0) > 1.0e-5, f.n(Y_0, z) / Y_0, 0.0) +
                       np.where(abs(A) > 1.0e-5, f.n(A, z) / A, 0.0))) * A_L
        Y_1 = f.linear(A_N1, 0.5 * h)

        A_N2 = np.exp(0.5 * h *
                     (np.where(abs(Y_1) > 1.0e-5, f.n(Y_1, z) / Y_1, 0.0) +
                      np.where(abs(A) > 1.0e-5, f.n(A, z) / A, 0.0))) * A_L
        Y_2 = f.linear(A_N2, 0.5 * h)

        return Y_2

    def ss_sym_midpoint(self, A, z, h, f):
        """ Symmetric split-step method (midpoint method for nonlinear) """
        A_L = f.linear(A, 0.5 * h)
        A_N = self.midpoint(A_L, z, h, f.n)

        return f.linear(A_N, 0.5 * h)

    def ss_sym_rk4(self, A, z, h, f):
        """
        Symmetric split-step method (classical Runge-Kutta for nonlinear)
        """
        A_L = f.linear(A, 0.5 * h)
        A_N = self.rk4(A_L, z, h, f.n)

        return f.linear(A_N, 0.5 * h)

    @staticmethod
    def rk4ip(A, z, h, f):
        """ Runge-Kutta in the interaction picture method """
        # Store half the step-size since it is used often:
        hh = 0.5 * h

        # Transform A into interaction picture:
        A_I = f.linear(A, hh)

        k0 = f.linear(h * f.n(A, z), hh)
        k1 = h * f.n(A_I + 0.5 * k0, z + hh)
        k2 = h * f.n(A_I + 0.5 * k1, z + hh)
        k3 = h * f.n(f.linear(A_I + k2, hh), z + h)

        # Transform back to normal picture (k3 already is) after the step:
        return (k3 / 6.0) + f.linear(A_I + (k0 + 2.0 * (k1 + k2)) / 6.0, hh)
