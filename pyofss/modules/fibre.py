
# ==============================================================================
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
# ==============================================================================

import numpy as np
# ==============================================================================
from linearity import Linearity
from nonlinearity import Nonlinearity
from stepper import Stepper
# ==============================================================================
class Fibre():
   """
   :param string name: Name of this module
   :param double length: Length of fibre
   :param object alpha: Attenuation of fibre
   :param object beta: Dispersion of fibre
   :param double gamma: Nonlinearity of fibre
   :param string sim_type: Type of simulation
   :param Uint traces: Number of field traces required
   :param double local_error: Relative local error used in adaptive stepper
   :param string method: Method to use in ODE solver
   :param Uint total_steps: Number of steps to use for ODE integration
   :param bool self_steepening: Toggles inclusion of self-steepening effects
   :param bool raman_scattering: Toggles inclusion of raman-scattering effects
   :param float rs_factor: Factor determining the amount of raman-scattering
   :param bool use_all: Toggles use of general expression for nonlinearity
   :param double centre_omega: Angular frequency used within dispersion class
   :param double tau_1: Constant used in Raman scattering calculation
   :param double tau_2: Constant used in Raman scattering calculation
   :param double f_R: Constant setting the fraction of Raman scattering used

   sim_type is either default or wdm.

   traces: If greater than 1, will save the field at uniformly-spaced points
   during fibre propagation. If zero, will output all saved points used.
   This is useful if using an adaptive stepper which will likely save
   points non-uniformly.

   method: simulation method such as RK4IP, ARK4IP.

   total_steps: If a non-adaptive stepper is used, this will be used to set the 
   step-size between successive points along the fibre.

   local_error: Relative local error to aim for between propagtion points.
   """
   # ===========================================================================
   def __init__( self, name = "fibre", length = 1.0, alpha = None,
                 beta = None, gamma = 0.0, sim_type = None, traces = 1, 
                 local_error = 1.0e-6, method = "RK4IP", total_steps = 100,
                 self_steepening = False, raman_scattering = False,
                 rs_factor = 0.003, use_all = False, centre_omega = None,
                 tau_1 = 12.2e-3, tau_2 = 32.0e-3, f_R = 0.18 ):

      use_cache = not( method.upper().startswith('A') )
      # ========================================================================
      self.name = name
      self.length = length
      self.linearity = Linearity( alpha, beta, sim_type,
                                  use_cache, centre_omega )
      self.nonlinearity = Nonlinearity( gamma, sim_type, self_steepening,
                                        raman_scattering, rs_factor, use_all,
                                        tau_1, tau_2, f_R )
      # ========================================================================
      # Generate a class to hold linear and nonlinear functions:
      class Function():
         def __init__( self, l, n, linear, nonlinear ):
            self.l = l
            self.n = n
            self.linear = linear
            self.nonlinear = nonlinear
         def __call__( self, A, z ):
            return self.l( A, z ) + self.n( A, z )

      self.function = Function( self.l, self.n, self.linear, self.nonlinear )
      # ========================================================================
      self.stepper = Stepper( traces, local_error, method, self.function,
                              self.length, total_steps )
   # ===========================================================================
   def __call__( self, domain, field ):

      self.linearity( domain )
      self.nonlinearity( domain )

      # Set temporal and spectral arrays for storage:
      self.stepper.storage.t = domain.t
      self.stepper.storage.nu = domain.nu

      # Propagate field through fibre:
      return self.stepper( field )
   # ===========================================================================
   # ===========================================================================
   def l( self, A, z ):
      return self.linearity.lin( A, z )
   def linear( self, A, h ):
      return self.linearity.exp_lin( A, h )
   def n( self, A, z ):
      return self.nonlinearity.non( A, z )
   def nonlinear( self, A, h, B ):
      return self.nonlinearity.exp_non( A, h, B )
# ==============================================================================
if __name__ == "__main__":
   """ 
   Plot the result of a Gaussian pulse propagating through optical fibre.
   Simulates both (third-order) dispersion and nonlinearity.
   Use five different methods: ss_simple, ss_symmetric, ss_sym_rk4, ss_sym_rkf,
   and rk4ip. Expect all five methods to produce similar results; plot traces 
   should all overlap. Separate traces should only be seen at a high zoom level.
   """
   # ===========================================================================
   import numpy as np

   from pyofss import *
   # ===========================================================================
   domain = Domain( bit_width = 200.0, samples_per_bit = 2048 )
   gaussian = Gaussian( peak_power = 1.0, width = 1.0 )
   
   P_ts = []
   methods = ['ss_simple', 'ss_symmetric', 'ss_sym_rk4', 'rk4ip']
   # ===========================================================================
   for m in methods:
      sys = System( domain )
      sys.add( gaussian )
      sys.add( Fibre(length = 5.0, method = m, total_steps = 50,
                     beta = [0.0, 0.0, 0.0, 1.0], gamma = 1.0) )
      sys.run()
      P_ts.append( temporal_power(sys.field) )
   # ===========================================================================
   multi_plot( sys.domain.t, P_ts, methods, labels["t"], labels["P_t"],
               methods, x_range = (80.0, 140.0) )
# ==============================================================================
