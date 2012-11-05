
# ==============================================================================
"""
    Copyright (C) 2012  David Bolt

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

from pyofss.modules.solver import Solver
from pyofss.modules.solver import NoDerivativeFunctionError

from pyofss.modules.stepper import Stepper

from numpy.testing.utils import assert_almost_equal # uses decimal places
from numpy.testing.utils import assert_equal

import numpy as np
from scipy.fftpack import fft, ifft, ifftshift

import unittest2
# ==============================================================================
def default_function( A, z ):
   return A
# ==============================================================================
class Default_parameters( unittest2.TestCase ):
   def test_no_derivative_function( self ):
      """ Should raise exception if derivative function is None (default) """
      self.assertRaises( NoDerivativeFunctionError, Solver )
   # ===========================================================================
   def test_default( self ):
      """ Should use default values (except for f) """
      solver = Solver( f = default_function )
      # ========================================================================
      self.assertEqual( solver.method, solver.rk4, \
                        "method should default to rk4" )
      self.assertFalse( solver.embedded, "RK4 is not an embedded method" )
   # ===========================================================================
   def test_attributes( self ):
      """ Should list ODE integration methods """
      self.assertEqual( Solver.explicit_solvers, ["euler", "midpoint", "rk4"] )
      self.assertEqual( Solver.embedded_solvers, ["bs", "rkf", "ck", "dp"] )
      self.assertEqual( Solver.ss_solvers, \
                        ["ss_simple", "ss_symmetric", "ss_reduced", \
                         "ss_agrawal", "ss_sym_midpoint", "ss_sym_rk4" ] )
      self.assertEqual( Solver.other_solvers, [ "rk4ip" ] )
# ==============================================================================
class Check_core_routines( unittest2.TestCase ):
   r"""
   ODE: $\frac{ dA }{ dz } + 2A - 3 \exp{-4z} = 0$
   Initial value: A(0) = 1.0
   Exact solution: $A(z) = \frac{5 \exp{-2z} - 3 \exp{-4z} }{ 2 }$
   A(0.5) = 0.71669567807368684
   """
   def simple( self, A, z ):
      """ Example of derivative to use for integration """
      return 3.0 * np.exp(-4.0 * z) - 2.0 * A
   # ===========================================================================
   def setUp( self ):
      """ Store common parameters for each test """
      self.parameters = { "traces" : 1, "local_error" : 1.0e-6,
                          "method" : None, "f" : self.simple, 
                          "length" : 0.5, "total_steps" : 50 }
      self.A_initial = 1.0
      self.A_analytical = 0.71669567807368684
   # ===========================================================================
   def test_euler( self ):
      """ A_numerical should be approximately equal to A_analytical """
      self.parameters[ "method" ] = "EULER"
      # ========================================================================
      self.parameters[ "total_steps" ] = 50
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_equal( A_numerical, 0.72152686293820223 )
      # ========================================================================
      self.parameters[ "total_steps" ] = 500
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_equal( A_numerical, 0.71717199224362427 )
   # ===========================================================================
   def test_midpoint( self ):
      """ A_numerical should be approximately equal to A_analytical """
      self.parameters[ "method" ] = "MIDPOINT"
      # ========================================================================
      self.parameters[ "total_steps" ] = 50
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_equal( A_numerical, 0.71663952968427402 )
      # ========================================================================
      self.parameters[ "total_steps" ] = 500
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_equal( A_numerical, 0.71669512731047913 )
   # ===========================================================================
   def test_rk4( self ):
      """ A_numerical should be approximately equal to A_analytical """
      self.parameters[ "method" ] = "RK4"
      # ========================================================================
      self.parameters[ "total_steps" ] = 50
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_almost_equal( A_numerical, self.A_analytical )
      assert_equal( A_numerical, 0.71669567757603803 )
      # ========================================================================
      self.parameters[ "total_steps" ] = 500
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_almost_equal( A_numerical, self.A_analytical )
      assert_equal( A_numerical, 0.71669567807363854 )
   # ===========================================================================
   def test_bs( self ):
      """ A_numerical should be approximately equal to A_analytical """
      self.parameters[ "method" ] = "BS"
      # ========================================================================
      self.parameters[ "total_steps" ] = 50
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_equal( A_numerical, 0.71669591977945657 )
      # ========================================================================
      self.parameters[ "total_steps" ] = 500
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_equal( A_numerical, 0.71669567831028291 )
   # ===========================================================================
   def test_rkf( self ):
      """ A_numerical should be approximately equal to A_analytical """
      self.parameters[ "method" ] = "RKF"
      # ========================================================================
      self.parameters[ "total_steps" ] = 50
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_almost_equal( A_numerical, self.A_analytical )
      assert_equal( A_numerical, 0.71669567807672185 )
      # ========================================================================
      self.parameters[ "total_steps" ] = 500
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_almost_equal( A_numerical, self.A_analytical )
      assert_equal( A_numerical, 0.71669567807368773 )
   # ===========================================================================
   def test_ck( self ):
      """ A_numerical should be approximately equal to A_analytical """
      self.parameters[ "method" ] = "CK"
      # ========================================================================
      self.parameters[ "total_steps" ] = 50
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_almost_equal( A_numerical, self.A_analytical )
      assert_equal( A_numerical, 0.71669567804234302 )
      # ========================================================================
      self.parameters[ "total_steps" ] = 500
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_almost_equal( A_numerical, self.A_analytical )
      assert_equal( A_numerical, 0.71669567807368351 )
   # ===========================================================================
   def test_dp( self ):
      """ A_numerical should be approximately equal to A_analytical """
      self.parameters[ "method" ] = "DP"
      # ========================================================================
      self.parameters[ "total_steps" ] = 50
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_almost_equal( A_numerical, self.A_analytical )
      assert_equal( A_numerical, 0.71669567807314427 )
      # ========================================================================
      self.parameters[ "total_steps" ] = 500
      stepper = Stepper( **self.parameters )
      A_numerical = stepper( self.A_initial )
      assert_almost_equal( A_numerical, self.A_analytical )
      assert_equal( A_numerical, 0.71669567807368462 )
# ==============================================================================
class Check_methods( unittest2.TestCase ):
   """
   Soliton example with periodic evolution. Fibre length is selected such that 
   the peak power is the same as that of the initial pulse.
   """
   class Function():
      def __init__( self, domain ):
         self.beta_2 = -1.0
         self.gamma = 1.0
         # =====================================================================
         self.Domega = domain.omega - domain.centre_omega
         self.factor = 1j * ifftshift( 0.5 * self.beta_2 * self.Domega**2 )
      # ========================================================================
      def l( self, A, z ):
         return fft( self.factor * ifft(A) )
      # ========================================================================
      def n( self, A, z ):
         return 1j * self.gamma * np.abs(A)**2 * A
      # ========================================================================
      def linear( self, A, h ):
         return fft( np.exp( h * self.factor ) * ifft(A) )
      # ========================================================================
      def nonlinear( self, A, h, B ):
         return np.exp( h * 1j * self.gamma * np.abs(A)**2 ) * B
   # ===========================================================================
   def setUp( self ):
      from pyofss.domain import Domain
      from pyofss.modules.sech import Sech
      # ========================================================================
      domain = Domain( bit_width = 200.0, samples_per_bit = 2048 )
      sech = Sech( peak_power = 4.0, width = 1.0 )
      # ========================================================================
      function = self.Function( domain )
      self.parameters = { "method" : None, "f" : function, 
                          "length" : 0.5 * np.pi, "total_steps" : 1000 }
      # ========================================================================
      self.A_analytical = 2.0
      self.P_analytical = np.abs( self.A_analytical )**2
      self.A_in = sech.generate( domain.t )
   # ===========================================================================
   def test_ss_simple( self ):
      self.parameters[ "method" ] = "SS_SIMPLE"
      stepper = Stepper( **self.parameters )
      A_out = stepper( self.A_in )
      # ========================================================================
      assert_almost_equal( max( np.abs(A_out)**2 ), self.P_analytical, 3 )
   # ===========================================================================
   def test_ss_symmetric( self ):
      self.parameters[ "method" ] = "SS_SYMMETRIC"
      stepper = Stepper( **self.parameters )
      A_out = stepper( self.A_in )
      # ========================================================================
      assert_almost_equal( max( np.abs(A_out)**2 ), self.P_analytical, 0 )
   # ===========================================================================
   def test_ss_reduced( self ):
      self.parameters[ "method" ] = "SS_REDUCED"
      stepper = Stepper( **self.parameters )
      A_out = stepper( self.A_in )
      # ========================================================================
      assert_almost_equal( max( np.abs(A_out)**2 ), self.P_analytical, 3 )
   # ===========================================================================
   def test_ss_agrawal( self ):
      self.parameters[ "method" ] = "SS_AGRAWAL"
      stepper = Stepper( **self.parameters )
      A_out = stepper( self.A_in )
      # ========================================================================
      assert_almost_equal( max( np.abs(A_out)**2 ), self.P_analytical, 4 )
   # ===========================================================================
   def test_ss_sym_midpoint( self ):
      self.parameters[ "method" ] = "SS_SYM_MIDPOINT"
      stepper = Stepper( **self.parameters )
      A_out = stepper( self.A_in )
      # ========================================================================
      assert_almost_equal( max( np.abs(A_out)**2 ), self.P_analytical, 2 )
   # ===========================================================================
   def test_ss_sym_rk4( self ):
      self.parameters[ "method" ] = "SS_SYM_RK4"
      stepper = Stepper( **self.parameters )
      A_out = stepper( self.A_in )
      # ========================================================================
      assert_almost_equal( max( np.abs(A_out)**2 ), self.P_analytical, 3 )
   # ===========================================================================
   def test_rk4ip( self ):
      self.parameters[ "method" ] = "RK4IP"
      stepper = Stepper( **self.parameters )
      A_out = stepper( self.A_in )
      # ========================================================================
      assert_almost_equal( max( np.abs(A_out)**2 ), self.P_analytical, 5 )
# ==============================================================================
if __name__ == "__main__":
   unittest2.main()
# ==============================================================================
