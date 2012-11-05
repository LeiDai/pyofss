
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

from pyofss.domain import Domain, OutOfRangeError, NotIntegerError
from pyofss.domain import nu_to_omega, nu_to_lambda
from pyofss.domain import omega_to_nu, omega_to_lambda
from pyofss.domain import lambda_to_nu, lambda_to_omega
from pyofss.domain import dnu_to_dlambda, dlambda_to_dnu

import unittest2
#from numpy.testing.utils import assert_approx_equal # uses significant figures
from numpy.testing.utils import assert_almost_equal # uses decimal places
#from numpy.testing.utils import assert_array_equal
from numpy.testing.utils import assert_array_almost_equal
# ==============================================================================
class Default_parameters( unittest2.TestCase ):
   def test_none( self ):
      """ Should use default value if no parameter given """
      domain = Domain()
      self.assertEqual( domain.total_bits, 1 )
      self.assertEqual( domain.samples_per_bit, 512 )
      self.assertEqual( domain.bit_width, 100.0 )
      self.assertEqual( domain.centre_nu, 193.1 )
      self.assertEqual( domain.channels, 1 )
# ==============================================================================
class Bad_parameters( unittest2.TestCase ):
   def test_too_low( self ):
      """ Should fail when parameters are too low """
      self.assertRaises( OutOfRangeError, Domain, total_bits = 0 )
      self.assertRaises( OutOfRangeError, Domain, samples_per_bit = 0 )
      self.assertRaises( OutOfRangeError, Domain, bit_width = 0.01 )
      self.assertRaises( OutOfRangeError, Domain, centre_nu = 185.0 )
      self.assertRaises( OutOfRangeError, Domain, channels = 0 )
   # ===========================================================================
   def test_too_high( self ):
      """ Should fail when parameters are too high """
      self.assertRaises( OutOfRangeError, Domain, total_bits = 4096 )
      self.assertRaises( OutOfRangeError, Domain, samples_per_bit = 262144 )
      self.assertRaises( OutOfRangeError, Domain, bit_width = 10000.0 )
      self.assertRaises( OutOfRangeError, Domain, centre_nu = 400.0 )
      self.assertRaises( OutOfRangeError, Domain, channels = 3 )
   # ===========================================================================
   def test_wrong_type( self ):
      """ Should fail if wrong type """
      self.assertRaises( NotIntegerError, Domain, total_bits = 32.5 )
      self.assertRaises( NotIntegerError, Domain, samples_per_bit = 8.3 )
      self.assertRaises( NotIntegerError, Domain, channels = 1.5 )
# ==============================================================================
class Check_conversions( unittest2.TestCase ):
   def test_derived( self ):
      """ Derived parameters should have expected values """
      domain = Domain( 2, 32, 200.0, 193.2 )
      self.assertEqual( domain.total_samples, 64 )
      self.assertEqual( domain.window_t, 400.0 )
      assert_almost_equal( domain.centre_omega, 1213.911401347, 9 )
      assert_almost_equal( domain.centre_lambda, 1551.720797, 6 )
   # ===========================================================================
   def test_deltas( self ):
      """ Delta parameters should have expected values """
      domain = Domain( 1, 64, 400.0, 193.3 )
      self.assertEqual( domain.dt, 400.0 / (1 * 64) )
      assert_almost_equal( domain.dnu, 1.0 / domain.window_t )
      assert_almost_equal( domain.dlambda, Domain.vacuum_light_speed * \
                           domain.dnu / (domain.centre_nu)**2 )
   # ===========================================================================
   def test_arrays( self ):
      """ Generated arrays should have expected values """
      domain = Domain( 4, 16, 50.0, 190.0 )
      self.assertEqual( domain.t[-1], domain.window_t - domain.dt )
      # Note division operator for integers is now //, / is for floats:
      self.assertEqual( domain.nu[domain.total_samples // 2], domain.centre_nu )
# ==============================================================================
class Check_utility_functions( unittest2.TestCase ):
   def test_from_nu( self ):
      """ Check conversions from nu to omega and nu to lambda """
      domain = Domain()
      self.assertEqual( nu_to_omega(domain.centre_nu), domain.centre_omega )
      self.assertEqual( nu_to_lambda(domain.centre_nu), domain.centre_lambda )
      assert_array_almost_equal( nu_to_omega(domain.nu), domain.omega )
      assert_array_almost_equal( nu_to_lambda(domain.nu), domain.Lambda )
   # ===========================================================================
   def test_from_omega( self ):
      """ Check conversions from omega to nu and omega to lambda """
      domain = Domain()
      self.assertEqual( omega_to_nu(domain.centre_omega), domain.centre_nu )
      self.assertEqual( omega_to_lambda(domain.centre_omega), \
                        domain.centre_lambda )
      assert_array_almost_equal( omega_to_nu(domain.omega), domain.nu )
      assert_array_almost_equal( omega_to_lambda(domain.omega), domain.Lambda )
   # ===========================================================================
   def test_from_lambda( self ):
      """ Check conversions from lambda to nu and lambda to omega """
      domain = Domain()
      self.assertEqual( lambda_to_nu(domain.centre_lambda), domain.centre_nu )
      self.assertEqual( lambda_to_omega(domain.centre_lambda), \
                        domain.centre_omega)
      assert_array_almost_equal( lambda_to_nu(domain.Lambda), domain.nu )
      assert_array_almost_equal( lambda_to_omega(domain.Lambda), domain.omega )
   # ===========================================================================
   def test_dnu_dlambda_conversion( self ):
      """ Check conversions from dnu to dlambda and dlambda to dnu """
      domain = Domain()
      self.assertEqual( dnu_to_dlambda(domain.dnu, domain.centre_nu), \
                        domain.dlambda )
      self.assertEqual( dlambda_to_dnu(domain.dlambda, domain.centre_lambda), \
                        domain.dnu )
# ==============================================================================
class Check_input_and_output( unittest2.TestCase ):
   def test_output( self ):
      """ Check Domain outputs its values correctly """
      domain = Domain( 2, 128, 25.0, 193.8, 2 )
      expected_string = ['Domain:', 'total_bits = 2', 'samples_per_bit = 128', \
         'bit_width = 25.0000 ps', 'centre_nu = 193.8000 THz', \
         'total_samples = 256', 'window_t = 50.0000 ps', \
         'centre_omega = 1217.6813 rad / ps', 'centre_lambda = 1546.9167 nm', \
         'dt = 0.1953 ps', 'window_nu = 5.1200 THz', 'dnu = 0.0200 THz', \
         'window_omega = 32.1699 rad / ps', 'domega = 0.1257 rad / ps', \
         'window_lambda = 40.8680 nm', 'dlambda = 0.1596 nm', 'channels = 2']
      self.assertEqual( str(domain), '\n\t'.join(expected_string) )
# ==============================================================================
if __name__ == "__main__":
   unittest2.main()
# ==============================================================================
