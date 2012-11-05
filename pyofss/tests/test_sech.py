
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

from pyofss.modules.sech import Sech, OutOfRangeError, NotIntegerError

from numpy.testing.utils import assert_almost_equal # uses decimal places

import unittest2
# ==============================================================================
class Default_parameters( unittest2.TestCase ):
   def test_none( self ):
      """ Should use default value if no parameter given """
      sech = Sech()
      self.assertEqual( sech.name, "sech" )
      self.assertEqual( sech.position, 0.5 )
      self.assertEqual( sech.width, 10.0 )
      self.assertEqual( sech.peak_power, 1e-3 )
      self.assertEqual( sech.offset_nu, 0.0 )
      self.assertEqual( sech.C, 0.0 )
      self.assertEqual( sech.initial_phase, 0.0 )
      self.assertEqual( sech.channel, 0 )
      self.assertIsNone( sech.fwhm )
# ==============================================================================
class Bad_parameters( unittest2.TestCase ):
   def test_too_low( self ):
      """ Should fail when parameters are too low """
      self.assertRaises( OutOfRangeError, Sech, position = -0.1 )
      self.assertRaises( OutOfRangeError, Sech, width = 1e-3 )
      self.assertRaises( OutOfRangeError, Sech, peak_power = -1e-9 )
      self.assertRaises( OutOfRangeError, Sech, offset_nu = -100.0 )
      self.assertRaises( OutOfRangeError, Sech, C = -1e3 )
      self.assertRaises( OutOfRangeError, Sech, initial_phase = -1.0 )
      self.assertRaises( OutOfRangeError, Sech, channel = -1 )
   # ===========================================================================
   def test_too_high( self ):
      """ Should fail when parameters are too high """
      from numpy import pi
      # ========================================================================
      self.assertRaises( OutOfRangeError, Sech, position = 1.1 )
      self.assertRaises( OutOfRangeError, Sech, width = 1e3 )
      self.assertRaises( OutOfRangeError, Sech, peak_power = 1e9 )
      self.assertRaises( OutOfRangeError, Sech, offset_nu = 100.0 )
      self.assertRaises( OutOfRangeError, Sech, C = 1e3 )
      self.assertRaises( OutOfRangeError, Sech, initial_phase = 2.0 * pi )
      self.assertRaises( OutOfRangeError, Sech, channel = 2 )
   # ===========================================================================
   def test_wrong_type( self ):
      """ Should fail if wrong type """
      self.assertRaises( NotIntegerError, Sech, channel = 0.5 )
# ==============================================================================
class Check_input_and_output( unittest2.TestCase ):
   def test_output( self ):
      """ Check Sech outputs its values correctly """
      sech = Sech( "sech", 0.2, 5.0, 1.4, 0.3, 0, 0.4, 0.0, 1 )
      expected_string = ['position = 0.200000', 'width = 5.000000 ps', \
         'fwhm = 8.813736 ps', 'peak_power = 1.400000 W', \
         'offset_nu = 0.300000 THz', 'C = 0.400000', \
         'initial_phase = 0.000000 rad', 'channel = 1']
      self.assertEqual( str(sech), '\n'.join(expected_string) )
# ==============================================================================
class Check_functions( unittest2.TestCase ):
   def test_conversion( self ):
      """ Should calculate a FWHM pulse width from the HWIeM value. """
      sech = Sech( width = 100.0 )
      fwhm = sech.calculate_fwhm()
      self.assertEqual( fwhm, 176.2747174039086 )
   # ===========================================================================
   def test_bad_t( self ):
      """ Should raise an exception when temporal array has too few values """
      from numpy import arange
      t = arange(0.0, 4.0)
      self.assertRaises( OutOfRangeError, Sech().generate, t )
   # ===========================================================================
   def test_call( self ):
      """ Check generated Sech function """
      from numpy import arange
      sech = Sech( "sech", 0.4, 5.0, 1.5, 0.3, 0, 1.5 )
      t = arange(0.0, 100.0)
      A = sech.generate( t )
      P = abs(A)**2
      assert_almost_equal( max( P ), 1.5 )
# ==============================================================================
if __name__ == "__main__":
   unittest2.main()
# ==============================================================================
