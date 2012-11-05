
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

from pyofss.modules.filter import Filter, OutOfRangeError, NotIntegerError

import unittest2
# ==============================================================================
class Default_parameters( unittest2.TestCase ):
   def test_none( self ):
      """ Should use default value if no parameter given """
      gfilter = Filter()
      self.assertEqual( gfilter.name, "filter" )
      self.assertEqual( gfilter.width_nu, 0.1 )
      self.assertEqual( gfilter.offset_nu, 0.0 )
      self.assertEqual( gfilter.m, 1 )
      self.assertEqual( gfilter.channel, 0 )
      self.assertIsNone( gfilter.fwhm_nu )
# ==============================================================================
class Bad_parameters( unittest2.TestCase ):
   def test_too_low( self ):
      """ Should fail when parameters are too low """
      self.assertRaises( OutOfRangeError, Filter, width_nu = 1e-6 )
      self.assertRaises( OutOfRangeError, Filter, offset_nu = -100.0 )
      self.assertRaises( OutOfRangeError, Filter, m = 0 )
      self.assertRaises( OutOfRangeError, Filter, channel = -1 )
   # ===========================================================================
   def test_too_high( self ):
      """ Should fail when parameters are too high """
      self.assertRaises( OutOfRangeError, Filter, width_nu = 1e3 )
      self.assertRaises( OutOfRangeError, Filter, offset_nu = 100.0 )
      self.assertRaises( OutOfRangeError, Filter, m = 50 )
      self.assertRaises( OutOfRangeError, Filter, channel = 2 )
# ==============================================================================
   def test_wrong_type( self ):
      """ Should fail if wrong type """
      self.assertRaises( NotIntegerError, Filter, m = 1.4 )
      self.assertRaises( NotIntegerError, Filter, channel = 0.5 )
# ==============================================================================
class Check_input_and_output( unittest2.TestCase ):
   def test_output( self ):
      """ Check filter outputs its values correctly """
      gfilter = Filter( "filter", 0.2, 0.5, 4, 1 )
      expected_string = ['width_nu = 0.200000 THz', 'fwhm_nu = 0.382088 THz', \
                         'offset_nu = 0.500000 THz', 'm = 4', 'channel = 1' ]
      self.assertEqual( str(gfilter), '\n'.join(expected_string) )
# ==============================================================================
class Check_functions( unittest2.TestCase ):
   def test_conversion( self ):
      """ Should calculate a FWHM bandwidth from the HWIeM value. """
      gfilter = Filter( width_nu = 0.2 )
      fwhm_nu = gfilter.calculate_fwhm()
      self.assertEqual( fwhm_nu, 0.33302184446307909 )
      # ========================================================================
      gfilter = Filter( width_nu = 0.2, m = 4 )
      fwhm_nu = gfilter.calculate_fwhm()
      self.assertEqual( fwhm_nu, 0.38208780263892828 )
# ==============================================================================
if __name__ == "__main__":
   unittest2.main()
# ==============================================================================
