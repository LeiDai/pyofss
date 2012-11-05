
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

from pyofss.modules.stepper import Stepper

import unittest2
# ==============================================================================
class Default_parameters( unittest2.TestCase ):
   def default_function( A, z ):
      return A
   # ===========================================================================
   def test_default( self ):
      """ Should use default values (except for f) """
      stepper = Stepper( f = self.default_function )
      # ========================================================================
      self.assertEqual( stepper.traces, 1 )
      self.assertEqual( stepper.local_error, 1e-6 )
      self.assertEqual( stepper.method, "RK4" )
      self.assertEqual( stepper.length, 1.0 )
      self.assertEqual( stepper.total_steps, 100 )
      self.assertEqual( stepper.adaptive, False )
# ==============================================================================
class Bad_parameters( unittest2.TestCase ):
   def test_too_low( self ):
      """ Should fail when parameters are too low """
      pass
   # ===========================================================================
   def test_too_high( self ):
      """ Should fail when parameters are too high """
      pass
   # ===========================================================================
   def test_wrong_type( self ):
      """ Should fail if wrong type """
      pass
# ==============================================================================
class Check_functions( unittest2.TestCase ):
   def test_standard_stepper( self ):
      """ Stepper should integrate function using fixed step-size """
      pass
   def test_adaptive_stepper( self ):
      """ Stepper should integrate function using adaptive step-size """
      pass
# ==============================================================================
if __name__ == "__main__":
   unittest2.main()
# ==============================================================================
