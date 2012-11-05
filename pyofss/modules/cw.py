
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

from numpy import pi
from scipy import sqrt, exp
# ==============================================================================
class Cw():
   """
   :param string name: Name of this module
   :param double peak_power: Peak power of the CW source
   :param double offset_nu: Offset frequency
   :param double initial_phase: Initial phase
   :param Uint channel: Channel sets the field to be modified

   Generate a continuous wave (CW) source. Add this to appropriate field.
   """
   # ===========================================================================
   def __init__( self, name = "cw", peak_power = 0.0, offset_nu = 0.0,
                 initial_phase = 0.0, channel = 0 ):

      self.name = name
      self.peak_power = peak_power
      self.offset_nu = offset_nu
      self.initial_phase = initial_phase
      self.channel = channel
   # ===========================================================================
   def __call__( self, domain, field ):

      self.field = field
      # ========================================================================
      phase = self.initial_phase
      phase += 2.0 * pi * self.offset_nu * domain.t
      shape = sqrt( self.peak_power ) * exp( -(1j * phase) )
      # ========================================================================
      if( domain.channels > 1 ):
         self.field[ self.channel ] += shape
      else:
         self.field += shape
      # ========================================================================
      return self.field
# ==============================================================================
