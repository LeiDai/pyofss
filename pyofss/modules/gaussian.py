
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
from scipy import sqrt, exp, log, power
# ==============================================================================
# Define exceptions
class GaussianError( Exception ):
   pass

class OutOfRangeError( GaussianError ):
   pass

class NotIntegerError( GaussianError ):
   pass
# ==============================================================================
class Gaussian():
   """
   :param string name: Name of this module
   :param double position: Relative position within time window
   :param double width: Half width at 1/e of maximum. *Unit: ps*
   :param double peak_power: Power at maximum. *Unit: W*
   :param double offset_nu: Offset frequency relative to domain centre 
                            frequency. *Unit: THz*
   :param Uint m: Order parameter. m > 1 describes a super-Gaussian pulse
   :param double C: Chirp parameter. *Unit: rad*
   :param double initial_phase: Control the initial phase. *Unit: rad*
   :param Uint channel: Channel of the field array to modify if multi-channel.
   :param bool using_fwhm: Determines whether the width parameter is a 
                           full-width at half-maximum measure (FWHM), or a 
                           half-width at 1/e-maximum measure (HWIeM).

   Generates a pulse with a Gaussian profile. A HWIeM pulse width is used 
   internally; a FWHM pulse width will be converted on initialisation.
   """
   # ===========================================================================
   def __init__( self, name = "gaussian", position = 0.5, width = 10.0,
                 peak_power = 1e-3, offset_nu = 0.0, m = 1, C = 0.0,
                 initial_phase = 0.0, channel = 0, using_fwhm = False ):

      if not (0.0 <= position <= 1.0):
         raise OutOfRangeError, \
            "position is out of range. Must be in [0.0, 1.0]"

      if not (1e-3 < width < 1e3):
         raise OutOfRangeError, \
            "width is out of range. Must be in (1e-3, 1e3)"

      if not (0.0 <= peak_power < 1e9):
         raise OutOfRangeError, \
            "peak_power is out of range. Must be in [0.0, 1e9)"

      if not (-100.0 < offset_nu < 100.0):
         raise OutOfRangeError, \
               "offset_nu is out of range. Must be in (-100.0, 100.0)"

      if not (0 < m < 50):
         raise OutOfRangeError, \
               "m is out of range. Must be in (0, 50)"

      if not (-1e3 < C < 1e3):
         raise OutOfRangeError, \
               "C is out of range. Must be in (-1e3, 1e3)"

      if not (0.0 <= initial_phase < 2.0 * pi):
         raise OutOfRangeError, \
            "initial_phase is out of range. Must be in [0.0, 2.0 * pi)"

      if not (0 <= channel < 2):
         raise OutOfRangeError, \
            "channel is out of range. Must be in [0, 2)"

      if int(m) <> m:
         raise NotIntegerError, "m must be an integer"

      if int(channel) <> channel:
         raise NotIntegerError, "channel must be an integer"
      # ========================================================================
      self.name = name
      self.position = position
      self.width = width # ps
      self.peak_power = peak_power # W
      self.offset_nu = offset_nu # THz
      self.m = m
      self.C = C # rad
      self.initial_phase = initial_phase # rad
      self.channel = channel
      self.fwhm = None
      # ========================================================================
      # For a FWHM pulse width, store then convert to a HWIeM pulse width:
      if( using_fwhm ):
         self.fwhm = width # store fwhm pulse width
         self.width *= 0.5 / power( log(2.0), 1.0 / (2 * m) )
   # ===========================================================================
   def calculate_fwhm( self ):
      if( self.fwhm is not None ):
         return self.fwhm
      else:
         return self.width * 2.0 * power( log(2.0), 1.0 / (2 * self.m) )
   # ===========================================================================
   def __call__( self, domain, field ):
      """
      :param object domain: A domain
      :param object field: Current field
      :return: Field after modification by Gaussian
      :rtype: Object
      """
      # ========================================================================
      self.field = field

      t_normalised = (domain.t - self.position * domain.window_t) / self.width
      time = power( t_normalised, (2 * self.m) )

      phase = self.initial_phase
      phase -= 2.0 * pi * self.offset_nu * domain.t + 0.5 * self.C * time

      magnitude = sqrt( self.peak_power ) * exp( -0.5 * time )

      if( domain.channels > 1 ):
         self.field[ self.channel ] += magnitude * exp( 1j * phase )
      else:
         self.field += magnitude * exp( 1j * phase )

      return self.field
   # ===========================================================================
   def generate( self, t ):
      """
      :param Dvector t: Temporal domain array
      :return: Array of complex values. *Unit:* :math:`\sqrt{W}`
      :rtype: Cvector
      
      Generate an array of complex values representing a Gaussian pulse.
      """
      # ========================================================================
      if len(t) < 8:
         raise OutOfRangeError, \
            "Require temporal array with at least 8 values"

      # Assume t[0] = t_0 and t[-1] = t_0 + t_range - dt, with dt = t[1] - t[0]
      t_range = t[-1] - t[0] + (t[1] - t[0])
      t_normalised = (t - self.position * t_range) / self.width
      time = power( t_normalised, (2 * self.m) )

      phase = self.initial_phase
      phase -= 2.0 * pi * self.offset_nu * t + 0.5 * self.C * time

      return sqrt(self.peak_power) * exp(-0.5 * time + 1j * phase)
   # ===========================================================================
   def __str__( self ):
      """
      :return: Information string
      :rtype: string

      Output information on Gaussian.
      """
      # ========================================================================
      output_string = ['position = {0:f}', 'width = {1:f} ps', \
         'fwhm = {2:f} ps', 'peak_power = {3:f} W', \
         'offset_nu = {4:f} THz', 'm = {5:d}', 'C = {6:f}', \
         'initial_phase = {7:f} rad', 'channel = {8:d}']

      return "\n".join( output_string ).format( self.position, self.width, \
         self.calculate_fwhm(), self.peak_power, self.offset_nu, \
         self.m, self.C, self.initial_phase, self.channel )
# ==============================================================================
if __name__ == "__main__":
   """ Plot a default Gaussian in temporal and spectral domain """
   # ===========================================================================
   from pyofss import *
   # ===========================================================================
   sys = System( Domain(bit_width = 500.0) )
   sys.add( Gaussian() )
   sys.run()

   double_plot( sys.domain.t, temporal_power( sys.field ),
                sys.domain.nu, spectral_power( sys.field, True ),
                labels["t"], labels["P_t"], labels["nu"], labels["P_nu"] )
# ==============================================================================
