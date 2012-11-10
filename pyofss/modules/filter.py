
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
from numpy import pi
from scipy import exp, power, log

from pyofss.field import fft, ifft, ifftshift
# ==============================================================================
# Define exceptions
class FilterError( Exception ):
   pass

class OutOfRangeError( FilterError ):
   pass

class NotIntegerError( FilterError ):
   pass
# ==============================================================================
class Filter( object ):
   """
   :param string name: Name of this module
   :param double width_nu: Spectral bandwidth (HWIeM). *Unit: THz*
   :param double offset_nu: Offset frequency relative to domain centre
                            frequency. *Unit: THz*
   :param Uint m: Order parameter. m > 1 describes a super-Gaussian filter
   :param Uint channel: Channel of the field array to modify if multi-channel.
   :param bool using_fwhm: Determines whether the width_nu parameter is a 
                           full-width at half-maximum measure (FWHM), or a 
                           half-width at 1/e-maximum measure (HWIeM).

   Generate a super-Gaussian filter. A HWIeM bandwidth is used internally; a
   FWHM bandwidth will be converted on initialisation.
   """
   # ===========================================================================
   def __init__( self, name = "filter", width_nu = 0.1, offset_nu = 0.0,
                 m = 1, channel = 0, using_fwhm = False ):

      if not (1e-6 < width_nu < 1e3):
         raise OutOfRangeError, \
            "width_nu is out of range. Must be in (1e-6, 1e3)"

      if not (-100.0 < offset_nu < 100.0):
         raise OutOfRangeError, \
            "offset_nu is out of range. Must be in (-100.0, 100.0)"

      if not (0 < m < 50):
         raise OutOfRangeError, \
            "m is out of range. Must be in (0, 50)"

      if not (0 <= channel < 2):
         raise OutOfRangeError, \
            "channel is out of range. Must be in [0, 2)"

      if int(m) <> m:
         raise NotIntegerError, "m must be an integer"

      if int(channel) <> channel:
         raise NotIntegerError, "channel must be an integer"
      # ========================================================================
      self.name = name
      self.width_nu = width_nu
      self.offset_nu = offset_nu
      self.m = m
      self.channel = channel
      self.fwhm_nu = None
      # ========================================================================
      # For a FWHM filter width, store then convert to a HWIeM filter width:
      if( using_fwhm ):
         self.fwhm_nu = width_nu # store fwhm filter width
         self.width_nu *= 0.5 / power( log(2.0), 1.0 / (2 * m) )
   # ===========================================================================
   def calculate_fwhm( self ):
      if( self.fwhm_nu is not None ):
         return self.fwhm_nu
      else:
         return self.width_nu * 2.0 * power( log(2.0), 1.0 / (2 * self.m) )
   # ===========================================================================
   def __call__( self, domain, field ):
      """
      :param object domain: A domain
      :param object field: Current field
      :return: Field after modification by Gaussian filter
      :rtype: Object
      """
      # ========================================================================
      # Convert field to spectral domain:
      self.field = fft( field )

      delta_nu = domain.nu - domain.centre_nu - self.offset_nu
      factor = power( delta_nu / self.width_nu, (2 * self.m) )
       # Frequency values are in order, inverse shift to put in fft order:
      self.shape = exp( -0.5 * ifftshift(factor) )

      if( domain.channels > 1 ):
         # Filter is applied only to one channel:
         self.field[ self.channel ] *= self.shape
      else:
         self.field *= self.shape

      # convert field back to temporal domain:
      return ifft( self.field )
   # ===========================================================================
   def transfer_function( self, nu, centre_nu ):
      """
      :param Dvector nu: Spectral domain array.
      :param double centre_nu: Centre frequency. 
      :return: Array of values.
      :rtype: Dvector
      
      Generate an array representing the filter power transfer function.
      """
      # ========================================================================
      if len(nu) < 8:
         raise OutOfRangeError, \
            "Require spectral array with at least 8 values"

      delta_nu = nu - centre_nu - self.offset_nu
      factor = power( delta_nu / self.width_nu, (2 * self.m) )
      self.shape = exp( -0.5 * factor )

      return np.abs(self.shape)**2
   # ===========================================================================
   def __str__( self ):
      """
      :return: Information string
      :rtype: string

      Output information on Filter.
      """
      # ========================================================================
      output_string = ['width_nu = {0:f} THz', 'fwhm_nu = {1:f} THz', \
         'offset_nu = {2:f} THz', 'm = {3:d}', 'channel = {4:d}']

      return "\n".join( output_string ).format( self.width_nu, \
         self.calculate_fwhm(), self.offset_nu, self.m, self.channel )
# ==============================================================================
if __name__ == "__main__":
   """ Plot the power transfer function of the filter """
   # ===========================================================================
   from pyofss import *
   # ===========================================================================
   domain = Domain( centre_nu = 193.0 )
   gauss_filter = Filter( offset_nu = 1.5 )
   filter_tf = gauss_filter.transfer_function( domain.nu, domain.centre_nu )

   # Expect the filter transfer function to be centred at 194.5 THz:
   single_plot( domain.nu, filter_tf )
# ==============================================================================
