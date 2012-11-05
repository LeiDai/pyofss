
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

from scipy.constants import constants
# ==============================================================================
# Define exceptions
class DomainError( Exception ):
   pass

class OutOfRangeError( DomainError ):
   pass

class NotIntegerError( DomainError ):
   pass
# ==============================================================================
def nu_to_omega( nu ):
   """
   :param double nu: Frequency to convert. *Unit: THz*
   :return: Angular frequency. *Unit: rad / ps*
   :rtype: double

   Convert from frequency to angular frequency
   """
   return 2.0 * pi * nu

def nu_to_lambda( nu ):
   """
   :param double nu: Frequency to convert. *Unit: THz*
   :return: Wavelength. *Unit: nm*
   :rtype: double

   Convert from frequency to wavelength
   """
   return Domain.vacuum_light_speed / nu

def omega_to_nu( omega ):
   """
   :param double omega: Angular frequency to convert. *Unit: rad / ps*
   :return: Frequency. *Unit: THz*
   :rtype: double

   Convert from angular frequency to frequency
   """
   return omega / (2.0 * pi)

def omega_to_lambda( omega ):
   """
   :param double omega: Angular frequency to convert. *Unit: rad / ps*
   :return: Wavelength. *Unit: nm*
   :rtype: double

   Convert from angular frequency to wavelength
   """
   return 2.0 * pi * Domain.vacuum_light_speed / omega

def lambda_to_nu( Lambda ):
   """
   :param double Lambda: Wavelength to convert. *Unit: nm*
   :return: Frequency. *Unit: THz*
   :rtype: double

   Convert from wavelength to frequency
   """
   return Domain.vacuum_light_speed / Lambda

def lambda_to_omega( Lambda ):
   """
   :param double Lambda: Wavelength to convert. *Unit: nm*
   :return: Angular frequency. *Unit: rad / ps*
   :rtype: double

   Convert from wavelength to angular frequency
   """
   return 2.0 * pi * Domain.vacuum_light_speed / Lambda

def dnu_to_dlambda( dnu, nu = 193.1 ):
   """
   :param double dnu: Small increment in frequency to convert. *Unit: THz*
   :param double nu: Reference frequency. *Unit: THz*
   :return: Small increment in wavelength. *Unit: nm*
   :rtype: double

   Convert a small change in frequency to a small change in wavelength,
   with reference to centre frequency
   """
   return Domain.vacuum_light_speed * dnu / (nu)**2

def dlambda_to_dnu( dlambda, Lambda = 1550.0 ):
   """
   :param double dlambda: Small increment in wavelength to convert. *Unit: nm*
   :param double Lambda: Reference wavelength. *Unit: nm*
   :return: Small increment in frequency. *Unit: THz*
   :rtype: double

   Convert a small change in wavelength to a small change in frequency,
   with reference to a centre wavelength
   """
   return Domain.vacuum_light_speed * dlambda / (Lambda)**2
# ==============================================================================
class Domain():
   """
   :param Uint total_bits: Total number of bits to generate
   :param Uint samples_per_bit: Number of samples to represent a bit
   :param double bit_width: Width of each bit. *Unit: ps*
   :param double centre_nu: Centre frequency. *Unit: THz*
   :param Uint channels: Number of channels to simulate.
                         Used for WDM simulations

   A domain consists of:
      **Bit data**:
         total_bits, bit_width
      **Samples data**:
         samples_per_bit, total_samples
      **Window size in each domain**:
         window_t, window_nu, window_omega, window_lambda
      **Increment of each domain between two adjacent samples**:
         dt, dnu, domega, dlambda
      **The generated domains**:
         t, nu, omega, Lambda
      **Centre of spectral domains**:
         centre_nu, centre_omega, centre_lambda

   .. note::
      Use of *Lambda*, and NOT the Python reserved word *lambda*
   """
   # ===========================================================================
   vacuum_light_speed = 1.0e-3 * constants.c # nm / ps
   """ Speed of light in a vacuum. *Unit: nm / ps* """
   # ===========================================================================
   def __init__( self, total_bits = 1, samples_per_bit = 512, \
                 bit_width = 100.0, centre_nu = 193.1, channels = 1 ):

      if not (0 < total_bits < 4096):
         raise OutOfRangeError, \
               "total_bits is out of range. Must be in (0, 4096)"

      if not (0 < samples_per_bit < 262144):
         raise OutOfRangeError, \
               "samples_per_bit is out of range. Must be in (0, 262144)"

      if not (0.01 < bit_width < 10000.0):
         raise OutOfRangeError, \
               "bit_width is out of range. Must be in (0.01, 10000.0)"

      if not (185.0 < centre_nu < 400.0):
         raise OutOfRangeError, \
               "centre_nu is out of range. Must be in (185.0, 400.0)"

      if not (0 < channels < 3):
         raise OutOfRangeError, \
               "channels is out of range. Must be in (0, 3)"

      if int( total_bits ) <> total_bits:
         raise NotIntegerError, "total_bits must be an integer"

      if int( samples_per_bit ) <> samples_per_bit:
         raise NotIntegerError, "samples_per_bit must be an integer"

      if int( channels ) <> channels:
         raise NotIntegerError, "channels must be an integer"
      # ========================================================================
      self.total_bits = total_bits
      self.samples_per_bit = samples_per_bit
      self.bit_width = bit_width # ps
      self.centre_nu = centre_nu # THz == ps^{-1}
      # ========================================================================
      self.total_samples = self.total_bits * self.samples_per_bit
      self.window_t = self.total_bits * self.bit_width
      # ========================================================================
      # Note the units used:
      # c = 3e8 m / s  = 3e5 nm / ps
      # centre_nu = 193.1e12 Hz = 193.1 THz (== 193.1 ps^{-1})
      # centre_lambda (nm) = c (nm / ps) / centre_nu (THz)

      self.centre_omega = 2.0 * pi * self.centre_nu # rad / ps
      self.centre_lambda = Domain.vacuum_light_speed / self.centre_nu # nm
      # ========================================================================
      # The last two parameters to linspace will set 'endpoint' to False and 
      # 'retstep' to True. The effect is to not use the final value, 
      # self.window_t, and to return a tuple (samples, step) instead of just 
      # samples; where step is the spacing between samples.

      (self.t, self.dt) = np.linspace( 0.0, self.window_t, \
                                       self.total_samples, False, True )
      # ========================================================================
      # Require nu = [_nu_min, _nu_min + window_nu)
      #            = [centre_nu - 0.5 * window_nu, centre_nu + 0.5 * window_nu)

      self.window_nu = 1.0 / self.dt

      # First value of nu (minimum value):
      _nu_min = self.centre_nu - 0.5 * self.window_nu 

      (self.nu, self.dnu) = np.linspace( _nu_min, _nu_min + self.window_nu, \
                                         self.total_samples, False, True )

      # Frequency values are in order. This convention is strict within pyofss.
      # If a calculation uses the frequency array then the values must be 
      # first transformed to fft order using:
      # fftshift( self.nu )
      # ========================================================================
      self.window_omega = 2.0 * pi * self.window_nu

      (self.omega, self.domega) = (2.0 * pi * self.nu, 2.0 * pi * self.dnu)
      # ========================================================================
      self.window_lambda = Domain.vacuum_light_speed * self.window_nu / \
                           (self.centre_nu)**2

      (self.Lambda, self.dlambda) = (Domain.vacuum_light_speed / self.nu, \
         Domain.vacuum_light_speed * self.dnu / (self.centre_nu)**2 )
      # ========================================================================
      self.channels = channels
   # ===========================================================================
   def __str__( self ):
      """
      :return: Information string
      :rtype: string

      Output information on Domain.
      """
      # ========================================================================
      output_string = ['Domain:', 'total_bits = {0:d}', \
         'samples_per_bit = {1:d}', 'bit_width = {2:.4f} ps', \
         'centre_nu = {3:.4f} THz', 'total_samples = {4:d}', \
         'window_t = {5:.4f} ps', 'centre_omega = {6:.4f} rad / ps', \
         'centre_lambda = {7:.4f} nm', 'dt = {8:.4f} ps', \
         'window_nu = {9:.4f} THz', 'dnu = {10:.4f} THz', \
         'window_omega = {11:.4f} rad / ps', 'domega = {12:.4f} rad / ps', \
         'window_lambda = {13:.4f} nm', 'dlambda = {14:.4f} nm', \
         'channels = {15:d}']

      return "\n\t".join( output_string ).format(self.total_bits, \
         self.samples_per_bit, self.bit_width, self.centre_nu, \
         self.total_samples, self.window_t, self.centre_omega, \
         self.centre_lambda, self.dt, self.window_nu, self.dnu, \
         self.window_omega, self.domega, self.window_lambda, \
         self.dlambda, self.channels)
# ==============================================================================
