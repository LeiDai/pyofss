
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

from scipy.constants import constants
from scipy.misc import factorial
from scipy import log10, exp
import numpy as np

from pyofss.field import fft, ifft, fftshift
# ==============================================================================
def convert_dispersion_to_physical( D = 0.0, S = 0.0, Lambda = 1550.0 ):
   """ 
   :param double D: Dispersion. *Unit:* :math:`ps / (nm \cdot km)`
   :param double S: Dispersion slope. *Unit:* :math:`ps / (nm^2 \cdot km)`
   :param double Lambda: Centre wavelength. *Unit: nm*
   :return: Second and third order dispersion
   :rtype: double, double

   Require D and S. Return a tuple containing beta_2 and beta_3.
   """
   # ===========================================================================
   if( (D == 0.0) and (S == 0.0) ):
      return (0.0, 0.0)

   # Constant 1.0e-3 modifies units of c to be nm / ps, see Domain.
   factor = Lambda**2 / (2.0 * np.pi * 1.0e-3 * constants.c)
   square_factor = factor**2

   beta_2 = -factor * D
   beta_3 = square_factor * ( S + (2.0 * D / Lambda) )

   return (beta_2, beta_3)
# ==============================================================================
def convert_dispersion_to_engineering( beta_2 = 0.0, beta_3 = 0.0,
                                       Lambda = 1550.0 ):
   """
   :param double beta_2: Second-order dispersion. *Unit:* :math:`ps^2 / km`
   :param double beta_3: Third-order dispersion. *Unit:* :math:`ps^3 / km`
   :param double Lambda: Centre wavelength. *Unit: nm*
   :return: Dispersion, dispersion slope
   :rtype: double, double

   Require beta_2 and beta_3. Return a tuple containing D and S.
   """
   # ===========================================================================
   if( (beta_2 == 0.0) and (beta_3 == 0.0) ):
      return (0.0, 0.0)

   factor_1 = (2.0 * np.pi * 1.0e-3 * constants.c) / (Lambda**2)
   square_factor_1 = factor_1 * factor_1

   factor_2 = (2.0 * factor_1) / Lambda

   D = -factor_1 * beta_2
   S = square_factor_1 * beta_3 + factor_2 * beta_2

   return (D, S)
# ==============================================================================
# ==============================================================================
def convert_alpha_to_linear( alpha_dB = 0.0 ):
   """
   :param double alpha_dB: Logarithmic attenuation factor
   :return: Linear attenuation factor
   :rtype: double

   Converts a logarithmic attenuation factor to a linear attenuation factor
   """
   # ===========================================================================
   factor = 10.0 * log10( exp(1.0) )
   return alpha_dB / factor
# ==============================================================================
def convert_alpha_to_dB( alpha_linear = 0.0 ):
   """
   :param double alpha_linear: Linear attenuation factor
   :return: Logarithmic attenuation factor
   :rtype: double

   Converts a linear attenuation factor to a logarithmic attenuation factor
   """
   # ===========================================================================
   factor = 10.0 * log10( exp(1.0) )
   return alpha_linear * factor
# ==============================================================================
# ==============================================================================
class Linearity():
   """
   :param double alpha: Attenuation factor
   :param array_like beta: Array of dispersion parameters
   :param string sim_type: Type of simulation, "default" or "wdm"
   :param bool use_cache: Cache calculated values if using fixed step-size
   :param double centre_omega: Angular frequency to use for dispersion array

   Dispersion is used by fibre to generate a fairly general dispersion array.
   """
   # ===========================================================================
   def __init__( self, alpha = None, beta = None, sim_type = None, 
                 use_cache = False, centre_omega = None ):

      self.alpha = alpha
      self.beta = beta
      self.centre_omega = centre_omega
      # ========================================================================
      self.generate_linearity = getattr( self, "%s_linearity" % sim_type,
                                          self.default_linearity )
      self.lin = getattr( self, "%s_f" % sim_type, self.default_f )
      # ========================================================================
      if( use_cache ):
         self.exp_lin = getattr( self, "%s_exp_f_cached" % sim_type,
                                 self.default_exp_f_cached )
      else:
         self.exp_lin = getattr( self, "%s_exp_f" % sim_type, 
                                 self.default_exp_f )
      # ========================================================================
      # Allows storing of calculation involving an exponential. Provides a 
      # significant speed increase if using a fixed step-size.
      self.cached_factor = None
   # ===========================================================================
   def __call__( self, domain ):

      return self.generate_linearity( domain )
   # ===========================================================================
   def default_linearity( self, domain ):
      # Calculate dispersive terms:
      if( self.beta is None ):
         self.factor = 0.0
      else:
         if( self.centre_omega is None ):
            self.Domega = domain.omega - domain.centre_omega
         else:
            self.Domega = domain.omega - self.centre_omega
         # =====================================================================
         # Allow general dispersion:
         terms = 0.0
         for n, beta in enumerate( self.beta ):
            terms += beta * np.power( self.Domega, n ) / factorial(n)
         self.factor = 1j * fftshift( terms )
      # ========================================================================
      # Include attenuation term if available:
      if( self.alpha is None ):
         return self.factor
      else:
         self.factor -= 0.5 * self.alpha
         return self.factor
   # ===========================================================================
   def wdm_linearity( self, domain ):
      # Calculate dispersive terms:
      if( self.beta is None ):
         self.factor = (0.0, 0.0)
      else:
         if( self.centre_omega is None ):
            self.Domega = (domain.omega - domain.centre_omega,
                           domain.omega - domain.centre_omega)
         else:
            self.Domega = (domain.omega - self.centre_omega[0],
                           domain.omega - self.centre_omega[1])
         # =====================================================================
         terms = [0.0, 0.0]
         for n, beta in enumerate( self.beta[0] ):
            terms[0] += beta * np.power( self.Domega[0], n ) / factorial(n)
         for n, beta in enumerate( self.beta[1] ):
            terms[1] += beta * np.power( self.Domega[1], n ) / factorial(n)
         self.factor = (1j * fftshift(terms[0]), 1j * fftshift(terms[1]) )
      # ========================================================================
      # Include attenuation terms if available:
      if( self.alpha is None ):
         return self.factor
      else:
         self.factor[0] -= 0.5 * self.alpha[0]
         self.factor[1] -= 0.5 * self.alpha[1]
         return factor
   # ===========================================================================
   # ===========================================================================
   def default_f( self, A, z ):
      return ifft( self.factor * fft(A) )
   # ===========================================================================
   def default_exp_f( self, A, h ):
      return ifft( np.exp(h * self.factor) * fft(A) )
   # ===========================================================================
   def default_exp_f_cached( self, A, h ):
      if( self.cached_factor is None ):
         print "Caching linear factor"
         self.cached_factor = np.exp( h * self.factor )
      # ========================================================================
      return ifft( self.cached_factor * fft(A) )
   # ===========================================================================
   def wdm_f( self, As, z ):
      return np.asarray([ ifft( self.factor[0] * fft(As[0]) ), \
                          ifft( self.factor[1] * fft(As[1]) ) ])
   # ===========================================================================
   def wdm_exp_f( self, As, h ):
      return np.asarray([ ifft( np.exp(h * self.factor[0]) * fft(As[0]) ), \
                          ifft( np.exp(h * self.factor[1]) * fft(As[1]) ) ])
   # ===========================================================================
   def wdm_exp_f_cached( self, As, h ):
      if( self.cached_factor is None ):
         print "Caching linear factor"
         self.cached_factor = [ np.exp(h * self.factor[0]),
                                np.exp(h * self.factor[1]) ]
      # ========================================================================
      return np.asarray([ ifft( self.cached_factor[0] * fft(As[0]) ), 
                          ifft( self.cached_factor[1] * fft(As[1]) ) ])
# ==============================================================================
