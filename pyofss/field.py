
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
import scipy.fftpack
# ==============================================================================
# Although use of global variables is generally a bad idea, in this case it is 
# a simple solution to recording the number of ffts used:
fft_counter = 0
# ==============================================================================
def temporal_power( A_t, normalise = False ):
   """
   :param array_like A_t: Input field array in the temporal domain
   :param bool normalise: Normalise returned array to that of maximum value
   :return: Array of power values
   :rtype: array_like

   Generate an array of temporal power values from complex amplitudes array.
   """
   # ===========================================================================
   P = np.abs( A_t )**2

   if( normalise ):
      P /= max(P)

   return P
# ==============================================================================
def spectral_power( A_t, normalise = False ):
   """
   :param array_like A_t: Input field array in the temporal domain
   :param bool normalise: Normalise returned array to that of maximum value
   :return: Array of power values
   :rtype: array_like

   Generate an array of spectral power values from complex amplitudes array.
   *Note: Expect input field to be in temporal domain. Never input A_nu!*
   """
   # ===========================================================================
   P = np.abs( fft(A_t) )**2

   if( normalise ):
      P /= max(P)

   return ifftshift( P )
# ==============================================================================
def phase( A_t, unwrap = True ):
   """
   :param array_like A_t: Input field array in the temporal domain
   :param bool unwrap: Whether to unwrap phase angles from fixed range
   :return: Array of phase angle values
   :rtype: array_like

   Generate an array of phase angles from complex amplitudes array.
   """
   # ===========================================================================
   if( unwrap ):
      return np.unwrap( np.angle(A_t) )
   else:
      return np.angle(A_t)
# ==============================================================================
def chirp( A_t, window_nu, unwrap = True ):
   """
   :param array_like A_t: Input field array in the temporal domain
   :param double window_nu: Spectral window of the simulation
   :param bool unwrap: Whether to unwrap phase angles from fixed range
   :return: Array of chirp values
   :rtype: array_like

   Generate an array of chirp values from complex amplitudes array.
   """
   # ===========================================================================
   if( unwrap ):
      return -np.gradient( phase(A_t, True) ) * window_nu
   else:
      return -np.gradient( phase(A_t, False) ) * window_nu
# ==============================================================================
def fft( A_t ):
   """
   :param array_like A_t: Input field array in the temporal domain
   :return: Output field array in the spectral domain
   :rtype: array_like

   Fourier transform field from temporal domain to spectral domain.
   *Note: Physics convention -- positive sign in exponential term.*
   """
   # ===========================================================================
   global fft_counter
   fft_counter += 1
   # ===========================================================================
   return scipy.fftpack.ifft( A_t )
# ==============================================================================
def ifft( A_nu ):
   """
   :param array_like A_nu: Input field array in the spectral domain
   :return: Output field array in the temporal domain
   :rtype: array_like

   Inverse Fourier transform field from spectral domain to temporal domain.
   *Note: Physics convention -- negative sign in exponential term.*
   """
   # ===========================================================================
   global fft_counter
   fft_counter += 1
   # ===========================================================================
   return scipy.fftpack.fft( A_nu )
# ==============================================================================
def ifftshift( A_nu ):
   """
   :param array_like A_nu: Input field array in the spectral domain
   :return: Shifted field array in the spectral domain
   :rtype: array_like

   Shift the field values from "FFT order" to "consecutive order".
   """
   # ===========================================================================
   return scipy.fftpack.fftshift( A_nu )
# ==============================================================================
def fftshift( A_nu ):
   """
   :param array_like A_nu: Input field array in the spectral domain
   :return: Shifted field array in the spectral domain
   :rtype: array_like

   Shift the field values from "consecutive order" to "FFT order".
   """
   # ===========================================================================
   return scipy.fftpack.ifftshift( A_nu )
# ==============================================================================
