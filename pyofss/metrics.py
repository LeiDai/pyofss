
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

from domain import Domain
from field import temporal_power
# ==============================================================================
def generate_eye_data( domain, field, scale = None ):
   """
   :param object domain: A domain
   :param object field: Current field
   :param double scale: Scale factor by which to divide power
   :return: Temporal array and power matrix
   :rtype: Tuple

   Generates a temporal array and matrix containing eye traces. Each eye trace 
   consists of the power array calculated for successive bits. 
   """
   # ===========================================================================
   # Use the first bit range for the temporal axis:
   t_eye = domain.t[0 : domain.samples_per_bit]

   # Calculate power array and normalise if requested:
   P_t = temporal_power( field )
   if( scale is not None ):
      P_t /= scale

   # Reshape into matrix with total_bits rows and samples_per_bit columns:
   P_eye = P_t.reshape( domain.total_bits, domain.samples_per_bit )

   # If passed a matrix, pyplot will plot each column as a separate data series.
   # Require each row to be plotted, so transpose the matrix:
   return t_eye, P_eye.T
# ==============================================================================
def calculate_regenerator_factor( alpha, D, gamma, length, peak_power,
                                  using_alpha_dB = False ):

   """
   :param double alpha: Attenuation factor
   :param double D: Dispersion
   :param double length: Fibre length
   :param double peak_power: Pulse peak power
   :param bool using_alpha_dB: Whether using a logarithmic attnuation factor
   :return: Maximum nonlinear phase and a factor indicating regeneration
   :rtype: double, double
   """
   # ===========================================================================
   if( using_alpha_dB ):
      factor = 10.0 * np.log10( np.exp(1.0) )
      alpha /= factor
   # ===========================================================================
   effective_length = ( 1.0 - np.exp(-alpha * length) ) / alpha
   # ===========================================================================
   max_nonlinear_phase = gamma * peak_power * effective_length
   regenerator_factor = abs( D * length ) / max_nonlinear_phase
   # ===========================================================================
   return (max_nonlinear_phase, regenerator_factor)
# ==============================================================================
class Metrics( object ):
   """
   Calculate useful metrics using a domain and a field. An example metric is Q.
   """
   # ===========================================================================
   def __init__( self, domain = None, field = None ):

      if( (domain is None) or (field is None) ):
         raise Exception( "Metrics require a domain AND a field" )
      # ========================================================================
      self.domain = domain
      self.field = field
      # ========================================================================
      self.max_Q_dB = None
      self.sample_time = None
      self.sample_threshold = None
      # ========================================================================
      self.ones = None
      self.zeros = None
      # ========================================================================
      self.mean_peak_power_ones = None
      self.mean_peak_power_zeros = None
      # ========================================================================
      self.amplitude_jitter = None
      self.extinction_ratio = None
   # ===========================================================================
   def __str__( self ):
      """
      :return: Information string
      :rtype: string

      Output information on metrics, including Q, extinction ratio, and jitter.
      """
      # ========================================================================
      output_string = ['max_Q = {0:f} dB', 'sample_time = {1:f} ps', \
         'threshold = {2:f} W', '<P_0,ones> = {3:f} W', \
         '<P_0,zeros> = {4:f} W', 'extinction_ratio = {5:f} dB', \
         'power_jitter = +/- {6:f} %']

      return "\n".join( output_string ).format( self.max_Q_dB, \
         self.sample_time, self.sample_threshold, self.mean_peak_power_ones,
         self.mean_peak_power_zeros, self.extinction_ratio,
         self.amplitude_jitter * 100.0 )
   # ===========================================================================
   def calculate( self ):

      maximum_absolute_difference = -1.0
      sample_threshold = -1
      maximum_Q = -1.0
      sample_time = -1.0
      # ========================================================================
      TB = self.domain.total_bits
      SPB = self.domain.samples_per_bit
      # ========================================================================
      for spb in range(SPB):
         data = [ temporal_power( self.field[tb * SPB + spb] ) \
                  for tb in range(TB) ]
         # =====================================================================
         threshold = np.sum( data ) / TB
         # =====================================================================
         zeros = []
         ones = []
         for datum in data:
            if( datum < threshold ):
               zeros.append( datum )
            else:
               ones.append( datum )
         # =====================================================================
         if( (len(zeros) < 4) or (len(ones) < 4) ):
            print "Not enough ones and zeros to calculate Q"
#            raise Exception( "Not enough ones and zeros to calculate Q" )
         # =====================================================================
         absolute_difference = np.abs( np.mean(ones) - np.mean(zeros) )
         if( absolute_difference > maximum_absolute_difference ):
            maximum_absolute_difference = absolute_difference
            maximum_Q = absolute_difference / \
               ( np.std(ones) + np.std(zeros) )
            sample_time = spb * self.domain.dt
            sample_threshold = threshold
            # ==================================================================
            # Store ones and zeros arrays:
            self.ones = ones
            self.zeros = zeros
      # ========================================================================
      if( maximum_Q < 0.0 ):
         raise Exception( "Unable to calculate maximum Q!" )
      # ========================================================================
      self.max_Q_dB = 10.0 * np.log10( maximum_Q )
      self.sample_time = sample_time
      self.sample_threshold = sample_threshold
      # ========================================================================
      self.mean_peak_power_ones = np.mean( self.ones )
      self.mean_peak_power_zeros = np.mean( self.zeros )
      # ========================================================================
      delta_ones = np.max( self.ones ) - np.min( self.ones )
      AJ = 0.5 * delta_ones / np.mean(self.ones)
      self.amplitude_jitter = AJ

      ER = np.max(self.zeros) / ( np.mean(self.ones) * (1.0 - AJ) )
      self.extinction_ratio = -10.0 * np.log10( ER )

      # Alternative definition for extinction ratio:
      #~self.extinction_ratio = -10.0 * np.log10( np.max( self.zeros ) /
                                                #~np.min( self.ones ) )
# ==============================================================================
