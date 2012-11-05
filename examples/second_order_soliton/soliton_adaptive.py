
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

import os.path
from itertools import cycle
import matplotlib.pyplot as plt

from pyofss import *
# ==============================================================================
def order_2_soliton( xi = 0.0, tau = 0.0 ):
   numerator = 4.0 * ( np.cosh(3.0 * tau) + \
                       3.0 * np.exp(4.0 * 1j * xi) * np.cosh(tau) ) * \
               np.exp(0.5 * 1j * xi)
   # ===========================================================================
   denominator = np.cosh(4.0 * tau) + 4.0 * np.cosh(2.0 * tau) + \
                 3.0 * np.cos(4.0 * xi)
   # ===========================================================================
   return numerator / denominator
# ==============================================================================
def generate_reference( domain ):

   offset = 0.5 * domain.bit_width
   A_true = np.array( order_2_soliton(0.5 * np.pi, domain.t - offset) )
   
   return A_true
# ==============================================================================
def save_simulations( domain, data_directory, methods, target_errors ):

   for method in methods:
      print "%s" % method

      for target_error in target_errors:
         print "\t%.1e" % target_error

         method_error = "-".join( [method, "%.0e" % (target_error)] )
         filename = os.path.join( data_directory, method_error )
         # =====================================================================
         system = System( domain )
         system.add( Sech(peak_power = 4.0, width = 1.0) )
         system.add( Fibre("fibre", length = 0.5 * np.pi, 
                           beta = [0.0, 0.0, -1.0, 0.0], gamma = 1.0, 
                           method = method, local_error = target_error) )
         system.run()
         # =====================================================================
         A_calc = system.fields['fibre']
         storage = system["fibre"].stepper.storage

         np.savez( filename, field = A_calc, ffts = storage.fft_total )
# ==============================================================================
def calculate_relative_error( data_directory, method, target_error, A_true ):

   method_error = "-".join( [method, "%.0e" % (target_error)] )
   filename = ".".join( [method_error, "npz"] )
   filename = os.path.join( data_directory, filename )

   try:
      npz_data = np.load( filename )
      A_calc = npz_data["field"]

      delta_power = temporal_power(A_calc) - temporal_power(A_true)
      mean_relative_error = np.mean( np.abs(delta_power) )
      mean_relative_error /= np.amax( temporal_power(A_true) )
      #~mean_relative_error /= np.mean( temporal_power(A_true) )
      #~mean_relative_error /= np.sum( temporal_power(A_true) )

      number_of_ffts = npz_data["ffts"]

      return (number_of_ffts, mean_relative_error)
   except IOError as e:
      print "Error opening file: %s" % filename
# ==============================================================================
def save_relative_errors( data_directory, methods, target_errors, A_true ):

   for method in methods:
      print "%s" % method

      filename = ".".join( (str(method),"dat") )
      filename = os.path.join( data_directory, filename )
      with open( filename, "w" ) as f:
         f.write( "Number of FFTs\t\tMean relative error" )

      results = []
      for target_error in target_errors:
         print "\t%.1e" % target_error

         result = calculate_relative_error( data_directory, method,
                                            target_error, A_true )
         results.append( result )

      with open( filename, "a" ) as f:
         for result in results:
            if( result is not None ):
               f.write( "\n%8d\t%.12e" % (result[0], result[1]) )
# ==============================================================================
def generate_plot( data_directory, methods ):

   labels = { "ass_simple" : "Simple split-step",
              "ass_symmetric" : "Symmetric split-step",
              "ass_reduced" : "Reduced split-step",
              "ass_agrawal" : "Agrawal split-step",
              "ass_sym_midpoint" : "Symmetric split-step midpoint",
              "ass_sym_rk4" : "Symmetric split-step RK4",
              "ark4ip" : "Runge-Kutta in the interaction picture" }

   # Generate a range of line and mark styles to cycle through:
   lines = [ "-", ":", "--" ]
   line_cycler = cycle( lines )

   marks = [ "*", "o", "s", "D", "^" ]
   mark_cycler = cycle( marks )
   # ===========================================================================
   for method in methods:
      filename = ".".join( (str(method),"dat") )
      filename = os.path.join( data_directory, filename )
      with open( filename, "r" ) as f:
         lines = f.readlines()
         # Skip first line, split on tab character, and unpack into two lists:
         (ffts, errors) = zip(*[ line.split("\t") for line in lines[1:] ])
      # ========================================================================
      plt.xlabel( "Number of FFTs" )
      plt.ylabel( "Mean relative error" )

      plt.plot( ffts, errors, label = labels[ method ], \
                linestyle = next(line_cycler), marker = next(mark_cycler) )
   # ===========================================================================
   plt.legend( loc = "lower left", frameon = False, prop = {"size":12} )

   plt.xlim( 50.0, 1.0e6 )

   plt.xscale( "log" )
   plt.yscale( "log" )

   plt.savefig( "soliton_error_vs_ffts" )
# ==============================================================================
if __name__ == "__main__":

   methods = [ "ass_simple", "ass_symmetric", "ass_reduced", \
               "ass_agrawal", "ass_sym_midpoint", "ass_sym_rk4", "ark4ip" ]
   target_errors = np.logspace( -1, -13, 13 )
   # ===========================================================================
   data_directory = "data_soliton_adaptive"
   # ===========================================================================
   domain = Domain( bit_width = 200.0, samples_per_bit = 4096 )
   # ===========================================================================
   # Step one: Save each simulation to indiviual file.
   save_simulations( domain, data_directory, methods, target_errors )
   # ===========================================================================
   # Step two: Using data from simulations, save relative error calculations.
   A_true = generate_reference( domain )
   save_relative_errors( data_directory, methods, target_errors, A_true )
   # ===========================================================================
   # Step three: Plot mean relative error data.
   generate_plot( data_directory, methods )
# ==============================================================================
