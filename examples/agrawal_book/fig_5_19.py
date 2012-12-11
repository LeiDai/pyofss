
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

import sys
from pyofss import *
# ==============================================================================
domain = Domain( bit_width = 0.4, samples_per_bit = 4096 )
s = 0.2
width = 1.0 / (s * domain.centre_omega)
beta_2 = width**2
# ==============================================================================
system = System( domain )
system.add( Sech(peak_power = 1.0, width = width) )
system.add( Fibre(length = 4.0, beta = [0.0, 0.0, -beta_2],
                  gamma = 4.0, self_steepening = True, 
                  traces = 100, method = 'ARK4IP') )
system.run()
# ==============================================================================
storage = system['fibre'].stepper.storage
(x, y, z) = storage.get_plot_data( False, (0.0, 600.0), normalised = True )
# ==============================================================================
map_plot( x, y, z, labels["nu"], labels["P_nu"], labels["z"],
          filename = "5-19_map_nu")
# ==============================================================================
waterfall_plot( x, y, z, labels["nu"], labels["z"], labels["P_nu"],
                filename = "5-19_waterfall_nu", y_range = (0.0, 1.1) )
# ==============================================================================
if( len(sys.argv) > 1 and sys.argv[1] == 'animate' ):
   animated_plot( x, y, z, labels["nu"], labels["P_nu"],
                  "$z = {0:7.3f} \, km$", (x[0], x[-1]), (0.0, 1.1), fps = 20,
                  frame_prefix = "nu_", filename = "5-19_animation_nu.avi" )
# ==============================================================================
(x, y, z) = storage.get_plot_data( reduced_range = (0.18, 0.23) )
# ==============================================================================
map_plot( x, y, z, labels["t"], labels["P_t"], labels["z"],
          filename = "5-19_map_t" )
# ==============================================================================
waterfall_plot( x, y, z, labels["t"], labels["z"], labels["P_t"],
                filename = "5-19_waterfall_t", y_range = (0.0, 3.6) )
# ==============================================================================
if( len(sys.argv) > 1 and sys.argv[1] == 'animate' ):
   animated_plot( x, y, z, labels["t"], labels["P_t"], "$z = {0:7.3f} \, km$", 
                  (x[0], x[-1]), (0.0, 3.6), fps = 20, frame_prefix = "t_",
                  filename = "5-19_animation_t.avi" )
# ==============================================================================
