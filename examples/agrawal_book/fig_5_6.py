
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

import sys
from pyofss import *
# ==============================================================================
system = System( Domain(bit_width = 100.0, samples_per_bit = 4096) )
system.add( Sech(peak_power = 1.0, width = 1.0) )
system.add( Fibre(length = 0.5 * np.pi, beta = [0.0, 0.0, -1.0, 0.0],
                  gamma = 9.0, traces = 100, method = 'ARK4IP') )
system.run()
# ==============================================================================
storage = system['fibre'].stepper.storage
(x, y, z) = storage.get_plot_data( False, (191.1, 195.1), normalised = True )
# ==============================================================================
map_plot( x, y, z, labels["nu"], labels["P_nu"], labels["z"],
          filename = "5-6_map_nu" )
# ==============================================================================
waterfall_plot( x, y, z, labels["nu"], labels["z"], labels["P_nu"],
                filename = "5-6_waterfall_nu", y_range = (0.0, 1.0) )
# ==============================================================================
if( len(sys.argv) > 1 and sys.argv[1] == 'animate' ):
   animated_plot( x, y, z, labels["nu"], labels["P_nu"],
                  "$z = {0:7.3f} \, km$", (x[0], x[-1]), (0.0, 1.1), fps = 20,
                  frame_prefix = "nu_", filename = "5-6_animation_nu.avi" )
# ==============================================================================
(x, y, z) = storage.get_plot_data( reduced_range = (46.0, 54.0) )
# ==============================================================================
map_plot( x, y, z, labels["t"], labels["P_t"], labels["z"],
          filename = "5-6_map_t" )
# ==============================================================================
waterfall_plot( x, y, z, labels["t"], labels["z"], labels["P_t"],
                filename = "5-6_waterfall_t", y_range = (0.0, 6.0) )
# ==============================================================================
if( len(sys.argv) > 1 and sys.argv[1] == 'animate' ):
   animated_plot( x, y, z, labels["t"], labels["P_t"], "$z = {0:7.3f} \, km$", 
                  (x[0], x[-1]), (0.0, 6.1), fps = 20, frame_prefix = "t_",
                  filename = "5-6_animation_t.avi" )
# ==============================================================================
