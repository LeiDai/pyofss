
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

"""
Note that although the temporal waterfall plot does not match that in Agrawal's 
"Nonlinear Fiber Optics" book (which uses the "NLSE" simulation software), it 
does match results using "LaserFOAM" simulation program.
"""

import sys
from pyofss import *
# ==============================================================================
system = System( Domain(bit_width = 100.0, samples_per_bit = 2048) )
system.add( Sech(peak_power = 1.0, width = 1.0, C = 0.5) )
system.add( Fibre(length = 20.0, beta = [0.0, 0.0, -1.0, 0.0],
                  gamma = 1.0, traces = 100, method = 'ARK4IP') )
system.run()
# ==============================================================================
storage = system['fibre'].stepper.storage
(x, y, z) = storage.get_plot_data( reduced_range = (40.0, 60.0) )
# ==============================================================================
map_plot( x, y, z, labels["t"], labels["P_t"], labels["z"],
          filename = "5-9_map_t" )
# ==============================================================================
waterfall_plot( x, y, z, labels["t"], labels["z"], labels["P_t"],
                filename = "5-9_waterfall_t", y_range = (0.0, 1.64) )
# ==============================================================================
if( len(sys.argv) > 1 and sys.argv[1] == 'animate' ):
   animated_plot( x, y, z, labels["t"], labels["P_t"], "$z = {0:7.3f} \, km$", 
                  (x[0], x[-1]), (0.0, 1.64), fps = 20, frame_prefix = "t_",
                  filename = "5-9_animation_t.avi" )
# ==============================================================================
