
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
for m in [1, 3]:
   system = System( Domain(bit_width = 200.0, samples_per_bit = 2048) )
   system.add( Gaussian(peak_power = 1.0, width = 1.0, m = m) )
   system.add( Fibre(length = 10.0, gamma = 1.0, traces = 50) )
   system.run()
   # ===========================================================================
   storage = system['fibre'].stepper.storage
   (x,y,z) = storage.get_plot_data( is_temporal = False, normalised = True )
   # ===========================================================================
   map_plot( x, y, z, labels["nu"], labels["P_nu"], labels["z"],
             filename = "4-4_map_m-{0:d}".format(m) )
   # ===========================================================================
   waterfall_plot( x, y, z, labels["nu"], labels["z"], labels["P_nu"],
                   filename = "4-4_waterfall_m-{0:d}".format(m) )
   # ===========================================================================
   if( len(sys.argv) > 1 and sys.argv[1] == 'animate' ):
      animated_plot( x, y, z, labels["nu"], labels["P_nu"],
                     "$z = {0:7.3f} \, km$", (x[0], x[-1]), (0.0, 1.1),
                     fps = 10, frame_prefix = "m{0:d}_".format(m), 
                     filename = "4-4_animation_m-{0:d}.avi".format(m) )
# ==============================================================================
