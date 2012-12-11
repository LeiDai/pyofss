
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
nu_0 = lambda_to_nu( 1060.0 )
nu_1 = lambda_to_nu( 1550.0 )

offset_nu = nu_0 - nu_1
# ==============================================================================
system = System( Domain(bit_width = 20.0, samples_per_bit = 8192, 
                        channels = 2, centre_nu = nu_0) )
system.add( Gaussian(width = 1.0, peak_power = 1000.0, channel = 0) )
system.add( Gaussian(width = 1.0, peak_power = 0.1, channel = 1,
                     offset_nu = -offset_nu) )
system.add( Fibre('fibre', length = 0.05, gamma = [0.9, 0.615483871],
            beta = [ [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, -1.0, 0.0] ],
            centre_omega = ( nu_to_omega(nu_0), nu_to_omega(nu_1) ),
            sim_type = 'wdm', method = 'ARK4IP', traces = 100) )
system.run()
# ==============================================================================
storage = system['fibre'].stepper.storage
(x, y, z_temp) = storage.get_plot_data( channel = 0 )
z_label = "Fibre length, $z \, (m)$"
z = z_temp * 1.0e3
# ==============================================================================
map_plot( x, y, z, labels["t"], labels["P_t"], z_label,
          filename = "7-9_map_t_pump" )
# ==============================================================================
waterfall_plot( x, y, z, labels["t"], z_label, labels["P_t"],
                filename = "7-9_waterfall_t_pump", y_range = (0.0, 1.0e3) )
# ==============================================================================
if( len(sys.argv) > 1 and sys.argv[1] == 'animate' ):
   animated_plot( x, y, z, labels["t"], labels["P_t"], "$z = {0:7.3f} \, m$", 
                  (x[0], x[-1]), (0.0, 1.0e3), fps = 20, frame_prefix = "pump_", filename = "7-9_animation_t_pump.avi" )
# ==============================================================================
(x, y, z_temp) = storage.get_plot_data( True, (6.0, 14.0), False, channel = 1 )
z = z_temp * 1.0e3
# ==============================================================================
map_plot( x, y, z, labels["t"], labels["P_t"], z_label,
          filename = "7-9_map_t_probe" )
# ==============================================================================
waterfall_plot( x, y, z, labels["t"], z_label, labels["P_t"],
                filename = "7-9_waterfall_t_probe", y_range = (0.0, 1.2) )
# ==============================================================================
if( len(sys.argv) > 1 and sys.argv[1] == 'animate' ):
   animated_plot( x, y, z, labels["t"], labels["P_t"], "$z = {0:7.3f} \, m$", 
                  (x[0], x[-1]), (0.0, 1.2), fps = 20, frame_prefix = "probe_",
                  filename = "7-9_animation_t_probe.avi" )
# ==============================================================================
