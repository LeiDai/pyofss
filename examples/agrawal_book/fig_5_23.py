
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
from pyofss.modules.dispersion import convert_dispersion_to_physical
# ==============================================================================
domain = Domain( bit_width = 2.0, samples_per_bit = 8192,
                 centre_nu = lambda_to_nu(1550.0) )

s = 0.05
width = 1.0 / (s * domain.centre_omega)

tau_R = 0.1
T_R = tau_R * width

D = 16.0
beta_2 = convert_dispersion_to_physical( D )[0]

delta_3 = 0.03
beta_3 = 6.0 * np.abs(beta_2) * width * delta_3

beta = [0.0, 0.0, beta_2, beta_3]

L_D = width**2 / np.abs(beta_2)
length = 4.0 * L_D

N = 2.0
P_0 = 1.0
gamma = N**2 / (L_D * P_0)
# ==============================================================================
system = System( domain )
system.add( Sech(peak_power = P_0, width = width) )
system.add( Fibre(length = length, gamma = gamma, beta = beta, rs_factor = T_R,
                  raman_scattering = True, self_steepening = True,
                  total_steps = 200, traces = 100, method = 'ARK4IP') )
system.run()
# ==============================================================================
storage = system['fibre'].stepper.storage
(x, y, z_temp) = storage.get_plot_data( False, (71.9, 314.9), True )
z_label = "Fibre length, $z \, (cm)$"
z = z_temp * 1.0e5
# ==============================================================================
map_plot( x, y, z, labels["nu"], labels["P_nu"], z_label,
          filename = "5-23_map_nu" )
# ==============================================================================
waterfall_plot( x, y, z, labels["nu"], z_label, labels["P_nu"],
                filename = "5-23_waterfall_nu", y_range = (0.0, 1.1) )
# ==============================================================================
if( len(sys.argv) > 1 and sys.argv[1] == 'animate' ):
   animated_plot( x, y, z, labels["nu"], labels["P_nu"],
                  "$z = {0:7.3f} \, cm$", (x[0], x[-1]), (0.0, 1.1), fps = 20,
                  frame_prefix = "nu_", filename = "5-23_animation_nu.avi" )
# ==============================================================================
(x, y, z_temp) = storage.get_plot_data( reduced_range = (0.9, 1.4) )
z = z_temp * 1.0e5
# ==============================================================================
map_plot( x, y, z, labels["t"], labels["P_t"], z_label,
          filename = "5-23_map_t" )
# ==============================================================================
waterfall_plot( x, y, z, labels["t"], z_label, labels["P_t"],
                filename = "5-23_waterfall_t", y_range = (0.0, 2.6) )
# ==============================================================================
if( len(sys.argv) > 1 and sys.argv[1] == 'animate' ):
   animated_plot( x, y, z, labels["t"], labels["P_t"], "$z = {0:7.3f} \, cm$", 
                  (x[0], x[-1]), (0.0, 2.6), fps = 20, frame_prefix = "t_",
                  filename = "5-23_animation_t.avi" )
# ==============================================================================
