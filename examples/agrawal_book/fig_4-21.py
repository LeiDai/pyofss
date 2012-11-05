
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

from pyofss import *
# ==============================================================================
domain = Domain( bit_width = 4.0, samples_per_bit = 4096 )

s = 0.01
width = 1.0 / (s * domain.centre_omega)
gamma = 100.0 / width**2

P_ts = []
P_nus = []
length = [20.0 / gamma, 40.0 / gamma]
# ==============================================================================
for l in length:
   system = System( domain )
   system.add( Gaussian(peak_power = 1.0, width = width) )
   system.add( Fibre(length = l, gamma = gamma, total_steps = 200,
                     self_steepening = True, beta = [0.0, 0.0, 1.0]) )
   system.run()
   # ===========================================================================
   field = system.fields['fibre']
   P_ts.append( temporal_power(field) )
   P_nus.append( spectral_power(field, True) )
# ==============================================================================
double_plot( system.domain.t, P_ts[0], system.domain.nu, P_nus[0],
             labels["t"], labels["P_t"], labels["nu"], labels["P_nu"],
             x_range = (1.5, 2.5), X_range = (146.1, 240.1),
             filename = "4-21a" )
# ==============================================================================
double_plot( system.domain.t, P_ts[1], system.domain.nu, P_nus[1],
             labels["t"], labels["P_t"], labels["nu"], labels["P_nu"],
             x_range = (1.0, 3.0), X_range = (146.1, 240.1),
             filename = "4-21b" )
# ==============================================================================
