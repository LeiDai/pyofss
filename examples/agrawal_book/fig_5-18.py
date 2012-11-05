
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
domain = Domain( bit_width = 0.4, samples_per_bit = 4096 )

s = 0.2
width = 1.0 / (s * domain.centre_omega)
beta_2 = width**2

P_ts = []
zs = [0.0, 5.0, 10.0]
# ==============================================================================
for z in zs:
   system = System( domain )
   system.add( Sech(peak_power = 1.0, width = width) )
   system.add( Fibre(length = z, gamma = 1.0, total_steps = 200,
                     self_steepening = True, beta = [0.0, 0.0, -beta_2]) )
   system.run()
   # ===========================================================================
   field = system.fields['fibre']
   P_ts.append( temporal_power(field) )
# ==============================================================================
multi_plot( system.domain.t, P_ts, zs, labels["t"], labels["P_t"],
            ["$z = {0:.0f} \, km$"], (0.175, 0.225), filename = "5-18" )
# ==============================================================================
