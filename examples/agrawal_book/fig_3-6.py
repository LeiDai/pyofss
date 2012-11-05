
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

from pyofss import *
# ==============================================================================
system = System( Domain(bit_width = 200.0, samples_per_bit = 2048) )
system.add( Gaussian("gaussian", peak_power = 1.0, width = 1.0) )
# ==============================================================================
system.run()
P_ts = [ temporal_power(system.fields['gaussian']) ]
# ==============================================================================
fibres = [ Fibre(length = 5.0, beta = [0.0, 0.0, 0.0, 1.0], total_steps = 100),
           Fibre(length = 5.0, beta = [0.0, 0.0, 1.0, 1.0], total_steps = 100) ]
# ==============================================================================
for fibre in fibres:
   system = System( Domain(bit_width = 200.0, samples_per_bit = 2048) )
   system.add( Gaussian(peak_power = 1.0, width = 1.0) )
   system.add( fibre )
   system.run()
   P_ts.append( temporal_power(system.fields['fibre']) )
# ==============================================================================
z_labels = [ r"$z = 0 \, km$", 
             r"$z = 5 \, km$, $\beta_2 = \, 0 \, ps / (nm \cdot km)$", \
             r"$z = 5 \, km$, $\beta_2 \neq \, 0 \, ps / (nm \cdot km)$" ]
multi_plot( system.domain.t, P_ts, z_labels, labels["t"], labels["P_t"],
            z_labels, [95.0, 115.0], filename = "3-6" )
# ==============================================================================
