
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
system = System( Domain(bit_width = 400.0, samples_per_bit = 2048) )
system.add( Gaussian(peak_power = 1.0, width = 30.0, C = -20.0) )
system.add( Fibre(length = 72.0, beta = [0.0, 0.0, 1.0, 0.0],
                  gamma = 1.0, total_steps = 500) )
system.run()
# ==============================================================================
P_t = temporal_power( system.fields['fibre'] )
P_nu_normalised = spectral_power( system.fields['fibre'], True )
# ==============================================================================
double_plot( system.domain.t, P_t, system.domain.nu, P_nu_normalised,
             labels['t'], labels['P_t'], labels['nu'], labels['P_nu'],
             x_range = (50.0, 350.0), X_range = (192.6, 193.6),
             filename = "4-13" )
# ==============================================================================
