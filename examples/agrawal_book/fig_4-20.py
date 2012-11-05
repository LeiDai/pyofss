
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
# ==============================================================================
system = System( domain )
system.add( Gaussian(peak_power = 1.0, width = width) )
system.add( Fibre(length = 20.0, gamma = 1.0, total_steps = 200,
                  self_steepening = True, method = "RK4IP") )
system.run()
# ==============================================================================
P_nu_normalised = spectral_power( system.fields['fibre'], True )
# ==============================================================================
single_plot( system.domain.nu, P_nu_normalised, labels["nu"], labels["P_nu"],
             filename = "4-20", x_range = (145.1, 256.1) )
# ==============================================================================
