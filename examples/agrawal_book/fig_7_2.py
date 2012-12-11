
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
system = System( Domain(bit_width = 200.0, samples_per_bit = 4096,
                        channels = 2) )
system.add( Gaussian(width = 1.0, peak_power = 1.0, channel = 0) )
system.add( Gaussian(width = 1.0, peak_power = 0.5, channel = 1) )
system.add( Fibre('fibre', length = 40.0, gamma = [1.0, 1.2],
            beta = [ [0.0, 0.0, 0.0, 0.0], [0.0, 0.125, 0.0, 0.0] ],
            sim_type = 'wdm', total_steps = 400, method = 'RK4IP') )
system.run()
# ==============================================================================
A_fs = system.fields['fibre']
# ==============================================================================
P_nu0 = spectral_power( A_fs[0], True )
P_nu1 = spectral_power( A_fs[1], True )
# ==============================================================================
double_plot( system.domain.nu, P_nu0, system.domain.nu, P_nu1,
             labels["nu"], labels["P_nu"], labels["nu"], labels["P_nu"],
             x_range = (181.1, 204.1), X_range = (181.1, 204.1),
             filename = "7-2" )
# ==============================================================================
