
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
nu_0 = 193.1
nu_1 = 1.2 * nu_0

offset_nu = 0.2 * 193.1
# ==============================================================================
system = System( Domain(bit_width = 30.0, samples_per_bit = 8192,
                        channels = 2) )
system.add( Gaussian(width = 1.0, peak_power = 100.0, channel = 0) )
system.add( Gaussian(width = 1.0, peak_power = 1.0, channel = 1,
                     offset_nu = offset_nu) )
system.add( Fibre('fibre', length = 0.4, gamma = [1.0, 1.2],
            beta = [ [0.0, 0.0, 0.0, 0.0], [0.0, 10.0, 0.0, 0.0] ],
            sim_type = 'wdm', method = 'ARK4IP') )
system.run()
# ==============================================================================
A_fs = system.fields['fibre']

P_nu0 = spectral_power( A_fs[0], True )
P_nu1 = spectral_power( A_fs[1], True )
# ==============================================================================
double_plot( system.domain.nu, P_nu0, system.domain.nu, P_nu1,
             labels["nu"], labels["P_nu"], labels["nu"], labels["P_nu"],
             x_range = (nu_0 - 8.0, nu_0 + 8.0), 
             X_range = (nu_1 - 8.0, nu_1 + 8.0),
             filename = "7-7" )
# ==============================================================================
