
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
A_fs = []
Cs = [0.0, 10.0, -10.0, -20.0]
for C in Cs:
   system = System( Domain(bit_width = 100.0, samples_per_bit = 2048) )
   system.add( Gaussian(width = 1.0, peak_power = 1.0, C = C) )
   system.add( Fibre('fibre', length = 1.0, gamma = 4.5 * np.pi) )
   system.run()
   A_fs.append( system.fields['fibre'] )
# ==============================================================================
P_nus = [ spectral_power(A_f, True) for A_f in A_fs ]
# ==============================================================================
quad_plot( system.domain.nu, P_nus, Cs, labels["nu"], labels["P_nu"],
           ["$C = {0:.0f}$"], (189.1, 197.1), filename = "4-5" )
# ==============================================================================
