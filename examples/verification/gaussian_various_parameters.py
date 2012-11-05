
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
system = System( Domain(bit_width = 10.0) )
t = system.domain.t
nu = system.domain.nu
window_nu = system.domain.window_nu
# ==============================================================================
system.add( Gaussian(initial_phase = 3.0, width = 1.0) )
system.run()

double_plot( t, phase(system.field), t, chirp(system.field, window_nu), 
             labels["t"], labels["phi"], labels["t"], labels["chirp"],
             filename = "1 - phase_offset" )
# ==============================================================================
system.clear( True )
system.add( Gaussian(offset_nu = 0.5, width = 1.0) )
system.run()

double_plot( t, phase(system.field), t, chirp(system.field, window_nu),
             labels["t"], labels["phi"], labels["t"], labels["chirp"],
             filename = "2 - frequency_offset" )
# ==============================================================================
system.clear( True )
system.add( Gaussian(width = 1.0, C = 0.5) )
system.run()

double_plot( t, phase(system.field), t, chirp(system.field, window_nu),
             labels["t"], labels["phi"], labels["t"], labels["chirp"],
             filename = "3 - chirp" )
# ==============================================================================
system.clear( True )
system.add( Gaussian(initial_phase = 3.0, offset_nu = 0.5, width = 1.0) )
system.run()

double_plot( t, phase(system.field), t, chirp(system.field, window_nu), 
             labels["t"], labels["phi"], labels["t"], labels["chirp"],
             filename = "4 - phase_and_frequency_offset" )
# ==============================================================================
system.clear( True )
system.add( Gaussian(initial_phase = 3.0, width = 1.0, C = 0.5) )
system.run()

double_plot( t, phase(system.field), t, chirp(system.field, window_nu),
             labels["t"], labels["phi"], labels["t"], labels["chirp"],
             filename = "5 - phase_and_chirp" )
# ==============================================================================
system.clear( True )
system.add( Gaussian(width = 1.0, offset_nu = 0.5, C = 0.5) )
system.run()

double_plot( t, phase(system.field), t, chirp(system.field, window_nu),
             labels["t"], labels["phi"], labels["t"], labels["chirp"],
             filename = "6 - frequency_offset_and_chirp" )
# ==============================================================================
system.clear( True )
system.add( Gaussian(width = 1.0, initial_phase = 3.0,
                     offset_nu = 0.5, C = 0.5) )
system.run()

double_plot( t, phase(system.field), t, chirp(system.field, window_nu),
             labels["t"], labels["phi"], labels["t"], labels["chirp"],
             filename = "7 - phase_frequency_offset_and_chirp" )
# ==============================================================================
