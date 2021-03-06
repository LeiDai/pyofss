
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

import sys
from pyofss import Domain, System, Gaussian, Fibre
from pyofss import map_plot, waterfall_plot, animated_plot, labels

system = System(Domain(bit_width=200.0, samples_per_bit=2048))
system.add(Gaussian(peak_power=1.0, width=1.0))
system.add(Fibre(length=5.0, beta=[0.0, 0.0, -1.0, 0.0],
                 gamma=1.0, traces=50))
system.run()

storage = system['fibre'].stepper.storage
(x, y, z) = storage.get_plot_data(False, (192.1, 194.1), normalised=True)

map_plot(x, y, z, labels["nu"], labels["P_nu"], labels["z"],
         filename="4-9_map_nu")

waterfall_plot(x, y, z, labels["nu"], labels["z"], labels["P_nu"],
               filename="4-9_waterfall_nu", y_range=(0.0, 1.9))

if (len(sys.argv) > 1) and (sys.argv[1] == 'animate'):
    animated_plot(x, y, z, labels["nu"], labels["P_nu"],
                  r"$z = {0:7.3f} \, km$", (x[0], x[-1]), (0.0, 1.9), fps=10,
                  frame_prefix="nu_", filename="4-9_animation_nu.avi")


(x, y, z) = storage.get_plot_data(reduced_range=(90.0, 110.0))

map_plot(x, y, z, labels["t"], labels["P_t"], labels["z"],
         filename="4-9_map_t")

waterfall_plot(x, y, z, labels["t"], labels["z"], labels["P_t"],
               filename="4-9_waterfall_t")

if (len(sys.argv) > 1) and (sys.argv[1] == 'animate'):
    animated_plot(x, y, z, labels["t"], labels["P_t"], r"$z = {0:7.3f} \, km$",
                  (x[0], x[-1]), (0.0, 1.1), fps=10, frame_prefix="t_",
                  filename="4-9_animation_t.avi")
