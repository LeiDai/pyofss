
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

import os.path
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt

from pyofss import Domain, System, Sech, Fibre
from pyofss import temporal_power, map_plot, waterfall_plot, labels


def order_2_soliton(xi=0.0, tau=0.0):
    """ Calculate a second order soliton. """
    numerator = \
        4.0 * np.exp(0.5 * 1j * xi) * \
        (np.cosh(3.0 * tau) + 3.0 * np.exp(4.0 * 1j * xi) * np.cosh(tau))

    denominator = \
        np.cosh(4.0 * tau) + 4.0 * np.cosh(2.0 * tau) + 3.0 * np.cos(4.0 * xi)

    return numerator / denominator


def generate_reference(domain):
    """ Generate a reference field, used as A_true. """
    offset = 0.5 * domain.bit_width
    A_true = np.array(order_2_soliton(0.5 * np.pi, domain.t - offset))

    return A_true


def save_simulations(domain, data_directory, methods, steps):
    """ Save data for each method and step to file. """
    for method in methods:
        print "%s" % method

        for step in steps:
            print "\t%d" % step

            method_step = "-".join([method, str(step)])
            filename = os.path.join(data_directory, method_step)

            system = System(domain)
            system.add(Sech(peak_power=4.0, width=1.0))
            system.add(Fibre("fibre", length=0.5 * np.pi,
                             beta=[0.0, 0.0, -1.0, 0.0], gamma=1.0,
                             method=method, total_steps=step))
            system.run()

            A_calc = system.fields['fibre']
            np.save(filename, A_calc)


def calculate_relative_error(data_directory, method, step, A_true):
    """ Calculate relative local error using A_true. """
    method_step = "-".join([method, str(step)])
    filename = ".".join([method_step, "npy"])
    filename = os.path.join(data_directory, filename)

    A_calc = np.load(filename)

    delta_power = temporal_power(A_calc) - temporal_power(A_true)
    mean_relative_error = np.mean(np.abs(delta_power))
    mean_relative_error /= np.amax(temporal_power(A_true))

    return mean_relative_error


def save_relative_errors(data_directory, methods, steps, A_true):
    """ Save calculated relative errors to file. """
    for method in methods:
        print "%s" % method

        filename = ".".join((str(method), "dat"))
        filename = os.path.join(data_directory, filename)
        with open(filename, "w") as f:
            f.write("Steps\t\tMean relative error")

        results = []
        for step in steps:
            print "\t%d" % step

            error = calculate_relative_error(
                data_directory, method, step, A_true)
            results.append((step, error))

        with open(filename, "a") as f:
            for result in results:
                f.write("\n%6d\t%.12e" % result)


def generate_plot(data_directory, methods):
    """ Plot steps against mean relative error. """
    labels = {"ss_simple": "Simple split-step",
              "ss_symmetric": "Symmetric split-step",
              "ss_reduced": "Reduced split-step",
              "ss_agrawal": "Agrawal split-step",
              "ss_sym_midpoint": "Symmetric split-step midpoint",
              "ss_sym_rk4": "Symmetric split-step RK4",
              "rk4ip": "Runge-Kutta in the interaction picture"}

    # Generate a range of line and mark styles to cycle through:
    lines = ["-", ":", "--"]
    line_cycler = cycle(lines)

    marks = ["*", "o", "s", "D", "^"]
    mark_cycler = cycle(marks)

    for method in methods:
        filename = ".".join((str(method), "dat"))
        filename = os.path.join(data_directory, filename)
        with open(filename, "r") as f:
            lines = f.readlines()
            # Skip first line, split on tab character,
            # and unpack into two lists:
            (steps, results) = zip(*[line.split("\t") for line in lines[1:]])

        plt.xlabel("Number of steps")
        plt.ylabel("Mean relative error")

        plt.plot(steps, results, label=labels[method],
                 linestyle=next(line_cycler), marker=next(mark_cycler))

    plt.legend(loc="lower left", frameon=False, prop={"size": 12})

    plt.xlim(20.0, 2.0e5)

    plt.xscale("log")
    plt.yscale("log")

    plt.savefig("soliton_error_vs_steps")


def generate_map_and_waterfall_plots(domain):
    """ Generate map and waterfall plots to visualise pulse propagation. """
    system = System(domain)
    system.add(Sech(peak_power=4.0, width=1.0))
    system.add(Fibre("fibre", length=0.5 * np.pi,
                     beta=[0.0, 0.0, -1.0, 0.0], gamma=1.0,
                     method="rk4ip", total_steps=1000, traces=50))
    system.run()

    storage = system['fibre'].stepper.storage
    (x, y, z) = storage.get_plot_data(reduced_range=(95.0, 105.0))

    map_plot(x, y, z, labels["t"], labels["P_t"], labels["z"],
             filename="soliton_map")

    waterfall_plot(x, y, z, labels["t"], labels["z"], labels["P_t"],
                   filename="soliton_waterfall", y_range=(0.0, 16.0))

if __name__ == "__main__":

    methods = ["ss_simple", "ss_symmetric", "ss_reduced",
               "ss_agrawal", "ss_sym_midpoint", "ss_sym_rk4", "rk4ip"]

    steps = [50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]

    data_directory = "data_soliton"

    domain = Domain(bit_width=200.0, samples_per_bit=4096)

    # Step one: Save each simulation to individual file.
    save_simulations(domain, data_directory, methods, steps)

    # Step two: Using data from simulations, save relative error calculations.
    A_true = generate_reference(domain)
    save_relative_errors(data_directory, methods, steps, A_true)

    # Step three: Plot mean relative error data.
    generate_plot(data_directory, methods)

    # Step four: Plot map and waterfall plots of the simulation.
    generate_map_and_waterfall_plots(domain)
