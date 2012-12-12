
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
from pyofss import temporal_power
from pyofss import single_plot, map_plot, waterfall_plot, labels


def generate_reference(domain, data_directory):
    """ Generate a reference field (to machine precision), used as A_true. """
    system = System(domain)
    system.add(Sech(peak_power=8.8e-3, width=(1.0 / 0.44),
                    position=0.625))
    system.add(Sech(peak_power=8.8e-3, width=(1.0 / 0.44),
                    position=0.375, offset_nu=-0.8))
    system.add(Fibre("fibre", length=400.0, beta=[0.0, 0.0, -0.1, 0.0],
                     gamma=2.2, method="ark4ip", local_error=1e-14))
    system.run()

    A_true = system.field

    filename = os.path.join(data_directory, "reference_field")
    np.save(filename, A_true)


def load_reference(data_directory):
    """ Load reference field data. """
    filename = ".".join(["reference_field", "npy"])
    filename = os.path.join(data_directory, filename)

    return np.load(filename)


def calculate_relative_error(data_directory, method, target_error, A_true):
    """ Calculate relative local error using reference field as A_true. """
    method_error = "-".join([method, "%.0e" % (target_error)])
    filename = ".".join([method_error, "npz"])
    filename = os.path.join(data_directory, filename)

    try:
        npz_data = np.load(filename)
        A_calc = npz_data["field"]

        delta_power = temporal_power(A_calc) - temporal_power(A_true)
        mean_relative_error = np.mean(np.abs(delta_power))
        mean_relative_error /= np.amax(temporal_power(A_true))
        #~mean_relative_error /= np.mean(temporal_power(A_true))
        #~mean_relative_error /= np.sum(temporal_power(A_true))

        number_of_ffts = npz_data["ffts"]

        return (number_of_ffts, mean_relative_error)
    except IOError:
        print "Error opening file: %s" % filename


def save_relative_errors(data_directory, methods, target_errors, A_true):
    """ Save calculate local relative errors to file. """
    for method in methods:
        print "%s" % method

        filename = ".".join((str(method), "dat"))
        filename = os.path.join(data_directory, filename)
        with open(filename, "w") as f:
            f.write("Number of FFTs\t\tMean relative error")

        results = []
        for target_error in target_errors:
            print "\t%.1e" % target_error

            result = calculate_relative_error(data_directory, method,
                                              target_error, A_true)
            results.append(result)

        with open(filename, "a") as f:
            for result in results:
                if result is not None:
                    f.write("\n%8d\t%.12e" % (result[0], result[1]))


def save_simulations(domain, data_directory, methods, target_errors):
    """ Save simulation data for each method to file. """
    for method in methods:
        print "%s" % method

        for target_error in target_errors:
            print "\t%.1e" % target_error

            method_error = "-".join([method, "%.0e" % (target_error)])
            filename = os.path.join(data_directory, method_error)

            system = System(domain)
            system.add(Sech(peak_power=8.8e-3, width=(1.0 / 0.44),
                            position=0.625))
            system.add(Sech(peak_power=8.8e-3, width=(1.0 / 0.44),
                            position=0.375, offset_nu=-0.8))
            system.add(Fibre("fibre", length=400.0,
                             beta=[0.0, 0.0, -0.1, 0.0], gamma=2.2,
                             method=method, local_error=target_error))
            system.run()

            A_calc = system.fields['fibre']
            storage = system["fibre"].stepper.storage

            np.savez(filename, field=A_calc, ffts=storage.fft_total)


def generate_plot(data_directory, methods):
    """ Plot number of FFTs against mean relative error for all methods. """
    labels = {"ass_simple": "Simple split-step",
              "ass_symmetric": "Symmetric split-step",
              "ass_reduced": "Reduced split-step",
              "ass_agrawal": "Agrawal split-step",
              "ass_sym_midpoint": "Symmetric split-step midpoint",
              "ass_sym_rk4": "Symmetric split-step RK4",
              "ark4ip": "Runge-Kutta in the interaction picture"}

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
            (ffts, errors) = zip(*[line.split("\t") for line in lines[1:]])

        plt.xlabel("Number of FFTs")
        plt.ylabel("Mean relative error")

        plt.plot(ffts, errors, label=labels[method],
                 linestyle=next(line_cycler), marker=next(mark_cycler))

    plt.legend(loc="lower left", frameon=False, prop={"size": 12})

    plt.xlim(50.0, 1.0e6)

    plt.xscale("log")
    plt.yscale("log")

    plt.savefig("soliton_collision_error_vs_ffts")


def generate_overview_plots(domain):
    """ Generate single, map, and waterfall plots of the soliton collision. """
    system = System(domain)
    system.add(Sech(peak_power=8.8e-3, width=(1.0 / 0.44),
                    position=0.625))
    system.add(Sech(peak_power=8.8e-3, width=(1.0 / 0.44),
                    position=0.375, offset_nu=-0.8))
    system.add(Fibre(length=400.0, beta=[0.0, 0.0, -0.1, 0.0],
                     gamma=2.2, total_steps=400, traces=100,
                     method='ARK4IP', local_error=1e-6))
    system.run()

    storage = system['fibre'].stepper.storage
    (x, y, z) = storage.get_plot_data(reduced_range=(140.0, 360.0))

    # Split step_sizes (list of tuples) into separate lists;
    # distances and steps:
    (distances, steps) = zip(*storage.step_sizes)

    print np.sum(steps)

    single_plot(distances, steps, labels["z"], "Step size, h (km)",
                filename="soliton_collision_steps")

    map_plot(x, y, z, labels["t"], labels["P_t"], labels["z"],
             filename="soliton_collision_map")

    waterfall_plot(x, y, z, labels["t"], labels["z"], labels["P_t"],
                   filename="soliton_collision_waterfall",
                   y_range=(0.0, 0.02))

if __name__ == "__main__":

    methods = ["ass_simple", "ass_symmetric", "ass_reduced",
               "ass_agrawal", "ass_sym_midpoint", "ass_sym_rk4", "ark4ip"]
    target_errors = np.logspace(-1, -12, 12)

    data_directory = "data_soliton_collision"

    domain = Domain(bit_width=400.0, samples_per_bit=4096)

    # Step one: Save each simulation to indiviual file.
    save_simulations(domain, data_directory, methods, target_errors)

    # Step two: Generate a reference field (to approximate A_true):
    generate_reference(domain, data_directory)

    # Step three: Using data from simulations,
    # save relative error calculations.
    A_true = load_reference(data_directory)
    save_relative_errors(data_directory, methods, target_errors, A_true)

    # Step four: Plot mean relative error data.
    generate_plot(data_directory, methods)

    # Step five: Generate overview plots of simulation (map, waterfall, steps).
    generate_overview_plots(domain)
