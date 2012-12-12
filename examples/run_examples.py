
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

import os
import fnmatch
import Image
from PIL import ImageChops


def find_files(directory, pattern, ignored=[], sort_results=False):
    """
    Input top-level directory and the pattern (e.g. "*.py") to search.
    Setting sort_results to True will sort the output list.
    Directories listed in ignored will not be included for pattern matching.
    Output a list of tuples in the form:
        [ (path_0, file_0), (path_1, file_1), ... ]
    """
    matches = []

    for path, dirs, files in os.walk(directory):
        for dir_ in ignored:
            if dir_ in dirs:
                dirs.remove(dir_)
        for file_ in fnmatch.filter(files, pattern):
            matches.append((path, file_))

    if(sort_results):
        matches.sort()

    return matches


def compare_images(first_image, second_image):
    """
    Compare the contents of two image files.
    Return True if both image files match exactly.
    """
    return ImageChops.difference(first_image, second_image).getbbox() is None


def compare_image_directory(reference_directory_of_images):
    """
    For each example output plot, check for a file with the same name as in
    the user-provided directory. If the reference plot exists, check that the
    contents of each are identical. Generate a list of those that do match
    and those that fail to match.
    """
    print reference_directory_of_images

    # Generate paths and figure names from examples:
    images = find_files(stored_path, "*.png", sort_results=True)
    figure_paths, figures = zip(*images)

    # Generate paths and figure names from reference plots:
    ref_images = find_files(reference_directory_of_images, "*.png")
    ref_paths, ref_figures = zip(*ref_images)

    match_pass = []
    match_fail = []

   # For each example figure, check if a reference figure exists
    for fig_index, figure in enumerate(figures):
        for ref_figure_index, ref_figure in enumerate(ref_figures):
            if(figure == ref_figure):
                # Join the path and figure name
                fig = os.path.join(figure_paths[fig_index], figure)
                ref_fig = os.path.join(ref_paths[ref_figure_index], ref_figure)

                # Open figure and reference figures
                im_fig = Image.open(fig)
                im_ref = Image.open(ref_fig)

                # Use compare_images() function which returns a boolean value:
                if(compare_images(im_fig, im_ref)):
                    match_pass.append(figure)
                    print "PASS - %s" % figure
                else:
                    match_fail.append(figure)
                    print "FAIL - %s" % figure

    print "\nThe following images matched their reference image:"
    print match_pass

    print "\nThe following images failed to match their reference image:"
    print match_fail


def run_scripts(stored_path):
    """ Provide a user interface for running Python scripts. """
    import subprocess

    matches = find_files(stored_path, "fig_*.py", sort_results=True)
    # Split a list of tuples into individual lists:
    paths, examples = zip(*matches)

    for index, script in enumerate(examples):
        print "%i - %s" % (index, script)

    user_input = raw_input("Choose script(s) to run: ")
    # Split string into whitespace separated components, return a list
    options = user_input.split()
    # Make sure something was entered:
    if(len(options) > 0):
        # If the second entry is '-', assume a range has been entered
        if(len(options) == 3 and options[1] == '-'):
            indices = range(int(options[0]), int(options[2]) + 1)
        else:
            indices = [int(option) for option in options]
        print "Selected scripts: ", indices

        use_animations = raw_input("Press A to allow animations: ")
        if(str.upper(use_animations) == 'A'):
            allow_animations = True
        else:
            allow_animations = False

        confirm = raw_input("Press Y to run chosen scripts: ")
        if(str.upper(confirm) == 'Y'):
            for index in indices:
                try:
                    print "\n-------------------------------------------------"
                    print "\tRunning script: %i - %s\n" % (index,
                                                           examples[index])
                    os.chdir(paths[index])
                    if(allow_animations):
                        subprocess.call(("python", examples[index], "animate"))
                    else:
                        subprocess.call(("python", examples[index]))
                    print "-------------------------------------------------\n"
                except:
                    print "Problem running selected script!"
                finally:
                    # Important to restore original directory:
                    os.chdir(stored_path)
        else:
            print "Not running any scripts. Exiting..."
    else:
        print("No scripts selected. Exiting...")

if __name__ == "__main__":

    # Store current directory
    stored_path = os.getcwd()
    print "Current path: %s" % os.getcwd()

    print "Select menu option number:"
    prompt = "\n1 - Run scripts; \n2 - Compare output plots to reference.\n-->"
    option = int(raw_input(prompt))

    if(option == 1):
        run_scripts(stored_path)
    elif(option == 2):
        input_dir = raw_input("Enter directory of reference plots:\n-->")
        compare_image_directory(input_dir)
    else:
        print "Option does not exist"
