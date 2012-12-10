
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

from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import subprocess

labels = {"t": "Time, $t \, (ps)$",
          "nu": r"Frequency, $\nu \, (THz)$",
          "Lambda": r"Wavelength, $\lambda \, (nm)$",
          "P_t": "Power, $|A(z, t)|^2 \, (W)$",
          "P_nu": r"Power, $|\tilde{A}(z, \nu)|^2 \, (a.u.)$",
          "P_lambda": r"Power, $|\tilde{A}(z, \lambda|^2 \, (a.u.)$",
          "z": "Fibre length, $z \, (km)$",
          "phi": "Phase, $\phi(t) \, (rad)$",
          "chirp": "Frequency chirp, $\delta \omega \, (rad / ps)$",
          "t_normal": r"Normalised time, $\frac{t}{T_0}$",
          "xi": r"Normalised distance, $\xi = \frac{z}{L_D}$",
          "xi_prime": r"Normalised distance, $\xi' = \frac{z}{L_D'}$"}


def map_plot(x, y, z, x_label="", y_label="", z_label="",
             interpolation='lanczos', use_colour=True,
             filename="", y_range=None):
    """
    :param Dvector x: First axis
    :param Dvector y: Second axis
    :param Dvector z: Third axis
    :param string x_label: Label for first axis
    :param string y_label: Label for second axis
    :param string z_label: Label for third axis
    :param string interpolation: Type of interpolation to use
    :param bool use_colour: Use colour plot (else black and white plot)
    :param string filename: Location to save file. If "", use interactive plot
    :param Dvector y_range: Change range of second axis if not None

    Generate a map plot.
    """
    print "\nGenerating map_plot..."
    plt.clf()

    plt.xlabel(x_label)
    plt.ylabel(z_label)

    if use_colour:
        cmap = cm.jet
    else:
        cmap = cm.gray

    im = plt.imshow(y, interpolation=interpolation, origin='lower',
                    aspect='auto', cmap=cmap,
                    extent=(x[0], x[-1], z[0], z[-1]))

    # Draw a colour bar showing mapping of colours to values for y-axis:
    cb = plt.colorbar(im, use_gridspec=True)
    cb.set_label(y_label)

    if y_range is not None:
        im.set_clim(y_range)

    # Avoid overlapping axis text:
    plt.tight_layout()

    if filename:
        plt.savefig(filename)
        print "Wrote file", filename
    else:
        plt.show()


def waterfall_plot(x, y, z, x_label="", y_label="", z_label="",
                   use_poly=True, alpha=0.2, filename="",
                   y_range=None, x_range=None):
    """
    :param Dvector x: First axis
    :param Dvector y: Second axis
    :param Dvector z: Third axis
    :param string x_label: Label for first axis
    :param string y_label: Label for second axis
    :param string z_label: Label for third axis
    :param bool use_poly: Whether to use filled polygons
    :param double alpha: Set transparency of filled polygons
    :param string filename: Location to save file. If "", use interactive plot
    :param Dvector y_range: Change range of second axis if not None
    :param Dvector x_range: Change range of first axis if not None

    Generate a waterfall plot.
    """
    print "\nGenerating waterfall_plot..."
    plt.clf()

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    data = [zip(x, y[n]) for n, z_n in enumerate(z)]

    if use_poly:
        from matplotlib.collections import PolyCollection
        # Parameter closed MUST be set False if axis does not start at zero:
        curves = PolyCollection(data, closed=False)
    else:
        from matplotlib.collections import LineCollection
        curves = LineCollection(data)

    curves.set_alpha(alpha)

    ax.add_collection3d(curves, zs=z, zdir='y')

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)

    if x_range is None:
        ax.set_xlim3d(x[0], x[-1])
    else:
        ax.set_xlim3d(x_range)

    # Note the use of z as a parameter for ylim and also the reverse:
    ax.set_ylim3d(z[0], z[-1])

    if y_range is not None:
        ax.set_zlim3d(y_range)

    # Avoid overlapping axis text:
    plt.tight_layout()

    if filename:
        plt.savefig(filename)
        print "Wrote file", filename
    else:
        plt.show()


def single_plot(x, y, x_label="", y_label="", label="",
                x_range=None, y_range=None, use_fill=True,
                alpha=0.2, filename="", style="b-", fill_colour="b"):
    """
    :param Dvector x: First axis
    :param Dvector y: Second axis
    :param string x_label: Label for first axis
    :param string y_label: Label for second axis
    :param string label: Label for plot area
    :param Dvector x_range: Change range of first axis if not None
    :param Dvector y_range: Change range of second axis if not None
    :param bool use_fill: Fill area between plot line and axis
    :param double alpha: Transparency of filled area
    :param string filename: Location to save file. If "", use interactive plot
    :param string style: Style of plot data. E.g. "g+" plots green plus signs
    :param string fill_colour: Filled region colour. E.g. "b" uses a blue fill

    Generate a single plot.
    """
    print "\nGenerating single_plot..."
    plt.clf()

    plt.plot(x, y, style, label=label)

    if use_fill:
        plt.fill_between(x, y, color=fill_colour, alpha=alpha)

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if x_range is not None:
        plt.xlim(x_range)
    if y_range is not None:
        plt.ylim(y_range)

    if label is not "":
        plt.legend()

    if filename:
        plt.savefig(filename)
        print "Wrote file", filename
    else:
        plt.show()


def double_plot(x, y, X, Y, x_label="", y_label="", X_label="",
                Y_label="", x_range=None, y_range=None, X_range=None,
                Y_range=None, use_fill=True, alpha=0.2, filename=""):
    """
    :param Dvector x: First axis of upper plot
    :param Dvector y: Second axis of upper plot
    :param Dvector X: First axis of lower plot
    :param Dvector Y: Second axis of lower plot
    :param string x_label: Label for first axis of upper plot
    :param string y_label: Label for second axis of upper plot
    :param string X_label: Label for first axis of lower plot
    :param string Y_label: Label for second axis of lower plot
    :param Dvector x_range: Change first axis range of upper plot if not None
    :param Dvector y_range: Change second axis range of upper plot if not None
    :param Dvector X_range: Change first axis range of lower plot if not None
    :param Dvector Y_range: Change second axis range of lower plot if not None
    :param bool use_fill: Fill area between plot line and axis
    :param double alpha: Transparency of filled area
    :param string filename: Location to save file. If "", use interactive plot

    Generate a double plot. The two plots will be arranged vertically,
    one above the other.
    """
    print "\nGenerating double_plot..."

    plt.clf()

    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    ax1.plot(x, y)
    ax2.plot(X, Y)

    if use_fill:
        ax1.fill_between(x, y, alpha=alpha)
        ax2.fill_between(X, Y, alpha=alpha)

    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)

    ax2.set_xlabel(X_label)
    ax2.set_ylabel(Y_label)

    if x_range is not None:
        ax1.set_xlim(x_range)
    if y_range is not None:
        ax1.set_ylim(y_range)

    if X_range is not None:
        ax2.set_xlim(X_range)
    if Y_range is not None:
        ax2.set_ylim(Y_range)

    # Avoid overlapping axis text:
    plt.tight_layout()

    if filename:
        plt.savefig(filename)
        print "Wrote file", filename
    else:
        plt.show()


def multi_plot(x, ys, zs=[], x_label="", y_label="",
               z_labels=[""], x_range=None, y_range=None,
               use_fill=True, alpha=0.2, filename=""):
    """
    :param Dvector x: First axis
    :param VDvector ys: Array of second axis
    :param VDvector zs: Array of third axis (used for legend, and colours)
    :param string x_label: Label for first axis
    :param string y_label: Label for second axis
    :param Vstring: List of strings for third axis
    :param Dvector x_range: Change range of first axis if not None
    :param Dvector y_range: Change range of second axis if not None
    :param bool use_fill: Fill area between plot line and axis
    :param double alpha: Transparency of filled area
    :param string filename: Location to save file. If "", use interactive plot

    Generate multiple overlapping plots.
    """
    print "\nGenerating multi_plot..."
    plt.clf()

    if len(zs) < 6:
        colours = ['blue', 'green', 'red', 'orange', 'purple']
    else:
        colours = ['blue' for z in zs]

    for n, y in enumerate(ys):
        if len(z_labels) > 1:
            plt.plot(x, y, label=z_labels[n])
        else:
            # If z_label does not contain {0} place-holder then format is
            # ignored:
            plt.plot(x, y, label=z_labels[0].format(zs[n]))

        if use_fill:
            plt.fill_between(x, y, alpha=alpha, color=colours[n])

        plt.xlabel(x_label)
        plt.ylabel(y_label)

    if x_range is not None:
        plt.xlim(x_range)
    if y_range is not None:
        plt.ylim(y_range)

    # Try to place legend in best location (avoiding plot):
    leg = plt.legend(loc="best", fancybox=True)
    # Use a semi-transparent legend box in case of overlap with plot:
    leg.get_frame().set_alpha(0.5)

    if filename:
        plt.savefig(filename)
        print "Wrote file", filename
    else:
        plt.show()


def quad_plot(x, ys, zs, x_label="", y_label="", z_labels=[""],
              x_range=None, y_range=None,
              use_fill=True, alpha=0.2, filename=""):
    """
    :param Dvector x: First axis
    :param VDvector ys: Array of second axis
    :param VDvector zs: Array of third axis
    :param string x_label: Label for first axis
    :param string y_label: Label for second axis
    :param Vstring: List of strings for third axis
    :param Dvector x_range: Change range of first axis if not None
    :param Dvector y_range: Change range of second axis if not None
    :param bool use_fill: Fill area between plot line and axis
    :param double alpha: Transparency of filled area
    :param string filename: Location to save file. If "", use interactive plot

    Generate four plots arranged in a two-by-two square.
    """
    print "\nGenerating quad_plot..."
    plt.clf()

    colours = ['blue', 'green', 'red', 'orange']

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True,
                                                 sharey=True)
    axs = [ax1, ax2, ax3, ax4]

    for n, y in enumerate(ys):
        if len(z_labels) > 1:
            axs[n].plot(x, y, label=z_labels[n], color=colours[n])
        else:
            axs[n].plot(x, y, label=z_labels[0].format(zs[n]),
                        color=colours[n])

        if use_fill:
            axs[n].fill_between(x, y, alpha=alpha, color=colours[n])

        if x_range is not None:
            axs[n].set_xlim(x_range)
        if y_range is not None:
            axs[n].set_ylim(y_range)

        # Try to place legend in best location (avoiding plot):
        leg = axs[n].legend(loc="best", fancybox=True)
        # Use a semi-transparent legend box in case of overlap with plot:
        leg.get_frame().set_alpha(0.5)

    ax1.set_ylabel(y_label)
    ax3.set_ylabel(y_label)
    ax3.set_xlabel(x_label)
    ax4.set_xlabel(x_label)

    # Avoid overlapping axis text:
    plt.tight_layout()

    if filename:
        plt.savefig(filename)
        print "Wrote file", filename
    else:
        plt.show()


def animated_plot(x, y, z, x_label="", y_label="", z_label="",
                  x_range=None, y_range=None, alpha=0.2, fps=5,
                  clear_temp=True, frame_prefix="_tmp", filename=""):
    """
    :param Dvector x: First axis
    :param VDvector y: Array of second axis
    :param Dvector z: Third axis
    :param string x_label: Label for first axis
    :param string y_label: Label for second axis
    :param string z_label: label for third axis
    :param Dvector x_range: Change range of first axis if not None
    :param Dvector y_range: Change range of second axis if not None
    :param double alpha: Transparency of filled area
    :param Uint fps: Number of frames to show every second (frames per second)
    :param bool clear_temp: Remove temporary image files after completion
    :param string frame_prefix: Prefix used for each temporary file
    :param string filename: Location to save file. If "" then interactive plot

    Generate an animated plot, either interactive or saved as a video.
    """
    print "\nGenerating animated_plot..."
    #~plt.clf()

    fig = plt.figure()

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if x_range is not None:
        plt.xlim(x_range)
    if y_range is not None:
        plt.ylim(y_range)

    ax = fig.add_subplot(111)

    images = []
    for i in range(len(y)):
        line, = plt.plot(x, y[i], 'b-')
        label = plt.legend([line], [z_label.format(z[i])])
        # The following line is necessary when using ani.save(), otherwise the
        # legend is only shown on the last frame!
        ax.add_artist(label)
        fill = plt.fill_between(x, y[i], facecolor='blue', alpha=alpha)
        images.append((line, label, fill))

    ani = animation.ArtistAnimation(fig, images, blit=True,
                                    interval=1e3 / fps)

    if filename:
        ani.save(filename, fps=fps, clear_temp=clear_temp,
                 frame_prefix=frame_prefix)
        print "Wrote file", filename
    else:
        plt.show()


def convert_video(filename, output="ogv"):
    """
    :param string filename: Location of movie file to convert
    :param string output: Convert movie to output video type

    Convert video to ogv, webm, or mp4 type.
    """
    import os
    out_file = os.path.splitext(filename)[0]

    if output == "ogv":
        # Use Theora/Vorbis codecs
        command = ('ffmpeg2theora', '{0}'.format(filename))
        #~command = ('ffmpeg', '-i', '{0}'.format( filename ), '-b', '1500k',
                   #~'-vcodec', 'libtheora', '-acodec', 'libvorbis', '-ab',
                   #~'160000', '-g', '30',
                   #~'{0}'.format( '.'.join((out_file, output)) ))
    elif output == "webm":
        # Use VP8/Vobis codecs
        command = ('ffmpeg', '-i', '{0}'.format(filename), '-b', '1500k',
                   '-vcodec', 'libvpx', '-acodec', 'libvorbis', '-ab',
                   '160000', '-f', 'webm', '-g', '30',
                   '{0}'.format('.'.join([out_file, output])))
    elif output == "mp4":
        # Use H.264/ACC codecs
        command = ('ffmpeg', '-i', '{0}'.format(filename), '-b', '1500k',
                   '-vcodec', 'libx264', '-vpre', 'slow', '-vpre', 'baseline',
                   '-g', '30', '{0}'.format('.'.join([out_file, output])))

    print "\n\nConverting video using command:\n%s\n\n" % ' '.join(command)
    subprocess.check_call(command)

    print "\nWrote file", '.'.join([out_file, output])
