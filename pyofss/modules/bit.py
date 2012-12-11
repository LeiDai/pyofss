
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

import random
import copy


class Bit(object):
    """
    :param double position: Position of pulse
    :param double width: Width of pulse
    :param double peak_power: Peak power of pulse
    :param double offset_nu: Offset frequency of pulse
    :param Uint m: Order parameter
    :param double C: Chirp parameter
    :param double initial_phase: Initial phase of the pulse
    :param bool using_fwhm: Determines whether the width parameter is a
                            full-width at half maximum measure, or a
                            half-width at 1/e maximum measure

    Each bit is represented by a pulse.
    """
    def __init__(self, position=0.5, width=10.0, peak_power=1e-3,
                 offset_nu=0.0, m=1, C=0.0, initial_phase=0.0,
                 channel=0, using_fwhm=False):

        self.data = {"position": position, "width": width,
                     "peak_power": peak_power, "offset_nu": offset_nu,
                     "m": m, "C": C, "initial_phase": initial_phase,
                     "channel": channel, "using_fwhm": using_fwhm}

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def __delitem__(self, key):
        del self.data[key]

    def __call__(self):
        return self.data

    #~def __repr__(self):
        #~return "Bit with data:\n", self.data


def generate_bitstream():
    """ Generate a bitstream from a single bit. """
    return [Bit()]


class Bit_stream(object):
    """
    A bit_stream consists of a number of bits.
    """
    def __init__(self):

        self.bits = []
        self.prbs = []

    def __getitem__(self, key):
        return self.bits[key]

    def add(self, bit):
        self.bits.append(bit)

    def __add__(self, bit):
        self.add(bit)

    #~def __repr__(self):
        #~return "Bit stream containing bits:\n", self.bits


def generate_prbs(domain, bit=Bit(), power_jitter=0, ghost_power=0):
    """
    :param object domain: Domain to use for calculations
    :param object bit: Bit parameters to use for generating a pulse
    :param double power_jitter: Variation as percentage of peak power
    :param double ghost_power: Power range as percentage of peak power
    :return: Bitstream
    :rtype: object

    Generate a list of bits, each of which describes a pulse.
    Each pulse is a variation of the pulse described by input parameter: 'bit'.
    """
    if not (0 <= power_jitter <= 100):
        raise Exception("power_jitter is out of range. Must be in [0, 100]")

    if not (0 <= ghost_power <= 100):
        raise Exception("ghost_power is out of range. Must be in [0, 100]")

    # The peak power of the bit will be stored as the mean peak power of ones:
    mean_peak_power = bit["peak_power"]

    # One bits will have peak power in the range given by the following tuple:
    one_range = (mean_peak_power * (1.0 - 0.01 * power_jitter),
                 mean_peak_power * (1.0 + 0.01 * power_jitter))

    # Zero bits will have peak_power in range given by the following tuple:
    ghost_range = (0.0, 0.01 * ghost_power * mean_peak_power)

    # Output details:
    #~print "mean_peak_power = %f W" % mean_peak_power
    #~print "one_range = [%.5f, %.5f] W" % one_range
    #~print "ghost_range = [%.5f, %.5f] W" % ghost_range

    # Store the relative position for each pulse within a bit width:
    relative_position = bit["position"]

    bit_stream = []
    for b in range(0, domain.total_bits):
        # Decide if one or zero bit:
        is_one = random.randint(0, 1)

        # Generate a pulse with peak_power in one_range or ghost_range:
        if is_one:
            peak_power = random.uniform(*one_range)
        else:
            peak_power = random.uniform(*ghost_range)

        # Add bit to bit_stream. Note that without using deepcopy,
        # peak_power would be the same for each bit!
        next_bit = copy.deepcopy(bit)
        next_bit["peak_power"] = peak_power
        next_bit["position"] = (b + relative_position) / domain.total_bits
        bit_stream.append(next_bit)

    return bit_stream
