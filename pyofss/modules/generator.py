
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

from scipy import sqrt, log, power

from bit import Bit
from gaussian import Gaussian
from sech import Sech


def convert_pulse_width(width, to_fwhm=True, shape="gaussian", m=1):
    """
    :param float width: Pulse width
    :param bool to_fwhm: Whether to convert to FWHM measure
    :param string shape: Shape of the pulse
    :param Uint m: Order parameter, used for super-Gaussian pulses

    Helper function to convert pulse widths between FWHM and HWIeM measures.
    """
    if shape.lower() == "gaussian":
        if to_fwhm:
            return width * 2.0 * power(log(2.0), 1.0 / (2 * m))
        else:
            return width * 0.5 / power(log(2.0), 1.0 / (2 * m))
    elif shape.lower() == "sech":
        if to_fwhm:
            return width * 2.0 * log(1.0 + sqrt(2.0))
        else:
            return width * 0.5 / log(1.0 + sqrt(2.0))
    else:
        print "Pulse shape not recognised: %s" % shape


class Generator(object):
    """
    :param string name: Name of this module
    :param object bit_stream: An array of Bit objects
    :param Uint channel: Channel to modify

    Generate a pulse with Gaussian or hyperbolic secant shape.
    Add this pulse to appropriate field, determined by channel.
    """
    def __init__(self, name="generator", bit_stream=None, channel=0):
        self.name = name
        self.bit_stream = bit_stream
        self.channel = channel

        # Always make a single bit as default, if no bitstream is passed:
        if self.bit_stream is None:
            self.bit_stream = [Bit()]

        self.field = None
        self.shape = None

    def __call__(self, domain, field):
        """
        :param object domain: A Domain
        :param object field: Current field
        :return: Field after modification by Generator
        :rtype: Object
        """
        self.field = field

        for b, bit in enumerate(self.bit_stream):
            if bit["m"] > 0:
                self.shape = Gaussian(**bit())
            else:
                self.shape = Sech(**bit())

            if domain.channels > 1:
                self.field[self.channel] += self.shape.generate(domain.t)
            else:
                self.field += self.shape.generate(domain.t)

            # Alternative: Only affect field of the current bit,
            # not the entire field:
            #~spb = domain.samples_per_bit
            #~bit_range = (b * spb, (b + 1) * spb)

            #~if domain.channels > 1:
                #~self.field[self.channel][bit_range[0]:bit_range[1]] += \
                    #~self.shape.generate(domain.t)[bit_range[0]:bit_range[1]]
            #~else:
                #~self.field[bit_range[0]:bit_range[1]] += \
                    #~self.shape.generate(domain.t)[bit_range[0]:bit_range[1]]

        return self.field
