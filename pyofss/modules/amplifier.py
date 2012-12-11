
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

from scipy import power, sqrt
from pyofss.field import fft, ifft


class Amplifier(object):
    """
    :param string name: Name of this module
    :param double gain: Amount of (logarithmic) gain. *Unit: dB*
    :param double power: Average power level to target

    Simple amplifier provides gain but no noise
    """
    def __init__(self, name="amplifier", gain=None, power=None):

        self.name = name

        if (gain is None) and (power is None):
            print "Amplifier has not been given gain or power parameter"
        if (gain is not None) and (power is not None):
            print "Amplifier given both gain and power parameters"

        self.gain = gain
        self.power = power

        self.field = None

    def __call__(self, domain, field):

        # Convert field to spectral domain:
        self.field = fft(field)

        if self.gain is not None:
            # Calculate linear gain from logarithmic gain (G_dB -> G_linear)
            G = power(10, 0.1 * self.gain)
            sqrt_G = sqrt(G)

            if domain.channels > 1:
                self.field[0] *= sqrt_G
                self.field[1] *= sqrt_G
            else:
                self.field *= sqrt_G

        if self.power is not None:
            pass

        # convert field back to temporal domain:
        return ifft(self.field)
