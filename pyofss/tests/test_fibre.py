
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

from pyofss.modules.fibre import Fibre

#from numpy.testing.utils import assert_almost_equal # uses decimal places

import unittest


class DefaultParameters(unittest.TestCase):
    """ Test default parameters. """
    def test_none(self):
        """ Should use default value if no parameter given """
        fibre = Fibre()
        self.assertEqual(fibre.name, "fibre")
        self.assertEqual(fibre.length, 1.0)


class BadParameters(unittest.TestCase):
    """ Test response to bad parameters. """
    def test_too_low(self):
        """ Should fail when parameters are too low """
        pass

    def test_too_high(self):
        """ Should fail when parameters are too high """
        pass

    def test_wrong_type(self):
        """ Should fail if wrong type """
        pass


class CheckInputAndOutput(unittest.TestCase):
    """ Test input and output of Fibre. """
    def test_output(self):
        """ Check Fibre outputs its values correctly """
        pass


class CheckFunctions(unittest.TestCase):
    """ Test class methods. """
    pass

if __name__ == "__main__":
    unittest.main()
