
Pyofss: Python-based optical fibre system simulator
===================================================

Pyofss allows construction of an optical fibre system from separate modules.
A typical system consists of a Gaussian pulse generator module and an optical 
fibre module.
The field generated is propagated through the fibre by numerical integration of 
an appropriate Schrödinger-type equation.
Simulated effects include that of dispersion (second, third, and higher order), attenuation, self-phase modulation, self-steepening, and Raman scattering.
Resulting field profiles (including multiple traces for the fibre module) may be visualised using a range of plot types.
These include standard x-y plots, top-down "map" plots, three-dimensional "waterfall" plots, and animation videos.

Installation
------------

Pyofss is available on Pypi and may be retrieved using the pip program::

   $ aptitude install pip
   $ pip install pyofss

Then import pyofss within scripts or in an interactive session:
   >>> from pyofss import *

.. note::
   If the required dependencies are not satisfied when installing pyofss, then manually install using either::

      $ aptitude install python-numpy python-scipy python-matplotlib

   or::

      $ pip install numpy scipy matplotlib

   The recommended versions are listed in the "requirements.txt" file within the pyofss package.
   Using this file, it is possible to automatically install all dependencies::

      $ pip install -r requirements.txt
