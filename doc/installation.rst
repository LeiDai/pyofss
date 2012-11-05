
Installation
============

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
