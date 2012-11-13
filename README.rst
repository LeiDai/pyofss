
Pyofss: Python-based optical fibre system simulator
===================================================

Pyofss allows construction of an optical fibre system from separate modules.
A typical system consists of a Gaussian pulse generator module and an optical 
fibre module.
The field generated is propagated through the fibre by numerical integration of 
an appropriate Schrödinger-type equation.

Simulated effects include that of dispersion (second, third, and higher order), attenuation, self-phase modulation, self-steepening, and Raman scattering.
Resulting field profiles (including multiple traces for the fibre module) can be visualised using a range of plot types.
These include standard x-y plots, top-down "map" plots, three-dimensional "waterfall" plots, and animation videos.

Installation
------------

Pyofss is available on PyPI and can be retrieved using the pip program:

.. code-block:: bash

   sudo aptitude install python-pip
   pip install pyofss

Then import pyofss within scripts or in an interactive session:

.. code-block:: pycon

   >>> import pyofss

Dependencies
------------

Pyofss depends on Numpy, Scipy, and Matplotlib.
They can be installed on Linux distributions using ``aptitude``:

.. code-block:: bash

   sudo aptitude install python-numpy python-scipy python-matplotlib

.. note::

   The recommended dependency versions are listed in the ``requirements.txt`` file within the pyofss package.
   Install each of these dependencies using ``pip``:

   .. code-block:: bash

      pip install -r requirements.txt

Development
-----------

It is recommended to install pyofss into a virtual environment, which can be initialised using:

.. code-block:: bash

   sudo aptitude install python-virtualenv
   sudo pip install virtualenvwrapper
   mkvirtualenv pyofss
   workon pyofss

Pyofss dependencies can then be satisfied using:

.. code-block:: bash

   pip install numpy
   pip install scipy
   sudo aptitude build-dep python-matplotlib
   pip install matplotlib

Install the latest development version of pyofss from GitHub:

.. code-block:: bash

   pip install git+https://github.com/daibo/pyofss.git

Tests
-----

Tests can be run within the pyofss package:

.. code-block:: bash

   python setup.py test
