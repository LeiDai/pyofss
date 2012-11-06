
Pyofss tutorial
===============

:Release: |version|
:Date: |today|

Pyofss has been developed to be fully interactive.
It allows the user to construct an optical fibre system, to get detailed information about the system, and to study the effects on an optical signal propagating through such a system.

Here is an example of a simple system consisting of a Gaussian pulse generator:

.. plot::
   :include-source:

   from pyofss import *
   sys = System()
   sys.add( Gaussian() )
   sys.run()
   single_plot( sys.domain.t, temporal_power(sys.field), 
                labels['t'], labels['P_t'] )

.. toctree::
   :numbered:

   system.rst
   domain.rst
   gaussian.rst
   fibre.rst
