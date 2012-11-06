
Gaussian
========

.. note::
   It is assumed that pyofss has been imported using

   >>> from pyofss import *

The Gaussian module allows a Gaussian shaped pulse to be generated.

A default Gaussian may be generated using:

   >>> gaussian = Gaussian()

An example of a Gaussian with a peak power of 0.5 W and a width of 20 ps:

   >>> gaussian = Gaussian( peak_power = 0.5, width = 20.0 )

The position parameter is relative to the temporal domain width (which may be found from sys.domain.window_t).

A simulation using a Gaussian pulse

.. plot::
   :include-source:

   from pyofss import *
   sys = System()
   sys.add( Gaussian(peak_power = 0.5, width = 20.0) )
   sys.run()

   single_plot( sys.domain.t, temporal_power(sys.field),
                labels['t'], labels['P_t'] )
