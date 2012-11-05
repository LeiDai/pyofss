
Complete reference
==================

Domain
------
.. autoclass:: pyofss.domain.Domain
   :members:
   :special-members:
.. autofunction:: pyofss.domain.nu_to_omega
.. autofunction:: pyofss.domain.nu_to_lambda
.. autofunction:: pyofss.domain.omega_to_nu
.. autofunction:: pyofss.domain.omega_to_lambda
.. autofunction:: pyofss.domain.lambda_to_nu
.. autofunction:: pyofss.domain.lambda_to_omega
.. autofunction:: pyofss.domain.dnu_to_dlambda
.. autofunction:: pyofss.domain.dlambda_to_dnu

Field
-----
.. autofunction:: pyofss.field.temporal_power
.. autofunction:: pyofss.field.spectral_power
.. autofunction:: pyofss.phase
.. autofunction:: pyofss.chirp
.. autofunction:: pyofss.fft
.. autofunction:: pyofss.ifft
.. autofunction:: pyofss.ifftshift
.. autofunction:: pyofss.fftshift

Metrics
-------
.. autoclass:: pyofss.metrics.Metrics
   :members:
   :undoc-members:
.. autofunction:: pyofss.metrics.calculate_regenerator_factor

System
------
.. autoclass:: pyofss.system.System
   :members:
   :undoc-members:

Amplifier
---------
.. autoclass:: pyofss.modules.amplifier.Amplifier
   :members:
   :undoc-members:

Bit
---
.. autoclass:: pyofss.modules.bit.Bit
   :members:
   :undoc-members:
.. autofunction:: pyofss.modules.bit.generate_bitstream
.. autoclass:: pyofss.modules.bit.Bit_stream
   :members:
   :undoc-members:
.. autofunction:: pyofss.modules.bit.generate_prbs

CW
--
.. autoclass:: pyofss.modules.cw.Cw
   :members:
   :special-members:

Linearity
---------
.. autoclass:: pyofss.modules.linearity.Linearity
   :members:
   :undoc-members:
.. autofunction:: pyofss.modules.linearity.convert_alpha_to_linear
.. autofunction:: pyofss.modules.linearity.convert_alpha_to_dB
.. autofunction:: pyofss.modules.linearity.convert_dispersion_to_physical
.. autofunction:: pyofss.modules.linearity.convert_dispersion_to_engineering

Fibre
-----
.. autoclass:: pyofss.modules.fibre.Fibre
   :members:
   :undoc-members:

Filter
------
.. autoclass:: pyofss.modules.filter.Filter
   :members:
   :undoc-members:

Gaussian
--------
.. autoclass:: pyofss.modules.gaussian.Gaussian
   :members:
   :special-members:

Generator
---------
.. autoclass:: pyofss.modules.generator.Generator
   :members:
   :special-members:

Nonlinearity
------------
.. autoclass:: pyofss.modules.nonlinearity.Nonlinearity
   :members:
   :special-members:
.. autofunction:: pyofss.modules.nonlinearity.calculate_gamma
.. autofunction:: pyofss.modules.nonlinearity.calculate_raman_term

Plotter
-------
.. autofunction:: pyofss.modules.plotter.map_plot
.. autofunction:: pyofss.modules.plotter.waterfall_plot
.. autofunction:: pyofss.modules.plotter.single_plot
.. autofunction:: pyofss.modules.plotter.double_plot
.. autofunction:: pyofss.modules.plotter.multi_plot
.. autofunction:: pyofss.modules.plotter.quad_plot
.. autofunction:: pyofss.modules.plotter.animated_plot
.. autofunction:: pyofss.modules.plotter.convert_video

Sech
----
.. autoclass:: pyofss.modules.sech.Sech
   :members:
   :special-members:

Solver
------
.. autoclass:: pyofss.modules.solver.Solver
   :members:
   :special-members:

Stepper
-------
.. autoclass:: pyofss.modules.stepper.Stepper
   :members:
   :special-members:

Storage
-------
.. autoclass:: pyofss.modules.stepper.Storage
   :members:
   :special-members:
.. autofunction:: pyofss.modules.storage.reduce_to_range
