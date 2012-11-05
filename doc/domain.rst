.. _tutorial-index:

Domain
======

.. note::
   It is assumed that pyofss has been imported using

   >>> from pyofss import *

The start of any system in pyofss is the domain. Generate a domain using

   >>> domain = Domain()

which will construct a default domain. Information on the domain may be viewed

   >>> print domain
   Domain:
	   total_bits = 1
	   samples_per_bit = 512
	   bit_width = 100.0000 ps
	   centre_nu = 193.1000 THz
	   total_samples = 512
	   window_t = 100.0000 ps
	   centre_omega = 1213.2831 rad / ps
	   centre_lambda = 1552.5244 nm
	   dt = 0.1953 ps
	   window_nu = 5.1200 THz
	   dnu = 0.0100 THz
	   window_omega = 32.1699 rad / ps
	   domega = 0.0628 rad / ps
	   window_lambda = 41.1648 nm
	   dlambda = 0.0804 nm
	   channels = 1

The default domain uses a centre frequency of 193.1 THz, though often a centre wavelength of 1550 nm is used.
There are helper functions to convert between frequency and wavelength.
To find the frequency corresponding to the desired wavelength:

   >>> nu_c = lambda_to_nu( 1550.0 )
   >>> print nu_c
   193.414489032

Now a domain may be constructed with this centre wavelength

   >>> domain = Domain( centre_nu = nu_c )
   >>> print domain
   Domain:
	   total_bits = 1
	   samples_per_bit = 512
	   bit_width = 100.0000 ps
	   centre_nu = 193.4145 THz
	   total_samples = 512
	   window_t = 100.0000 ps
	   centre_omega = 1215.2591 rad / ps
	   centre_lambda = 1550.0000 nm
	   dt = 0.1953 ps
	   window_nu = 5.1200 THz
	   dnu = 0.0100 THz
	   window_omega = 32.1699 rad / ps
	   domega = 0.0628 rad / ps
	   window_lambda = 41.0311 nm
	   dlambda = 0.0801 nm
	   channels = 1

The temporal and spectral arrays may be accessed using

   >>> domain.t
   >>> domain.nu
   >>> domain.Lambda
