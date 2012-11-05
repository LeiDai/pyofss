
Overview
========

Domain
------

A domain in pyofss contains the temporal and spectral axis (arrays) which 
are used by the system modules during simulation.
There are helper functions which convert between frequency, angular frequency, 
and wavelength.

   * TB -- total bits
   * SPB -- samples per bit
   * BW -- bit width

   - TS -- total samples (TS = TB x SPB)
   - TW -- time window (TW = TB x BW)

Using :math:`\Delta t_{max} \equiv TW` and :math:`N \equiv TS`, the temporal 
array corresponding to :math:`t \in [0, \Delta t_{max})` is generated from

.. math::
   \Delta t &= \Delta t_{max} / N \\
   t_n &= n \Delta t \qquad n \in [0, N)

The first and last elements of the temporal array will have the values:

.. math::
   t[0] &= 0.0 \\
   t[N-1] &= \Delta t_{max} - \Delta t

For the spectral array there will be a fixed shift such that the middle of the 
array corresponds to a central frequency :math:`\nu_c`.
Then the spectral array corresponding to 
:math:`\nu \in \left[ \nu_c - \frac{ \Delta \nu_{max} }{ 2 }, \nu_c + \frac{ \Delta \nu_{max} }{ 2 } \right)`
is generated from

.. math::
   \Delta \nu &= \frac{ 1 }{ \Delta t_{max} } \\
   \Delta \nu_{max} &= N \Delta \nu \\
   \nu_n &= \nu_c - \frac{ \Delta \nu_{max} }{ 2 } + n \Delta \nu \qquad n \in [0, N)

The first and last elements of the spectral array will have the values:

.. math::
   \nu[0] &= \nu_c - \frac{ \Delta \nu_{max} }{ 2 } \\
   \nu[N-1] &= \nu_c + \frac{ \Delta \nu_{max} }{ 2 } - \Delta \nu

An angular frequency array and a wavelength array are also generated using the 
following two relations

.. math::
   \omega &= 2 \pi \nu \\
   \lambda &= c / \nu

Field
-----

The field in pyofss represents the complex valued slowly varying envelope of 
the electric field that will be propagated through the optical fibre system. 
A fixed polarisation is used for the field.
It is possible to use separate field arrays for each frequency/wavelength 
(channel) to be simulated.


Modules
-------

A module in pyofss is a class that may be called with a domain and a field as 
parameters. The following modules may be generated:


Gaussian
^^^^^^^^

.. math::
   A(0,t) = \sqrt{P_0} \, \exp \left[ \frac{ -(1+iC) }{ 2 } \left( \
                                    \frac{ t-t_0 }{ \Delta t_{1/e} } \
                                    \right)^{2m} + i \left( \phi_0 - \
                                    2 \pi (\nu_c + \nu_\text{offset})t \
                                    \right) \right]

The gaussian module will generate a (super-)Gaussian shaped pulse using the 
following parameters:

   :math:`P_0` -- Peak power (W)

   :math:`C` -- Chirp parameter (rad)

   :math:`t_0` -- pulse position (ps)

   :math:`\Delta t_{1/e}` -- pulse HWIeM width (ps) [**note: not FWHM**]

   :math:`m` -- Order paramater

   :math:`\phi_0` -- initial phase (rad)

   :math:`\nu_c` -- centre frequency (THz)

   :math:`\nu_\text{offset}` -- offset frequency (THz)

Since the centre frequency :math:`\nu_c` is fixed by the domain, only the 
offset frequency :math:`\nu_\text{offset}` is provided by the user.

The pulse position :math:`t_0` is given as a factor of the time window 
:math:`\Delta t_{max}` such that :math:`t_0 = f \Delta t_{max}`, and the user 
supplies the factor :math:`f`. 


Sech
^^^^

.. math::
   A(0,t) = \sqrt{P_0} \, \text{sech} \left( \
                                \frac{ t-t_0 }{ \Delta t_{1/e} } \right) \
                                \, \exp \left[ \frac{ -iC }{ 2 } \left( \
                                \frac{ t-t_0 }{ \Delta t_{1/e} } \
                                \right)^2 + i \left( \phi_0 - \
                                2 \pi (\nu_c + \nu_\text{offset}) t \
                                \right) \right]

The sech module will generate a hyperbolic-Secant shaped pulse using the 
following parameters:

   :math:`P_0` -- Peak power (W)

   :math:`C` -- Chirp parameter (rad)

   :math:`t_0` -- pulse position (ps)

   :math:`\Delta t_{1/e}` -- pulse HWIeM width (ps) [**note: not FWHM**]

   :math:`\phi_0` -- initial phase (rad)

   :math:`\nu_c` -- centre frequency (THz)

   :math:`\nu_\text{offset}` -- offset frequency (THz)

Since the centre frequency :math:`\nu_c` is fixed by the domain, only the 
offset frequency :math:`\nu_\text{offset}` is provided by the user.

The pulse position :math:`t_0` is given as a factor of the time window 
:math:`\Delta t_{max}` such that :math:`t_0 = f \Delta t_{max}`, and the user 
supplies the factor :math:`f`. 


Generator
^^^^^^^^^

A generator module can generate Gaussian or hyperbolic-secant (Sech) shaped 
pulses.
The main difference to directly calling a Gaussian or Sech module is that the 
pulse position paramter :math:`t_0` is given as a factor of the bit width, and 
not as a factor of the time window.
The user provides the factor :math:`f` where :math:`t_0 = f \Delta t_{bit}`.

Fibre
^^^^^

The fibre module propagates the input field incrementally using the generalised 
non-linear Schrödinger equation:

.. math::
   \frac{ \partial{A} }{ \partial{z} } = \left[ \hat{L} + \
   \hat{N}\left(A\right) \right] A 

where :math:`A` is the complex field envelope of the pulse and :math:`z` is the 
dimension along the fibre length. The linear operator :math:`\hat{L}` and 
non-linear operator :math:`\hat{N}` are usually written:

.. math::
   \hat{L} = -\frac{\alpha}{2} - \frac{i\beta_2}{2} \
   \frac{\partial^2}{\partial t^2} + \frac{\beta_3}{6} \
   \frac{\partial^3}{\partial t^3} + \ldots

.. math::
   \hat{N} = i \gamma \left( |A|^2 + \frac{i}{\omega_0} \frac{1}{A} \
   \frac{\partial |A|^2 A}{\partial t} - \
   t_R \frac{\partial |A|^2}{\partial t} \right)

The linear operator contains terms for attenuation and (second order and 
higher) dispersion. The nonlinear operator contains terms for self-phase 
modulation (SPM), self-steepening, and Raman scattering.

It is useful to apply the linear operator to the field in the frequency domain 
using the property :math:`\frac{\partial A}{\partial t} \leftrightarrow -i 
\omega \tilde{A}`:

.. math::
   \hat{L} = - \frac{ \alpha }{ 2 } + i \left( \frac{ \beta_2 }{ 2 } \
   \omega^2 + \frac{ \beta_3 }{ 6 } \omega^3 + \ldots \right)

Filter
^^^^^^

.. math::
   H(\nu) = \exp \left[ \left(\frac{- 2 \pi \
   \nu_\text{filter}}{2 \Delta \nu} \right)^{2m} \right]

The filter module will slice a specified band of the current field in the 
spectral domain.

Amplifier
^^^^^^^^^

.. math::
   \tilde{A}_\text{out} = \sqrt{G} \tilde{A}_\text{in}

The amplifier module applies a gain to the input field.
