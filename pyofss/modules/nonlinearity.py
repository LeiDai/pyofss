
# ==============================================================================
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
# ==============================================================================

import numpy as np
from numpy import pi

from pyofss.field import fft, ifft, fftshift
from pyofss.domain import Domain
# ==============================================================================
def calculate_gamma( nonlinear_index, effective_area,
                     centre_omega = 2.0 * pi * 193.1 ):
   """
   :param double nonlinear_index: :math:`n_2`. *Unit:* :math:`m^2 / W`
   :param double effective_area: :math:`A_{eff}`. *Unit:* :math:`\mu m^2`
   :param double centre_omega: Angular frequency of carrier. *Unit: rad / ps*
   :return: Nonlinear parameter. :math:`\gamma`. *Unit:* :math:`rad / (W km)`
   :rtype: double

   Calculate nonlinear parameter from the nonlinear index and effective area.
   """
   # ===========================================================================
   gamma = 1.0e24 * nonlinear_index * centre_omega / \
           (effective_area * Domain.vacuum_light_speed)

   return gamma
# ==============================================================================
def calculate_raman_term( domain, tau_1 = 12.2e-3, tau_2 = 32.0e-3 ):
   """
   :param object domain: A domain
   :param double tau_1: First adjustable parameter. *Unit: ps*
   :param double tau_2: Second adjustable parameter. *Unit: ps*
   :return: Raman response function
   :rtype: double

   Calculate raman response function from tau_1 and tau_2.
   """
   # ===========================================================================
   h_R = (tau_2**2 + tau_1**2) / (tau_1 * tau_2**2)
   h_R *= np.exp( -domain.t / tau_2 ) * np.sin( domain.t / tau_1 )

   return h_R
# ==============================================================================
class Nonlinearity():
   """
   Nonlinearity is used by fibre to generate a nonlinear factor.
   """
   # ===========================================================================
   def __init__( self, gamma = None, sim_type = None, self_steepening = False, 
                 raman_scattering = False, rs_factor = 3e-3, use_all = False,
                 tau_1 = 12.2e-3, tau_2 = 32.0e-3, f_R = 0.18 ):

      self.gamma = gamma
      self.self_steepening = self_steepening
      self.raman_scattering = raman_scattering
      self.rs_factor = rs_factor
      self.use_all = use_all

      self.tau_1 = tau_1
      self.tau_2 = tau_2
      self.f_R = f_R
      # ========================================================================
      self.generate_nonlinearity = getattr( self, "%s_nonlinearity" % sim_type,
                                            self.default_nonlinearity )

      if( use_all ):
         print "Using general expression for nonlinearity"
         self.non = getattr( self, "%s_f_all" % sim_type, self.default_f_all )
         self.exp_non = getattr( self, "%s_exp_f_all" % sim_type,
                                 self.default_exp_f_all )
      else:
         if( self.self_steepening and self.raman_scattering ):
            print "Using self_steepening and raman_scattering"
            self.non = getattr( self, "%s_f_with_ss_and_rs" % sim_type,
                                self.default_f_with_ss_and_rs )
            self.exp_non = getattr( self, "%s_exp_f_with_ss_and_rs" % sim_type,
                                    self.default_exp_f_with_ss_and_rs )
         elif( self.self_steepening ):
            print "Using self_steepening"
            self.non = getattr( self, "%s_f_with_ss" % sim_type,
                                self.default_f_with_ss )
            self.exp_non = getattr( self, "%s_exp_f_with_ss" % sim_type,
                                    self.default_exp_f_with_ss )
         elif( self.raman_scattering ):
            print "Using raman_scattering"
            self.non = getattr( self, "%s_f_with_rs" % sim_type,
                                self.default_f_with_rs )
            self.exp_non = getattr( self, "%s_exp_f_with_rs" % sim_type,
                                    self.default_exp_f_with_rs )
         else:
            self.non = getattr( self, "%s_f" % sim_type,
                                self.default_f )
            self.exp_non = getattr( self, "%s_exp_f" % sim_type,
                                    self.default_exp_f )
   # ===========================================================================
   def __call__( self, domain ):

      self.centre_omega = domain.centre_omega
      self.omega = fftshift( domain.omega - domain.centre_omega )

      if( self.self_steepening ):
         self.ss_factor = 1.0 / self.centre_omega
      else:
         self.ss_factor = 0.0

      if( self.use_all ):
         # Require h_R in spectral domain, so take FFT of returned value:
         self.h_R = fft( calculate_raman_term(domain, self.tau_1, self.tau_2) )
      else:
         self.h_R = 0.0

      self.generate_nonlinearity()
   # ===========================================================================
   def default_nonlinearity( self ):
      if( self.gamma is None ):
         self.factor = 0.0
      else:
         self.factor = 1j * self.gamma
   # ===========================================================================
   def wdm_nonlinearity( self ):
      if( self.gamma is None ):
         self.factor = (0.0, 0.0)
      else:
         self.factor = (1j * self.gamma[0], 1j * self.gamma[1])
   # ===========================================================================
   # ===========================================================================
   def default_f_all( self, A, z ):
      term_spm = np.abs(A)**2
      convolution = ifft( self.h_R * fft(term_spm) )
      p = fft( A * (1.0 - self.f_R) * term_spm + self.f_R * A * convolution )

      return ifft( self.factor * (1.0  + self.omega * self.ss_factor) * p )
   # ===========================================================================
   def default_exp_f_all( self, A, h, B ):
      term_spm = np.abs(A)**2
      convolution = ifft( self.h_R * fft(term_spm) )
      p = fft( (1.0 - self.f_R) * term_spm + self.f_R * convolution )

      return np.exp( h * ifft(self.factor * \
         (1.0 + self.omega * self.ss_factor) * p) ) * B
   # ===========================================================================
   # ===========================================================================
   def default_f_with_ss( self, A, z ):
      term_spm = np.abs(A)**2 * A
      term_ss = self.ss_factor * ifft( self.omega * fft(term_spm) )
      
      return self.factor * (term_spm + term_ss)
   # ===========================================================================
   def default_exp_f_with_ss( self, A, h, B):
      term_spm = np.abs(A)**2
      term_ss = (self.ss_factor / B) * ifft( self.omega * fft(term_spm * A) )
      return np.exp( h * self.factor * (term_spm + term_ss) ) * B
   # ===========================================================================
   # ===========================================================================
   def default_f_with_rs( self, A, z ):
      term_spm = np.abs(A)**2
      term_rs = self.rs_factor * ifft( 1j * self.omega * fft(term_spm) )

      return self.factor * (term_spm + term_rs) * A
   # ===========================================================================
   def default_exp_f_with_rs( self, A, h, B):
      term_spm = np.abs(A)**2
      term_rs = self.rs_factor * ifft( 1j * self.omega * fft(term_spm) )

      return np.exp( h * self.factor * (term_spm + term_rs) ) * B
   # ===========================================================================
   # ===========================================================================
   def default_f_with_ss_and_rs( self, A, z ):
      term_spm = np.abs(A)**2
      term_ss = self.ss_factor * ifft( self.omega * fft(term_spm * A) )
      term_rs = self.rs_factor * ifft( 1j * self.omega * fft(term_spm) )

      return self.factor * (term_spm + term_rs) * A + term_ss
   # ===========================================================================
   def default_exp_f_with_ss_and_rs( self, A, h, B ):
      term_spm = np.abs(A)**2
      term_ss = (self.ss_factor / B) * ifft( self.omega * fft(term_spm * A) )
      term_rs = self.rs_factor * ifft( 1j * self.omega * fft(term_spm) )

      return np.exp( h * self.factor * (term_spm + term_ss + term_rs) ) * B
   # ===========================================================================
   def default_f( self, A, z ):
      return self.factor * np.abs(A)**2 * A
   # ===========================================================================
   def default_exp_f( self, A, h, B ):
      return np.exp( h * self.factor * np.abs(A)**2 ) * B
   # ===========================================================================
   # ===========================================================================
   # TODO: implement higher order nonlinear effects for wdm.
   def wdm_f_with_ss( self, As, z ):
      self.wdm_f( As, z )
   def wdm_exp_f_with_ss( self, As, h, Bs ):
      self.wdm_exp_f( As, h, Bs )
   # ===========================================================================
   def wdm_f_with_rs( self, As, z ):
      self.wdm_f( As, z )
   def wdm_exp_f_with_rs( self, As, h, Bs ):
      self.wdm_exp_f( As, h, Bs )
   # ===========================================================================
   def wdm_f_with_ss_and_rs( self, As, z ):
      self.wdm_f( As, z )
   def wdm_exp_f_with_ss_and_rs( self, As, h, Bs ):
      self.wdm_exp_f( As, h, Bs )
   # ===========================================================================
   def wdm_f( self, As, z ):
      return np.asarray( [ \
         self.factor[0] * (np.abs(As[0])**2 + 2.0 * np.abs(As[1])**2) * As[0], \
         self.factor[1] * (np.abs(As[1])**2 + 2.0 * np.abs(As[0])**2) * As[1] ])
   # ===========================================================================
   def wdm_exp_f( self, As, h, Bs ):
      return np.asarray([ np.exp( h * self.factor[0] * \
         (np.abs(As[0])**2 + 2.0 * np.abs(As[1])**2) ) * Bs[0], \
         np.exp( h * self.factor[1] * \
         (np.abs(As[1])**2 + 2.0 * np.abs(As[0])**2) ) * Bs[1] ])
# ==============================================================================
