
# ==============================================================================
"""
    Copyright (C) 2011, 2012  David Bolt

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

from domain import Domain
# ==============================================================================
class System():
   """
   :param object domain: A domain to be used with contained modules

   A system consists of a list of modules, each of which may be called with a 
   domain and field as parameters. The result of each module call is stored in 
   a dictionary.
   """
   # ===========================================================================
   def __init__( self, domain = Domain() ):
      self.domain = domain
      self.clear( remove_modules = True )
   # ===========================================================================
   def clear( self, remove_modules = False ):
      """
      Clear contents of all fields.
      Clear (remove) all modules if requested.
      """
      # ========================================================================
      if( self.domain.channels > 1 ):
         self.field = [ np.zeros([self.domain.total_samples], complex)
                        for n in range(self.domain.channels) ]
      else:
         self.field = np.zeros( [self.domain.total_samples], complex )

      self.fields = {}

      if( remove_modules ):
         self.modules = []
   # ===========================================================================
   def add( self, module ):
      self.modules.append( module )
   # ===========================================================================
   def __getitem__( self, module_name ):
      for n, module in enumerate( self.modules ):
         if( module.name == module_name ):
            return self.modules[n]
   # ===========================================================================
   def __setitem__( self, module_name, new_module ):
      for n, module in enumerate( self.modules ):
         if( module.name == module_name ):
            self.modules[n] = new_module
            return

      raise Exception( "Tried to modify non-existing module in system" )
   # ===========================================================================
   def run( self ):
   # For each module add result of call to dictionary using module name as key:
      for module in self.modules:
         self.field = module( self.domain, self.field )
         self.fields[module.name] = self.field
# ==============================================================================
