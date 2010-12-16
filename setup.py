# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#  Last Changed Rev:  $Rev$
#  Last Changed Date: $Date$
#  Last Changed By:   $Author$

try:
    import ez_setup
    ez_setup.use_setuptools()
    from setuptools import setup, Extension
    
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension


def name():
    return 'STEPS'
    
def version():
    return '1.1.4'
    
def author():
    return 'STEPS Development Team'
    
def email():
    return 'steps.dev@gmail.com'
    
def url():
    return 'http://steps.sourceforge.net'
    
def desc():
    return 'STochastic Engine for Pathway Simulation'

def download():
    return 'http://sourceforge.net/projects/steps/files/src'

def platforms():
    return ['Mac OS X', 'Windows XP', 'Windows Vista', 'Linux', 'Unix']
    
def license():
    return 'GNU General Public License Version 3.0'
    
def packages():
    return ['steps', 'steps/utilities']
  
def steps_ext():
    return dict(
        name='_steps_swig',
        
        sources=['cpp/error.cpp',
        'cpp/model/model.cpp', 'cpp/model/diff.cpp',
        'cpp/model/reac.cpp','cpp/model/spec.cpp','cpp/model/sreac.cpp',
        'cpp/model/surfsys.cpp','cpp/model/volsys.cpp',
        
        'cpp/math/tetrahedron.cpp','cpp/geom/tetmesh.cpp',
        'cpp/math/linsolve.cpp','cpp/math/triangle.cpp',
        'cpp/geom/comp.cpp','cpp/geom/geom.cpp','cpp/geom/patch.cpp',
        'cpp/geom/tet.cpp', 'cpp/geom/tetmesh_rw.cpp',
        'cpp/geom/tmcomp.cpp','cpp/geom/tmpatch.cpp', 'cpp/geom/diffboundary.cpp',
        'cpp/geom/tri.cpp',

        
        'cpp/math/tools.cpp',
        'cpp/rng/rng.cpp', 'cpp/rng/mt19937.cpp',

        'cpp/solver/api_comp.cpp','cpp/solver/api_main.cpp',
        'cpp/solver/api_patch.cpp','cpp/solver/api_tet.cpp',
        'cpp/solver/api_tri.cpp', 'cpp/solver/api_diffboundary.cpp',
        'cpp/solver/compdef.cpp',
        'cpp/solver/diffdef.cpp','cpp/solver/patchdef.cpp',
        'cpp/solver/reacdef.cpp','cpp/solver/specdef.cpp',
        'cpp/solver/sreacdef.cpp', 'cpp/solver/diffboundarydef.cpp',
        'cpp/solver/statedef.cpp',
        
        'cpp/tetexact/comp.cpp','cpp/tetexact/diff.cpp',
        'cpp/tetexact/kproc.cpp','cpp/tetexact/patch.cpp',
        'cpp/tetexact/reac.cpp','cpp/tetexact/sreac.cpp',
        'cpp/tetexact/tet.cpp','cpp/tetexact/tetexact.cpp',
        'cpp/tetexact/tri.cpp', 'cpp/tetexact/diffboundary.cpp',
        
        'cpp/wmdirect/comp.cpp','cpp/wmdirect/kproc.cpp',
        'cpp/wmdirect/patch.cpp','cpp/wmdirect/reac.cpp',
        'cpp/wmdirect/sreac.cpp','cpp/wmdirect/wmdirect.cpp',
        
        'cpp/wmrk4/wmrk4.cpp',
        
        'swig/steps_wrap.cpp'],
        
        undef_macros=['NDEBUG']
    )
        
        
def ext_modules():
    modules = [steps_ext()]
    return modules

ExtModule = lambda extension:  Extension(**extension)

setup(name = name(),
      version = version(),
      author = author(),
      author_email = email(),
      url = url(),
      description = desc(),
      download_url = download(),
      platforms = platforms(),
      license = license(),
      
      packages = packages(), 
      
      ext_package = 'steps',
      
      ext_modules  = [ExtModule(ext) for ext in ext_modules()]
      )


