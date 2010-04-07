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

from distutils.core import setup
from distutils.extension import Extension


def name():
    return 'STEPS'
    
def version():
    return '1.1.1'
    
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
    
def model_ext():
    return dict(
        name='_model_swig',
        
        sources=['cpp/error.cpp',
        'cpp/model/model.cpp', 'cpp/model/diff.cpp',
        'cpp/model/reac.cpp','cpp/model/spec.cpp','cpp/model/sreac.cpp',
        'cpp/model/surfsys.cpp','cpp/model/volsys.cpp',
        'swig/model_wrap.cpp'],
        
        depends=['cpp/error.hpp',
        'cpp/model/model.hpp', 'cpp/model/diff.hpp',
        'cpp/model/reac.hpp','cpp/model/spec.hpp','cpp/model/sreac.hpp',
        'cpp/model/surfsys.hpp','cpp/model/volsys.hpp'],
        )
        
def geom_ext():
    return dict(
        name='_geom_swig',
        
        sources=['cpp/error.cpp','cpp/math/tetrahedron.cpp','cpp/geom/tetmesh.cpp',
        'cpp/math/linsolve.cpp','cpp/math/triangle.cpp',
        'cpp/geom/comp.cpp','cpp/geom/geom.cpp','cpp/geom/patch.cpp',
        'cpp/geom/tet.cpp', 'cpp/geom/tetmesh_rw.cpp',
        'cpp/geom/tmcomp.cpp','cpp/geom/tmpatch.cpp','cpp/geom/tri.cpp',
        'swig/geom_wrap.cpp'],
        
        depends=['cpp/error.hpp','cpp/math/tetrahedron.hpp','cpp/geom/tetmesh.hpp',
        'cpp/math/linsolve.hpp','cpp/math/triangle.hpp',
        'cpp/geom/comp.hpp','cpp/geom/geom.hpp','cpp/geom/patch.hpp',
        'cpp/geom/tet.hpp', 'cpp/geom/tetmesh_rw.hpp',
        'cpp/geom/tmcomp.hpp','cpp/geom/tmpatch.hpp','cpp/geom/tri.hpp'],
        )
        
def rng_ext():
    return dict(
        name='_rng',
        
        sources=['cpp/math/tools.cpp','cpp/error.cpp',
        'cpp/rng/rng.cpp', 'cpp/rng/mt19937.cpp',
        'swig/rng_wrap.cpp'],
        
        depends=['cpp/math/tools.hpp','cpp/error.hpp',
        'cpp/rng/rng.hpp', 'cpp/rng/mt19937.hpp'],
        )
        
def solver_ext_Mac():
    return dict(
        name='_solver_swig',
        
        sources=['cpp/error.cpp',

            'cpp/solver/api_comp.cpp','cpp/solver/api_main.cpp',
            'cpp/solver/api_patch.cpp','cpp/solver/api_tet.cpp',
            'cpp/solver/api_tri.cpp','cpp/solver/compdef.cpp',
            'cpp/solver/diffdef.cpp','cpp/solver/patchdef.cpp',
            'cpp/solver/reacdef.cpp','cpp/solver/specdef.cpp',
            'cpp/solver/sreacdef.cpp','cpp/solver/statedef.cpp',
        
            'cpp/tetexact/comp.cpp','cpp/tetexact/diff.cpp',
            'cpp/tetexact/kproc.cpp','cpp/tetexact/patch.cpp',
            'cpp/tetexact/reac.cpp','cpp/tetexact/sreac.cpp',
            'cpp/tetexact/tet.cpp','cpp/tetexact/tetexact.cpp',
            'cpp/tetexact/tri.cpp',
        
            'cpp/wmdirect/comp.cpp','cpp/wmdirect/kproc.cpp',
            'cpp/wmdirect/patch.cpp','cpp/wmdirect/reac.cpp',
            'cpp/wmdirect/sreac.cpp','cpp/wmdirect/wmdirect.cpp',
        
            'cpp/wmrk4/wmrk4.cpp',
        
            'swig/solver_wrap.cpp'],
        
        depends=['cpp/error.hpp',

            'cpp/solver/api.hpp','cpp/solver/compdef.hpp',
            'cpp/solver/diffdef.hpp','cpp/solver/patchdef.hpp',
            'cpp/solver/reacdef.hpp','cpp/solver/specdef.hpp',
            'cpp/solver/sreacdef.hpp','cpp/solver/statedef.hpp',
        
            'cpp/tetexact/comp.hpp','cpp/tetexact/diff.hpp',
            'cpp/tetexact/kproc.hpp','cpp/tetexact/patch.hpp',
            'cpp/tetexact/reac.hpp','cpp/tetexact/sreac.hpp',
            'cpp/tetexact/tet.hpp','cpp/tetexact/tetexact.hpp',
            'cpp/tetexact/tri.hpp',
        
            'cpp/wmdirect/comp.hpp','cpp/wmdirect/kproc.hpp',
            'cpp/wmdirect/patch.hpp','cpp/wmdirect/reac.hpp',
            'cpp/wmdirect/sreac.hpp','cpp/wmdirect/wmdirect.hpp',
        
            'cpp/wmrk4/wmrk4.hpp'],
        )
    
def solver_ext():
    return dict(
        name='_solver_swig',
        
        sources=['cpp/error.cpp',

            'cpp/solver/api_comp.cpp','cpp/solver/api_main.cpp',
            'cpp/solver/api_patch.cpp','cpp/solver/api_tet.cpp',
            'cpp/solver/api_tri.cpp','cpp/solver/compdef.cpp',
            'cpp/solver/diffdef.cpp','cpp/solver/patchdef.cpp',
            'cpp/solver/reacdef.cpp','cpp/solver/specdef.cpp',
            'cpp/solver/sreacdef.cpp','cpp/solver/statedef.cpp',
        
            'cpp/tetexact/comp.cpp','cpp/tetexact/diff.cpp',
            'cpp/tetexact/kproc.cpp','cpp/tetexact/patch.cpp',
            'cpp/tetexact/reac.cpp','cpp/tetexact/sreac.cpp',
            'cpp/tetexact/tet.cpp','cpp/tetexact/tetexact.cpp',
            'cpp/tetexact/tri.cpp',
        
            'cpp/wmdirect/comp.cpp','cpp/wmdirect/kproc.cpp',
            'cpp/wmdirect/patch.cpp','cpp/wmdirect/reac.cpp',
            'cpp/wmdirect/sreac.cpp','cpp/wmdirect/wmdirect.cpp',
        
            'cpp/wmrk4/wmrk4.cpp',
        
            'cpp/model/model.cpp', 'cpp/model/diff.cpp',
            'cpp/model/reac.cpp','cpp/model/spec.cpp','cpp/model/sreac.cpp',
            'cpp/model/surfsys.cpp','cpp/model/volsys.cpp'
            
            'cpp/math/tetrahedron.cpp','cpp/geom/tetmesh.cpp',
            'cpp/math/linsolve.cpp','cpp/math/triangle.cpp',
            'cpp/geom/comp.cpp','cpp/geom/geom.cpp','cpp/geom/patch.cpp',
            'cpp/geom/tet.cpp', 'cpp/geom/tetmesh_rw.cpp',
            'cpp/geom/tmcomp.cpp','cpp/geom/tmpatch.cpp','cpp/geom/tri.cpp',
        
            'cpp/math/tools.cpp','cpp/rng/rng.cpp', 'cpp/rng/mt19937.cpp',
        
            'swig/solver_wrap.cpp'],
        
        depends=['cpp/error.hpp',

            'cpp/solver/api.hpp','cpp/solver/compdef.hpp',
            'cpp/solver/diffdef.hpp','cpp/solver/patchdef.hpp',
            'cpp/solver/reacdef.hpp','cpp/solver/specdef.hpp',
            'cpp/solver/sreacdef.hpp','cpp/solver/statedef.hpp',
        
            'cpp/tetexact/comp.hpp','cpp/tetexact/diff.hpp',
            'cpp/tetexact/kproc.hpp','cpp/tetexact/patch.hpp',
            'cpp/tetexact/reac.hpp','cpp/tetexact/sreac.hpp',
            'cpp/tetexact/tet.hpp','cpp/tetexact/tetexact.hpp',
            'cpp/tetexact/tri.hpp',
        
            'cpp/wmdirect/comp.hpp','cpp/wmdirect/kproc.hpp',
            'cpp/wmdirect/patch.hpp','cpp/wmdirect/reac.hpp',
            'cpp/wmdirect/sreac.hpp','cpp/wmdirect/wmdirect.hpp',
        
            'cpp/wmrk4/wmrk4.hpp'],
        )
        
def ext_modules():
    import sys
    if sys.platform == 'darwin':
        modules = [model_ext(), geom_ext(), rng_ext(), solver_ext_Mac()]
    else:
        modules = [model_ext(), geom_ext(), rng_ext(), solver_ext()]
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