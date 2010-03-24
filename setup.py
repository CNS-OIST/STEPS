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
import sys

if sys.platform != 'win32':
    solver_src = ['cpp/error.cpp',

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
        
            'swig/solver_wrap.cpp']
else:
    solver_src = ['cpp/error.cpp',
            'cpp/rng/rng.cpp', 'cpp/rng/mt19937.cpp',
            'cpp/math/tools.cpp',
            
            'cpp/model/model.cpp', 'cpp/model/diff.cpp',
            'cpp/model/reac.cpp','cpp/model/spec.cpp','cpp/model/sreac.cpp',
            'cpp/model/surfsys.cpp','cpp/model/volsys.cpp',
            
            'cpp/math/tetrahedron.cpp','cpp/geom/tetmesh.cpp',
            'cpp/math/linsolve.cpp','cpp/math/triangle.cpp',
            'cpp/geom/comp.cpp','cpp/geom/geom.cpp','cpp/geom/patch.cpp',
            'cpp/geom/tet.cpp', 'cpp/geom/tetmesh_rw.cpp',
            'cpp/geom/tmcomp.cpp','cpp/geom/tmpatch.cpp','cpp/geom/tri.cpp',
        
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
        
            'swig/solver_wrap.cpp']

setup(name = 'STEPS',
      version = '1.1.0',
      author = 'STEPS Development Team',
      author_email = 'steps.dev@gmail.com',
      url = 'http://steps.sourceforge.net',
      description = 'STochastic Engine for Pathway Simulation',
      download_url = 'http://sourceforge.net/projects/steps/files/',
      platforms = ['Mac OS X', 'Windows XP', 'Windows Vista', 'Linux', 'Unix'],
      license = 'GNU General Public License Version 3.0',
      
      py_modules=['steps/__init__','steps/rng','steps/model','steps/geom', 'steps/solver',
        'steps/model_swig','steps/geom_swig', 'steps/solver_swig',
        'steps/utilities/meshio', 'steps/utilities/visual', 'steps/utilities/__init__'],
        
      ext_package='steps',
      
      ext_modules=[Extension('_rng',
        sources = ['cpp/math/tools.cpp','cpp/error.cpp',
        'cpp/rng/rng.cpp', 'cpp/rng/mt19937.cpp',
        'swig/rng_wrap.cpp']),
        
        Extension('_model_swig',
        sources = ['cpp/error.cpp',
        'cpp/model/model.cpp', 'cpp/model/diff.cpp',
        'cpp/model/reac.cpp','cpp/model/spec.cpp','cpp/model/sreac.cpp',
        'cpp/model/surfsys.cpp','cpp/model/volsys.cpp',
        
        'swig/model_wrap.cpp']),
        
        Extension('_geom_swig',
        sources = ['cpp/error.cpp','cpp/math/tetrahedron.cpp','cpp/geom/tetmesh.cpp',
        'cpp/math/linsolve.cpp','cpp/math/triangle.cpp',
        'cpp/geom/comp.cpp','cpp/geom/geom.cpp','cpp/geom/patch.cpp',
        'cpp/geom/tet.cpp', 'cpp/geom/tetmesh_rw.cpp',
        'cpp/geom/tmcomp.cpp','cpp/geom/tmpatch.cpp','cpp/geom/tri.cpp',
        
        'swig/geom_wrap.cpp']),
        
        Extension('_solver_swig',
        sources = solver_src),
      ]
      )