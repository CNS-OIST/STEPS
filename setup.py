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

setup(name='steps',
      version='1.1.0',
      py_modules=['steps/__init__','steps/rng','steps/model','steps/geom', 'steps/solver',
        'steps/model_swig','steps/geom_swig', 'steps/solver_swig',
        'steps/utilities/meshio', 'steps/utilities/visual', 'steps/utilities/__init__'],
      ext_package='steps',
      ext_modules=[Extension('_rng',
        sources = ['cpp/steps/math/tools.cpp','cpp/steps/error.cpp',
        'cpp/steps/rng/rng.cpp', 'cpp/steps/rng/mt19937.cpp',
        'steps/rng_wrap.cpp']),
        Extension('_model_swig',
        sources = ['cpp/steps/error.cpp',
        'cpp/steps/model/model.cpp', 'cpp/steps/model/diff.cpp',
        'cpp/steps/model/reac.cpp','cpp/steps/model/spec.cpp','cpp/steps/model/sreac.cpp',
        'cpp/steps/model/surfsys.cpp','cpp/steps/model/volsys.cpp',
        'steps/model_wrap.cpp']),
        Extension('_geom_swig',
        sources = ['cpp/steps/error.cpp','cpp/steps/math/tetrahedron.cpp','cpp/steps/geom/tetmesh.cpp',
        'cpp/steps/math/linsolve.cpp','cpp/steps/math/triangle.cpp',
        'cpp/steps/geom/comp.cpp','cpp/steps/geom/geom.cpp','cpp/steps/geom/patch.cpp',
        'cpp/steps/geom/tet.cpp', 'cpp/steps/geom/tetmesh_rw.cpp',
        'cpp/steps/geom/tmcomp.cpp','cpp/steps/geom/tmpatch.cpp','cpp/steps/geom/tri.cpp',
        'steps/geom_wrap.cpp']),
        Extension('_solver_swig',
        sources = ['cpp/steps/error.cpp',
        'cpp/steps/rng/rng.cpp', 'cpp/steps/rng/mt19937.cpp',
        'cpp/steps/model/model.cpp', 'cpp/steps/model/diff.cpp',
        'cpp/steps/model/reac.cpp','cpp/steps/model/spec.cpp','cpp/steps/model/sreac.cpp',
        'cpp/steps/model/surfsys.cpp','cpp/steps/model/volsys.cpp',
        'cpp/steps/math/tetrahedron.cpp','cpp/steps/geom/tetmesh.cpp',
        'cpp/steps/math/linsolve.cpp','cpp/steps/math/triangle.cpp',
        'cpp/steps/geom/comp.cpp','cpp/steps/geom/geom.cpp','cpp/steps/geom/patch.cpp',
        'cpp/steps/geom/tet.cpp', 'cpp/steps/geom/tetmesh_rw.cpp',
        'cpp/steps/geom/tmcomp.cpp','cpp/steps/geom/tmpatch.cpp','cpp/steps/geom/tri.cpp',
        'cpp/steps/solver/api_comp.cpp','cpp/steps/solver/api_main.cpp',
        'cpp/steps/solver/api_patch.cpp','cpp/steps/solver/api_tet.cpp',
        'cpp/steps/solver/api_tri.cpp','cpp/steps/solver/compdef.cpp',
        'cpp/steps/solver/diffdef.cpp','cpp/steps/solver/patchdef.cpp',
        'cpp/steps/solver/reacdef.cpp','cpp/steps/solver/specdef.cpp',
        'cpp/steps/solver/sreacdef.cpp','cpp/steps/solver/statedef.cpp',
        'cpp/steps/tetexact/comp.cpp','cpp/steps/tetexact/diff.cpp',
        'cpp/steps/tetexact/kproc.cpp','cpp/steps/tetexact/patch.cpp',
        'cpp/steps/tetexact/reac.cpp','cpp/steps/tetexact/sreac.cpp',
        'cpp/steps/tetexact/tet.cpp','cpp/steps/tetexact/tetexact.cpp',
        'cpp/steps/tetexact/tri.cpp',
        'cpp/steps/wmdirect/comp.cpp','cpp/steps/wmdirect/kproc.cpp',
        'cpp/steps/wmdirect/patch.cpp','cpp/steps/wmdirect/reac.cpp',
        'cpp/steps/wmdirect/sreac.cpp','cpp/steps/wmdirect/wmdirect.cpp',
        'cpp/steps/wmrk4/wmrk4.cpp',
        'steps/solver_wrap.cpp']),
      ]
      )