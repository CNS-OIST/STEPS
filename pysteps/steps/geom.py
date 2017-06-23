####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   
###

from steps import stepslib
castToTmComp = stepslib.castToTmComp
castToTmPatch = stepslib.castToTmPatch

class Geom(stepslib._py_Geom): pass
class Comp(stepslib._py_Comp): pass
class Tetmesh(stepslib._py_Tetmesh): pass
class Patch(stepslib._py_Patch): pass
class TmComp(stepslib._py_TmComp): pass
class TmPatch(stepslib._py_TmPatch): pass
class DiffBoundary(stepslib._py_DiffBoundary): pass
class SDiffBoundary(stepslib._py_SDiffBoundary): pass
class Memb(stepslib._py_Memb): pass

#Do we want to extract them module level?
ELEM_VERTEX    = stepslib._py_ElementType.ELEM_VERTEX
ELEM_TRI       = stepslib._py_ElementType.ELEM_TRI
ELEM_TET       = stepslib._py_ElementType.ELEM_TET
ELEM_UNDEFINED = stepslib._py_ElementType.ELEM_UNDEFINED
