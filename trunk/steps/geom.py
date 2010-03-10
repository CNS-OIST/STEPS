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


"""
This file is the user-interface file for all geom objects in STEPS.
All objects are directly derived from the corresponding swig objects.
Geom and Tetmesh(derived from Geom) container objects are owned by Python
All other objects are owned by c++ and container is responsible for 
all the cleaning-up of these objects (see cpp/steps/geom/geom.cpp class 
destructor).
"""


from . import geom_swig
import _geom_swig

### Now defunct mesh saving/loading tool ###
# from geom_swig import loadASCII, saveASCII


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

class Geom(geom_swig.Geom):
    
    def __init__(self, *args): 
        """__init__(self) -> Geom"""
        this = _geom_swig.new_Geom(*args)
        try: self.this.append(this)
        except: self.this = this
        # set Geom object to do all the cleaning-up
        self.thisown = True


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


class Patch(geom_swig.Patch):
    
    def __init__(self, *args, **kwargs): 
        """
        __init__(self, string id, Geom container, Comp icomp, Comp ocomp=0, 
            double area=0.0) -> Patch
        """
        this = _geom_swig.new_Patch(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _geom_swig.Patch_setID
        self.__swig_getmethods__["id"] = _geom_swig.Patch_getID
        self.__swig_getmethods__["container"] = _geom_swig.Patch_getContainer
        self.__swig_getmethods__["surfsys"] = _geom_swig.Patch_getSurfsys
        self.__swig_setmethods__["area"] = _geom_swig.Patch_setArea
        self.__swig_getmethods__["area"] = _geom_swig.Patch_getArea
        self.__swig_getmethods__["icomp"] = _geom_swig.Patch_getIComp
        self.__swig_getmethods__["ocomp"] = _geom_swig.Patch_getOComp
        
    id = geom_swig._swig_property(_geom_swig.Patch_getID, _geom_swig.Patch_setID)
    container = geom_swig._swig_property(_geom_swig.Patch_getContainer)
    surfsys = geom_swig._swig_property(_geom_swig.Patch_getSurfsys)
    area = geom_swig._swig_property(_geom_swig.Patch_getArea, _geom_swig.Patch_setArea)
    icomp = geom_swig._swig_property(_geom_swig.Patch_getIComp)
    ocomp = geom_swig._swig_property(_geom_swig.Patch_getOComp)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


class Comp(geom_swig.Comp):
    
    def __init__(self, *args, **kwargs): 
        """__init__(self, string id, Geom container, double vol=0.0) -> Comp"""
        this = _geom_swig.new_Comp(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _geom_swig.Comp_setID
        self.__swig_getmethods__["id"] = _geom_swig.Comp_getID
        self.__swig_getmethods__["container"] = _geom_swig.Comp_getContainer
        self.__swig_getmethods__["volsys"] = _geom_swig.Comp_getVolsys
        self.__swig_setmethods__["vol"] = _geom_swig.Comp_setVol
        self.__swig_getmethods__["vol"] = _geom_swig.Comp_getVol
        self.__swig_getmethods__["ipatches"] = _geom_swig.Comp_getIPatches
        self.__swig_getmethods__["opatches"] = _geom_swig.Comp_getOPatches
    
    id = geom_swig._swig_property(_geom_swig.Comp_getID, _geom_swig.Comp_setID)
    container = geom_swig._swig_property(_geom_swig.Comp_getContainer)
    volsys = geom_swig._swig_property(_geom_swig.Comp_getVolsys)
    vol = geom_swig._swig_property(_geom_swig.Comp_getVol, _geom_swig.Comp_setVol)
    ipatches = geom_swig._swig_property(_geom_swig.Comp_getIPatches)
    opatches = geom_swig._swig_property(_geom_swig.Comp_getOPatches)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


class Tetmesh(geom_swig.Tetmesh):
    
    def __init__(self, *args): 
        """
        __init__(self, unsigned int nverts, unsigned int ntets, unsigned int ntris) -> Tetmesh
        __init__(self, vector_dbl verts, vector_uint tets, vector_uint tris=std::vector< unsigned int >()) -> Tetmesh
        __init__(self, vector_dbl verts, vector_uint tets) -> Tetmesh
        """
        this = _geom_swig.new_Tetmesh(*args)
        try: self.this.append(this)
        except: self.this = this
        # set Tetmesh object to do all the cleaning up
        self.thisown = True
        self.__swig_getmethods__["nverts"] = _geom_swig.Tetmesh_countVertices
        self.__swig_getmethods__["ntris"] = _geom_swig.Tetmesh_countTris
        self.__swig_getmethods__["ntets"] = _geom_swig.Tetmesh_countTets
    
    nverts = geom_swig._swig_property(_geom_swig.Tetmesh_countVertices)
    ntris = geom_swig._swig_property(_geom_swig.Tetmesh_countTris)
    ntets = geom_swig._swig_property(_geom_swig.Tetmesh_countTets)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

"""
////////////////////////////////////////////////////////////////////////////////
//////// OBJECT REMOVED BECAUSE OF MEMORY ISSUES. SEE TODO NOTE IN C++ /////////
///////////////////////// CONSTRUCTOR FOR DETAILS //////////////////////////////
////////////////////////////////////////////////////////////////////////////////
"""

# class Tet(geom_swig.Tet):
#    def __init__(self, *args): 
#        """__init__(self, Tetmesh mesh, unsigned int tidx) -> Tet"""
#        this = _geom_swig.new_Tet(*args)
#        try: self.this.append(this)
#        except: self.this = this
#        # this object owned by python
#        self.thisown = True


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

"""
////////////////////////////////////////////////////////////////////////////////
//////// OBJECT REMOVED BECAUSE OF MEMORY ISSUES. SEE TODO NOTE IN C++ /////////
///////////////////////// CONSTRUCTOR FOR DETAILS //////////////////////////////
////////////////////////////////////////////////////////////////////////////////
"""

# class Tri(geom_swig.Tri):
#
#    def __init__(self, *args): 
#        """__init__(self, Tetmesh mesh, unsigned int tidx) -> Tri"""
#        this = _geom_swig.new_Tri(*args)
#        try: self.this.append(this)
#        except: self.this = this
#        # this object is owned by python
#        self.thisown = True


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


class TmComp(geom_swig.TmComp):
    
    def __init__(self, *args, **kwargs): 
        """__init__(self, string id, Tetmesh container, vector_uint tets) -> TmComp"""
        this = _geom_swig.new_TmComp(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_getmethods__["tets"] = _geom_swig.TmComp_getAllTetIndices
    
    tets = geom_swig._swig_property(_geom_swig.TmComp_getAllTetIndices)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


class TmPatch(geom_swig.TmPatch):
    
    def __init__(self, *args, **kwargs): 
        """
        __init__(self, string id, Tetmesh container, vector_uint tris, TmComp icomp, 
            TmComp ocomp=0) -> TmPatch
        """
        this = _geom_swig.new_TmPatch(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_getmethods__["tris"] = _geom_swig.TmPatch_getAllTriIndices
        
    tris = geom_swig._swig_property(_geom_swig.TmPatch_getAllTriIndices)
    

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
loadASCII.__doc__ = 
Reads a tetrahedral mesh from a simple ASCII format. Please refer to 
the documentation of steps.geom.saveASCII for more information on this
file format.  

PARAMETERS
    pathname
        The name of the file and its location.

RETURNS
    A steps.geom.Tetmesh object.

EXCEPTIONS
    steps.ArgErr
        Something is wrong with the pathname.
    steps.IOErr
        Cannot open the file or the file is badly formatted.

SEE ALSO
    steps.geom.saveASCII


saveASCII.__doc__ = 
saveASCII() saves a tetmesh to a very simple STEPS-centered ASCII 
(i.e. text based) format. The format is straightforward and consists 
of five parts, following the five parts of the tetmesh itself.

1/ A list of vertices: first the number of vertices is given. Then
   for each vertex there is a line with its x, y and z coordinates
   in scientific notation.

2/ A list of triangles. First the total number of triangles is given.  
   After this, we have 1 line per triangle giving its three corner 
   points as indices into the vertex list.

3/ A list of tetrahedrons. On a first line, the number of tets is 
   specified. After this, for each tet, follows a line with its four
   corner points, given as indices into the vertex list.

4/ The compartments in the mesh. First the number of compartments is
   given on a single line. Then each of these compartments is 
   specified using the following structure:
   
    a. Its name (ID string) on a single line.
    b. The number of volume systems that have been added to the
       compartment.
    c. For each volume system, its name (ID string) is given on a 
       separate line.
    d. The number of tetrahedrons associated with the compartment.
    e. The indices of these tetrahedrons (8 per line).

5/ A similar approach is taken to list all patches in the meshes. 
   first, we specify the number of patches attached to the mesh.
   Then, for each patch, we describe like this:
 
    a. The patch's name on a single line.
    b. If the patch has NO inner compartment, we print a "0"
       on a single line. Otherwise the line will start with a "1",
       followed by the name (ID string) of this inner compartment.
    c. The same for the outer compartment.
    d. The number of surface systems that have been added to the
       patch.
    e. For each of these surface systems, we specify its name on
       a separate line.
    f. The number of triangles in the patch.
    g. The indices of these triangles (8 per line).

All indices start from zero in this format, as they do in 
STEPS (and C++) internally.

PARAMETERS
    pathname
        The file and optionally the location that the mesh 
        will be written to.
    m
        A reference to a steps.geom.Tetmesh object.

RETURNS
    /

EXCEPTIONS
    steps.ArgErr
        The pathname is malformed or no mesh was specified.
    steps.IOErr
        Error while writing the mesh file (disk full; unable to 
        open file).

SEE ALSO
    steps.geom.loadASCII
"""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
