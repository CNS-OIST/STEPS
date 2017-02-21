# -*- coding: utf-8 -*-

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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This file is the user-interface file for all geom objects in STEPS.
# All objects are directly derived from the corresponding swig objects.
# Geom and Tetmesh(derived from Geom) container objects are owned by Python
# All other objects are owned by c++ and container is responsible for 
# all the cleaning-up of these objects (see cpp/geom/geom.cpp class 
# destructor).
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

try:
    from . import steps_swig_numpy as steps_swig
    import _steps_swig_numpy as _steps_swig
except:
    from . import steps_swig
    import _steps_swig

import steps

### Now defunct mesh saving/loading tool ###
# from steps_swig import loadASCII, saveASCII

###
# add for mesh object casting from wm objects
# Weiliang 20130521
def castToTmComp(c):
    """
        Construction::
        tmcomp = steps.geom.castToTmComp(c)
        
        Arguments:
        * steps.geom.Comp c
        
        Return:
        steps.geom.TmComp
        
        Try to cast a steps.geom.Comp object to steps.geom.TmComp.
        """
    return steps_swig.castToTmComp(c)

def castToTmPatch(p):
    """
        Construction::
        tmpatch = steps.geom.castToTmPatch(p)
        
        Arguments:
        * steps.geom.Patch p
        
        Return:
        steps.geom.TmPatch
        
        Try to cast a steps.geom.Patch object to steps.geom.TmPatch.
        """
    return steps_swig.castToTmPatch(p)

###

ELEM_VERTEX = steps_swig.ELEM_VERTEX
ELEM_TRI = steps_swig.ELEM_TRI
ELEM_TET = steps_swig.ELEM_TET
ELEM_UNDEFINED = steps_swig.ELEM_UNDEFINED

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Geom(steps_swig.Geom):
    """
    Top-level geometry container to which a number of compartment objects 
    and patches objects may be grouped. 
    
    A steps.geom.Geom object is parent to the following objects:

    * steps.geom.Comp
    * steps.geom.Patch
    """
    def __init__(self, *args): 
        """
        Construction::
        
            g = steps.geom.Geom()
            
        Create a geometry container object.
            
        Arguments: 
            None
        """
        this = _steps_swig.new_Geom(*args)
        try: self.this.append(this)
        except: self.this = this
        # set Geom object to do all the cleaning-up
        self.thisown = 1


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


class Patch(steps_swig.Patch):
    """
    Base class for patch objects. A patch is a piece of 2D surface surrounding 
    (part of) a 3D compartment, which may be connected to another compartment. 
    It provides basic functionality and descriptive data that is shared by the derived 
    class steps.geom.TmPatch (that is used to describe a surface comprised of 
    triangles in a tetrahedral mesh):

    * Getting and setting a valid patch identifier string, and handling 
      the interaction with the container object.
    * Getting (and at least in this base class also setting) the total area of the patch.
    * Adding surface systems to the patch.
    * References to inside/outside compartments to which the patch is adjoined.

    This base class can be used directly with well-mixed solvers. 
    """
    def __init__(self, *args, **kwargs): 
        """
        Construction::
        
            patch = steps.geom.Patch(id, container, icomp, ocomp = None, area = 0.0)
            
        Construct a patch object with identifier string id, assign container 
        as the parent geometry container and assign icomp as the "inner" 
        compartment. Assign also ocomp as the "outer" compartment (if required) 
        and optionally set the area to area (in m^2). 
                    
        Arguments: 
            * string id
            * steps.geom.Geom container
            * steps.geom.Comp icomp
            * steps.geom.Comp ocomp (default = None)
            * float area (default = 0.0)
            
        .. note: "Inner" compartment and "outer" compartment are purely defined 
           by their order to the class constructor.
        
        """
        this = _steps_swig.new_Patch(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_setmethods__["id"] = _steps_swig.Patch_setID
        self.__swig_getmethods__["id"] = _steps_swig.Patch_getID
        self.__swig_getmethods__["container"] = _steps_swig.Patch_getContainer
        self.__swig_getmethods__["surfsys"] = _steps_swig.Patch_getSurfsys
        self.__swig_setmethods__["area"] = _steps_swig.Patch_setArea
        self.__swig_getmethods__["area"] = _steps_swig.Patch_getArea
        self.__swig_getmethods__["icomp"] = _steps_swig.Patch_getIComp
        self.__swig_getmethods__["ocomp"] = _steps_swig.Patch_getOComp
        
    id = steps_swig._swig_property(_steps_swig.Patch_getID, _steps_swig.Patch_setID)
    """Identifier string of the patch."""
    container = steps_swig._swig_property(_steps_swig.Patch_getContainer)
    """Reference to parent steps.geom.Geom container."""
    surfsys = steps_swig._swig_property(_steps_swig.Patch_getSurfsys)
    """Reference to assocated surface system."""
    area = steps_swig._swig_property(_steps_swig.Patch_getArea, _steps_swig.Patch_setArea)
    """Area of the patch."""
    icomp = steps_swig._swig_property(_steps_swig.Patch_getIComp)
    """Reference to the inner compartment."""
    ocomp = steps_swig._swig_property(_steps_swig.Patch_getOComp)
    """reference to the outer compartment."""


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


class Comp(steps_swig.Comp):
    """
    Base class for compartment objects. It provides basic functionality and data 
    that is shared by the derived class steps.geom.TmComp (that is used to 
    describe a compartment in a tetrahedral mesh):
    
    * Getting and setting a valid compartment identifier string, and handling 
      the interaction with the container object. 
    * Getting (and at least in this base class also setting) the total volume 
      of the compartment. 
    * Adding volume systems to the compartment. 
    * References to steps.geom.Patch objects to which the compartment is adjoined.
    
    This base class can be used directly with well-mixed solvers.
    
    **Relationship between Compartments and Patches**
    
    It is necessary to explain the inner/outer relationship between compartments 
    and patches. When a patch object is created (a surface of a compartment, which 
    may be shared with another compartment). it is necessary to arbitrarily label the 
    compartment(s) "inner" and "outer" (if a patch is connected to only one compartment
    then the compartment must be labelled "inner" by convention). This is necessary 
    in order to fully describe the surface reaction rules. Accordingly, compartments 
    also store a list of connections, "inner" patches and "outer" patches. So if a 
    patch1 is created with comp1 as it's "inner" compartment, comp1 knows patch 1 as 
    an "outer" patch. The labelling is purely defined when creating the Patch objects, 
    bearing in mind the stoichiometry defined in the surface reaction objects. This may 
    seem a little confusing at first, but will become clearer when experience is gained 
    with these objects.
     
    .. seealso:: :ref:`ip3`
     
    """     
    def __init__(self, *args, **kwargs): 
        """
        Construction::
        
            comp = steps.geom.Comp(id, container, vol = 0.0)
            
        Construct a compartment object with identifier string id and assign 
        container as the parent geometry container. Optionally set volume 
        to vol (in m^3).
            
        Arguments: 
            * string id
            * steps.geom.Geom container
            * float vol (default = 0.0)
        """
        this = _steps_swig.new_Comp(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_setmethods__["id"] = _steps_swig.Comp_setID
        self.__swig_getmethods__["id"] = _steps_swig.Comp_getID
        self.__swig_getmethods__["container"] = _steps_swig.Comp_getContainer
        self.__swig_getmethods__["volsys"] = _steps_swig.Comp_getVolsys
        self.__swig_setmethods__["vol"] = _steps_swig.Comp_setVol
        self.__swig_getmethods__["vol"] = _steps_swig.Comp_getVol
        self.__swig_getmethods__["ipatches"] = _steps_swig.Comp_getIPatches
        self.__swig_getmethods__["opatches"] = _steps_swig.Comp_getOPatches
    
    id = steps_swig._swig_property(_steps_swig.Comp_getID, _steps_swig.Comp_setID)
    """Identifier string of the compartment."""
    container = steps_swig._swig_property(_steps_swig.Comp_getContainer)
    """Reference to parent steps.geom.Geom container."""
    volsys = steps_swig._swig_property(_steps_swig.Comp_getVolsys)
    """Reference to assocated volume system."""
    vol = steps_swig._swig_property(_steps_swig.Comp_getVol, _steps_swig.Comp_setVol)
    """Volume of the compartment."""
    ipatches = steps_swig._swig_property(_steps_swig.Comp_getIPatches)
    """List of reference to inner patches."""
    opatches = steps_swig._swig_property(_steps_swig.Comp_getOPatches)
    """List of reference to outer patches."""


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


class Tetmesh(steps_swig.Tetmesh):
    """
    Main container class for static tetrahedral meshes. This class stores the 
    vertices points, 3D tetrahedral and 2D triangular elements that comprise 
    the mesh. The indices of the elements will be stored as unsigned integers 
    (a positive integer or zero) beginning at zero and incremented by 1. For 
    example, if there are ntets number of tetrahedrons in the mesh, the indices 
    of the tetrahedrons will be [0,1,2,..., (ntets-1)]. In addition, this class 
    computes and contains some auxiliary data from the mesh:

    * Rectangular, axis-aligned bounding box.
    * Overall volume.
    * The total number of tetrahedrons in the mesh.
    * The total number of triangles in the mesh.
    * The total number of vertices in the mesh.

    Auxiliary data is also stored for each tetrahedron:
    
    * Volume of the tetrahedron.
    * Indices of the 4 neighbouring tetrahedrons. If there is no neighbour 
      (i.e. if the tetrahedron lies on the border), this index will be -1. 
      The sequence of neighbours is determined by the following common boundary 
      triangles: (0,1,2), (0,1,3), (0,2,3), (1,2,3).
    * Indices of the 4 neighbouring boundary triangles. The sequence of neighbours 
      is also determined by (0,1,2), (0,1,3), (0,2,3), (1,2,3).
    * Compartment (steps.geom.TmComp object) that the tetrahedron belongs to. 
      Stores zero pointer if the tetrahedron has not been assigned to a compartment.

    And for each triangle:
    
    * Area of the triangle.
    * Normal vector of the triangle, normalized to length 1.0.
    * Indices of the 2 neighbouring tetrahedrons. If one tetrahedron does not exist 
      (because the triangle lies on the outer boundary), this index will be -1.
    * Patch (steps.geom.TmPatch object) that a triangle belongs to. Stores zero pointer 
      if triangle has not been assigned to a patch.
    """
    def __init__(self, *args): 
        """
        Construction1::
        
            mesh = steps.geom.Tetmesh(verts, tets, tris)
            
        Construct a Tetmesh container by the “first” method: Supply a list of all 
        vertices verts (by Cartesian coordinates), supply a list of all tetrahedrons 
        tets (by indices of the 4 vertices) and supply a full or partial list of 
        triangles tris (by indices of the 3 vertices). Indexing in STEPS begins at 
        zero, so the first 3 coordinates in verts will describe the zeroth vertex, t
        he next 3 coordinates will describe the 1st vertex and so on. Labelling of 
        the vertices in tets and tris should follow this indexing. Lists must be 
        one-dimensional. Length of verts = nverts*3 where nverts is the total 
        number of vertices; length of tets = ntets*4 where ntets is the total 
        number of tetrahedrons; maximum length of tris ntris*3 where ntris is 
        the total number of triangles. For example, if we have just three tetrahedrons; 
        tet0=[0,1,2,3], tet1=[0,1,3,4] and tet2=[1,3,4,5] then the required 
        one-dimensional list tets=[0,1,2,3,0,1,3,4,1,3,4,5]. 
            
        Arguments: 
            * list<float> verts
            * list<uint> tets
            * list<unit> tris
            
        Construction2::
            mesh = steps.geom.Tetmesh(nverts, ntets, ntris)
            
        Construct a Tetmesh container by the “second” method: Supply only the 
        number of vertices nverts, the number of tetrahedrons ntets and the number 
        of triangles ntris to the initializer, then use set methods to supply the 
        vertex, tetrahedron and triangle information one by one. It is up to the 
        user to make sure all information is supplied and then call setup() explicitly. 
        It is highly recommended to use the first constructor wherever possible due to t
        he scope for user error when using this method.
        
        Arguments: 
            * uint nverts
            * uint ntets
            * uint ntris
        """        
        this = _steps_swig.new_Tetmesh(*args)
        try: self.this.append(this)
        except: self.this = this
        # set Tetmesh object to do all the cleaning up
        self.thisown = 1
        self.__swig_getmethods__["nverts"] = _steps_swig.Tetmesh_countVertices
        self.__swig_getmethods__["ntris"] = _steps_swig.Tetmesh_countTris
        self.__swig_getmethods__["ntets"] = _steps_swig.Tetmesh_countTets
    
    nverts = steps_swig._swig_property(_steps_swig.Tetmesh_countVertices)
    """Number of vertices in the mesh."""
    ntris = steps_swig._swig_property(_steps_swig.Tetmesh_countTris)
    """Number of triangles in the mesh."""
    ntets = steps_swig._swig_property(_steps_swig.Tetmesh_countTets)
    """Number of tetrahedrons in the mesh."""


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


class TmComp(steps_swig.TmComp):
    """
    Derived class from base steps.geom.Comp class. It provides the same 
    functionality as the steps.geom.Comp class extended for annotation of 
    a group of tetrahedrons in a Tetmesh. The volume is the total volume of 
    the encapsulated tetrahedrons.
    """
    def __init__(self, *args, **kwargs): 
        """
        Construction::
        
            tmcomp = steps.geom.Comp(id, container, tets)
            
        Construct a TmComp object with identifier string id and assign container 
        as the parent Tetmesh container. Set the group of tetrahedrons that describe 
        this compartment with tets.
            
        Arguments: 
            * string id
            * steps.geom.Tetmesh container
            * list<uint> tets
        """
        this = _steps_swig.new_TmComp(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_getmethods__["tets"] = _steps_swig.TmComp_getAllTetIndices
    
    tets = steps_swig._swig_property(_steps_swig.TmComp_getAllTetIndices)
    """List of indices of tetrahedrons associated to the compartment."""


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


class TmPatch(steps_swig.TmPatch):
    """
    Derived class from base steps.geom.Patch class. It provides the same 
    functionality as the steps.geom.Patch class extended for annotation of 
    a group of triangles in a Tetmesh. The area is the total area of the 
    encapsulated triangles.
    """    
    def __init__(self, *args, **kwargs): 
        """
        Construction::
        
            tmpatch = steps.geom.Comp(id, container, tris, icomp, ocomp = None)
            
            Construct a TmPatch object with identifier string id and assign container 
            as the parent geometry container. Set the collection of triangles in 
            the patch to tris and assign icomp as the “inner” compartment 
            and assign also ocomp as the “outer” compartment (if required). 
            
        Arguments: 
            * string id
            * steps.geom.Tetmesh container
            * list<uint> tris
            * steps.geom.TmComp icomp
            * steps.geom.TmComp ocomp (default = None)
        """
        this = _steps_swig.new_TmPatch(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_getmethods__["tris"] = _steps_swig.TmPatch_getAllTriIndices
        
    tris = steps_swig._swig_property(_steps_swig.TmPatch_getAllTriIndices)
    """List of indices of triangles associated to the patch."""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class Memb(steps_swig.Memb):
    """
        This class provides annotation for a group of triangles that comprise  
        a surface describing a membrane in a Tetmesh. This may be the same 
        description of one or several TmPatches in order for voltage-dependent 
        transitions, currents and so on to be inserted in the membrane. 
        A Memb object must be available if membrane potential calculation is to be 
        performed.
        """
    
    def __init__(self, *args, **kwargs):
        """
            Construction:
            
            memb = steps.geom.Memb(id, container, patches, verify = False, opt_method=1, opt_file_name = '')			
            
            Construct a Memb object with identifier string id and assign container
            as the parent geometry container. Set the collection of triangles in 
            the membrane to those that belong to all TmPatches in patches. Perform some checks on 
            suitability of membrane if verify is True: these checks will print
            warnings if membrane forms an open surface or if any triangle is found to have 
            more than 3 neighbours. Specify optimization method with opt_method (default = 1): 
            1 = principal axis ordering (quick to set up but usually results in slower simulation than method 2). 
            2 = breadth first search (can be time-consuming to set up, but usually faster simulation), 
            If a filename (with full path) is given in optional argument opt_file_name the membrane optimization will be loaded from file,
            which was saved previously for this membrane with solver method steps.solver.Tetexact.saveMembOpt()
            
            Arguments:
            * string id
            * steps.geom.Tetmesh container
            * list<steps.geom.TmPatch> patches
            * bool verify (default = False)
            * uint opt_method (default = 1)
            * string opt_file_name (default = '')
            """
        this = _steps_swig.new_Memb(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_getmethods__["tris"] = _steps_swig.Memb_getAllTriIndices
    
    tris = steps_swig._swig_property(_steps_swig.Memb_getAllTriIndices)
    """List of indices of triangles associated to the membrane."""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class DiffBoundary(steps_swig.DiffBoundary):
    """
    Annotation of a group of triangles in a Tetmesh. The triangles form 
    a boundary between two compartments, that may allow diffusion of some
    specified species.
    """    
    def __init__(self, *args, **kwargs): 
        """
        Construction::
        
            diffb = steps.geom.DiffBoundary(id, container, tris)
            
        Construct a DiffBoundary object with identifier string id and assign container 
        as the parent geometry container, described by group of triangles tris.
            
        Arguments: 
            * string id
            * steps.geom.Tetmesh container
            * list<uint> tris
        """
        this = _steps_swig.new_DiffBoundary(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = 0
        self.__swig_getmethods__["tris"] = _steps_swig.DiffBoundary_getAllTriIndices
        self.__swig_setmethods__["id"] = _steps_swig.DiffBoundary_setID
        self.__swig_getmethods__["id"] = _steps_swig.DiffBoundary_getID
        self.__swig_getmethods__["container"] = _steps_swig.DiffBoundary_getContainer
        
    tris = steps_swig._swig_property(_steps_swig.DiffBoundary_getAllTriIndices)
    """List of indices of triangles associated to the patch."""
    id = steps_swig._swig_property(_steps_swig.DiffBoundary_getID, _steps_swig.DiffBoundary_setID)
    """Identifier string of the diffusion boundary."""
    container = steps_swig._swig_property(_steps_swig.DiffBoundary_getContainer)
    """Reference to parent steps.tetmesh.Tetmesh container."""
    
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
