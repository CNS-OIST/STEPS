####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

import numpy as np

from steps import stepslib
castToTmComp = stepslib.castToTmComp
castToTmPatch = stepslib.castToTmPatch

class Geom(stepslib._py_Geom): 
    """
    Top-level geometry container to which a number of compartment objects
    and patches objects may be grouped.
    
    A steps.geom.Geom object is parent to the following objects:
        * steps.geom.Comp
        * steps.geom.Patch
    """
    pass


class Comp(stepslib._py_Comp): 
    """
    Base class for compartment objects. It provides basic functionality and data
    that is shared by the derived class steps.geom.TmComp (that is used to
    describe a compartment in a tetrahedral mesh):

    * Getting and setting a valid compartment identifier string, and handling the interaction with the container object.
    * Getting (and at least in this base class also setting) the total volume of the compartment.
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
    
    .. seealso:: :ref:`/ip3.ipynb`
    
    """
    pass


class Tetmesh(stepslib._py_Tetmesh): 
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
    * Indices of the 4 neighbouring tetrahedrons. If there is no neighbour (i.e. if the tetrahedron lies on the border), this index will be UNKNOWN_TET. The sequence of neighbours is determined by the following common boundary triangles: (0,1,2), (0,1,3), (0,2,3), (1,2,3).
    * Indices of the 4 neighbouring boundary triangles. The sequence of neighbours is also determined by (0,1,2), (0,1,3), (0,2,3), (1,2,3).
    * Compartment (steps.geom.TmComp object) that the tetrahedron belongs to. Stores zero pointer if the tetrahedron has not been assigned to a compartment.
    
    And for each triangle:

    * Area of the triangle.
    * Normal vector of the triangle, normalized to length 1.0.
    * Indices of the 2 neighbouring tetrahedrons. If one tetrahedron does not exist (because the triangle lies on the outer boundary), this index will be UNKNOWN_TET.
    * Patch (steps.geom.TmPatch object) that a triangle belongs to. Stores zero pointer if triangle has not been assigned to a patch.
    """
    pass


try:
    class Library(stepslib._py_Library):
        def __init__(self, comm=None):
            if comm is None:
                import mpi4py
                comm = mpi4py.MPI.COMM_WORLD
            super().__init__(comm)

    class DistMesh(stepslib._py_DistMesh):
        pass

    class DistComp(stepslib._py_DistComp):
        pass

    class DistPatch(stepslib._py_DistPatch):
        pass

    class DistMemb(stepslib._py_DistMemb):
        pass
except AttributeError:
    pass

class Patch(stepslib._py_Patch):
    """
    Base class for patch objects. A patch is a piece of 2D surface surrounding
    (part of) a 3D compartment, which may be connected to another compartment.
    It provides basic functionality and descriptive data that is shared by the derived
    class steps.geom.TmPatch (that is used to describe a surface comprised of
    triangles in a tetrahedral mesh):
    
    * Getting and setting a valid patch identifier string, and handling the interaction with the container object.
    * Getting (and at least in this base class also setting) the total area of the patch.
    * Adding surface systems to the patch.
    * References to inside/outside compartments to which the patch is adjoined.
    
    This base class can be used directly with well-mixed solvers.
    """
    pass


class TmComp(stepslib._py_TmComp): 
    """
    Derived class from base steps.geom.Comp class. It provides the same
    functionality as the steps.geom.Comp class extended for annotation of
    a group of tetrahedrons in a Tetmesh. The volume is the total volume of
    the encapsulated tetrahedrons.
    """
    pass


class TmPatch(stepslib._py_TmPatch): 
    """
    Derived class from base steps.geom.Patch class. It provides the same
    functionality as the steps.geom.Patch class extended for annotation of
    a group of triangles in a Tetmesh. The area is the total area of the
    encapsulated triangles.
    """
    pass


class DiffBoundary(stepslib._py_DiffBoundary): 
    """
    Annotation of a group of triangles in a Tetmesh. The triangles form
    a boundary between two compartments, that may allow diffusion of some
    specified species.
    """
    pass


class SDiffBoundary(stepslib._py_SDiffBoundary): 
    """
    Annotation of a group of bars in a Tetmesh. The bars form
    a boundary between two patches, that may allow diffusion of some
    specified species.
    """
    pass


class Memb(stepslib._py_Memb): 
    """
    This class provides annotation for a group of triangles that comprise
    a surface describing a membrane in a Tetmesh. This may be the same
    description of one or several TmPatches in order for voltage-dependent
    transitions, currents and so on to be inserted in the membrane.
    A Memb object must be available if membrane potential calculation is to be
    performed.
    """
    pass



#Do we want to extract them module level?
ELEM_VERTEX     = stepslib._py_ElementType.ELEM_VERTEX
ELEM_TRI        = stepslib._py_ElementType.ELEM_TRI
ELEM_TET        = stepslib._py_ElementType.ELEM_TET
ELEM_UNDEFINED  = stepslib._py_ElementType.ELEM_UNDEFINED
UNKNOWN_TET     = stepslib.UNKNOWN_TET
UNKNOWN_TRI     = stepslib.UNKNOWN_TRI
INDEX_NUM_BYTES = stepslib.INDEX_NUM_BYTES
INDEX_DTYPE     = np.uint32 if INDEX_NUM_BYTES == 4 else np.uint64
DIST_INDEX_DTYPE= np.int64
