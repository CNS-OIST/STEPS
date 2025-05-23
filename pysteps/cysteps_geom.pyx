# cython:language_level=3str
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

"""
This file is the user-interface file for all geom objects in STEPS.  All
objects are directly derived from the corresponding swig objects.  Geom and
Tetmesh(derived from Geom) container objects are owned by Python All other
objects are owned by c++ and container is responsible for all the cleaning-up
of these objects (see cpp/geom/geom.cpp class destructor).
"""

from steps_wm cimport *
from steps_tetmesh cimport *
import warnings

# ======================================================================================================================
# Python Wrappers to namespace steps::wm
# ======================================================================================================================

cdef extern from "util/vocabulary.hpp":
    cdef index_t _UNKNOWN_TET "steps::tetrahedron_id_t::unknown_value()"
    cdef index_t _UNKNOWN_TRI "steps::triangle_id_t::unknown_value()"
    cdef index_t _UNKNOWN_VERT "steps::vertex_id_t::unknown_value()"
    cdef int _INDEX_NUM_BYTES "sizeof(steps::index_t)";

UNKNOWN_TET = _UNKNOWN_TET
UNKNOWN_TRI = _UNKNOWN_TRI
UNKNOWN_VERT = _UNKNOWN_VERT
INDEX_NUM_BYTES = _INDEX_NUM_BYTES

#Functions previously defined in the .i Swig files(!!)
def castToTmComp(_py_Comp base):
    """
    Construction::
        tmcomp = steps.geom.castToTmComp(c)

    Arguments:
        * steps.geom.Comp c

    Return:
    steps.geom.TmComp

    Try to cast a steps.geom.Comp object to steps.geom.TmComp.
    """
    return _py_TmComp.from_ptr( <TmComp*>(base.ptr()) )

def castToTmPatch(_py_Patch base):
    """
    Construction::
        tmpatch = steps.geom.castToTmPatch(p)

    Arguments:
        * steps.geom.Patch p

    Return:
    steps.geom.TmPatch

    Try to cast a steps.geom.Patch object to steps.geom.TmPatch.
    """
    return _py_TmPatch.from_ptr( <TmPatch*>(base.ptr()) )


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Geom(_py__base):
    "Python wrapper class for Geom"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Geom] _autodealoc
    cdef Geom *ptr(self):
        return <Geom*> self._ptr

    def __init__(self):
        """
        Construction::

            g = steps.geom.Geom()

        Create a geometry container object.

        Arguments:
        None
        """
        self._ptr = new Geom()      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getComp(self, str id):
        """
        Returns a reference to the steps.model.Comp compartment object with
        identifier string comp_id (if defined).

        Syntax::

            getComp(comp_id)

        Arguments:
        string comp_id

        Return:
        steps.model.Comp

        """
        return _py_Comp.from_ptr(&self.ptr().getComp(to_std_string(id)))

    def delComp(self, str id):
        """
        Removes the steps.geom.Comp object with identifier string comp_id (if defined)
        from the geometry container.

        Syntax::

            delComp(comp_id)

        Arguments:
        string comp_id

        Return:
        None

        """
        self.ptr().delComp(to_std_string(id))

    def getAllComps(self, ):
        """
        Returns a list of references to all steps.geom.Comp compartment objects in the
        geometry container.

        Syntax::

            getAllComps()

        Arguments:
        None

        Return:
        list<steps.geom.Comp>

        """
        return _py_Comp.vector2list(self.ptr().getAllComps())

    def getPatch(self, str id):
        """
        Removes the steps.geom.Patch object with identifier string patch_id (if defined)
        from the geometry container.

        Syntax::

            getPatch(patch_id)

        Arguments:
        string patch_id

        Return:
        steps.geom.Patch

        """
        return _py_Patch.from_ptr(&self.ptr().getPatch(to_std_string(id)))

    def delPatch(self, str id):
        """
        Removes the steps.geom.Patch object with identifier string patch_id (if defined)
        from the geometry container.

        Syntax::

            delPatch(patch_id)

        Arguments:
        string patch_id

        Return:
        None

        """
        self.ptr().delPatch(to_std_string(id))

    def getAllPatches(self, ):
        """
        Returns a list of references to all steps.geom.Patch patch objects in the
        geometry container.

        Syntax::

            getAllPatches()

        Arguments:
        None

        Return:
        list<steps.geom.Patch>

        """
        return _py_Patch.vector2list(self.ptr().getAllPatches())

    @staticmethod
    cdef _py_Geom from_ptr(Geom *ptr):
        if (ptr == NULL):
            return None
        cdef _py_Geom obj = _py_Geom.__new__(_py_Geom)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Geom from_ref(const Geom &ref):
        return _py_Geom.from_ptr(<Geom*>&ref)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Patch(_py__base):
    "Python wrapper class for Patch"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Patch] _autodealoc
    cdef Patch *ptr(self):
        return <Patch*> self._ptr

    def __init__(self, str id, _py_Geom container, _py_Comp icomp, _py_Comp ocomp=None, double area=0):
        """
        Construction::

            patch = steps.geom.Patch(id, container, icomp, ocomp = None, area = 0.0)

        Construct a patch object with identifier string id, assign container
        as the parent geometry container and assign icomp as the "inner"
        compartment. Assign also ocomp as the "outer" compartment (if required)
        and optionally set the area to area (in m^2).

        Arguments:
        string id
        steps.geom.Geom container
        steps.geom.Comp icomp
        steps.geom.Comp ocomp (default = None)
        float area (default = 0.0)

        .. note::
            "Inner" compartment and "outer" compartment are purely defined
            by their order to the class constructor.

        """
        self._ptr = new Patch(to_std_string(id), deref(container.ptr()), deref(icomp.ptr()), ocomp.ptr() if ocomp else NULL, area )

    def getID(self, ):
        """
        Get the identifier string of the patch.

        Syntax::

            getID()

        Arguments:
        None

        Return:
        string

        """
        return from_std_string(self.ptr().getID())

    def setID(self, str id):
        """
        Set the identifier string of the patch.

        Syntax::

            setID(name)

        Arguments:
        string name

        Return:
        None

        """
        self.ptr().setID(to_std_string(id))

    def getContainer(self, ):
        """
        Returns a reference to the parent steps.geom.Geom container object.

        Syntax::

            getContainer()

        Arguments:
        None

        Return:
        steps.geom.Geom

        """
        return _py_Geom.from_ptr(&self.ptr().getContainer())

    def getArea(self, ):
        """
        Get the area of the patch (in m^2).

        Syntax::

            getArea()

        Arguments:
        None

        Return:
        float

        """
        return self.ptr().getArea()

    def setArea(self, double vol):
        """
        Set the area of the patch (in m^2).

        Syntax::

            setArea(area)

        Arguments:
        float area

        Return:
        None

        """
        self.ptr().setArea(vol)

    def addSurfsys(self, str id):
        """
        Add surface system identifier string surfsys_id to the patch object.

        Syntax::

            addSurfsys(surfsys_id)

        Arguments:
        string surfsys_id

        Return:
        None

        """
        self.ptr().addSurfsys(to_std_string(id))

    def getSurfsys(self, ):
        """
        Returns a list of the surface system identifier strings which have
        been added to the patch.

        Syntax::

            getSurfsys()

        Arguments:
        None

        Return:
        list<string>

        """
        return string_flat_set_to_list(self.ptr().getSurfsys())

    def delSurfsys(self, str id):
        """
        Removes surface system identifier string surfsys_id from this patch.

        Syntax::

            delSurfsys(surfsys_id)

        Arguments:
        string surfsys_id

        Return:
        None

        """
        self.ptr().delSurfsys(to_std_string(id))

    def getAllSpecs(self, _py_Model model):
        """
        Given a steps.model.Model, return all species in the patch. That
        is, all 'surface' species referenced in surface system surface
        reactions, voltage-dependent surface reactions and surface diffusion.
        This also includes ohmic current and ghk current 'channel states'.

        Syntax::

            getAllSpecs(model)

        Arguments:
        steps.model.Model model

        Return:
        list<steps.model.Spec>

        """
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs(deref(model.ptr())))

    def getAllSReacs(self, _py_Model model):
        """
        Given a steps.model.Model, return all surface reactions in the compartment.

        Syntax::

            getAllSReacs(model)

        Arguments:
        steps.model.Model model

        Return:
        list<steps.model.SReac>

        """
        return _py_SReac.flat_set2list(self.ptr().getAllSReacs(deref(model.ptr())))

    def getIComp(self, ):
        """
        Returns a reference to the steps.geom.Comp compartment object representing
        the inner compartment.

        Syntax::

            getIComp()

        Arguments:
        None

        Return:
        steps.geom.Comp

        """
        return _py_Comp.from_ptr(&self.ptr().getIComp())

    def getOComp(self, ):
        """
        Returns a reference to the steps.geom.Comp compartment object representing
        the outer compartment.

        Syntax::

            getOComp()

        Arguments:
        None

        Return:
        steps.geom.Comp

        """
        return _py_Comp.from_ptr(self.ptr().getOComp())

    @staticmethod
    cdef _py_Patch from_ptr(Patch *ptr):
        if (ptr == NULL):
            return None
        cdef _py_Patch obj = _py_Patch.__new__(_py_Patch)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Patch from_ref(const Patch &ref):
        return _py_Patch.from_ptr(<Patch*>&ref)

    @staticmethod
    cdef std.vector[Patch*] *list2vector(list patchList, std.vector[Patch*] *dstVec):
        for item in patchList:
            assert isinstance(item, _py_Patch), "Wrong type. Expected _py_Patch, given: " + str(type(item))
            dstVec.push_back( (<_py_Patch>item).ptr())
        return dstVec

    @staticmethod
    cdef list vector2list(std.vector[Patch*] patches):
        return [ _py_Patch.from_ptr(p) for p in patches ]

    @staticmethod
    cdef set stdset2set(std.set[Patch*] patches):
        return { _py_Patch.from_ptr(p) for p in patches }

    @staticmethod
    cdef set flatset2set(flat_set[Patch*] patches):
        return { _py_Patch.from_ptr(p) for p in patches }

    def __hash__(self):
        return hash(self.id)

    ## properties ##
    id = property(getID, setID, doc="Identifier string of the patch.")
    container = property(getContainer, doc="Reference to parent steps.geom.Geom container.")
    surfsys = property(getSurfsys, doc="Reference to assocated surface system.")
    area = property(getArea, setArea, doc="Area of the patch.")
    icomp = property(getIComp, doc="Reference to the inner compartment.")
    ocomp = property(getOComp, doc="Reference to the outer compartment.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Comp(_py__base):
    "Python wrapper class for Comp"
# ----------------------------------------------------------------------------------------------------------------------
    cdef Comp *ptr(self):
        return <Comp*> self._ptr

    def __init__(self, str id, _py_Geom container, double vol=0):
        """
        Construction::

            comp = steps.geom.Comp(id, container, vol = 0.0)

        Construct a compartment object with identifier string id and assign
        container as the parent geometry container. Optionally set volume
        to vol (in m^3).

        Arguments:
        string id
        steps.geom.Geom container
        float vol (default = 0.0)
        """
        self._ptr = new Comp(to_std_string(id), deref(container.ptr()), vol)     # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        """
        Get the identifier string of the compartment.

        Syntax::

            getID()

        Arguments:
        None

        Return:
        string

        """
        return from_std_string(self.ptr().getID())

    def setID(self, str id):
        """
        Set the identifier string of the compartment.

        Syntax::

            setID(name)

        Arguments:
        string name

        Return:
        None

        """
        self.ptr().setID(to_std_string(id))

    def getContainer(self, ):
        """
        Returns a reference to the parent steps.geom.Geom container object.

        Syntax::

            getContainer()

        Arguments:
        None

        Return:
        steps.geom.Geom

        """
        return _py_Geom.from_ptr(&self.ptr().getContainer())

    def getVol(self, ):
        """
        Get the volume of the compartment (in m^3).

        Syntax::

            getVol()

        Arguments:
        None

        Return:
        float

        """
        return self.ptr().getVol()

    def setVol(self, double vol):
        """
        Set the volume of the compartment (in m^3).

        Syntax::

            setVol(vol)

        Arguments:
        float vol

        Return:
        None

        """
        self.ptr().setVol(vol)

    def addVolsys(self, str id):
        """
        Add volume system identifier string volsys_id to the compartment object.

        Syntax::

            addVolsys(volsys_id)

        Arguments:
        string volsys_id

        Return:
        None

        """
        return self.ptr().addVolsys(to_std_string(id))

    def getVolsys(self, ):
        """
        Returns a list of the volume system identifier strings which have been
        added to the compartment.

        Syntax::

            getVolsys()

        Arguments:
        None

        Return:
        list<string>

        """
        return string_flat_set_to_list(self.ptr().getVolsys())

    def delVolsys(self, str id):
        """
        Removes volume system identifier string volsys_id from this compartment.

        Syntax::

            delVolsys(volsys_id)

        Arguments:
        string volsys_id

        Return:
        None

        """
        self.ptr().delVolsys(to_std_string(id))

    def getAllSpecs(self, _py_Model model):
        """
        Given a steps.model.Model, return all species in the compartment. That
        is, all species that appear in the compartment due to volume system
        reaction and diffusion rules, plus surface system surface reactions,
        voltage-dependent surface reactions etc where reactants or products
        reference this compartment as 'inner' or 'outer' compartment.

        Syntax::

            getAllSpecs(model)

        Arguments:
        steps.model.Model model

        Return:
        list<steps.model.Spec>

        """
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs(deref(model.ptr())))

    def getAllReacs(self, _py_Model model):
        """
        Given a steps.model.Model, return all reactions in the compartment.

        Syntax::

            getAllReacs(model)

        Arguments:
        steps.model.Model model

        Return:
        list<steps.model.Reac>

        """
        return _py_Reac.flat_set2list(self.ptr().getAllReacs(deref(model.ptr())))

    def getAllDiffs(self, _py_Model model):
        """
        Given a steps.model.Model, return all diffusions in the compartment.

        Syntax::

            getAllDiffs(model)

        Arguments:
        steps.model.Model model

        Return:
        list<steps.model.Diff>

        """
        return _py_Diff.flat_set2list(self.ptr().getAllDiffs(deref(model.ptr())))

    # Sets, not vector!?
    def getIPatches(self, ):
        """
        Returns a list of references to steps.geom.Patch patch objects:
        the 'inner' patches.

        Syntax::

            getIPatches()

        Arguments:
        None

        Return:
        list<steps.geom.Patch>

        """
        return _py_Patch.flatset2set(self.ptr().getIPatches())

    def getOPatches(self, ):
        """
        Returns a list of references to steps.geom.Patch patch objects:
        the 'outer' patches.

        Syntax::

            getOPatches()

        Arguments:
        None

        Return:
        list<steps.geom.Patch>

        """
        return _py_Patch.flatset2set(self.ptr().getOPatches())

    @staticmethod
    cdef _py_Comp from_ptr(Comp *ptr):
        if (ptr == NULL):
            return None
        cdef _py_Comp obj = _py_Comp.__new__(_py_Comp)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Comp from_ref(const Comp &ref):
        return _py_Comp.from_ptr(<Comp*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[Comp*] vec):
        return [ _py_Comp.from_ptr(elem) for elem in vec ]

    ## properties ##
    id = property(getID, setID, doc="Identifier string of the compartment.")
    container = property(getContainer, doc="Reference to parent steps.geom.Geom container.")
    volsys = property(getVolsys, doc="Reference to assocated volume system.")
    vol = property(getVol, setVol, doc="Volume of the compartment.")
    ipatches = property(getIPatches, doc="List of reference to inner patches.")
    opatches = property(getOPatches, doc="List of reference to outer patches.")

# ======================================================================================================================
# Python bindings to namespace steps::tetmesh
# ======================================================================================================================

## Enums ##
cimport steps_tetmesh
cdef class _py_ElementType:
    ELEM_VERTEX = steps_tetmesh.ELEM_VERTEX
    ELEM_TRI = steps_tetmesh.ELEM_TRI
    ELEM_TET = steps_tetmesh.ELEM_TET
    ELEM_UNDEFINED = steps_tetmesh.ELEM_UNDEFINED

cdef class _py_ROISet:
    cdef readonly ElementType type
    cdef readonly std.vector[index_t] indices
    def __init__(self, ElementType t=ELEM_UNDEFINED, std.vector[index_t] indices=[]):
        self.type = t
        self.indices = indices

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Tetmesh(_py_Geom):
    "Python wrapper class for Tetmesh"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Tetmesh] _autodealoc
    cdef Tetmesh *ptrx(self):
        return <Tetmesh*> self._ptr

    def __init__(self, std.vector[double] verts, std.vector[index_t] tets, std.vector[index_t] tris=[]):
        """
        Syntax::

            mesh = steps.geom.Tetmesh(verts, tets, tris)

        Construct a Tetmesh container: Supply a list of all
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
        list<double> verts
        list<index_t> tets
        list<index_t> tris

        """
        self._ptr = new Tetmesh( verts, tets, tris )

    def getAllComps(self, ):
        """
        Returns a list of references to all steps.geom.TmComp compartment objects in the
        tetmesh container.

        Syntax::

            getAllComps()

        Arguments:
        None

        Return:
        list<steps.geom.TmComp>

        """
        return _py_TmComp.vector2list(self.ptr().getAllComps())

    def getAllPatches(self, ):
        """
        Returns a list of references to all steps.geom.TmPatch patch objects in the
        tetmesh container.

        Syntax::

            getAllPatches()

        Arguments:
        None

        Return:
        list<steps.geom.TmPatch>

        """
        return _py_TmPatch.vector2list(self.ptr().getAllPatches())

    def getVertex(self, index_t vidx):
        """
        Returns the coordinates of vertex with index vidx in the container.

        Syntax::

            getVertex(vidx)

        Arguments:
        index_t vidx

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getVertex(vertex_id_t(vidx))

    def countVertices(self, ):
        """
        Returns the total number of vertices in the mesh.

        Syntax::

            countVertices()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().countVertices()

    def getBar(self, index_t bidx):
        """
        Returns the vertices of bar with index bidx in the container.

        Syntax::

            getBar(bidx)

        Arguments:
        index_t bidx

        Return:
        list<index_t, length = 2>

        """
        return self.ptrx().getBar(bar_id_t(bidx))

    def countBars(self, ):
        """
        Returns the total number of bars in the mesh.

        Syntax::

            countBars()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().countBars()

    def getTri(self, index_t tidx):
        """
        Returns the triangle with index tidx in the container by its three vertex indices.

        Syntax::

            getTri(tidx)

        Arguments:
        index_t tidx

        Return:
        list<index_t, length = 3>

        """
        return self.ptrx().getTri(triangle_global_id(tidx))

    def countTris(self, ):
        """
        Returns the total number of triangles in the mesh.

        Syntax::

            countTris()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().countTris()

    def getTriArea(self, index_t tidx):
        """
        Returns the area of the triangle with index tidx.

        Syntax::

            getTriArea(tidx)

        Arguments:
        index_t tidx

        Return:
        float

        """
        return self.ptrx().getTriArea(triangle_global_id(tidx))

    def getTriBars(self, index_t tidx):
        """
        Returns the index of the bars that comprise the triangle.

        Syntax::

            getTriBars(tidx)

        Arguments:
        index_t tidx

        Return:
        list<int, length = 3>

        """
        return self.ptrx().getTriBars(triangle_global_id(tidx))

    def getTriBarycenter(self, index_t tidx):
        """
        Returns the Cartesian coordinates of the barycenter of triangle with index tidx.

        Syntax::

            getTriBarycenter(tidx)

        Arguments:
        index_t tidx

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getTriBarycenter(triangle_global_id(tidx))

    def getTriNorm(self, index_t tidx):
        """
        Returns the normal vector of the triangle with index tidx.

        Syntax::

            getTriNorm(tidx)

        Arguments:
        index_t tidx

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getTriNorm(triangle_global_id(tidx))

    def getTriPatch(self, index_t tidx):
        """
        Returns a reference to a step.geom.TmPatch object: the patch which triangle
        with index tidx belongs to. Returns None if triangle not assigned to a patch.

        Syntax::

            getTriPatch(tidx)

        Arguments:
        index_t tidx

        Return:
        steps.geom.TmPatch

        """
        return _py_TmPatch.from_ptr(self.ptrx().getTriPatch(triangle_global_id(tidx)))

    def setTriPatch(self, index_t tidx, _py_TmPatch patch):
        """
        Set the patch which triangle with index tidx belongs to..

        Syntax::

            setTriPatch(tidx, patch)

        Arguments:
        index_t tidx
        steps.geom.TmPatch patch

        Return:
        None

        """
        self.ptrx().setTriPatch(triangle_global_id(tidx), patch.ptrx())

    def setTriDiffBoundary(self, index_t tidx, _py_DiffBoundary diffb):
        """
        Set the diffusion boundary which triangle with index tidx belongs to..

        Syntax::

            setTriDiffBoundary(triangle_global_id(tidx), diffb)

        Arguments:
        index_t tidx
        steps.geom.DiffBoundary diffb

        Return:
        None

        """
        self.ptrx().setTriDiffBoundary(triangle_global_id(tidx), diffb.ptr())

    def getTriDiffBoundary(self, index_t tidx):
        """
        Returns a reference to a step.geom.Diffboundary object: the diffusion boundary triangle
        with index tidx belongs to. Returns None if triangle not assigned to a diffusion boundary.

        Syntax::

            getTriDiffBoundary(tidx)

        Arguments:
        index_t tidx

        Return:
        steps.geom.DiffBoundary

        """
        return _py_DiffBoundary.from_ptr(self.ptrx().getTriDiffBoundary(triangle_global_id(tidx)))

    def getTriTetNeighb(self, index_t tidx):
        """
        Returns the indices of the two neighbouring tetrahedrons of triangle with
        index tidx. An index of UNKNOWN_TET indicates no neighbour (triangle is on the mesh border).

        Syntax::

            getTriTetNeighb(tidx)

        Arguments:
        index_t tidx

        Return:
        list<index_t, length = 2>

        """
        return self.ptrx().getTriTetNeighb(triangle_global_id(tidx))

    def getTriTriNeighb(self, index_t tidx, _py_TmPatch tmpatch):
        """
        Returns the indices of the neighbouring triangles (that is all triangles that
        share a 'bar') of triangle with index tidx within the same patch.

        Syntax::

            getTriTriNeighb(int)

        Arguments:
        index_t tidx

        Returns:
        list<index_t>

        """
        return self.ptrx().getTriTriNeighb(triangle_global_id(tidx), tmpatch.ptrx())

    def getTriTriNeighbs(self, index_t tidx):
        """
        Returns the indices of the neighbouring triangles (that is all triangles that
        share a 'bar') of triangle with index tidx.

        Syntax::

            getTriTriNeighbs(int)

        Arguments:
        index_t tidx

        Returns:
        list<index_t>

        """
        return self.ptrx().getTriTriNeighbs(triangle_global_id(tidx))

    def getTet(self, index_t tidx):
        """
        Returns the tetrahedron with index tidx in the container by its four vertex indices.

        Syntax::
            getTet(tidx)

        Arguments:
        index_t tidx

        Return:
        list<index_t, length = 4>

        """
        return self.ptrx().getTet(tetrahedron_global_id(tidx))

    def countTets(self, ):
        """
        Returns the total number of tetrahedrons in the mesh.

        Syntax::

            countTets()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().countTets()

    def getTetVol(self, index_t tidx):
        """
        Returns the volume of the tetrahedron with index tidx.

        Syntax::

            getTetVol(tidx)

        Arguments:
        index_t tidx

        Return:
        float

        """
        return self.ptrx().getTetVol(tetrahedron_global_id(tidx))

    def getTetQualityRER(self, index_t tidx):
        """
        Returns the radius-edge-ratio (a quality measurement) of tetrahedron with index tidx.

        Syntax::

            getTetQualityRER(tidx)

        Arguments:
        index_t tidx

        Return:
        float

        """
        return self.ptrx().getTetQualityRER(tetrahedron_global_id(tidx))

    def getTetBarycenter(self, index_t tidx):
        """
        Returns the barycenter of the tetrahedron with index tidx.

        Syntax::

            getTetBarycenter(tidx)

        Arguments:
        index_t tidx

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getTetBarycenter(tetrahedron_global_id(tidx))

    def getTetComp(self, index_t tidx):
        """
        Returns a reference to a steps.geom.Comp object: the compartment which
        tetrahedron with index tidx belongs to. Returns None if tetrahedron not
        assigned to a compartment.

        Syntax::

            getTetComp(tidx)

        Arguments:
        index_t tidx

        Return:
        steps.geom.TmComp

        """
        return _py_TmComp.from_ptr(self.ptrx().getTetComp(tetrahedron_global_id(tidx)))

    def setTetComp(self, index_t tidx, _py_TmComp comp):
        """
        Set the compartment which tetrahedron with index tidx belongs to..

        Syntax::

            setTetComp(tidx, comp)

        Arguments:
        index_t tidx
        steps.geom.TmComp comp

        Return:
        None

        """
        self.ptrx().setTetComp(tetrahedron_global_id(tidx), comp.ptrx())

    def getTetTriNeighb(self, index_t tidx):
        """
        Returns the indices of the four neighbouring triangles of tetrahedron with index tidx.

        Syntax::

            getTetTriNeighb(tidx)

        Arguments:
        index_t tidx

        Return:
        list<index_t, length = 4>

        """
        return self.ptrx().getTetTriNeighb(tetrahedron_global_id(tidx))

    def getTetTetNeighb(self, index_t tidx):
        """
        Returns the indices of the four neighbouring tetrahedrons of tetrahedron with index tidx.
        An index of -1 indicates no neighbour (tetrahedron is on the mesh border).

        Syntax::

            getTetTetNeighb(tidx)

        Arguments:
        index_t tidx

        Return:
        list<index_t, length = 4>

        """
        return self.ptrx().getTetTetNeighb(tetrahedron_global_id(tidx))

    def findTetByPoint(self, std.vector[double] p):
        """
        Returns the index of the tetrahedron which encompasses a given point
        p (given in Cartesian coordinates x,y,z). Returns -1 if p is a position
        outside the mesh. It uses either findTetByPointLinear or findTetByPointWalk
        depending on the mesh size.

        Syntax::

            findTetByPoint(p)

        Arguments:
        list<float, length = 3> p

        Return:
        index_t

        """
        return self.ptrx().findTetByPoint(p).get()

    def findTetByPointLinear(self, std.vector[double] p):
        """
        Returns the index of the tetrahedron which encompasses a given point
        p (given in Cartesian coordinates x,y,z). Returns -1 if p is a position
        outside the mesh. Linear search. Only suitable for small meshes.

        Syntax::

            findTetByPointLinear(p)

        Arguments:
        list<double, length = 3> p

        Return:
        index_t

        """
        return self.ptrx().findTetByPointLinear(p).get()

    def findTetByPointWalk(self, std.vector[double] p):
        """
        Returns the index of the tetrahedron which encompasses a given point
        p (given in Cartesian coordinates x,y,z). Returns -1 if p is a position
        outside the mesh. A* search. After a seeding round of random points 
        the algorithm walks towards the closest tetrahedron. It works for disconnected
        meshes too. Usually faster than a linear search for normal-sized meshes.

        Syntax::

            findTetByPointWalk(p)

        Arguments:
        list<double, length = 3> p

        Return:
        index_t

        """
        return self.ptrx().findTetByPointWalk(p).get()

    def isPointInTet(self, std.vector[double] p, index_t tidx):
        """
        Check if point belongs to the tetrahedron or not

        Syntax::

            isPointInTet(p, tidx)

        Arguments:
        list<float, length = 3> p
        int tetrahedron tidx

        Return:
        bool
        """
        return self.ptrx().isPointInTet(p, tetrahedron_global_id(tidx))

    def getBoundMin(self, ):
        """
        Returns the minimal Cartesian coordinate of the rectangular bounding box of the mesh.

        Syntax::

            getBoundMin()

        Arguments:
        None

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMin()

    def getBoundMax(self, ):
        """
        Returns the maximal Cartesian coordinate of the rectangular bounding box of the mesh.

        Syntax::

            getBoundMax()

        Arguments:
        None

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMax()

    def getMeshVolume(self, ):
        """
        Returns the total volume of the mesh.

        Syntax::

            getMeshVolume()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getMeshVolume()

    def getSurfTris(self, ):
        """
        Returns a list of triangles that form the mesh boundary.
        Support function for steps.utilities.visual.

        Syntax::

            getTriBoundary()

        Arguments:
        None

        Return:
        list<index_t>

        """
        return self.ptrx().getSurfTris()

    ## Batch related
    def getBatchVertices(self, std.vector[index_t] verts):
        """
        Get coordinates of a list of vertices.

        Syntax::

            getBatchVertices(verts)

        Arguments:
        list<index_t> verts

        Return:
        list<float, length = len(verts) * 3>

        """
        return self.ptrx().getBatchVertices(verts)

    def getBatchVerticesNP(self, index_t[:] indices, double[:] coordinates):
        """
        Get coordinates of a list of vertices.

        Syntax::

            import numpy as np
            indices = np.array([0, 1, 2], dtype= np.uint64)
            coordinates = np.zeros(len(indices) * 3)
            getBatchVertices(indices, coordinates)

        Arguments:
        numpy.array<index_t> indices
        numpy.array<float, length = len(indices) * 3> coordinates

        Return:
        None

        """
        self.ptrx().getBatchVerticesNP(&indices[0], indices.shape[0], &coordinates[0], coordinates.shape[0])

    def getBatchTris(self, std.vector[index_t] tris):
        """
        Get vertex indices of a list of triangles.

        Syntax::

            getBatchTris(tris)

        Arguments:
        list<index_t> tris

        Return:
        list<index_t, length = len(tris) * 3>

        """
        return self.ptrx().getBatchTris(tris)

    def getBatchTrisNP(self, index_t[:] t_indices, index_t[:] v_indices):
        """
        Get vertex indices of a list of triangles.

        Syntax::

            getBatchTrisNP(t_indices, v_indices)

        Arguments:
        numpy.array<index_t> t_indices
        numpy.array<index_t, length = len(t_indices) * 3> v_indices

        Return:
        None

        """
        self.ptrx().getBatchTrisNP(&t_indices[0], t_indices.shape[0], &v_indices[0], v_indices.shape[0])

    def getBatchTets(self, std.vector[index_t] tets):
        """
        Get vertex indices of a list of tetrahedrons.

        Syntax::

            getBatchTets(tets)

        Arguments:
        list<index_t> tets

        Return:
        list<index_t, length = len(tets) * 4>

        """
        return self.ptrx().getBatchTets(tets)

    def getBatchTetsNP(self, index_t[:] t_indices, index_t[:] v_indices):
        """
        Get vertex indices of a list of triangles.

        Syntax::

            getBatchTetsNP(t_indices, v_indices)

        Arguments:
        numpy.array<index_t> t_indices
        numpy.array<index_t, length = len(t_indices) * 4> v_indices

        Return:
        None

        """
        self.ptrx().getBatchTetsNP( & t_indices[0], t_indices.shape[0], & v_indices[0], v_indices.shape[0])

    def getTriVerticesSetSizeNP(self, index_t[:] t_indices):
        """
        Return the size of a set with unique vertex indices of a list of triangles,
        preparation function for furture numpy data access.

        Syntax::

            getTriVerticesSetSizeNP(t_indices)

        Arguments:
        numpy.array<index_t> t_indices

        Return:
        int

        """
        return self.ptrx().getTriVerticesSetSizeNP( & t_indices[0], t_indices.shape[0])

    def getTetVerticesSetSizeNP(self, index_t[:] t_indices):
        """
        Return the size of a set with unique vertex indices of a list of tetrahedrons,
        preparation function for furture numpy data access.

        Syntax::

            getTetVerticesSetSizeNP(t_indices)

        Arguments:
        numpy.array<index_t> t_indices

        Return:
        index_t

        """
        return self.ptrx().getTetVerticesSetSizeNP( & t_indices[0], t_indices.shape[0])

    def getTriVerticesMappingSetNP(self, index_t[:] t_indices, index_t[:] t_vertices, index_t[:] v_set):
        """
        Get the vertex indices of a list of triangles.
        The vertex indices are reindexed, with their oringinal STEPS indices stored in a given array,
        whose size is provided by getTriVerticesSetSizeNP().

        Syntax::

            getTriVerticesMappingSetNP(t_indices, t_vertices, v_set)

        Arguments:
        numpy.array<index_t> t_indices
        numpy.array<index_t, length = length(t_indices) * 3> t_vertices
        numpy.array<index_t, length = getTriVerticesSetSizeNP(t_indices)> v_set

        Return:
        None

        """
        self.ptrx().getTriVerticesMappingSetNP( & t_indices[0], t_indices.shape[0], & t_vertices[0], t_vertices.shape[0], & v_set[0], v_set.shape[0])

    def getTetVerticesMappingSetNP(self, index_t[:] t_indices, index_t[:] t_vertices, index_t[:] v_set):
        """
        Get the vertex indices of a list of tetrahedrons.
        The vertex indices are reindexed, with their oringinal STEPS indices stored in a given array,
        whose size is provided by getTriVerticesSetSizeNP().

        Syntax::

            getTetVerticesMappingSetNP(t_indices, t_vertices, v_set)

        Arguments:
        numpy.array<index_t> t_indices
        numpy.array<index_t, length = length(t_indices) * 4> t_vertices
        numpy.array<index_t, length = getTriVerticesSetSizeNP(t_indices)> v_set

        Return:
        None

        """
        self.ptrx().getTetVerticesMappingSetNP( & t_indices[0], t_indices.shape[0], & t_vertices[0], t_vertices.shape[0], & v_set[0], v_set.shape[0])

    def genPointsInTet(self, index_t tidx, uint npnts, double[:] coords):
        """
        Generate npnts random point coordinates x,y,z within a tetraedron with index tidx, export it to NumPy array cords

        Syntax::

            genPointsInTet(t_idx, npnts, coords)

        Arguments:
        index_t tidx
        int npnts
        numpy.array<float, length = npnts * 3> coords

        Return:
        None

        """
        if not len(coords): return False
        return self.ptrx().genPointsInTet(tetrahedron_global_id(tidx), npnts, &coords[0], coords.shape[0])

    def genPointsInTri(self, index_t tidx, uint npnts, double[:] coords):
        """
        Generate npnts random point coordinates x,y,z within a triangle with index tidx, export it to NumPy array cords

        Syntax::

            genPointsInTri(t_idx, npnts, coords)

        Arguments:
        index_t tidx
        int npnts
        numpy.array<float, length = npnts * 3> coords

        Return:
        None

        """
        if not len(coords): return False
        return self.ptrx().genPointsInTri(triangle_global_id(tidx), npnts, &coords[0], coords.shape[0])

    def genTetVisualPointsNP(self, index_t[:] indices, uint[:] point_counts, double[:] coords):
        """
        For each tetrahedron index in indices, randomly generate a set of point coordinates x,y,z within the tetrahedron, where n is
        stored in point_counts. The number of points required to be generated for tetrahedron indices[i] is point_counts[i].
        All generated points are stored in cords.

        Syntax::

            genTetVisualPointsNP(indices, point_counts, coords)

        Arguments:
        numpy.array<index_t> indices
        numpy.array<uint, length = length(indices)> point_counts
        numpy.array<double, length = sum(point_counts) * 3> coords

        Return:
        None

        """
        if not len(coords): return False
        return self.ptrx().genTetVisualPointsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], &coords[0], coords.shape[0])

    def genTriVisualPointsNP(self, index_t[:] indices, uint[:] point_counts, double[:] coords):
        """
        For each triangle index in indices, randomly generate a set of point coordinates x,y,z within the triangle, where n is
        stored in point_counts. The number of points required to be generated for triangle indices[i] is point_counts[i].
        All generated points are stored in cords.

        Syntax::

            genTriVisualPointsNP(indices, point_counts, coords)

        Arguments:
        numpy.array<index_t> indices
        numpy.array<uint, length = length(indices)> point_counts
        numpy.array<double, length = sum(point_counts) * 3> coords

        Return:
        None

        """
        if not len(coords): return False
        return self.ptrx().genTriVisualPointsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], &coords[0], coords.shape[0])

    def getBatchTetVolsNP(self, index_t[:] indices, double[:] volumes):
        """
        Get the volumes of a list of tetrahedrons in indices and stored in volumes.

        Syntax::

            getBatchTetVolsNP(indices, volumes)

        Arguments:
        numpy.array<index_t> indices
        numpy.array<double, length = length(indices)> volumes

        Return:
        None

        """
        return self.ptrx().getBatchTetVolsNP(&indices[0], indices.shape[0], &volumes[0], volumes.shape[0])

    def getBatchTriAreasNP(self, index_t[:] indices, double[:] areas):
        """
        Get the areas of a list of triangles in indices and stored in areas.

        Syntax::

            getBatchTriAreasNP(indices, areas)

        Arguments:
        numpy.array<index_t> indices
        numpy.array<float, length = length(indices)> areas

        Return:
        None

        """
        return self.ptrx().getBatchTriAreasNP(&indices[0], indices.shape[0], &areas[0], areas.shape[0])

    def reduceBatchTetPointCountsNP(self, index_t[:] indices, uint[:] point_counts, double max_density):
        """
        Reduce the number of random point coordinates generated for each tetrahedron in indices so that the point density of the tetrahedron is below max_density. If the density is already below max_density for that tetrahedron, the count stored in point_counts is intacted.

        Syntax::

            reduceBatchTetPointCountsNP(indices, point_counts, max_density)

        Arguments:
        numpy.array<index_t> indices
        numpy.array<uint, length = length(indices)> point_counts
        float max_density

        Return:
        None

        """
        return self.ptrx().reduceBatchTetPointCountsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], max_density)

    def reduceBatchTriPointCountsNP(self, index_t[:] indices, uint[:] point_counts, double max_density):
        """
        Reduce the number of random point coordinates generated for each triangle in indices so that the point density of the triangle is below max_density. If the density is already below max_density for that triangle, the count stored in point_counts is intacted.

        Syntax::

            reduceBatchTriPointCountsNP(indices, point_counts, max_density)

        Arguments:
        numpy.array<index_t> indices
        numpy.array<uint, length = length(indices)> point_counts
        double max_density

        Return:
        None

        """
        return self.ptrx().reduceBatchTriPointCountsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], max_density)

    ## ROI related ##
    def addROI(self, str id, ElementType type, std.set[index_t] indices):
        """
        Add a Region of Interest data record with name id to the ROI dataset.
        The type of elements stored in the ROI data can be one of the follows:
        steps.geom.ELEM_VERTEX, steps.geom.ELEM_TET, steps.geom.ELEM_TRI, steps.geom.ELEM_UNDEFINED.

        Syntax::

            addROI(id, type, indices)

        Arguments:
        string id
        ElementType type
        set<index_t> indices

        Return:
        None

        """
        self.ptrx().addROI(to_std_string(id), type, indices)

    def removeROI(self, str id):
        """
        Remove a Region of Interest data record with name id.

        Syntax::

            removeROI(id)

        Arguments:
        string id

        Return:
        None

        """
        self.ptrx().removeROI(to_std_string(id))

    def replaceROI(self, str id, ElementType type, std.set[index_t] indices):
        """
        Replace a Region of Interest data record with name id with new data.

        Syntax::

            replaceROI(id, type, indices)

        Arguments:
        string id
        ElementType type
        set<index_t> indices

        Return:
        None

        """
        self.ptrx().replaceROI(to_std_string(id), type, indices)

    def getROIType(self, str id):
        """
        Get the element type of a Region of Interest data record with name id.

        Syntax::

            getROIType(id)

        Arguments:
        string id

        Return:
        ElementType

        """
        return self.ptrx().getROIType(to_std_string(id))

    def getROIData(self, str id):
        """
        Get the stored data of a Region of Interest data record with name id.

        Syntax::

            getROIData(id)

        Arguments:
        string id

        Return:
        list<index_t>

        """
        return self.ptrx().getROIData(to_std_string(id))

    def getROIDataSize(self, str id):
        """
        Get the number of elements stored in a Region of Interest data record with name id.

        Syntax::

            getROIDataSize(id)

        Arguments:
        string id

        Return:
        int

        """
        return self.ptrx().getROIDataSize(to_std_string(id))

    def getNROIs(self):
        """
        Get the number of Region of Interest data stored in the ROI dataset.

        Syntax::

            getNROIs()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().getNROIs()

    def getROI(self, str id):
        """
        Get a Region of Interest data record with name id.

        Syntax::

            getROI(id)

        Arguments:
        string id

        Return:
        ROISet

        """
        cdef ROISet rois = self.ptrx().getROI(to_std_string(id))
        return _py_ROISet(rois.type, rois.indices)

    def getROIArea(self, str ROI_id):
        """
        Get summed area of all triangles stored in a triangular ROI.

        Syntax::

            getROIArea(id)

        Arguments:
        string id

        Return:
        float

        """
        return self.ptrx().getROIArea(to_std_string(ROI_id))

    def getROIVol(self, str ROI_id):
        """
        Get summed area of all tetrahedrons stored in a tetrahedral ROI.

        Syntax::

            getROIVol(id)

        Arguments:
        string id

        Return:
        float

        """
        return self.ptrx().getROIVol(to_std_string(ROI_id))

    def getAllROINames(self):
        """
        Get a list of the names of all Region of Interest data stored in ROI dataset.

        Syntax::

            getAllROINames()

        Arguments:
        None

        Return:
        list<string>

        """
        all_names = self.ptrx().getAllROINames()
        return [from_std_string(s) for s in all_names]

    def checkROI(self, str id, ElementType type, uint count=0, bool warning=True):
        """
        Check if an ROI enquire is valid.

        Syntax::

            checkROI(id, type, count=0, warning=true)

        Arguments:
        string id
        ElementType type
        int count
        bool warning

        Return:
        bool

        """
        return self.ptrx().checkROI(to_std_string(id), type, count, warning )

    def getROITetBarycentres(self, str ROI_id):
        """
        Get barycenters of elements stored in a tetrahedral ROI.

        This function has been renamed by replacing "centre" with "center".
        It will be deprecated in the next release.

        Syntax::

            getROITetBarycenters(ROI_id)

        Arguments:
        string ROI_id

        Return:
        list<float>

        """
        warnings.warn('This function has been renamed by replacing "centre" with "center".\n')
        warnings.warn('The original function will be deprecated in the next release.\n')
        return self.ptrx().getROITetBarycenters(to_std_string(ROI_id))

    def getROITetBarycenters(self, str ROI_id):
        """
        Get barycenters of elements stored in a tetrahedral ROI.

        Syntax::

            getROITetBarycenters(ROI_id)

        Arguments:
        string ROI_id

        Return:
        list<float>

        """
        return self.ptrx().getROITetBarycenters(to_std_string(ROI_id))

    def getROITetBarycentresNP(self, str ROI_id, double[:] centers):
        """
        Get barycenters of elements stored in a tetrahedral ROI and write to a NumPy array centers.
        The size of centers should be the same as the number of elements stored in the ROI.

        This function has been renamed by replacing "centre" with "center".
        It will be deprecated in the next release.

        Syntax::

            getROITetBarycentersNP(ROI_id, centers)

        Arguments:
        string ROI_id
        numpy.array<float> centers

        Return:
        None

        """
        warnings.warn('This function has been renamed by replacing "centre" with "center".\n')
        warnings.warn('The original function will be deprecated in the next release.\n')
        return self.ptrx().getROITetBarycentersNP(to_std_string(ROI_id), &centers[0], centers.shape[0])

    def getROITetBarycentersNP(self, str ROI_id, double[:] centers):
        """
        Get barycenters of elements stored in a tetrahedral ROI and write to a NumPy array centers.
        The size of centers should be the same as the number of elements stored in the ROI.

        Syntax::

            getROITetBarycentersNP(ROI_id, centers)

        Arguments:
        string ROI_id
        numpy.array<float> centers

        Return:
        None

        """

        return self.ptrx().getROITetBarycentersNP(to_std_string(ROI_id), &centers[0], centers.shape[0])

    def getROITriBarycentres(self, str ROI_id):
        """
        Get barycenters of elements stored in a triangular ROI.

        This function has been renamed by replacing "centre" with "center".
        It will be deprecated in the next release.

        Syntax::

            getROITriBarycenters(ROI_id)

        Arguments:
        string ROI_id

        Return:
        list<float>

        """
        warnings.warn('This function has been renamed by replacing "centre" with "center".\n')
        warnings.warn('The original function will be deprecated in the next release.\n')
        return self.ptrx().getROITriBarycenters(to_std_string(ROI_id))

    def getROITriBarycenters(self, str ROI_id):
        """
        Get barycenters of elements stored in a triangular ROI.

        Syntax::

            getROITriBarycenters(ROI_id)

        Arguments:
        string ROI_id

        Return:
        list<float>

        """

        return self.ptrx().getROITriBarycenters(to_std_string(ROI_id))

    def getROITriBarycentresNP(self, str ROI_id, double[:] centers):
        """
        Get barycenters of elements stored in a triangular ROI and write to a NumPy array centers.
        The size of centers should be the same as the number of elements stored in the ROI.

        This function has been renamed by replacing "centre" with "center".
        It will be deprecated in the next release.

        Syntax::

            getROITriBarycentersNP(ROI_id, centers)

        Arguments:
        string ROI_id
        numpy.array<float> centers

        Return:
        None

        """
        warnings.warn('This function has been renamed by replacing "centre" with "center".\n')
        warnings.warn('The original function will be deprecated in the next release.\n')
        return self.ptrx().getROITriBarycentersNP(to_std_string(ROI_id), &centers[0], centers.shape[0])

    def getROITriBarycentersNP(self, str ROI_id, double[:] centers):
        """
        Get barycenters of elements stored in a triangular ROI and write to a NumPy array centers.
        The size of centers should be the same as the number of elements stored in the ROI.

        Syntax::

            getROITriBarycenters(ROI_id)

        Arguments:
        string ROI_id

        Return:
        list<float>

        """
        return self.ptrx().getROITriBarycenters(to_std_string(ROI_id))

    def getROIVertices(self, str ROI_id):
        """
        Get coordinates of elements stored in a vertices ROI.

        Syntax::

            getROIVertices(ROI_id)

        Arguments:
        string ROI_id

        Return:
        list<float>

        """
        return self.ptrx().getROIVertices(to_std_string(ROI_id))

    def getROIVerticesNP(self, str ROI_id, double[:] coordinates):
        """
        Get coordinates of elements stored in a vertices ROI and write to a NumPy array coordinates.
        The size of coordinates should be the same as the number of elements stored in the ROI.

        Syntax::

            getROIVerticesNP(ROI_id, coordinates)

        Arguments:
        string ROI_id
        numpy.array<float> coordinates

        Return:
        None

        """
        return self.ptrx().getROIVerticesNP(to_std_string(ROI_id), &coordinates[0], coordinates.shape[0])

    def getROITris(self, str ROI_id):
        """
        Get vertices of elements stored in a triangular ROI.

        Syntax::

            getROITris(ROI_id)

        Arguments:
        string ROI_id

        Return:
        list<index_t>

        """
        return self.ptrx().getROITris(to_std_string(ROI_id))

    def getROITrisNP(self, str ROI_id, index_t[:] v_indices):
        """
        Get vertices of elements stored in a triangular ROI and write to a NumPy array v_indices.
        The size of v_indices should be 3 * the number of elements stored in the ROI.

        Syntax::

            getROITrisNP(ROI_id, v_indices)

        Arguments:
        string ROI_id
        numpy.array<index_t> v_indices

        Return:
        None

        """
        return self.ptrx().getROITrisNP(to_std_string(ROI_id), &v_indices[0], v_indices.shape[0])

    def getROITets(self, str ROI_id):
        """
        Get vertices of elements stored in a tetrahedral ROI.

        Syntax::

            getROITets(ROI_id)

        Arguments:
        string ROI_id

        Return:
        list<index_t>

        """
        return self.ptrx().getROITets(to_std_string(ROI_id))

    def getROITetsNP(self, str ROI_id, index_t[:] v_indices):
        """
        Get vertices of elements stored in a tetrahedral ROI and write to a NumPy array v_indices.
        The size of v_indices should be 3 * the number of elements stored in the ROI.

        Syntax::

            getROITetsNP(ROI_id, v_indices)

        Arguments:
        string ROI_id
        numpy.array<index_t> v_indices

        Return:
        None

        """
        return self.ptrx().getROITetsNP(to_std_string(ROI_id), &v_indices[0], v_indices.shape[0])

    def getROITriVerticesSetSizeNP(self, str ROI_id):
        """
        Add all vertex indices of a list of triangles in a ROI to a set and return its size.

        Syntax::

            getROITriVerticesSetSizeNP(ROI_id)

        Arguments:
        string ROI_id

        Return:
        int

        """
        return self.ptrx().getROITriVerticesSetSizeNP(to_std_string(ROI_id))

    def getROITetVerticesSetSizeNP(self, str ROI_id):
        """
        Add all vertex indices of a list of tetrahedrons in a ROI to a set and return its size.

        Syntax::

            getROITetVerticesSetSizeNP(ROI_id)

        Arguments:
        string ROI_id

        Return:
        int

        """
        return self.ptrx().getROITetVerticesSetSizeNP(to_std_string(ROI_id))

    def getROITriVerticesMappingSetNP(self, str ROI_id, index_t[:] t_vertices, index_t[:] v_set):
        """
        Add all vertex indices of a list of triangles in a ROI to a set and write it to a NumPy array v_set.
        For each of the triangle, t_vertices records the positions of its vertices in v_set.
        i.e. For the i triangle in the ROI, the STEPS indices of its vertices are
        v_set[t_vertices[3*i]], v_set[t_vertices[3*i + 1]], v_set[t_vertices[3*i + 2]]

        Syntax::

            getROITriVerticesMappingSetNP(ROI_id, t_vertices, v_set)

        Arguments:
        string ROI_id
        numpy.array<index_t> v_indices
        numpy.array<index_t> v_set

        Return:
        None

        """
        return self.ptrx().getROITriVerticesMappingSetNP(to_std_string(ROI_id), &t_vertices[0], t_vertices.shape[0], &v_set[0], v_set.shape[0])

    def getROITetVerticesMappingSetNP(self, str ROI_id, index_t[:] t_vertices, index_t[:] v_set):
        """
        Add all vertex indices of a list of tetrahedrons in a ROI to a set and write it to a NumPy array v_set.
        For each of the tetrahedron, t_vertices records the positions of its vertices in v_set.
        i.e. For the i tetrahedron in the ROI, the STEPS indices of its vertices are
        v_set[t_vertices[4*i]], v_set[t_vertices[4*i + 1]], v_set[t_vertices[4*i + 2]], v_set[t_vertices[4*i + 3]]

        Syntax::

            getROITetVerticesMappingSetNP(ROI_id, t_vertices, v_set)

        Arguments:
        string ROI_id
        numpy.array<index_t> v_indices
        numpy.array<index_t> v_set

        Return:
        None

        """
        return self.ptrx().getROITetVerticesMappingSetNP(to_std_string(ROI_id), &t_vertices[0], t_vertices.shape[0], &v_set[0], v_set.shape[0])

    def genROITetVisualPointsNP(self, str ROI_id, uint[:] point_counts, double[:] coords):
        """
        For each tetrahedron index in a ROI, randomly generate a set of point coordinates x,y,z within the tetrahedron, where n is
        stored in point_counts. The number of points required to be generated for tetrahedron i in the ROI is point_counts[i].
        All generated points are stored in cords.

        Syntax::

            genTetVisualPointsNP(ROI_id, point_counts, coords)

        Arguments:
        string ROI_id
        numpy.array<uint> point_counts
        numpy.array<double> coords

        Return:
        None

        """
        return self.ptrx().genROITetVisualPointsNP(to_std_string(ROI_id), &point_counts[0], point_counts.shape[0], &coords[0], coords.shape[0])

    def genROITriVisualPointsNP(self, str ROI_id, uint[:] point_counts, double[:] coords):
        """
        For each triangle index in a ROI, randomly generate a set of point coordinates x,y,z within the triangle, where n is
        stored in point_counts. The number of points required to be generated for triangle i in the ROI is point_counts[i].
        All generated points are stored in cords.

        Syntax::

            genROITriVisualPointsNP(ROI_id, point_counts, coords)

        Arguments:
        string ROI_id
        numpy.array<uint> point_counts
        numpy.array<double> coords

        Return:
        None

        """
        return self.ptrx().genROITriVisualPointsNP(to_std_string(ROI_id), &point_counts[0], point_counts.shape[0], &coords[0], coords.shape[0])

    def getROITetVolsNP(self, str ROI_id, double[:] volumes):
        """
        Get the volumes of a list of tetrahedrons in a ROI and stored in volumes.

        Syntax::

            getROITetVolsNP(indices, volumes)

        Arguments:
        string ROI_id
        numpy.array<float> volumes

        Return:
        None

        """
        return self.ptrx().getROITetVolsNP(to_std_string(ROI_id), &volumes[0], volumes.shape[0])

    def getROITriAreasNP(self, str ROI_id, double[:] areas):
        """
        Get the areas of a list of triangles in a ROI and stored in areas.

        Syntax::

            getROITriAreasNP(indices, areas)

        Arguments:
        string ROI_id
        numpy.array<float> areas

        Return:
        None

        """
        return self.ptrx().getROITriAreasNP(to_std_string(ROI_id), &areas[0], areas.shape[0])

    def reduceROITetPointCountsNP(self, str ROI_id, uint[:] point_counts, double max_density):
        """
        Reduce the number of random point coordinates generated for each tetrahedron in a ROI so that the point density of the tetrahedron is below max_density. If the density is already below max_density for that tetrahedron, the count stored in point_counts is intacted.

        Syntax::

            reduceROITetPointCountsNP(indices, point_counts, max_density)

        Arguments:
        string ROI_id
        numpy.array<uint> point_counts
        double max_density

        Return:
        None

        """
        return self.ptrx().reduceROITetPointCountsNP(to_std_string(ROI_id), &point_counts[0], point_counts.shape[0], max_density)

    def reduceROITriPointCountsNP(self, str ROI_id, uint[:] point_counts, double max_density):
        """
        Reduce the number of random point coordinates generated for each triangle in a ROI so that the point density of the triangle is below max_density. If the density is already below max_density for that triangle, the count stored in point_counts is intacted.

        Syntax::

            reduceROITriPointCountsNP(ROI_id, point_counts, max_density)

        Arguments:
        string ROI_id
        numpy.array<uint, length = length(indices)> point_counts
        double max_density

        Return:
        None

        """
        return self.ptrx().reduceROITriPointCountsNP(to_std_string(ROI_id), &point_counts[0], point_counts.shape[0], max_density)

    def intersect(self, double[:, :] point_coords, int sampling=-1):
        """
        Computes the intersection of line segment(s) given the vertices with the current mesh

        Args:
            points: A 2-D numpy array (/memview), where each position contains the 3 point
                    coordinates
            int sampling: not specified or sampling < 1 --> use deterministic method
                          sampling > 0 --> use montecarlo method with sampling points

        Returns:
            A list where each position contains the list of intersected tets (and respective
            intersection ratio) of each line segment.
        """
        if (point_coords.strides[0] != 24 or point_coords.strides[1] != 8):
            raise Exception("Wrong memory layout for point_coords, np array should be [pts,3] and row major")
        return self.ptrx().intersect(&point_coords[0][0], point_coords.shape[0], sampling)

    def intersectIndependentSegments(self, double[:, :] point_coords, int sampling=-1):
        """
        Similar to the intersect method but here we deal with independent segments, i.e.
        every two points we have a segment not related to previous or following ones.
        E.g. seg0 = (points[0], points[1]), seg1 = (points[2], points[3]), etc.

        Args:
            points: A 2-D numpy array (/memview), where each position contains the 3 point
                    coordinates
            int sampling: not specified or sampling < 1 --> use deterministic method
                          sampling > 0 --> use montecarlo method with sampling points

        Returns:
            A list where each position contains the list of intersected tets (and respective
            intersection ratio) of each line segment.
        """
        if (point_coords.strides[0] != 24 or point_coords.strides[1] != 8):
            raise Exception("Wrong memory layout for point_coords, np array should be [pts,3] and row major")
        return self.ptrx().intersectIndependentSegments(&point_coords[0][0], point_coords.shape[0], sampling)

    @staticmethod
    cdef _py_Tetmesh from_ptr(Tetmesh *ptr):
        if (ptr == NULL):
            return None
        cdef _py_Tetmesh obj = _py_Tetmesh.__new__(_py_Tetmesh)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Tetmesh from_ref(const Tetmesh &ref):
        return _py_Tetmesh.from_ptr(<Tetmesh*>&ref)

    ## properties ##
    nverts = property(countVertices, doc="Number of vertices in the mesh.")
    ntris  = property(countTris, doc="Number of triangles in the mesh.")
    ntets  = property(countTets, doc="Number of tetrahedrons in the mesh.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TmComp(_py_Comp):
    "Python wrapper class for TmComp"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[TmComp] _autodealoc
    cdef TmComp *ptrx(self):
        return <TmComp*> self._ptr

    def __init__(self, str id, _py_Tetmesh container, std.vector[index_t] tets):
        """
        Construction::

            tmcomp = steps.geom.Comp(id, container, tets)

        Construct a TmComp object with identifier string id and assign container
        as the parent Tetmesh container. Set the group of tetrahedrons that describe
        this compartment with tets.

        Arguments:
        string id
        steps.geom.Tetmesh container
        list<index_t> tets
        """
        self._ptr = new TmComp(to_std_string(id), deref(container.ptrx()), tets)

    def setVol(self, double vol):
        """Obsolete"""
        self.ptrx().setVol(vol)

    def getAllTetIndices(self, ):
        """
        Returns a list of all tetrahedrons assigned to the compartment.

        Syntax::

            getAllTetIndices()

        Arguments:
        None

        Return:
        list<int>

        """
        return self.ptrx().getAllTetIndices()

    def countTets(self, ):
        """
        Returns the number of tetrahedrons assigned to the compartment.

        Syntax::

            countTets()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().countTets()

    def isTetInside(self, std.vector[index_t] tets):
        """
        Returns a list of Booleans describing if tetrahedrons tets are
        assigned to the compartment.

        Syntax::

            isTetInside(tets)

        Arguments:
        list<index_t> tets

        Return:
        list<bool, length = length(tets)>

        """
        return self.ptrx().isTetInside(tets)

    def getBoundMin(self, ):
        """
        Returns the minimal Cartesian coordinate of the rectangular bounding box
        of the compartment.

        Syntax::

            getBoundMin()

        Arguments:
        None

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMin()

    def getBoundMax(self, ):
        """
        Returns the maximal Cartesian coordinate of the rectangular bounding box
        of the compartment.

        Syntax::

            getBoundMax()

        Arguments:
        None

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMax()

    @staticmethod
    cdef _py_TmComp from_ptr(TmComp *ptr):
        if (ptr == NULL):
            return None
        cdef _py_TmComp obj = _py_TmComp.__new__(_py_TmComp)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_TmComp from_ref(const TmComp &ref):
        return _py_TmComp.from_ptr(<TmComp*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[Comp*] vec):
        return [ _py_TmComp.from_ptr(<TmComp*?>elem) for elem in vec ]

    ## properties ##
    tets = property(getAllTetIndices, doc="List of indices of tetrahedrons associated to the compartment.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TmPatch(_py_Patch):
    "Python wrapper class for TmPatch"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[TmPatch] _autodealoc
    cdef TmPatch *ptrx(self):
        return <TmPatch*> self._ptr

    def __init__(self, str id, _py_Tetmesh container, std.vector[index_t] tris, _py_Comp icomp, _py_Comp ocomp=None):
        """
        Construction::

            tmpatch = steps.geom.Comp(id, container, tris, icomp, ocomp = None)

        Construct a TmPatch object with identifier string id and assign container
        as the parent geometry container. Set the collection of triangles in
        the patch to tris and assign icomp as the "inner" compartment
        and assign also ocomp as the "outer" compartment (if required).

        Arguments:
        string id
        steps.geom.Tetmesh container
        list<index_t> tris
        steps.geom.TmComp icomp
        steps.geom.TmComp ocomp (default = None)
        """
        self._ptr = new TmPatch(to_std_string(id), deref(container.ptrx()), tris, deref(icomp.ptr()), ocomp.ptr() if ocomp else NULL)

    def isTriInside(self, std.vector[index_t] tris):
        """
        Returns a list of Booleans describing if triangles tris are
        assigned to the patch.

        Syntax::

            isTriInside(tris)

        Arguments:
        list<index_t> tris

        Return:
        list<bool, length = length(tris)>

        """
        return self.ptrx().isTriInside(tris)

    def getAllTriIndices(self, ):
        """
        Returns a list of indices of all triangles assigned to the patch.

        Syntax::

            getAllTriIndices()

        Arguments:
        None

        Return:
        list<index_t>

        """
        return self.ptrx().getAllTriIndices()

    def getBoundMin(self, ):
        """
        Returns the minimal Cartesian coordinate of the rectangular bounding box
        of the compartment.

        Syntax::

            getBoundMin()

        Arguments:
        None

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMin()

    def getBoundMax(self, ):
        """
        Returns the maximal Cartesian coordinate of the rectangular bounding box
        of the compartment.

        Syntax::

            getBoundMax()

        Arguments:
        None

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMax()

    def getAllEndocyticZones(self, ):
        """
        Get all endocytic zones declared in the patch


        :rtype: List[EndocyticZone]
        """
        cdef std.vector[EndocyticZone*] vec = self.ptrx().getAllEndocyticZones()
        return [_py_EndocyticZone.from_ptr(ptr) for ptr in vec]

    @staticmethod
    cdef _py_TmPatch from_ptr(TmPatch *ptr):
        if (ptr == NULL):
            return None
        cdef _py_TmPatch obj = _py_TmPatch.__new__(_py_TmPatch)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_TmPatch from_ref(const TmPatch &ref):
        return _py_TmPatch.from_ptr(<TmPatch*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[Patch*] patches):
        return [ _py_TmPatch.from_ptr(<TmPatch*?>p) for p in patches ]

    @staticmethod
    cdef list flat_set2list(flat_set[Patch*] patches):
        return [ _py_TmPatch.from_ptr(<TmPatch*?>p) for p in patches ]

    ## properties ##
    tris = property(getAllTriIndices, doc="List of indices of triangles associated to the patch.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_EndocyticZone(_py__base):
    "Python wrapper class for EndocyticZone"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[EndocyticZone] _autodealoc
    cdef EndocyticZone *ptr(self):
        return <EndocyticZone*> self._ptr

    def __init__(self, str id, _py_Patch patch, std.vector[index_t] tris):
        """
        Constructor

        :param id: ID of the endocytic zone
        :type id: str
        :param patch: patch in which the zone is defined
        :type patch: TmPatch
        :param tris: A sequence of triangles (by index)
        :type tris: List[index_t]

        """
        self._ptr = new EndocyticZone(to_std_string(id), deref(<TmPatch *>(patch.ptr())), tris)

    def getID(self, ):
        """
        Return the endocytic zone id.

        :returns: ID of the endocytic zone.


        :rtype: str
        """
        return from_std_string(self.ptr().getID())

    def getPatch(self, ):
        """
        Return a pointer to the patch container object.

        :returns: The parent patch


        :rtype: TmPatch
        """
        return _py_TmPatch.from_ptr(&self.ptr().getPatch())

    def getAllTriIndices(self, ):
        """
        Return all triangles (by index) in the endocytic zone.

        :returns: List of indices of triangles.


        :rtype: List[triangle_global_id]
        """
        return [_index.get() for _index in self.ptr().getAllTriIndices()]

    @staticmethod
    cdef _py_EndocyticZone from_ptr(EndocyticZone *ptr):
        if (ptr == NULL):
            return None
        cdef _py_EndocyticZone obj = _py_EndocyticZone.__new__(_py_EndocyticZone)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_EndocyticZone from_ref(const EndocyticZone &ref):
        return _py_EndocyticZone.from_ptr(<EndocyticZone*>&ref)



# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Memb(_py__base):
    "Python wrapper class for Memb"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Memb] _autodealoc
    cdef Memb *ptr(self):
        return <Memb*> self._ptr

    def __init__(self, str id, _py_Tetmesh container, list patches, bool verify=False, uint opt_method=1, double search_percent=100.0, str opt_file_name="", list supplementary_comps=[]):
        """
        Construction:

        memb = steps.geom.Memb(id, container, patches, verify = False, opt_method=1, opt_file_name = '', supplementary_comps=[])

        Construct a Memb object with identifier string id and assign container
        as the parent geometry container. Set the collection of triangles in
        the membrane to those that belong to all TmPatches in patches. Perform some checks on
        suitability of membrane if verify is True: these checks will print
        warnings if membrane forms an open surface or if any triangle is found to have
        more than 3 neighbours. Specify optimization method with opt_method (default = 1):
        1 = principal axis ordering (quick to set up but usually results in slower simulation than method 2).
        2 = breadth first search (can be time-consuming to set up, but usually faster simulation.
        If 2:breadth first search is chosen then argument search_percent can specify the number of starting points to search for the
        lowest bandwidth.
        If a filename (with full path) is given in optional argument opt_file_name the membrane optimization will be loaded from file,
        which was saved previously for this membrane with solver method steps.solver.Tetexact.saveMembOpt()

        By default STEPS only adds the inner compartments of the patches to the conduction volume, if other compartments should also
        be added to the conduction volume, they should be supplied through the supplementary_comps keyword parameter. 

        Arguments:
        string id
        steps.geom.Tetmesh container
        list<steps.geom.TmPatch> patches
        bool verify (default = False)
        int opt_method (default = 1)
        float search_percent (default=100)
        string opt_file_name (default = '')
        list<steps.geom.TmComp> supplementary_comps
        """
        cdef std.vector[TmPatch*] _patches
        _patches.reserve(len(patches))
        for elem in patches:
            _patches.push_back( (<_py_TmPatch>elem).ptrx() )
        cdef std.vector[TmComp*] _compartments
        _compartments.reserve(len(supplementary_comps))
        for elem in supplementary_comps:
            _compartments.push_back( (<_py_TmComp>elem).ptrx() )
        self._ptr = new Memb(to_std_string(id), deref(container.ptrx()), _patches, _compartments, verify, opt_method, search_percent, to_std_string(opt_file_name))

    def getContainer(self, ):
        """
        Returns a reference to the parent steps.geom.Tetmesh container object.

        Syntax::

            getTetmesh()

        Arguments:
        None

        Return:
        steps.tetmesh.Tetmesh

        """
        return _py_Tetmesh.from_ptr(&self.ptr().getContainer())

    def getID(self, ):
        """
        Get the identifier string of the membrane.

        Syntax::

            getID()

        Arguments:
        None

        Return:
        string

        """
        return from_std_string(self.ptr().getID())

    def isTriInside(self, std.vector[index_t] tri):
        """
        Returns a list of Booleans describing if triangles tris are
        assigned to the membrane.

        Syntax::

            isTriInside(tris)

        Arguments:
        list<index_t> tris

        Return:
        list<bool, length = length(tris)>

        """
        return self.ptr().isTriInside(tri)

    def getAllTriIndices(self, ):
        """
        Returns a list of indices of all triangles assigned to the membrane.

        Syntax::

            getAllTriIndices()

        Arguments:
        None

        Return:
        list<index_t>

        """
        return self.ptr().getAllTriIndices()

    def countTris(self, ):
        """
        Returns the number of triangles assigned to the membrane.

        Syntax:

        countTris()

        Arguments:
        None

        Return:
        int


        """
        return self.ptr().countTris()

    def getAllVolTetIndices(self, ):
        """
        Returns a list of indices of all tetrahedrons assigned to the conduction volume.

        Syntax::

            getAllVolTetIndices()

        Arguments:
        None

        Return:
        list<index_t>

        """
        return self.ptr().getAllVolTetIndices()

    def countVolTets(self, ):
        """
        Returns the number of tetrahedrons assigned to the conduction volume.

        Syntax::

            countVolTets()

        Arguments:
        None

        Return:
        uint


        """
        return self.ptr().countVolTets()

    def getAllVirtTriIndices(self, ):
        """
        Returns a list of all virtual triangles for the membrane forming a closed surface.

        Syntax:

        getAllVirtTris()

        Arguments:
        None

        Return:
        list<index_t>


        """
        return self.ptr().getAllVirtTriIndices()

    def countVirtTris(self, ):
        """
        Returns the number of virtual triangles for the membrane forming a closed surface.

        Syntax:

        countTris()

        Arguments:
        None

        Return:
        uint


        """
        return self.ptr().countVirtTris()

    def getAllVertIndices(self, ):
        """
        Returns a list of all vertices in the conduction volume.

        Syntax:

        getAllVertices()

        Arguments:
        None

        Return:
        list<index_t>

        """
        return self.ptr().getAllVertIndices()

    def countVerts(self, ):
        """
        Returns the number of vertices in the conduction volume and membrane surface.

        Syntax:

        countVertices()

        Arguments:
        None

        Returns:
        uint


        """
        return self.ptr().countVerts()

    def open(self, ):
        """
        Returns whether a membrane is open or not, that is whether it contains any holes or not.

        Syntax:

        open()

        Arguments:
        None

        Returns:
        bool


        """
        return self.ptr().open()

    @staticmethod
    cdef _py_Memb from_ptr(Memb *ptr):
        if (ptr == NULL):
            return None
        cdef _py_Memb obj = _py_Memb.__new__(_py_Memb)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Memb from_ref(const Memb &ref):
        return _py_Memb.from_ptr(<Memb*>&ref)

    ## properties ##
    tris = property(getAllTriIndices, doc="List of indices of triangles associated to the membrane.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_DiffBoundary(_py__base):
    "Python wrapper class for DiffBoundary"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[DiffBoundary] _autodealoc
    cdef DiffBoundary *ptr(self):
        return <DiffBoundary*> self._ptr

    def __init__(self, str id, _py_Tetmesh container, std.vector[index_t] tris):
        """
        Construction::

            diffb = steps.geom.DiffBoundary(id, container, tris)

        Construct a DiffBoundary object with identifier string id and assign container
        as the parent geometry container, described by group of triangles tris.

        Arguments:
        string id
        steps.geom.Tetmesh container
        list<index_t> tris
        """
        self._ptr = new DiffBoundary(to_std_string(id), deref(container.ptrx()), tris)

    def getID(self, ):
        """
        Get the identifier string of the diffusion boundary.

        Syntax::

            getID()

        Arguments:
        None

        Return:
        string

        """
        return from_std_string(self.ptr().getID())

    def setID(self, str id):
        """
        Set the identifier string of the diffusion boundary.

        Syntax::

            setID(name)

        Arguments:
        string name

        Return:
        None

        """
        self.ptr().setID(to_std_string(id))

    def getContainer(self, ):
        """
        Returns a reference to the parent steps.tetmesh.Tetmesh container object.

        Syntax::

            getContainer()

        Arguments:
        None

        Return:
        steps.tetmesh.Tetmesh

        """
        return _py_Tetmesh.from_ptr(&self.ptr().getContainer())

    def isTriInside(self, std.vector[index_t] tri):
        """
        Returns a list of Booleans describing if triangles tris are
        assigned to the diffusion boundary.

        Syntax::

            isTriInside(tris)

        Arguments:
        list<index_t> tris

        Return:
        list<bool, length = length(tris)>

        """
        return (self.ptr().isTriInside(tri))

    def getAllTriIndices(self, ):
        """
        Returns a list of indices of all triangles assigned to the diffusion boundary.

        Syntax::

            getAllTriIndices()

        Arguments:
        None

        Return:
        list<index_t>

        """
        return self.ptr().getAllTriIndices()

    def getComps(self, ):
        """
        Returns a list of the two compartments this diffusion boundary connects.

        Syntax::

            getComps()

        Arguments:
        None

        Return:
        list<steps::wm::Comp, length = 2>

        """
        return _py_Comp.vector2list(self.ptr().getComps())

    @staticmethod
    cdef _py_DiffBoundary from_ptr(DiffBoundary *ptr):
        if (ptr == NULL):
            return None
        cdef _py_DiffBoundary obj = _py_DiffBoundary.__new__(_py_DiffBoundary)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_DiffBoundary from_ref(const DiffBoundary &ref):
        return _py_DiffBoundary.from_ptr(<DiffBoundary*>&ref)

    ## properties ##
    tris      = property(getAllTriIndices, doc="List of indices of triangles associated to the diffusion boundary.")
    id        = property(getID, setID, doc="Identifier string of the diffusion boundary.")
    container = property(getContainer, doc="Reference to parent steps.tetmesh.Tetmesh container.")
    comps     = property(getComps, doc="Reference to two steps.tetmesh.Comp compartments connected by this diffusion boundary.")

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_SDiffBoundary(_py__base):
    "Python wrapper class for SDiffBoundary"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[SDiffBoundary] _autodealoc
    cdef SDiffBoundary *ptr(self):
        return <SDiffBoundary*> self._ptr

    def __init__(self, str id, _py_Tetmesh container, std.vector[index_t] bars, list patches):
        """
        Construction::

            sdiffb = steps.geom.SDiffBoundary(id, container, bars, patches)

        Construct a SDiffBoundary object with identifier string id and assign container
        as the parent geometry container, described by group of bars. Specify the patches to be
        connected by this surface diffusion boundary to avoid potential ambiguity.

        Arguments:
        string id
        steps.geom.Tetmesh container
        list<index_t> bars
        list<steps.geom.TmPatch> (length 2) patches
        """
        cdef std.vector[TmPatch*] _patches
        for elem in patches:
            _patches.push_back( (<_py_TmPatch>elem).ptrx() )
        self._ptr = new SDiffBoundary(to_std_string(id), deref(container.ptrx()), bars, _patches)

    def getID(self, ):
        """
        Get the identifier string of the surface diffusion boundary.

        Syntax::

            getID()

        Arguments:
        None

        Return:
        string

        """
        return from_std_string(self.ptr().getID())

    def setID(self, str id):

        """
        Set the identifier string of the surface diffusion boundary.

        Syntax::

            setID(name)

        Arguments:
        string name

        Return:
        None

        """

        self.ptr().setID(to_std_string(id))

    def getContainer(self, ):
        """
        Get a reference to the parent steps.tetmesh.Tetmesh container object.

        Syntax::

            getContainer()

        Arguments:
        None

        Return:
        steps.tetmesh.Tetmesh

        """
        return _py_Tetmesh.from_ptr(&self.ptr().getContainer())

    def isBarInside(self, std.vector[index_t] bars):
        """
        Returns a list of Booleans describing if bars are
        assigned to the surface diffusion boundary.

        Syntax::

            isBarInside(bars)

        Arguments:
        list<index_t> bars

        Return:
        list<bool, length = length(bars)>

        """
        return (self.ptr().isBarInside(bars))

    def getAllBarIndices(self, ):
        """
        Returns a list of indices of all bars assigned to the surface diffusion boundary.

        Syntax::

            getAllBarIndices()

        Arguments:
        None

        Return:
        list<int>

        """
        return self.ptr().getAllBarIndices()

    def getPatches(self, ):
        """
        Returns a list of the two patches this surface diffusion boundary connects.

        Syntax::

            getPatches()

        Arguments:
        None

        Return:
        list<steps::wm::Patch, length = 2>

        """
        return _py_Patch.vector2list(self.ptr().getPatches())

    @staticmethod
    cdef _py_SDiffBoundary from_ptr(SDiffBoundary *ptr):
        if (ptr == NULL):
            return None
        cdef _py_SDiffBoundary obj = _py_SDiffBoundary.__new__(_py_SDiffBoundary)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_SDiffBoundary from_ref(const SDiffBoundary &ref):
        return _py_SDiffBoundary.from_ptr(<SDiffBoundary*>&ref)

    ## properties ##
    bars      = property(getAllBarIndices, doc="List of indices of bars associated to the surface diffusion boundary.")
    id        = property(getID, setID, doc="Identifier string of the surface diffusion boundary.")
    container = property(getContainer, doc="Reference to parent steps.tetmesh.Tetmesh container.")
    patches   = property(getPatches, doc="Reference to two steps.tetmesh.Patch patches connected by this surface diffusion boundary.")
