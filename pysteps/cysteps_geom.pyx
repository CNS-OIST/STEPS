###___license_placeholder___###

"""
This file is the user-interface file for all geom objects in STEPS.  All
objects are directly derived from the corresponding swig objects.  Geom and
Tetmesh(derived from Geom) container objects are owned by Python All other
objects are owned by c++ and container is responsible for all the cleaning-up
of these objects (see cpp/geom/geom.cpp class destructor).
"""

from steps_wm cimport *
from steps_tetmesh cimport *

# ======================================================================================================================
# Python Wrappers to namespace steps::wm
# ======================================================================================================================

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
        return _py_Comp.from_ptr(self.ptr().getComp(to_std_string(id)))

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
        return _py_Patch.from_ptr(self.ptr().getPatch(to_std_string(id)))

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

        .. note: "Inner" compartment and "outer" compartment are purely defined
        by their order to the class constructor.

        """
        self._ptr = new Patch(to_std_string(id), container.ptr(), icomp.ptr(), ocomp.ptr() if ocomp else NULL, area )

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
        return _py_Geom.from_ptr(self.ptr().getContainer())

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
        all_names = self.ptr().getSurfsys()
        return [from_std_string(s) for s in all_names]

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
        Giving a steps.model.Model, return all species in the compartment.

        Syntax::

            getAllSpecs(model)

        Arguments:
        steps.model.Model model

        Return:
        list<steps.model.Spec>

        """
        return _py_Spec.vector2list(self.ptr().getAllSpecs(model.ptr()))

    def getAllSReacs(self, _py_Model model):
        """
        Giving a steps.model.Model, return all surface reactions in the compartment.

        Syntax::

            getAllSReacs(model)

        Arguments:
        steps.model.Model model

        Return:
        list<steps.model.SReac>

        """
        return _py_SReac.vector2list(self.ptr().getAllSReacs(model.ptr()))

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
        return _py_Comp.from_ptr(self.ptr().getIComp())

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
    #cdef unique_ptr[Comp] _autodealoc
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
        self._ptr = new Comp(to_std_string(id), container.ptr(), vol)     # We create an object
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
        return _py_Geom.from_ptr(self.ptr().getContainer())

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
        all_names = self.ptr().getVolsys()
        return [from_std_string(s) for s in all_names]

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
        Giving a steps.model.Model, return all species in the compartment.

        Syntax::

            getAllSpecs(model)

        Arguments:
        steps.model.Model model

        Return:
        list<steps.model.Spec>

        """
        return _py_Spec.vector2list(self.ptr().getAllSpecs(model.ptr()))

    def getAllReacs(self, _py_Model model):
        """
        Giving a steps.model.Model, return all reactions in the compartment.

        Syntax::

            getAllReacs(model)

        Arguments:
        steps.model.Model model

        Return:
        list<steps.model.Reac>

        """
        return _py_Reac.vector2list(self.ptr().getAllReacs(model.ptr()))

    def getAllDiffs(self, _py_Model model):
        """
        Giving a steps.model.Model, return all diffusions in the compartment.

        Syntax::

            getAllDiffs(model)

        Arguments:
        steps.model.Model model

        Return:
        list<steps.model.Diff>

        """
        return _py_Diff.vector2list(self.ptr().getAllDiffs(model.ptr()))

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
        return _py_Patch.stdset2set(self.ptr().getIPatches())

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
        return _py_Patch.stdset2set(self.ptr().getOPatches())

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
    cdef readonly std.vector[uint] indices
    def __init__(self, ElementType t=ELEM_UNDEFINED, std.vector[uint] indices=[]):
        self.type = t
        self.indices = indices

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Tetmesh(_py_Geom):
    "Python wrapper class for Tetmesh"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Tetmesh] _autodealoc
    cdef Tetmesh *ptrx(self):
        return <Tetmesh*> self._ptr

    def __init__(self, std.vector[double] verts, std.vector[unsigned int] tets, std.vector[unsigned int] tris=[]):
        """
        Construction1::

        mesh = steps.geom.Tetmesh(verts, tets, tris)

        Construct a Tetmesh container by the "first" method: Supply a list of all
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
        list<float> verts
        list<int> tets
        list<int> tris

        Construction2::
        mesh = steps.geom.Tetmesh(nverts, ntets, ntris)

        Construct a Tetmesh container by the "second" method: Supply only the
        number of vertices nverts, the number of tetrahedrons ntets and the number
        of triangles ntris to the initializer, then use set methods to supply the
        vertex, tetrahedron and triangle information one by one. It is up to the
        user to make sure all information is supplied and then call setup() explicitly.
        It is highly recommended to use the first constructor wherever possible due to t
        he scope for user error when using this method.

        Arguments:
        int nverts
        int ntets
        int ntris
        """
        self._ptr = new Tetmesh( verts, tets, tris )

    def getVertex(self, unsigned int vidx):
        """
        Returns the coordinates of vertex with index vidx in the container.

        Syntax::

            getVertex(vidx)

        Arguments:
        int vidx

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getVertex(vidx)

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

    def getBar(self, unsigned int bidx):
        """
        Returns the vertices of bar with index bidx in the container.

        Syntax::

            getBar(bidx)

        Arguments:
        int bidx

        Return:
        list<int, length = 2>

        """
        return self.ptrx().getBar(bidx)

    def countBars(self, ):
        """
        Returns the total nubmer of bars in the mesh.

        Syntax::

            countBars()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().countBars()

    def getTri(self, unsigned int tidx):
        """
        Returns the triangle with index tidx in the container by its three vertex indices.

        Syntax::

            getTri(tidx)

        Arguments:
        int tidx

        Return:
        list<int, length = 3>

        """
        return self.ptrx().getTri(tidx)

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

    def getTriArea(self, unsigned int tidx):
        """
        Returns the area of the triangle with index tidx.

        Syntax::

            getTriArea(tidx)

        Arguments:
        int tidx

        Return:
        float

        """
        return self.ptrx().getTriArea(tidx)

    def getTriBars(self, unsigned int tidx):
        """
        Returns the index of the bars that comprise the triangle.

        Syntax::

            getTriBars(tidx)

        Arguments:
        int tidx

        Return:
        list<int, length = 3>

        """
        return self.ptrx().getTriBars(tidx)

    def getTriBarycenter(self, unsigned int tidx):
        """
        Returns the Cartesian coordinates of the barycenter of triangle with index tidx.

        Syntax::

            getTriBarycenter(tidx)

        Arguments:
        int tidx

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getTriBarycenter(tidx)

    def getTriNorm(self, unsigned int tidx):
        """
        Returns the normal vector of the triangle with index tidx.

        Syntax::

            getTriNorm(tidx)

        Arguments:
        int tidx

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getTriNorm(tidx)

    def getTriPatch(self, unsigned int tidx):
        """
        Returns a reference to a step.geom.TmPatch object: the patch which triangle
        with index tidx belongs to. Returns None if triangle not assigned to a patch.

        Syntax::

            getTriPatch(tidx)

        Arguments:
        int tidx

        Return:
        steps.geom.TmPatch

        """
        return _py_TmPatch.from_ptr(self.ptrx().getTriPatch(tidx))

    def setTriPatch(self, unsigned int tidx, _py_TmPatch patch):
        """
        Set the patch which triangle with index tidx belongs to..

        Syntax::

            setTriPatch(tidx, patch)

        Arguments:
        int tidx
        steps.geom.TmPatch patch

        Return:
        None

        """
        self.ptrx().setTriPatch(tidx, patch.ptrx())

    def setTriDiffBoundary(self, unsigned int tidx, _py_DiffBoundary diffb):
        """
        Set the diffusion boundary which triangle with index tidx belongs to..

        Syntax::

            setTriDiffBoundary(tidx, diffb)

        Arguments:
        int tidx
        steps.geom.DiffBoundary diffb

        Return:
        None

        """
        self.ptrx().setTriDiffBoundary(tidx, diffb.ptr())

    def getTriDiffBoundary(self, unsigned int tidx):
        """
        Returns a reference to a step.geom.Diffboundary object: the diffusion boundary triangle
        with index tidx belongs to. Returns None if triangle not assigned to a diffusion boundary.

        Syntax::

            getTriDiffBoundary(tidx)

        Arguments:
        int tidx

        Return:
        steps.geom.DiffBoundary

        """
        return _py_DiffBoundary.from_ptr(self.ptrx().getTriDiffBoundary(tidx))

    def getTriTetNeighb(self, unsigned int tidx):
        """
        Returns the indices of the two neighbouring tetrahedrons of triangle with
        index tidx. An index of -1 indicates no neighbour (triangle is on the mesh border).

        Syntax::

            getTriTetNeighb(tidx)

        Arguments:
        int tidx

        Return:
        list<int, length = 2>

        """
        return self.ptrx().getTriTetNeighb(tidx)

    def getTriTriNeighb(self, unsigned int tidx, _py_TmPatch tmpatch):
        """
        Returns the indices of the neighbouring triangles (that is all triangles that
        share a 'bar') of triangle with index tidx within the same patch.

        Syntax::

            getTriTriNeighb(int)

        Arguments:
        int tidx

        Returns:
        list<int>

        """
        return self.ptrx().getTriTriNeighb(tidx, tmpatch.ptrx())

    def getTriTriNeighbs(self, unsigned int tidx):
        """
        Returns the indices of the neighbouring triangles (that is all triangles that
        share a 'bar') of triangle with index tidx.

        Syntax::

            getTriTriNeighbs(int)

        Arguments:
        int tidx

        Returns:
        list<int>

        """
        return self.ptrx().getTriTriNeighbs(tidx)

    def getTet(self, unsigned int tidx):
        """
        Returns the tetrahedron with index tidx in the container by its four vertex indices.

        Syntax::
            getTet(tidx)

        Arguments:
        int tidx

        Return:
        list<int, length = 4>

        """
        return self.ptrx().getTet(tidx)

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

    def getTetVol(self, unsigned int tidx):
        """
        Returns the volume of the tetrahedron with index tidx.

        Syntax::

            getTetVol(tidx)

        Arguments:
        int tidx

        Return:
        float

        """
        return self.ptrx().getTetVol(tidx)

    def getTetQualityRER(self, unsigned int tidx):
        """
        Returns the radius-edge-ratio (a quality measurement) of tetrahedron with index tidx.

        Syntax::

            getTetQualityRER(tidx)

        Arguments:
        int tidx

        Return:
        float

        """
        return self.ptrx().getTetQualityRER(tidx)

    def getTetBarycenter(self, unsigned int tidx):
        """
        Returns the barycenter of the tetrahedron with index tidx.

        Syntax::

            getTetBarycenter(tidx)

        Arguments:
        int tidx

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getTetBarycenter(tidx)

    def getTetComp(self, unsigned int tidx):
        """
        Returns a reference to a steps.geom.Comp object: the compartment which
        tetrahedron with index tidx belongs to. Returns None if tetrahedron not
        assigned to a compartment.

        Syntax::

            getTetComp(tidx)

        Arguments:
        int tidx

        Return:
        steps.geom.TmComp

        """
        return _py_TmComp.from_ptr(self.ptrx().getTetComp(tidx))

    def setTetComp(self, unsigned int tidx, _py_TmComp comp):
        """
        Set the compartment which tetrahedron with index tidx belongs to..

        Syntax::

            setTetComp(tidx, comp)

        Arguments:
        int tidx
        steps.geom.TmComp comp

        Return:
        None

        """
        self.ptrx().setTetComp(tidx, comp.ptrx())

    def getTetTriNeighb(self, unsigned int tidx):
        """
        Returns the indices of the four neighbouring triangles of tetrahedron with index tidx.

        Syntax::

            getTetTriNeighb(tidx)

        Arguments:
        int tidx

        Return:
        list<int, length = 4>

        """
        return self.ptrx().getTetTriNeighb(tidx)

    def getTetTetNeighb(self, unsigned int tidx):
        """
        Returns the indices of the four neighbouring tetrahedrons of tetrahedron with index tidx.
        An index of -1 indicates no neighbour (tetrahedron is on the mesh border).

        Syntax::

            getTetTetNeighb(tidx)

        Arguments:
        int tidx

        Return:
        list<int, length = 4>

        """
        return self.ptrx().getTetTetNeighb(tidx)

    def findTetByPoint(self, std.vector[double] p):
        """
        Returns the index of the tetrahedron which encompasses a given point
        p (given in Cartesian coordinates x,y,z). Returns -1 if p is a position
        outside the mesh.

        Syntax::

            findTetByPoint(p)

        Arguments:
        list<float, length = 3> p

        Return:
        int

        """
        return self.ptrx().findTetByPoint(p)

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
        list<int>

        """
        return self.ptrx().getSurfTris()

    ## Batch related
    def getBatchVertices(self, std.vector[uint] verts):
        """
        Get coordinates of a list of vertices.

        Syntax::

            getBatchVertices(verts)

        Arguments:
        list<int> verts

        Return:
        list<float, length = len(verts) * 3>

        """
        return self.ptrx().getBatchVertices(verts)

    def getBatchVerticesNP(self, uint[:] indices, double[:] coordinates):
        """
        Get coordinates of a list of vertices.

        Syntax::

            import numpy as np
            indices = np.array([0, 1, 2], dtype= np.uint32)
            coordinates = np.zeros(len(indices) * 3)
            getBatchVertices(indices, coordinates)

        Arguments:
        numpy.array<uint> indices
        numpy.array<float, length = len(indices) * 3> coordinates

        Return:
        None

        """
        self.ptrx().getBatchVerticesNP(&indices[0], indices.shape[0], &coordinates[0], coordinates.shape[0])

    def getBatchTris(self, std.vector[uint] tris):
        """
        Get vertex indices of a list of triangles.

        Syntax::

            getBatchTris(tris)

        Arguments:
        list<int> tris

        Return:
        list<int, length = len(tris) * 3>

        """
        return self.ptrx().getBatchTris(tris)

    def getBatchTrisNP(self, uint[:] t_indices, uint[:] v_indices):
        """
        Get vertex indices of a list of triangles.

        Syntax::

            getBatchTrisNP(t_indices, v_indices)

        Arguments:
        numpy.array<uint> t_indices
        numpy.array<uint, length = len(t_indices) * 3> v_indices

        Return:
        None

        """
        self.ptrx().getBatchTrisNP(&t_indices[0], t_indices.shape[0], &v_indices[0], v_indices.shape[0])

    def getBatchTets(self, std.vector[uint] tets):
        """
        Get vertex indices of a list of tetrahedrons.

        Syntax::

            getBatchTets(tets)

        Arguments:
        list<int> tets

        Return:
        list<int, length = len(tets) * 4>

        """
        return self.ptrx().getBatchTets(tets)

    def getBatchTetsNP(self, uint[:] t_indices, uint[:] v_indices):
        """
        Get vertex indices of a list of triangles.

        Syntax::

            getBatchTetsNP(t_indices, v_indices)

        Arguments:
        numpy.array<uint> t_indices
        numpy.array<uint, length = len(t_indices) * 4> v_indices

        Return:
        None

        """
        self.ptrx().getBatchTetsNP( & t_indices[0], t_indices.shape[0], & v_indices[0], v_indices.shape[0])

    def getTriVerticesSetSizeNP(self, uint[:] t_indices):
        """
        Return the size of a set with unique vertex indices of a list of triangles,
        preparation function for furture numpy data access.

        Syntax::

            getTriVerticesSetSizeNP(t_indices)

        Arguments:
        numpy.array<uint> t_indices

        Return:
        int

        """
        return self.ptrx().getTriVerticesSetSizeNP( & t_indices[0], t_indices.shape[0])

    def getTetVerticesSetSizeNP(self, uint[:] t_indices):
        """
        Return the size of a set with unique vertex indices of a list of tetrahedrons,
        preparation function for furture numpy data access.

        Syntax::

            getTetVerticesSetSizeNP(t_indices)

        Arguments:
        numpy.array<uint> t_indices

        Return:
        int

        """
        return self.ptrx().getTetVerticesSetSizeNP( & t_indices[0], t_indices.shape[0])

    def getTriVerticesMappingSetNP(self, uint[:] t_indices, uint[:] t_vertices, uint[:] v_set):
        """
        Get the vertex indices of a list of triangles.
        The vertex indices are reindexed, with their oringinal STEPS indices stored in a given array,
        whose size is provided by getTriVerticesSetSizeNP().

        Syntax::

            getTriVerticesMappingSetNP(t_indices, t_vertices, v_set)

        Arguments:
        numpy.array<uint> t_indices
        numpy.array<uint, length = length(t_indices) * 3> t_vertices
        numpy.array<uint, length = getTriVerticesSetSizeNP(t_indices)> v_set

        Return:
        None

        """
        self.ptrx().getTriVerticesMappingSetNP( & t_indices[0], t_indices.shape[0], & t_vertices[0], t_vertices.shape[0], & v_set[0], v_set.shape[0])

    def getTetVerticesMappingSetNP(self, uint[:] t_indices, uint[:] t_vertices, uint[:] v_set):
        """
        Get the vertex indices of a list of tetrahedrons.
        The vertex indices are reindexed, with their oringinal STEPS indices stored in a given array,
        whose size is provided by getTriVerticesSetSizeNP().

        Syntax::

            getTetVerticesMappingSetNP(t_indices, t_vertices, v_set)

        Arguments:
        numpy.array<uint> t_indices
        numpy.array<uint, length = length(t_indices) * 4> t_vertices
        numpy.array<uint, length = getTriVerticesSetSizeNP(t_indices)> v_set

        Return:
        None

        """
        self.ptrx().getTetVerticesMappingSetNP( & t_indices[0], t_indices.shape[0], & t_vertices[0], t_vertices.shape[0], & v_set[0], v_set.shape[0])

    def genPointsInTet(self, uint tidx, uint npnts, double[:] coords):
        """
        Generate npnts random point coordinates x,y,z within a tetraedron with index tidx, export it to NumPy array cords

        Syntax::

            genPointsInTet(t_idx, npnts, coords)

        Arguments:
        int tidx
        int npnts
        numpy.array<float, length = npnts * 3> coords

        Return:
        None

        """
        if not len(coords): return False
        return self.ptrx().genPointsInTet(tidx, npnts, &coords[0], coords.shape[0])

    def genPointsInTri(self, uint tidx, uint npnts, double[:] coords):
        """
        Generate npnts random point coordinates x,y,z within a triangle with index tidx, export it to NumPy array cords

        Syntax::

            genPointsInTri(t_idx, npnts, coords)

        Arguments:
        int tidx
        int npnts
        numpy.array<float, length = npnts * 3> coords

        Return:
        None

        """
        if not len(coords): return False
        return self.ptrx().genPointsInTri(tidx, npnts, &coords[0], coords.shape[0])

    def genTetVisualPointsNP(self, uint[:] indices, uint[:] point_counts, double[:] coords):
        """
        For each tetrahedron index in indices, randomly generate a set of point coordinates x,y,z within the tetrahedron, where n is
        stored in point_counts. The number of points required to be generated for tetrahedron indices[i] is point_counts[i].
        All generated points are stored in cords.

        Syntax::

            genTetVisualPointsNP(indices, point_counts, coords)

        Arguments:
        numpy.array<uint> indices
        numpy.array<uint, length = length(indices)> point_counts
        numpy.array<float, length = sum(point_counts) * 3> coords

        Return:
        None

        """
        if not len(coords): return False
        return self.ptrx().genTetVisualPointsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], &coords[0], coords.shape[0])

    def genTriVisualPointsNP(self, uint[:] indices, uint[:] point_counts, double[:] coords):
        """
        For each triangle index in indices, randomly generate a set of point coordinates x,y,z within the triangle, where n is
        stored in point_counts. The number of points required to be generated for triangle indices[i] is point_counts[i].
        All generated points are stored in cords.

        Syntax::

            genTriVisualPointsNP(indices, point_counts, coords)

        Arguments:
        numpy.array<uint> indices
        numpy.array<uint, length = length(indices)> point_counts
        numpy.array<float, length = sum(point_counts) * 3> coords

        Return:
        None

        """
        if not len(coords): return False
        return self.ptrx().genTriVisualPointsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], &coords[0], coords.shape[0])

    def getBatchTetVolsNP(self, uint[:] indices, double[:] volumes):
        """
        Get the volumes of a list of tetrahedrons in indices and stored in volumes.

        Syntax::

            getBatchTetVolsNP(indices, volumes)

        Arguments:
        numpy.array<uint> indices
        numpy.array<float, length = length(indices)> volumes

        Return:
        None

        """
        return self.ptrx().getBatchTetVolsNP(&indices[0], indices.shape[0], &volumes[0], volumes.shape[0])

    def getBatchTriAreasNP(self, uint[:] indices, double[:] areas):
        """
        Get the areas of a list of triangles in indices and stored in areas.

        Syntax::

            getBatchTriAreasNP(indices, areas)

        Arguments:
        numpy.array<uint> indices
        numpy.array<float, length = length(indices)> areas

        Return:
        None

        """
        return self.ptrx().getBatchTriAreasNP(&indices[0], indices.shape[0], &areas[0], areas.shape[0])

    def reduceBatchTetPointCountsNP(self, uint[:] indices, uint[:] point_counts, double max_density):
        """
        Reduce the number of random point coordinates generated for each tetrahedron in indices so that the point density of the tetrahedron is below max_density. If the density is already below max_density for that tetrahedron, the count stored in point_counts is intacted.

        Syntax::

            reduceBatchTetPointCountsNP(indices, point_counts, max_density)

        Arguments:
        numpy.array<uint> indices
        numpy.array<uint, length = length(indices)> point_counts
        float max_density

        Return:
        None

        """
        return self.ptrx().reduceBatchTetPointCountsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], max_density)

    def reduceBatchTriPointCountsNP(self, uint[:] indices, uint[:] point_counts, double max_density):
        """
        Reduce the number of random point coordinates generated for each triangle in indices so that the point density of the triangle is below max_density. If the density is already below max_density for that triangle, the count stored in point_counts is intacted.

        Syntax::

            reduceBatchTriPointCountsNP(indices, point_counts, max_density)

        Arguments:
        numpy.array<uint> indices
        numpy.array<uint, length = length(indices)> point_counts
        float max_density

        Return:
        None

        """
        return self.ptrx().reduceBatchTriPointCountsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], max_density)

    ## ROI related ##
    def addROI(self, str id, ElementType type, std.set[unsigned int] indices):
        """
        Add a Region of Interest data record with name id to the ROI dataset.
        The type of elements stored in the ROI data can be one of the follows:
        steps.geom.ELEM_VERTEX, steps.geom.ELEM_TET, steps.geom.ELEM_TRI, steps.geom.ELEM_UNDEFINED.

        Syntax::

            addROI(id, type, indices)

        Arguments:
        string id
        ElementType type
        list<int> indices

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

    def replaceROI(self, str id, ElementType type, std.set[unsigned int] indices):
        """
        Replace a Region of Interest data record with name id with new data.

        Syntax::

            replaceROI(id, type, indices)

        Arguments:
        string id
        ElementType type
        list<int> indices

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
        list<int>

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
        cdef ROISet roi = self.ptrx().getROI(to_std_string(id))
        return _py_ROISet(roi.type, roi.indices)

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
        Get barycentres of elements stored in a tetrahedral ROI.

        Syntax::

            getROITetBarycentres(ROI_id)

        Arguments:
        string ROI_id

        Return:
        list<float>

        """
        return self.ptrx().getROITetBarycentres(to_std_string(ROI_id))

    def getROITetBarycentresNP(self, str ROI_id, double[:] centres):
        """
        Get barycentres of elements stored in a tetrahedral ROI and write to a NumPy array centres.
        The size of centres should be the same as the number of elements stored in the ROI.

        Syntax::

            getROITetBarycentresNP(ROI_id, centres)

        Arguments:
        string ROI_id
        numpy.array<float> centres

        Return:
        None

        """
        return self.ptrx().getROITetBarycentresNP(to_std_string(ROI_id), &centres[0], centres.shape[0])

    def getROITriBarycentres(self, str ROI_id):
        """
        Get barycentres of elements stored in a triangular ROI.

        Syntax::

            getROITriBarycentres(ROI_id)

        Arguments:
        string ROI_id

        Return:
        list<float>

        """
        return self.ptrx().getROITriBarycentres(to_std_string(ROI_id))

    def getROITriBarycentresNP(self, str ROI_id, double[:] centres):
        """
        Get barycentres of elements stored in a triangular ROI and write to a NumPy array centres.
        The size of centres should be the same as the number of elements stored in the ROI.

        Syntax::

            getROITriBarycentresNP(ROI_id, centres)

        Arguments:
        string ROI_id
        numpy.array<float> centres

        Return:
        None

        """
        return self.ptrx().getROITriBarycentresNP(to_std_string(ROI_id), &centres[0], centres.shape[0])

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
        list<int>

        """
        return self.ptrx().getROITris(to_std_string(ROI_id))

    def getROITrisNP(self, str ROI_id, uint[:] v_indices):
        """
        Get vertices of elements stored in a triangular ROI and write to a NumPy array v_indices.
        The size of v_indices should be 3 * the number of elements stored in the ROI.

        Syntax::

            getROITrisNP(ROI_id, v_indices)

        Arguments:
        string ROI_id
        numpy.array<uint> v_indices

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
        list<int>

        """
        return self.ptrx().getROITets(to_std_string(ROI_id))

    def getROITetsNP(self, str ROI_id, uint[:] v_indices):
        """
        Get vertices of elements stored in a tetrahedral ROI and write to a NumPy array v_indices.
        The size of v_indices should be 3 * the number of elements stored in the ROI.

        Syntax::

            getROITetsNP(ROI_id, v_indices)

        Arguments:
        string ROI_id
        numpy.array<uint> v_indices

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

    def getROITriVerticesMappingSetNP(self, str ROI_id, uint[:] t_vertices, uint[:] v_set):
        """
        Add all vertex indices of a list of triangles in a ROI to a set and write it to a NumPy array v_set.
        For each of the triangle, t_vertices records the positions of its vertices in v_set.
        i.e. For the i triangle in the ROI, the STEPS indices of its vertices are
        v_set[t_vertices[3*i]], v_set[t_vertices[3*i + 1]], v_set[t_vertices[3*i + 2]]

        Syntax::

            getROITriVerticesMappingSetNP(ROI_id, t_vertices, v_set)

        Arguments:
        string ROI_id
        numpy.array<uint> v_indices
        numpy.array<uint> v_set

        Return:
        None

        """
        return self.ptrx().getROITriVerticesMappingSetNP(to_std_string(ROI_id), &t_vertices[0], t_vertices.shape[0], &v_set[0], v_set.shape[0])

    def getROITetVerticesMappingSetNP(self, str ROI_id, uint[:] t_vertices, uint[:] v_set):
        """
        Add all vertex indices of a list of tetrahedrons in a ROI to a set and write it to a NumPy array v_set.
        For each of the tetrahedron, t_vertices records the positions of its vertices in v_set.
        i.e. For the i tetrahedron in the ROI, the STEPS indices of its vertices are
        v_set[t_vertices[4*i]], v_set[t_vertices[4*i + 1]], v_set[t_vertices[4*i + 2]], v_set[t_vertices[4*i + 3]]

        Syntax::

            getROITetVerticesMappingSetNP(ROI_id, t_vertices, v_set)

        Arguments:
        string ROI_id
        numpy.array<uint> v_indices
        numpy.array<uint> v_set

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
        numpy.array<float> coords

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
        numpy.array<float> coords

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
        float max_density

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
        float max_density

        Return:
        None

        """
        return self.ptrx().reduceROITriPointCountsNP(to_std_string(ROI_id), &point_counts[0], point_counts.shape[0], max_density)

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

    def __init__(self, str id, _py_Tetmesh container, std.vector[unsigned int] tets):
        """
        Construction::

            tmcomp = steps.geom.Comp(id, container, tets)

        Construct a TmComp object with identifier string id and assign container
        as the parent Tetmesh container. Set the group of tetrahedrons that describe
        this compartment with tets.

        Arguments:
        string id
        steps.geom.Tetmesh container
        list<int> tets
        """
        self._ptr = new TmComp(to_std_string(id), container.ptrx(), tets)

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

    def isTetInside(self, std.vector[unsigned int] tets):
        """
        Returns a list of Booleans describing if tetrahedrons tets are
        assigned to the compartment.

        Syntax::

            isTetInside(tets)

        Arguments:
        list<int> tets

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

    ## properties ##
    tets = property(getAllTetIndices, doc="List of indices of tetrahedrons associated to the compartment.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TmPatch(_py_Patch):
    "Python wrapper class for TmPatch"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[TmPatch] _autodealoc
    cdef TmPatch *ptrx(self):
        return <TmPatch*> self._ptr

    def __init__(self, str id, _py_Tetmesh container, std.vector[unsigned int] tris, _py_Comp icomp, _py_Comp ocomp=None):
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
        list<int> tris
        steps.geom.TmComp icomp
        steps.geom.TmComp ocomp (default = None)
        """
        self._ptr = new TmPatch(to_std_string(id), container.ptrx(), tris, icomp.ptr(), ocomp.ptr() if ocomp else NULL)

    def isTriInside(self, std.vector[unsigned int] tris):
        """
        Returns a list of Booleans describing if triangles tris are
        assigned to the patch.

        Syntax::

            isTriInside(tris)

        Arguments:
        list<int> tris

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
        list<int>

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

    ## properties ##
    tris = property(getAllTriIndices, doc="List of indices of triangles associated to the patch.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Memb(_py__base):
    "Python wrapper class for Memb"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Memb] _autodealoc
    cdef Memb *ptr(self):
        return <Memb*> self._ptr

    def __init__(self, str id, _py_Tetmesh container, list patches, bool verify=False, unsigned int opt_method=1, double search_percent=100.0, str opt_file_name=""):
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
        2 = breadth first search (can be time-consuming to set up, but usually faster simulation.
        If 2:breadth first search is chosen then argument search_percent can specify the number of starting points to search for the 
        lowest bandwidth. 
        If a filename (with full path) is given in optional argument opt_file_name the membrane optimization will be loaded from file,
        which was saved previously for this membrane with solver method steps.solver.Tetexact.saveMembOpt()

        Arguments:
        string id
        steps.geom.Tetmesh container
        list<steps.geom.TmPatch> patches
        bool verify (default = False)
        int opt_method (default = 1)
        float search_percent (default=100)
        string opt_file_name (default = '')
        """
        cdef std.vector[TmPatch*] _patches
        for elem in patches:
            _patches.push_back( (<_py_TmPatch>elem).ptrx() )
        self._ptr = new Memb(to_std_string(id), container.ptrx(), _patches, verify, opt_method, search_percent, to_std_string(opt_file_name))

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
        return _py_Tetmesh.from_ptr(self.ptr().getContainer())

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

    def isTriInside(self, std.vector[unsigned int] tri):
        """
        Returns a list of Booleans describing if triangles tris are
        assigned to the membrane.

        Syntax::

            isTriInside(tris)

        Arguments:
        list<int> tris

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
        list<int>

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
        list<int>

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
        int


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
        list<int>


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
        int


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
        list<int>

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
        int


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

    def __init__(self, str id, _py_Tetmesh container, std.vector[unsigned int] tris):
        """
        Construction::

            diffb = steps.geom.DiffBoundary(id, container, tris)

        Construct a DiffBoundary object with identifier string id and assign container
        as the parent geometry container, described by group of triangles tris.

        Arguments:
        string id
        steps.geom.Tetmesh container
        list<int> tris
        """
        self._ptr = new DiffBoundary(to_std_string(id), container.ptrx(), tris)

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
        return _py_Tetmesh.from_ptr(self.ptr().getContainer())

    def isTriInside(self, std.vector[unsigned int] tri):
        """
        Returns a list of Booleans describing if triangles tris are
        assigned to the diffusion boundary.

        Syntax::

            isTriInside(tris)

        Arguments:
        list<int> tris

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
        list<int>

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

    def __init__(self, str id, _py_Tetmesh container, std.vector[unsigned int] bars, list patches):
        """
        Construction::

            sdiffb = steps.geom.SDiffBoundary(id, container, bars, patches)

        Construct a SDiffBoundary object with identifier string id and assign container
        as the parent geometry container, described by group of bars. Specify the patches to be
        connected by this surface diffusion boundary to avoid potential ambiguity.

        Arguments:
        string id
        steps.geom.Tetmesh container
        list<int> bars
        list<steps.geom.TmPatch> (length 2) patches
        """
        cdef std.vector[TmPatch*] _patches
        for elem in patches:
            _patches.push_back( (<_py_TmPatch>elem).ptrx() )
        self._ptr = new SDiffBoundary(to_std_string(id), container.ptrx(), bars, _patches)

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
        
        if not isinstance(id, bytes):
            id = id.encode()

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
        return _py_Tetmesh.from_ptr(self.ptr().getContainer())

    def isBarInside(self, std.vector[unsigned int] bars):
        """
        Returns a list of Booleans describing if bars are
        assigned to the surface diffusion boundary.

        Syntax::

            isBarInside(bars)

        Arguments:
        list<int> bars

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
