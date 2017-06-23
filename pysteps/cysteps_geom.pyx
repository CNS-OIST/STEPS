###___license_placeholder___###

from steps_wm cimport *
from steps_tetmesh cimport *

# ======================================================================================================================
# Python Wrappers to namespace steps::wm
# ======================================================================================================================

#Functions previously defined in the .i Swig files(!!)
def castToTmComp(_py_Comp base):
    return _py_TmComp.from_ptr( <TmComp*>(base.ptr()) )

def castToTmPatch(_py_Patch base):
    return _py_TmPatch.from_ptr( <TmPatch*>(base.ptr()) )


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Geom(_py__base):
    "Python wrapper class for Geom"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Geom] _autodealoc
    cdef Geom *ptr(self):
        return <Geom*> self._ptr

    def __init__(self):
        self._ptr = new Geom()      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getComp(self, std.string id):
        return _py_Comp.from_ptr(self.ptr().getComp(id))

    def delComp(self, std.string id):
        self.ptr().delComp(id)

    def getAllComps(self, ):
       return _py_Comp.vector2list(self.ptr().getAllComps())

    def getPatch(self, std.string id):
        return _py_Patch.from_ptr(self.ptr().getPatch(id))

    def delPatch(self, std.string id):
        self.ptr().delPatch(id)

    def getAllPatches(self, ):
       return _py_Patch.vector2list(self.ptr().getAllPatches())

    @staticmethod
    cdef _py_Geom from_ptr(Geom *ptr):
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

    def __init__(self, std.string id, _py_Geom container, _py_Comp icomp, _py_Comp ocomp=None, double area=0):
        self._ptr = new Patch( id, container.ptr(), icomp.ptr(), ocomp.ptr() if ocomp else NULL, area )

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getContainer(self, ):
        return _py_Geom.from_ptr(self.ptr().getContainer())

    def getArea(self, ):
        return self.ptr().getArea()

    def setArea(self, double vol):
        self.ptr().setArea(vol)

    def addSurfsys(self, std.string id):
        self.ptr().addSurfsys(id)

    def getSurfsys(self, ):
        return self.ptr().getSurfsys()

    def delSurfsys(self, std.string id):
        self.ptr().delSurfsys(id)

    def getAllSpecs(self, _py_Model model):
        return _py_Spec.vector2list(self.ptr().getAllSpecs(model.ptr()))

    def getAllSReacs(self, _py_Model model):
        return _py_SReac.vector2list(self.ptr().getAllSReacs(model.ptr()))

    def getIComp(self, ):
        return _py_Comp.from_ptr(self.ptr().getIComp())

    def getOComp(self, ):
        return _py_Comp.from_ptr(self.ptr().getOComp())

    @staticmethod
    cdef _py_Patch from_ptr(Patch *ptr):
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
    id = property(getID, setID)
    container = property(getContainer)
    surfsys = property(getSurfsys)
    area = property(getArea, setArea)
    icomp = property(getIComp)
    ocomp = property(getOComp)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Comp(_py__base):
    "Python wrapper class for Comp"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Comp] _autodealoc
    cdef Comp *ptr(self):
        return <Comp*> self._ptr

    def __init__(self, std.string id, _py_Geom container, double vol=0):
        self._ptr = new Comp(id, container.ptr(), vol)     # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getContainer(self, ):
        return _py_Geom.from_ptr(self.ptr().getContainer())

    def getVol(self, ):
        return self.ptr().getVol()

    def setVol(self, double vol):
        self.ptr().setVol(vol)

    def addVolsys(self, std.string id):
        return self.ptr().addVolsys(id)

    def getVolsys(self, ):
        return self.ptr().getVolsys()

    def delVolsys(self, std.string id):
        self.ptr().delVolsys(id)

    def getAllSpecs(self, _py_Model model):
        return _py_Spec.vector2list(self.ptr().getAllSpecs(model.ptr()))

    def getAllReacs(self, _py_Model model):
        return _py_Reac.vector2list(self.ptr().getAllReacs(model.ptr()))

    def getAllDiffs(self, _py_Model model):
        return _py_Diff.vector2list(self.ptr().getAllDiffs(model.ptr()))

    # Sets, not vector!?
    def getIPatches(self, ):
        return _py_Patch.stdset2set(self.ptr().getIPatches())

    def getOPatches(self, ):
        return _py_Patch.stdset2set(self.ptr().getOPatches())

    @staticmethod
    cdef _py_Comp from_ptr(Comp *ptr):
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
    id = property(getID, setID)
    container = property(getContainer)
    volsys = property(getVolsys)
    vol = property(getVol)
    ipatches = property(getIPatches)
    opatches = property(getOPatches)


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
        self._ptr = new Tetmesh( verts, tets, tris )

    def getVertex(self, unsigned int vidx):
        return self.ptrx().getVertex(vidx)

    def countVertices(self, ):
        return self.ptrx().countVertices()

    def getBar(self, unsigned int bidx):
        return self.ptrx().getBar(bidx)

    def countBars(self, ):
        return self.ptrx().countBars()

    def getTri(self, unsigned int tidx):
        return self.ptrx().getTri(tidx)

    def countTris(self, ):
        return self.ptrx().countTris()

    def getTriArea(self, unsigned int tidx):
        return self.ptrx().getTriArea(tidx)

    def getTriBars(self, unsigned int tidx):
        return self.ptrx().getTriBars(tidx)

    def getTriBarycenter(self, unsigned int tidx):
        return self.ptrx().getTriBarycenter(tidx)

    def getTriNorm(self, unsigned int tidx):
        return self.ptrx().getTriNorm(tidx)

    def getTriPatch(self, unsigned int tidx):
       return _py_TmPatch.from_ptr(self.ptrx().getTriPatch(tidx))

    def setTriPatch(self, unsigned int tidx, _py_TmPatch patch):
       self.ptrx().setTriPatch(tidx, patch.ptrx())

    def setTriDiffBoundary(self, unsigned int tidx, _py_DiffBoundary diffb):
       self.ptrx().setTriDiffBoundary(tidx, diffb.ptr())

    def getTriDiffBoundary(self, unsigned int tidx):
       return _py_DiffBoundary.from_ptr(self.ptrx().getTriDiffBoundary(tidx))

    def getTriTetNeighb(self, unsigned int tidx):
        return self.ptrx().getTriTetNeighb(tidx)

    def getTriTriNeighb(self, unsigned int tidx, _py_TmPatch tmpatch):
        return self.ptrx().getTriTriNeighb(tidx, tmpatch.ptrx())

    def getTriTriNeighbs(self, unsigned int tidx):
        return self.ptrx().getTriTriNeighbs(tidx)

    def getTet(self, unsigned int tidx):
        return self.ptrx().getTet(tidx)

    def countTets(self, ):
        return self.ptrx().countTets()

    def getTetVol(self, unsigned int tidx):
        return self.ptrx().getTetVol(tidx)

    def getTetQualityRER(self, unsigned int tidx):
        return self.ptrx().getTetQualityRER(tidx)

    def getTetBarycenter(self, unsigned int tidx):
        return self.ptrx().getTetBarycenter(tidx)

    def getTetComp(self, unsigned int tidx):
        return _py_TmComp.from_ptr(self.ptrx().getTetComp(tidx))

    def setTetComp(self, unsigned int tidx, _py_TmComp comp):
        self.ptrx().setTetComp(tidx, comp.ptrx())

    def getTetTriNeighb(self, unsigned int tidx):
        return self.ptrx().getTetTriNeighb(tidx)

    def getTetTetNeighb(self, unsigned int tidx):
        return self.ptrx().getTetTetNeighb(tidx)

    def findTetByPoint(self, std.vector[double] p):
        return self.ptrx().findTetByPoint(p)

    def getBoundMin(self, ):
        return self.ptrx().getBoundMin()

    def getBoundMax(self, ):
        return self.ptrx().getBoundMax()

    def getMeshVolume(self, ):
        return self.ptrx().getMeshVolume()

    def getSurfTris(self, ):
        return self.ptrx().getSurfTris()

    ## Batch related
    def getBatchVertices(self, std.vector[uint] verts):
        return self.ptrx().getBatchVertices(verts)

    def getBatchVerticesNP(self, uint[:] indices, double[:] coordinates):
        self.ptrx().getBatchVerticesNP(&indices[0], indices.shape[0], &coordinates[0], coordinates.shape[0])

    def getBatchTris(self, std.vector[uint] tris):
        return self.ptrx().getBatchTris(tris)

    def getBatchTrisNP(self, uint[:] t_indices, uint[:] v_indices):
        self.ptrx().getBatchTrisNP(&t_indices[0], t_indices.shape[0], &v_indices[0], v_indices.shape[0])

    def getBatchTets(self, std.vector[uint] tets):
        return self.ptrx().getBatchTets(tets)

    def getBatchTetsNP(self, uint[:] t_indices, uint[:] v_indices):
        self.ptrx().getBatchTetsNP( & t_indices[0], t_indices.shape[0], & v_indices[0], v_indices.shape[0])

    def getTriVerticesSetSizeNP(self, uint[:] t_indices):
        return self.ptrx().getTriVerticesSetSizeNP( & t_indices[0], t_indices.shape[0])

    def getTetVerticesSetSizeNP(self, uint[:] t_indices):
        return self.ptrx().getTetVerticesSetSizeNP( & t_indices[0], t_indices.shape[0])

    def getTriVerticesMappingSetNP(self, uint[:] t_indices, uint[:] t_vertices, uint[:] v_set):
        self.ptrx().getTriVerticesMappingSetNP( & t_indices[0], t_indices.shape[0], & t_vertices[0], t_vertices.shape[0], & v_set[0], v_set.shape[0])

    def getTetVerticesMappingSetNP(self, uint[:] t_indices, uint[:] t_vertices, uint[:] v_set):
        self.ptrx().getTetVerticesMappingSetNP( & t_indices[0], t_indices.shape[0], & t_vertices[0], t_vertices.shape[0], & v_set[0], v_set.shape[0])

    def genPointsInTet(self, uint tidx, uint npnts, double[:] coords):
        if not len(coords): return False
        return self.ptrx().genPointsInTet(tidx, npnts, &coords[0], coords.shape[0])

    def genPointsInTri(self, uint tidx, uint npnts, double[:] coords):
        if not len(coords): return False
        return self.ptrx().genPointsInTri(tidx, npnts, &coords[0], coords.shape[0])

    def genTetVisualPointsNP(self, uint[:] indices, uint[:] point_counts, double[:] coords):
        if not len(coords): return False
        return self.ptrx().genTetVisualPointsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], &coords[0], coords.shape[0])

    def genTriVisualPointsNP(self, uint[:] indices, uint[:] point_counts, double[:] coords):
        if not len(coords): return False
        return self.ptrx().genTriVisualPointsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], &coords[0], coords.shape[0])

    def getBatchTetVolsNP(self, uint[:] indices, double[:] volumes):
        return self.ptrx().getBatchTetVolsNP(&indices[0], indices.shape[0], &volumes[0], volumes.shape[0])

    def getBatchTriAreasNP(self, uint[:] indices, double[:] areas):
        return self.ptrx().getBatchTriAreasNP(&indices[0], indices.shape[0], &areas[0], areas.shape[0])

    def reduceBatchTetPointCountsNP(self, uint[:] indices, uint[:] point_counts, double max_density):
        return self.ptrx().reduceBatchTetPointCountsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], max_density)

    def reduceBatchTriPointCountsNP(self, uint[:] indices, uint[:] point_counts, double max_density):
        return self.ptrx().reduceBatchTriPointCountsNP(&indices[0], indices.shape[0], &point_counts[0], point_counts.shape[0], max_density)

    ## ROI related ##
    def addROI(self, std.string id, ElementType type, std.set[unsigned int] indices):
        self.ptrx().addROI(id, type, indices)

    def removeROI(self, std.string id):
        self.ptrx().removeROI(id)

    def replaceROI(self, std.string id, ElementType type, std.set[unsigned int] indices):
        self.ptrx().replaceROI(id, type, indices)

    def getROIType(self, std.string id):
        return self.ptrx().getROIType(id)

    def getROIData(self, std.string id):
        return self.ptrx().getROIData(id)

    def getROIDataSize(self, std.string id):
        return self.ptrx().getROIDataSize(id)

    def getNROIs(self):
        return self.ptrx().getNROIs()

    def getROI(self, std.string id):
        cdef ROISet roi = self.ptrx().getROI(id)
        return _py_ROISet(roi.type, roi.indices)

    def getROIArea(self, std.string ROI_id):
        return self.ptrx().getROIArea(ROI_id)

    def getROIVol(self, std.string ROI_id):
        return self.ptrx().getROIVol(ROI_id)

    def getAllROINames(self):
        return self.ptrx().getAllROINames()

    def checkROI(self, std.string id, ElementType type, uint count=0, bool warning=True):
        return self.ptrx().checkROI(id, type, count, warning )

    def getROITetBarycentres(self, std.string ROI_id):
        return self.ptrx().getROITetBarycentres(ROI_id)

    def getROITetBarycentresNP(self, std.string ROI_id, double[:] centres):
        return self.ptrx().getROITetBarycentresNP(ROI_id, &centres[0], centres.shape[0])

    def getROITriBarycentres(self, std.string ROI_id):
        return self.ptrx().getROITriBarycentres(ROI_id)

    def getROITriBarycentresNP(self, std.string ROI_id, double[:] centres):
        return self.ptrx().getROITriBarycentresNP(ROI_id, &centres[0], centres.shape[0])

    def getROIVertices(self, std.string ROI_id):
        return self.ptrx().getROIVertices(ROI_id)

    def getROIVerticesNP(self, std.string ROI_id, double[:] coordinates):
        return self.ptrx().getROIVerticesNP(ROI_id, &coordinates[0], coordinates.shape[0])

    def getROITris(self, std.string ROI_id):
        return self.ptrx().getROITris(ROI_id)

    def getROITrisNP(self, std.string ROI_id, uint[:] v_indices):
        return self.ptrx().getROITrisNP(ROI_id, &v_indices[0], v_indices.shape[0])

    def getROITets(self, std.string ROI_id):
        return self.ptrx().getROITets(ROI_id)

    def getROITetsNP(self, std.string ROI_id, uint[:] v_indices):
        return self.ptrx().getROITetsNP(ROI_id, &v_indices[0], v_indices.shape[0])

    def getROITriVerticesSetSizeNP(self, std.string ROI_id):
        return self.ptrx().getROITriVerticesSetSizeNP(ROI_id)

    def getROITetVerticesSetSizeNP(self, std.string ROI_id):
        return self.ptrx().getROITetVerticesSetSizeNP(ROI_id)

    def getROITriVerticesMappingSetNP(self, std.string ROI_id, uint[:] t_vertices, uint[:] v_set):
        return self.ptrx().getROITriVerticesMappingSetNP(ROI_id, &t_vertices[0], t_vertices.shape[0], &v_set[0], v_set.shape[0])

    def getROITetVerticesMappingSetNP(self, std.string ROI_id, uint[:] t_vertices, uint[:] v_set):
        return self.ptrx().getROITetVerticesMappingSetNP(ROI_id, &t_vertices[0], t_vertices.shape[0], &v_set[0], v_set.shape[0])

    def genROITetVisualPointsNP(self, std.string ROI_id, uint[:] point_counts, double[:] coords):
        return self.ptrx().genROITetVisualPointsNP(ROI_id, &point_counts[0], point_counts.shape[0], &coords[0], coords.shape[0])

    def genROITriVisualPointsNP(self, std.string ROI_id, uint[:] point_counts, double[:] coords):
        return self.ptrx().genROITriVisualPointsNP(ROI_id, &point_counts[0], point_counts.shape[0], &coords[0], coords.shape[0])

    def getROITetVolsNP(self, std.string ROI_id, double[:] volumes):
        return self.ptrx().getROITetVolsNP(ROI_id, &volumes[0], volumes.shape[0])

    def getROITriAreasNP(self, std.string ROI_id, double[:] areas):
        return self.ptrx().getROITriAreasNP(ROI_id, &areas[0], areas.shape[0])

    def reduceROITetPointCountsNP(self, std.string ROI_id, uint[:] point_counts, double max_density):
        return self.ptrx().reduceROITetPointCountsNP(ROI_id, &point_counts[0], point_counts.shape[0], max_density)

    def reduceROITriPointCountsNP(self, std.string ROI_id, uint[:] point_counts, double max_density):
        return self.ptrx().reduceROITriPointCountsNP(ROI_id, &point_counts[0], point_counts.shape[0], max_density)

    @staticmethod
    cdef _py_Tetmesh from_ptr(Tetmesh *ptr):
        cdef _py_Tetmesh obj = _py_Tetmesh.__new__(_py_Tetmesh)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Tetmesh from_ref(const Tetmesh &ref):
        return _py_Tetmesh.from_ptr(<Tetmesh*>&ref)

    ## properties ##
    nverts = property(countVertices)
    ntris  = property(countTris)
    ntets  = property(countTets)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TmComp(_py_Comp):
    "Python wrapper class for TmComp"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[TmComp] _autodealoc
    cdef TmComp *ptrx(self):
        return <TmComp*> self._ptr

    def __init__(self, std.string id, _py_Tetmesh container, std.vector[unsigned int] tets):
        self._ptr = new TmComp(id, container.ptrx(), tets)

    def setVol(self, double vol):
        self.ptrx().setVol(vol)

    def getAllTetIndices(self, ):
        return self.ptrx().getAllTetIndices()

    def countTets(self, ):
        return self.ptrx().countTets()

    def isTetInside(self, std.vector[unsigned int] tets):
        return self.ptrx().isTetInside(tets)

    def getBoundMin(self, ):
        return self.ptrx().getBoundMin()

    def getBoundMax(self, ):
        return self.ptrx().getBoundMax()

    @staticmethod
    cdef _py_TmComp from_ptr(TmComp *ptr):
        cdef _py_TmComp obj = _py_TmComp.__new__(_py_TmComp)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_TmComp from_ref(const TmComp &ref):
        return _py_TmComp.from_ptr(<TmComp*>&ref)

    ## properties ##
    tets = property(getAllTetIndices)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_TmPatch(_py_Patch):
    "Python wrapper class for TmPatch"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[TmPatch] _autodealoc
    cdef TmPatch *ptrx(self):
        return <TmPatch*> self._ptr

    def __init__(self, std.string id, _py_Tetmesh container, std.vector[unsigned int] tris, _py_Comp icomp, _py_Comp ocomp=None):
        self._ptr = new TmPatch(id, container.ptrx(), tris, icomp.ptr(), ocomp.ptr() if ocomp else NULL)

    def isTriInside(self, std.vector[unsigned int] tris):
        return self.ptrx().isTriInside(tris)

    def getAllTriIndices(self, ):
        return self.ptrx().getAllTriIndices()

    def getBoundMin(self, ):
        return self.ptrx().getBoundMin()

    def getBoundMax(self, ):
        return self.ptrx().getBoundMax()

    @staticmethod
    cdef _py_TmPatch from_ptr(TmPatch *ptr):
        cdef _py_TmPatch obj = _py_TmPatch.__new__(_py_TmPatch)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_TmPatch from_ref(const TmPatch &ref):
        return _py_TmPatch.from_ptr(<TmPatch*>&ref)

    ## properties ##
    tris = property(getAllTriIndices)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Memb(_py__base):
    "Python wrapper class for Memb"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Memb] _autodealoc
    cdef Memb *ptr(self):
        return <Memb*> self._ptr

    def __init__(self, std.string id, _py_Tetmesh container, list patches, bool verify=False, unsigned int opt_method=1, double search_percent=100.0, std.string opt_file_name=""):
        cdef std.vector[TmPatch*] _patches
        for elem in patches:
            _patches.push_back( (<_py_TmPatch>elem).ptrx() )
        self._ptr = new Memb(id, container.ptrx(), _patches, verify, opt_method, search_percent, opt_file_name)

    def getContainer(self, ):
        return _py_Tetmesh.from_ptr(self.ptr().getContainer())

    def getID(self, ):
        return self.ptr().getID()

    def isTriInside(self, std.vector[unsigned int] tri):
        return self.ptr().isTriInside(tri)

    def getAllTriIndices(self, ):
        return self.ptr().getAllTriIndices()

    def countTris(self, ):
        return self.ptr().countTris()

    def getAllVolTetIndices(self, ):
        return self.ptr().getAllVolTetIndices()

    def countVolTets(self, ):
        return self.ptr().countVolTets()

    def getAllVirtTriIndices(self, ):
        return self.ptr().getAllVirtTriIndices()

    def countVirtTris(self, ):
        return self.ptr().countVirtTris()

    def getAllVertIndices(self, ):
        return self.ptr().getAllVertIndices()

    def countVerts(self, ):
        return self.ptr().countVerts()

    def open(self, ):
        return self.ptr().open()

    @staticmethod
    cdef _py_Memb from_ptr(Memb *ptr):
        cdef _py_Memb obj = _py_Memb.__new__(_py_Memb)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Memb from_ref(const Memb &ref):
        return _py_Memb.from_ptr(<Memb*>&ref)

    ## properties ##
    tris = property(getAllTriIndices)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_DiffBoundary(_py__base):
    "Python wrapper class for DiffBoundary"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[DiffBoundary] _autodealoc
    cdef DiffBoundary *ptr(self):
        return <DiffBoundary*> self._ptr

    def __init__(self, std.string id, _py_Tetmesh container, std.vector[unsigned int] tris):
        self._ptr = new DiffBoundary(id, container.ptrx(), tris)

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getContainer(self, ):
        return _py_Tetmesh.from_ptr(self.ptr().getContainer())

    def isTriInside(self, std.vector[unsigned int] tri):
        return (self.ptr().isTriInside(tri))

    def getAllTriIndices(self, ):
        return self.ptr().getAllTriIndices()

    def getComps(self, ):
        return _py_Comp.vector2list(self.ptr().getComps())

    @staticmethod
    cdef _py_DiffBoundary from_ptr(DiffBoundary *ptr):
        cdef _py_DiffBoundary obj = _py_DiffBoundary.__new__(_py_DiffBoundary)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_DiffBoundary from_ref(const DiffBoundary &ref):
        return _py_DiffBoundary.from_ptr(<DiffBoundary*>&ref)

    ## properties ##
    tris      = property(getAllTriIndices)
    id        = property(getID, setID)
    container = property(getContainer)

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_SDiffBoundary(_py__base):
    "Python wrapper class for SDiffBoundary"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[SDiffBoundary] _autodealoc
    cdef SDiffBoundary *ptr(self):
        return <SDiffBoundary*> self._ptr

    def __init__(self, std.string id, _py_Tetmesh container, std.vector[unsigned int] bars, list patches):
        cdef std.vector[TmPatch*] _patches
        for elem in patches:
            _patches.push_back( (<_py_TmPatch>elem).ptrx() )
        self._ptr = new SDiffBoundary(id, container.ptrx(), bars, _patches)

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getContainer(self, ):
        return _py_Tetmesh.from_ptr(self.ptr().getContainer())

    def isBarInside(self, std.vector[unsigned int] bar):
        return (self.ptr().isBarInside(bar))

    def getAllBarIndices(self, ):
        return self.ptr().getAllBarIndices()

    def getPatches(self, ):
        return _py_Patch.vector2list(self.ptr().getPatches())

    @staticmethod
    cdef _py_SDiffBoundary from_ptr(SDiffBoundary *ptr):
        cdef _py_SDiffBoundary obj = _py_SDiffBoundary.__new__(_py_SDiffBoundary)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_SDiffBoundary from_ref(const SDiffBoundary &ref):
        return _py_SDiffBoundary.from_ptr(<SDiffBoundary*>&ref)

    ## properties ##
    bars      = property(getAllBarIndices)
    id        = property(getID, setID)
    container = property(getContainer)
