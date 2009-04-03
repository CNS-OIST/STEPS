# geom.py

"""
This file is the user-interface file for all geom objects in steps.
All objects are directly derived from the corresponding swig objects.
Geom and Tetmesh(derived from Geom) container objects are owned by Python
All other objects are owned by c++ and container is responsible for 
all the cleaning-up of these objects (see cpp/steps/geom/geom.cpp class destructor).

"""


import geom_swig
import _geom_swig

class Geom(geom_swig.Geom) :
    def __init__(self, *args): 
        """__init__(self) -> Geom"""
        this = _geom_swig.new_Geom(*args)
        try: self.this.append(this)
        except: self.this = this
        # set Geom object to do all the cleaning-up
        self.thisown = True

class Patch(geom_swig.Patch) :
    def __init__(self, *args, **kwargs): 
        """
        __init__(self, string id, Geom container, Comp icomp, Comp ocomp=0, 
            Surfsys surfsys=0, double area=0.0) -> Patch
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

class Comp(geom_swig.Comp) :
    def __init__(self, *args, **kwargs): 
        """__init__(self, string id, Geom container, Volsys volsys=0, double vol=0.0) -> Comp"""
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

class Tetmesh(geom_swig.Tetmesh) :
    def __init__(self, *args): 
        """
        __init__(self, unsigned int nverts, unsigned int ntris, unsigned int ntets) -> Tetmesh
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

class Tet(geom_swig.Tet) :
    def __init__(self, *args): 
        """__init__(self, Tetmesh mesh, unsigned int tidx) -> Tet"""
        this = _geom_swig.new_Tet(*args)
        try: self.this.append(this)
        except: self.this = this
        # this object owned by python
        self.thisown = True

class Tri(geom_swig.Tri) :
    def __init__(self, *args): 
        """__init__(self, Tetmesh mesh, unsigned int tidx) -> Tri"""
        this = _geom_swig.new_Tri(*args)
        try: self.this.append(this)
        except: self.this = this
        # this object is owned by python
        self.thisown = True

class TmComp(geom_swig.TmComp) :
    def __init__(self, *args, **kwargs): 
        """__init__(self, string id, Tetmesh container, vector_uint tets, Volsys volsys=0) -> TmComp"""
        this = _geom_swig.new_TmComp(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_getmethods__["tets"] = _geom_swig.TmComp_getAllTetIndices
    tets = geom_swig._swig_property(_geom_swig.TmComp_getAllTetIndices)

class TmPatch(geom_swig.TmPatch) :
    def __init__(self, *args, **kwargs): 
        """
        __init__(self, string id, Tetmesh container, vector_uint tris, TmComp icomp, 
            TmComp ocomp=0, Surfsys surfsys=0) -> TmPatch
        """
        this = _geom_swig.new_TmPatch(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_getmethods__["tris"] = _geom_swig.TmPatch_getAllTriIndices
    tris = geom_swig._swig_property(_geom_swig.TmPatch_getAllTriIndices)


