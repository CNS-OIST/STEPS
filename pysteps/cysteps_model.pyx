###___license_placeholder___###

from steps_model cimport *

# ======================================================================================================================
# Python bindings to namespace steps::model
# ======================================================================================================================

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Model(_py__base):
    "Python wrapper class for Model"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef shared_ptr[Model] _autodealoc
    cdef Model *ptr(self):
        return <Model*> self._ptr

    def __init__(self):
        self._ptr = new Model()      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getSpec(self, std.string id):
        return _py_Spec.from_ptr(self.ptr().getSpec(id))

    def delSpec(self, std.string id):
        self.ptr().delSpec(id)

    def getAllSpecs(self, ):
        return _py_Spec.vector2list(self.ptr().getAllSpecs())

    def getChan(self, std.string id):
        return _py_Chan.from_ptr(self.ptr().getChan(id))

    def getAllChans(self, ):
        return _py_Chan.vector2list(self.ptr().getAllChans())

    def getVolsys(self, std.string id):
        return _py_Volsys.from_ptr(self.ptr().getVolsys(id))

    def delVolsys(self, std.string id):
        self.ptr().delVolsys(id)

    def getAllVolsyss(self, ):
        return _py_Volsys.vector2list(self.ptr().getAllVolsyss())

    def getSurfsys(self, std.string id):
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys(id))

    def delSurfsys(self, std.string id):
        self.ptr().delSurfsys(id)

    def getAllSurfsyss(self, ):
        return _py_Surfsys.vector2list(self.ptr().getAllSurfsyss())

    @staticmethod
    cdef _py_Model from_ptr(Model *ptr):
        cdef _py_Model obj = _py_Model.__new__(_py_Model)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Model from_ref(const Model &ref):
        return _py_Model.from_ptr(<Model*>&ref)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Spec(_py__base):
    "Python wrapper class for Spec"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Spec] _autodealoc
    cdef Spec *ptr(self):
        return <Spec*> self._ptr

    def __init__(self, std.string id, _py_Model model, int valence=0):
        self._ptr = new Spec(id, model.ptr(), valence)      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getModel(self, ):
        return _py_Model.from_ptr(self.ptr().getModel())

    def setValence(self, int valence):
        self.ptr().setValence(valence)

    def getValence(self, ):
        return self.ptr().getValence()

    @staticmethod
    cdef _py_Spec from_ptr(Spec *ptr):
        cdef _py_Spec obj = _py_Spec.__new__(_py_Spec)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Spec from_ref(const Spec &ref):
        return _py_Spec.from_ptr(<Spec*>&ref)

    ## Converting lists ##
    @staticmethod
    cdef std.vector[Spec*] *list2vector(list specList, std.vector[Spec*] *dstVec):
        for item in specList:
            assert isinstance(item, _py_Spec), "Wrong type of spec: " + str(type(item))
            dstVec.push_back( (<_py_Spec>item).ptr())
        return dstVec

    @staticmethod
    cdef list vector2list(std.vector[Spec*] specVec):
        return [ _py_Spec.from_ptr(elem) for elem in specVec ]

    ## properties ##
    id = property(getID, setID)
    model = property(getModel)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Chan(_py__base):
    "Python wrapper class for Chan"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Chan] _autodealoc
    cdef Chan *ptr(self):
        return <Chan*> self._ptr

    def __init__(self, std.string id='', _py_Model model=None):
        self._ptr = new Chan(id, model.ptr())      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getModel(self, ):
        return _py_Model.from_ptr(self.ptr().getModel())

    def getChanState(self, std.string id):
        return _py_ChanState.from_ptr(self.ptr().getChanState(id))

    def getAllChanStates(self, ):
        return _py_ChanState.vector2list2(self.ptr().getAllChanStates())

    @staticmethod
    cdef _py_Chan from_ptr(Chan *ptr):
        cdef _py_Chan obj = _py_Chan.__new__(_py_Chan)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Chan from_ref(const Chan &ref):
        return _py_Chan.from_ptr(<Chan*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[Chan*] vec):
        return [ _py_Chan.from_ptr(elem) for elem in vec ]

    ## properties ##
    id = property(getID, setID)
    model = property(getModel)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_ChanState(_py_Spec):
    "Python wrapper class for ChanState"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[ChanState] _autodealoc
    cdef ChanState *ptrx(self):
        return <ChanState*> self._ptr

    def __init__(self, std.string id, _py_Model model, _py_Chan chan):
        self._ptr = new ChanState(id, model.ptr(), chan.ptr())      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getChan(self, ):
        return _py_Chan.from_ptr(self.ptrx().getChan())

    def setID(self, std.string id):
        self.ptrx().setID(id)

    @staticmethod
    cdef _py_ChanState from_ptr(ChanState *ptr):
        cdef _py_ChanState obj = _py_ChanState.__new__(_py_ChanState)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_ChanState from_ref(const ChanState &ref):
        return _py_ChanState.from_ptr(<ChanState*>&ref)

    @staticmethod
    cdef list vector2list2(std.vector[ChanState*] vec):
        return [ _py_ChanState.from_ptr(elem) for elem in vec ]

    ## properties ##
    chan = property(getChan)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Surfsys(_py__base):
    "Python wrapper class for Surfsys"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Surfsys] _autodealoc
    cdef Surfsys *ptr(self):
        return <Surfsys*> self._ptr

    def __init__(self, std.string id, _py_Model model):
        self._ptr = new Surfsys( id, model.ptr() )

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getModel(self, ):
        return _py_Model.from_ptr(self.ptr().getModel())

    def getSReac(self, std.string id):
        return _py_SReac.from_ptr(self.ptr().getSReac(id))

    def delSReac(self, std.string id):
        self.ptr().delSReac(id)

    def getAllSReacs(self, ):
        return _py_SReac.vector2list(self.ptr().getAllSReacs())

    def getDiff(self, std.string id):
        return _py_Diff.from_ptr(self.ptr().getDiff(id))

    def delDiff(self, std.string id):
        self.ptr().delDiff(id)

    def getAllDiffs(self, ):
        return _py_Diff.vector2list(self.ptr().getAllDiffs())

    def getVDepTrans(self, std.string id):
        return _py_VDepTrans.from_ptr(self.ptr().getVDepTrans(id))

    def delVDepTrans(self, std.string id):
        self.ptr().delVDepTrans(id)

    def getAllVDepTrans(self, ):
        return _py_VDepTrans.vector2list(self.ptr().getAllVDepTrans())

    def getVDepSReac(self, std.string id):
        return _py_VDepSReac.from_ptr(self.ptr().getVDepSReac(id))

    def delVDepSReac(self, std.string id):
        self.ptr().delVDepSReac(id)

    def getAllVDepSReacs(self, ):
        return _py_VDepSReac.vector2list(self.ptr().getAllVDepSReacs())

    def getOhmicCurr(self, std.string id):
        return _py_OhmicCurr.from_ptr(self.ptr().getOhmicCurr(id))

    def delOhmicCurr(self, std.string id):
        self.ptr().delOhmicCurr(id)

    def getAllOhmicCurrs(self, ):
        return _py_OhmicCurr.vector2list(self.ptr().getAllOhmicCurrs())

    def getGHKcurr(self, std.string id):
        return _py_GHKcurr.from_ptr(self.ptr().getGHKcurr(id))

    def delGHKcurr(self, std.string id):
        self.ptr().delGHKcurr(id)

    def getAllGHKcurrs(self, ):
        return _py_GHKcurr.vector2list(self.ptr().getAllGHKcurrs())

    def getAllSpecs(self, ):
        return _py_Spec.vector2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_Surfsys from_ptr(Surfsys *ptr):
        cdef _py_Surfsys obj = _py_Surfsys.__new__(_py_Surfsys)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Surfsys from_ref(const Surfsys &ref):
        return _py_Surfsys.from_ptr(<Surfsys*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[Surfsys*] vec):
        return [ _py_Surfsys.from_ptr(elem) for elem in vec ]

    ## properties ##
    id = property(getID, setID)
    model = property(getModel)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Volsys(_py__base):
    "Python wrapper class for Volsys"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef shared_ptr[Volsys] _autodealoc
    cdef Volsys *ptr(self):
        return <Volsys*> self._ptr

    def __init__(self, std.string id, _py_Model model):
        self._ptr = new Volsys(id, model.ptr())
        #self._autodealoc.reset(self.ptr())

    def getID(self, ):
        return deref(self.ptr()).getID()

    def setID(self, str newid):
        cdef std.string id = newid
        self.ptr().setID(id)

    def getModel(self, ):
        return _py_Model.from_ptr(deref(self.ptr()).getModel())

    def getReac(self, std.string id):
        return _py_Reac.from_ptr(self.ptr().getReac(id))

    def delReac(self, std.string id):
        return self.ptr().delReac(id)

    def getAllReacs(self, ):
       return _py_Reac.vector2list(self.ptr().getAllReacs())

    def getDiff(self, std.string id):
        return _py_Diff.from_ptr(self.ptr().getDiff(id))

    def delDiff(self, std.string id):
        return self.ptr().delDiff(id)

    def getAllDiffs(self, ):
       return _py_Diff.vector2list(self.ptr().getAllDiffs())

    def getAllSpecs(self, ):
       return _py_Spec.vector2list(self.ptr().getAllSpecs())


    @staticmethod
    cdef _py_Volsys from_ptr(Volsys *ptr):
        cdef _py_Volsys obj = _py_Volsys.__new__(_py_Volsys)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Volsys from_ref(const Volsys &ref):
        return _py_Volsys.from_ptr(<Volsys*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[Volsys*] vec):
        return [ _py_Volsys.from_ptr(elem) for elem in vec ]

    ## properties ##
    id = property(getID, setID)
    model = property(getModel)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Diff(_py__base):
    "Python wrapper class for Diff"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Diff] _autodealoc
    cdef Diff *ptr(self):
        return <Diff*> self._ptr

    def __init__(self, std.string id, _py__base volsys_or_surfsys, _py_Spec lig, double dcst=0 ):
        if isinstance(volsys_or_surfsys, _py_Volsys):
            self._ptr = new Diff(id, <Volsys*>(volsys_or_surfsys._ptr), lig.ptr(), dcst)
        elif isinstance(volsys_or_surfsys, _py_Surfsys):
            self._ptr = new Diff(id, <Surfsys*>(volsys_or_surfsys._ptr), lig.ptr(), dcst)
        else:
            raise Exception("Wrong argument: model.Diff(std.string id, Volsys_or_Surfsys, Spec lig [, double dcst])")

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getVolsys(self, ):
        return _py_Volsys.from_ptr(self.ptr().getVolsys())

    def getSurfsys(self, ):
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        return _py_Model.from_ptr(self.ptr().getModel())

    def getLig(self, ):
        return _py_Spec.from_ptr(self.ptr().getLig())

    def setLig(self, _py_Spec lig):
        self.ptr().setLig(lig.ptr())

    def getDcst(self, ):
        return self.ptr().getDcst()

    def setDcst(self, double dcst):
        self.ptr().setDcst(dcst)

    def getAllSpecs(self, ):
       return _py_Spec.vector2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_Diff from_ptr(Diff *ptr):
        cdef _py_Diff obj = _py_Diff.__new__(_py_Diff)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Diff from_ref(const Diff &ref):
        return _py_Diff.from_ptr(<Diff*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[Diff*] vec):
        return [ _py_Diff.from_ptr(elem) for elem in vec ]

    ## properties ##
    id      = property(getID, setID)
    model   = property(getModel)
    volsys  = property(getVolsys)
    surfsys = property(getSurfsys)
    dcst    = property(getDcst, setDcst)
    lig     = property(getLig)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Reac(_py__base):
    "Python wrapper class for Reac"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Reac] _autodealoc
    cdef Reac *ptr(self):
        return <Reac*> self._ptr

    def __init__(self, std.string id, _py_Volsys volsys, list lhs=[], list rhs=[], double kcst=0):
        if id=="" or not volsys:
            raise Exception("React requires at least two arguments: std.string id and Volsys volsys")
        #Convert vectors
        cdef std.vector[Spec*] cpp_lhs
        cdef std.vector[Spec*] cpp_rhs
        for elem in lhs:
            cpp_lhs.push_back( (<_py_Spec>elem).ptr())
        for elem in rhs:
            cpp_rhs.push_back( (<_py_Spec>elem).ptr())

        self._ptr = new Reac(id, volsys.ptr(), cpp_lhs, cpp_rhs, kcst)

    def getID(self):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getModel(self):
        return _py_Model.from_ptr(self.ptr().getModel())

    def getVolsys(self):
        return _py_Volsys.from_ptr(self.ptr().getVolsys())

    def getLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getLHS())

    def setLHS(self, list lhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(lhs, &vec)
        self.ptr().setLHS(vec)

    def getRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getRHS())

    def setRHS(self, list rhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(rhs, &vec)
        self.ptr().setRHS(vec)

    def getAllSpecs(self, ):
        return _py_Spec.vector2list(self.ptr().getAllSpecs())

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, double kcst):
        self.ptr().setKcst(kcst)

    def getOrder(self):
        return self.ptr().getOrder()

    @staticmethod
    cdef _py_Reac from_ptr(Reac *ptr):
        cdef _py_Reac obj = _py_Reac.__new__(_py_Reac)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Reac from_ref(const Reac &ref):
        return _py_Reac.from_ptr(<Reac*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[Reac*] vec):
        return [ _py_Reac.from_ptr(elem) for elem in vec ]

    ## properties ##
    id = property(getID, setID)
    model = property(getModel)
    volsys = property(getVolsys)
    kcst = property(getKcst, setKcst)
    lhs = property(getLHS, setLHS)
    rhs = property(getRHS, setRHS)
    order = property(getOrder)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_SReac(_py__base):
    "Python wrapper class for SReac"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[SReac] _autodealoc
    cdef SReac *ptr(self):
        return <SReac*> self._ptr

    def __init__(self, std.string id, _py_Surfsys surfsys, list olhs=[], list ilhs=[], list slhs=[], list irhs=[], list srhs=[], list orhs=[], std.vector[double] ktab=[], double kcst=0):
        cdef std.vector[Spec*] _olhs, _ilhs, _slhs, _irhs, _srhs, _orhs
        _py_Spec.list2vector(olhs, &_olhs)
        _py_Spec.list2vector(ilhs, &_ilhs)
        _py_Spec.list2vector(slhs, &_slhs)
        _py_Spec.list2vector(irhs, &_irhs)
        _py_Spec.list2vector(srhs, &_srhs)
        _py_Spec.list2vector(orhs, &_orhs)
        # Instantiate
        self._ptr = new SReac( id, surfsys.ptr(), _olhs, _ilhs, _slhs, _irhs, _srhs, _orhs, kcst )

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        return self.ptr().setID(id)

    def getSurfsys(self, ):
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        return _py_Model.from_ptr(self.ptr().getModel())

    def getInner(self, ):
        return self.ptr().getInner()

    def getOuter(self, ):
        return self.ptr().getOuter()

    def getOLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getOLHS())

    def setOLHS(self, list olhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(olhs, &vec)
        self.ptr().setOLHS(vec)

    def getILHS(self, ):
        return _py_Spec.vector2list(self.ptr().getILHS())

    def setILHS(self, list ilhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(ilhs, &vec)
        self.ptr().setILHS(vec)

    def getSLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getSLHS())

    def setSLHS(self, list slhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(slhs, &vec)
        self.ptr().setSLHS(vec)

    def getIRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getIRHS())

    def setIRHS(self, list irhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(irhs, &vec)
        self.ptr().setIRHS(vec)

    def getSRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getSRHS())

    def setSRHS(self, list srhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(srhs, &vec)
        self.ptr().setSRHS(vec)

    def getORHS(self, ):
        return _py_Spec.vector2list(self.ptr().getORHS())

    def setORHS(self, list orhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(orhs, &vec)
        self.ptr().setORHS(vec)

    def getOrder(self, ):
        return self.ptr().getOrder()

    def getKcst(self, ):
        return self.ptr().getKcst()

    def setKcst(self, double kcst):
        self.ptr().setKcst(kcst)

    def getAllSpecs(self, ):
        return _py_Spec.vector2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_SReac from_ptr(SReac *ptr):
        cdef _py_SReac obj = _py_SReac.__new__(_py_SReac)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_SReac from_ref(const SReac &ref):
        return _py_SReac.from_ptr(<SReac*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[SReac*] vec):
        return [ _py_SReac.from_ptr(elem) for elem in vec ]

    ## properties ##
    id      = property(getID, setID)
    model   = property(getModel)
    surfsys = property(getSurfsys)
    kcst    = property(getKcst, setKcst)
    outer   = property(getOuter)
    olhs    = property(getOLHS, setOLHS)
    ilhs    = property(getILHS, setILHS)
    slhs    = property(getSLHS, setSLHS)
    orhs    = property(getORHS, setORHS)
    irhs    = property(getIRHS, setIRHS)
    srhs    = property(getSRHS, setSRHS)
    order   = property(getOrder)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_VDepTrans(_py__base):
    "Python wrapper class for VDepTrans"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[VDepTrans] _autodealoc
    cdef VDepTrans *ptr(self):
        return <VDepTrans*> self._ptr

    def __init__(self, std.string id, _py_Surfsys surfsys, _py_ChanState src, _py_ChanState dst, std.vector[double] ktab, double vmin, double vmax, double dv, uint tablesize):
            self._ptr = new VDepTrans(id, surfsys.ptr(), src.ptrx(), dst.ptrx(), ktab, vmin, vmax, dv, tablesize)

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getSurfsys(self, ):
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        return _py_Model.from_ptr(self.ptr().getModel())

    def getChan(self, ):
        return _py_Chan.from_ptr(self.ptr().getChan())

    def getSrc(self, ):
        return _py_ChanState.from_ptr(self.ptr().getSrc())

    def setSrc(self, _py_ChanState src):
        self.ptr().setSrc(src.ptrx())

    def getDst(self, ):
        return _py_ChanState.from_ptr(self.ptr().getDst())

    def setDst(self, _py_ChanState dst):
        self.ptr().setDst(dst.ptrx())

    def getRate(self, ):
        return self.ptr().getRate()

    @staticmethod
    cdef _py_VDepTrans from_ptr(VDepTrans *ptr):
        cdef _py_VDepTrans obj = _py_VDepTrans.__new__(_py_VDepTrans)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_VDepTrans from_ref(const VDepTrans &ref):
        return _py_VDepTrans.from_ptr(<VDepTrans*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[VDepTrans*] vec):
        return [ _py_VDepTrans.from_ptr(elem) for elem in vec ]

    ## properties ##
    id      = property(getID, setID)
    model   = property(getModel)
    surfsys = property(getSurfsys)
    src     = property(getSrc, setSrc)
    dst     = property(getDst, setDst)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_VDepSReac(_py__base):
    "Python wrapper class for VDepSReac"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[VDepSReac] _autodealoc
    cdef VDepSReac *ptr(self):
        return <VDepSReac*> self._ptr

    def __init__(self, std.string id, _py_Surfsys surfsys, list olhs=[], list ilhs=[], list slhs=[], list irhs=[], list srhs=[], list orhs=[], std.vector[double] ktab=[], double vmin=0, double vmax=0, double dv=0, unsigned int tablesize=0):
        cdef std.vector[Spec*] _olhs, _ilhs, _slhs, _irhs, _srhs, _orhs
        _py_Spec.list2vector(olhs, &_olhs)
        _py_Spec.list2vector(ilhs, &_ilhs)
        _py_Spec.list2vector(slhs, &_slhs)
        _py_Spec.list2vector(irhs, &_irhs)
        _py_Spec.list2vector(srhs, &_srhs)
        _py_Spec.list2vector(orhs, &_orhs)
        # Instantiate
        self._ptr = new VDepSReac(  id, surfsys.ptr(), _olhs, _ilhs, _slhs, _irhs, _srhs, _orhs, ktab, vmin, vmax, dv, tablesize )      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getSurfsys(self, ):
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        return _py_Model.from_ptr(self.ptr().getModel())

    def getInner(self, ):
        return self.ptr().getInner()

    def getOuter(self, ):
        return self.ptr().getOuter()

    def getOLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getOLHS())

    def setOLHS(self, list olhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(olhs, &vec)
        self.ptr().setOLHS(vec)

    def getILHS(self, ):
        return _py_Spec.vector2list(self.ptr().getILHS())

    def setILHS(self, list ilhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(ilhs, &vec)
        self.ptr().setILHS(vec)

    def getSLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getSLHS())

    def setSLHS(self, list slhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(slhs, &vec)
        self.ptr().setSLHS(vec)

    def getIRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getIRHS())

    def setIRHS(self, list irhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(irhs, &vec)
        self.ptr().setIRHS(vec)

    def getSRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getSRHS())

    def setSRHS(self, list srhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(srhs, &vec)
        self.ptr().setSRHS(vec)

    def getORHS(self, ):
        return _py_Spec.vector2list(self.ptr().getORHS())

    def setORHS(self, list orhs):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(orhs, &vec)
        self.ptr().setORHS(vec)

    def getOrder(self, ):
        return self.ptr().getOrder()

    def getK(self, ):
        return self.ptr().getK()

    def getAllSpecs(self, ):
        return _py_Spec.vector2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_VDepSReac from_ptr(VDepSReac *ptr):
        cdef _py_VDepSReac obj = _py_VDepSReac.__new__(_py_VDepSReac)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_VDepSReac from_ref(const VDepSReac &ref):
        return _py_VDepSReac.from_ptr(<VDepSReac*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[VDepSReac*] vec):
        return [ _py_VDepSReac.from_ptr(elem) for elem in vec ]

    ## properties ##
    id      = property(getID, setID)
    model   = property(getModel)
    surfsys = property(getSurfsys)
    olhs    = property(getOLHS, setOLHS)
    ilhs    = property(getILHS, setILHS)
    slhs    = property(getSLHS, setSLHS)
    orhs    = property(getORHS, setORHS)
    irhs    = property(getIRHS, setIRHS)
    srhs    = property(getSRHS, setSRHS)
    order   = property(getOrder)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_OhmicCurr(_py__base):
    "Python wrapper class for OhmicCurr"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[OhmicCurr] _autodealoc
    cdef OhmicCurr *ptr(self):
        return <OhmicCurr*> self._ptr

    def __init__(self, std.string id, _py_Surfsys surfsys, _py_ChanState chanstate, double erev, double g):
        self._ptr = new OhmicCurr(id, surfsys.ptr(), chanstate.ptrx(), erev, g)
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getSurfsys(self, ):
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        return _py_Model.from_ptr(self.ptr().getModel())

    def getChanState(self, ):
        return _py_ChanState.from_ptr(self.ptr().getChanState())

    def setChanState(self, _py_ChanState chanstate):
        self.ptr().setChanState(chanstate.ptrx())

    def getERev(self, ):
        return self.ptr().getERev()

    def setERev(self, double erev):
        self.ptr().setERev(erev)

    def getG(self, ):
        return self.ptr().getG()

    def setG(self, double g):
        self.ptr().setG(g)

    @staticmethod
    cdef _py_OhmicCurr from_ptr(OhmicCurr *ptr):
        cdef _py_OhmicCurr obj = _py_OhmicCurr.__new__(_py_OhmicCurr)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_OhmicCurr from_ref(const OhmicCurr &ref):
        return _py_OhmicCurr.from_ptr(<OhmicCurr*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[OhmicCurr*] vec):
        return [ _py_OhmicCurr.from_ptr(elem) for elem in vec ]

    ## properties ##
    id      = property(getID, setID)
    model   = property(getModel)
    surfsys = property(getSurfsys)
    erev    = property(getERev, setERev)
    g       = property(getG, setG)


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_GHKcurr(_py__base):
    "Python wrapper class for GHKcurr"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[GHKcurr] _autodealoc
    cdef GHKcurr *ptr(self):
        return <GHKcurr*> self._ptr

    def __init__(self, std.string id, _py_Surfsys surfsys, _py_ChanState chanstate, _py_Spec ion, bool computeflux=True, double virtual_oconc=-1, double vshift=0.0):
        self._ptr = new GHKcurr(id, surfsys.ptr(), chanstate.ptrx(), ion.ptr(), computeflux, virtual_oconc, vshift)
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        return self.ptr().getID()

    def setID(self, std.string id):
        self.ptr().setID(id)

    def getSurfsys(self, ):
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        return _py_Model.from_ptr(self.ptr().getModel())

    def getChanState(self, ):
        return _py_ChanState.from_ptr(self.ptr().getChanState())

    def setChanState(self, _py_ChanState chanstate):
        self.ptr().setChanState(chanstate.ptrx())

    def getIon(self, ):
        return _py_Spec.from_ptr(self.ptr().getIon())

    def setIon(self, _py_Spec ion):
        self.ptr().setIon(ion.ptr())

    def setPInfo(self, double g, double V, double T, double oconc, double iconc):
        self.ptr().setPInfo(g, V, T, oconc, iconc)

    def setP(self, double p):
        self.ptr().setP(p)

    @staticmethod
    cdef _py_GHKcurr from_ptr(GHKcurr *ptr):
        cdef _py_GHKcurr obj = _py_GHKcurr.__new__(_py_GHKcurr)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_GHKcurr from_ref(const GHKcurr &ref):
        return _py_GHKcurr.from_ptr(<GHKcurr*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[GHKcurr*] vec):
        return [ _py_GHKcurr.from_ptr(elem) for elem in vec ]

    ## properties ##
    id      = property(getID, setID)  #Identifier string of the ghk current
    model   = property(getModel)      #Reference to parent model
    surfsys = property(getSurfsys)    #Reference to parent surface system
    ion     = property(getIon, setIon)#The current ion.
    chanState = property(getChanState, setChanState)
