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
        """
        Construction::
        
            m = steps.model.Model()
        
        Create a model container object.
        
        Arguments:
        None
        """
        self._ptr = new Model()      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getSpec(self, str id):
        """
        Returns a reference to the steps.model.Spec species object with
        identifier string spec_id (if defined).
        
        Syntax::
        
            getSpec(spec_id)
        
        Arguments:
        string spec_id
        
        Return:
        steps.model.Spec
        
        """
        return _py_Spec.from_ptr(self.ptr().getSpec(to_std_string(id)))

    def delSpec(self, str id):
        """
        Remove the steps.model.Spec species object with identifier
        string spec_id (if defined) from the model.
        
        Syntax::
        
            delSpec(spec_id)
        
        Arguments:
        string spec_id
        
        Return:
        None
        
        """
        self.ptr().delSpec(to_std_string(id))

    def getAllSpecs(self, ):
        """
        Returns a list of steps.model.Spec object references of all species in the model.
        
        Syntax::
        
            getAllSpecs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getAllSpecs())

    def getChan(self, str id):
        """
        Returns a reference to the steps.model.Chan channel object with
        identifier string chan_id (if defined).
        
        Syntax::
        
            getSpec(chan_id)
        
        Arguments:
        string chan_id
        
        Return:
        steps.model.Chan
        
        """
        return _py_Chan.from_ptr(self.ptr().getChan(to_std_string(id))) 

    def getAllChans(self, ):
        """
        Returns a list of steps.model.Chan object references of all channels in the model.
        
        Syntax::
        
            getAllChans()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Chan>
        
        """
        return _py_Chan.vector2list(self.ptr().getAllChans())

    def getVolsys(self, str id):
        """
        Returns a reference to the steps.model.Volsys volume system object with
        identifier string vsys_id (if defined).
        
        Syntax::
        
            getVolsys(vsys_id)
        
        Arguments:
        string vsys_id
        
        Return:
        steps.model.Volsys
        
        """
        return _py_Volsys.from_ptr(self.ptr().getVolsys(to_std_string(id)))

    def delVolsys(self, str id):
        """
        Remove the steps.model.Volsys volume system object with identifier string
        vsys_id (if defined) from the model.
        
        Syntax::
        
            delVolsys(vsys_id)
        
        Arguments:
        string vsys_id
        
        Return:
        None
        
        """
        self.ptr().delVolsys(to_std_string(id))

    def getAllVolsyss(self, ):
        """
        Returns a list of steps.model.Volsys object references of all volume systems in the model.
        
        Syntax::
        
            getAllVolsyss()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Volsys>
        
        """
        return _py_Volsys.vector2list(self.ptr().getAllVolsyss())

    def getSurfsys(self, str id):
        """
        Returns a reference to the steps.model.Surfsys surface system object with
        identifier string ssys_id (if defined).
        
        Syntax::
        
            getSurfsys(ssys_id)
        
        Arguments:
        string ssys_id
        
        Return:
        steps.model.Surfsys
        
        """
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys(to_std_string(id)))

    def delSurfsys(self, str id):
        """
        Remove the steps.model.Surfsys surface system object with identifier string
        ssys_id (if defined) from the model.
        
        Syntax::
        
            delSurfsys(ssys_id)
        
        Arguments:
        string ssys_id
        
        Return:
        None
        
        """
        self.ptr().delSurfsys(to_std_string(id))

    def getAllSurfsyss(self, ):
        """
        Returns a list of steps.model.Surfsys object references of all surface systems in the model.
        
        Syntax::
        
            getAllSurfsyss()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Surfsys>
        
        """
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

    def __init__(self, str id, _py_Model model, int valence=0):
        """
        Construction::
        
            s = steps.model.Spec(id, mdl)
        
        Create a species object with identifier string id and assign the
        object mdl as its parent model.
        
        Arguments:
        string id
        steps.model.Model mdl
        int valence (default=0)
        """
        self._ptr = new Spec(to_std_string(id), model.ptr(), valence)      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        """
        Get the identifier string of the species.
        
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
        Set the identifier string of the species.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        self.ptr().setID(to_std_string(id))

    def getModel(self, ):
        """
        Returns a reference to the parent steps.model.Model container object.
        
        Syntax::
        
            getModel()
        
        Arguments:
        None
        
        Return:
        steps.model.Model
        
        Attribute:
        model
        
        """
        return _py_Model.from_ptr(self.ptr().getModel())

    def setValence(self, int valence):
        """
        Set the valence of the species.
        
        Syntax::
        
            setValence(valence)
        
        Arguments:
        int valence
        
        Return:
        None
        
        """
        self.ptr().setValence(valence)

    def getValence(self, ):
        """
        Returns the valence of the species.
        
        Syntax::
        
            getValence()
        
        Arguments:
        None
        
        Return:
        int
        
        """
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
    id = property(getID, setID, doc="Identifier string of the species.")
    model = property(getModel, doc="Reference to parent model.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Chan(_py__base):
    "Python wrapper class for Chan"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Chan] _autodealoc
    cdef Chan *ptr(self):
        return <Chan*> self._ptr

    def __init__(self, str id='', _py_Model model=None):
        """
        Construction:
        
        chan = steps.model.Chan(id, mdl)
        
        Create a channel object with identifier string id and assign the object
        mdl as its parent model.
        
        Arguments:
        string id
        steps.model.Model mdl
        """
        self._ptr = new Chan(to_std_string(id), model.ptr())      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        """
        Get the identifier string of the channel.
        
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
        Set the identifier string of the channel.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        self.ptr().setID(to_std_string(id))

    def getModel(self, ):
        """
        Returns a reference to the parent steps.model.Model container object.
        
        Syntax::
        
            getModel()
        
        Arguments:
        None
        
        Return:
        steps.model.Model
        
        Attribute:
        model
        
        """
        return _py_Model.from_ptr(self.ptr().getModel())

    def getChanState(self, str id):
        """
        Returns a reference to channel state of the channel with string identifier id.
        
        Syntax::
        
            getChanState(id)
        
        Arguments:
        string id
        
        Return:
        steps.model.Chanstate
        
        """
        return _py_ChanState.from_ptr(self.ptr().getChanState(to_std_string(id)))

    def getAllChanStates(self, ):
        """
        Returns a list of steps.model.Chanstate object references of all channel states in the channel.
        
        Syntax::
        
            getAllChanStates()
        
        Arguments:
        None
        
        Return:
        list<steps.model.ChanState>
        
        """
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
    id = property(getID, setID, doc=" Identifier string of the channel. ")
    model = property(getModel, doc=" Refernece to parent model. ")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_ChanState(_py_Spec):
    "Python wrapper class for ChanState"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[ChanState] _autodealoc
    cdef ChanState *ptrx(self):
        return <ChanState*> self._ptr

    def __init__(self, str id, _py_Model model, _py_Chan chan):
        """
        Construction:
        
        chanstate = steps.model.ChanState(id, mdl, chan)
        
        Create a channel state object with identifier string id, assign the object
        mdl as its parent model. This object describes one of the possible states
        of channel object chan.
        
        Arguments:
        string id
        steps.model.Model mdl
        steps.model.Chan chan
        """
        self._ptr = new ChanState(to_std_string(id), model.ptr(), chan.ptr())      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getChan(self, ):
        """
        Returns a reference to the parent steps.model.Chan container object.
        
        Syntax::
        
            getChan()
        
        Arguments:
        None
        
        Return:
        steps.model.Chan
        
        """
        return _py_Chan.from_ptr(self.ptrx().getChan())

    def setID(self, str id):
        """
        Set the identifier string of the channel state.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        self.ptrx().setID(to_std_string(id))

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
    chan = property(getChan, doc=" Reference to parent channel. ")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Surfsys(_py__base):
    "Python wrapper class for Surfsys"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Surfsys] _autodealoc
    cdef Surfsys *ptr(self):
        return <Surfsys*> self._ptr

    def __init__(self, str id, _py_Model model):
        """
        Construction::
        
            s = steps.model.Surfsys(id, mdl)
        
        Construct a surface system object with identifier string id and assign
        the object mdl as its parent model.
        
        Arguments:
        string id
        steps.model.Model mdl
        """
        self._ptr = new Surfsys(to_std_string(id), model.ptr() )

    def getID(self, ):
        """
        Get the identifier string of the surface system.
        
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
        Set the identifier string of the surface system.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        self.ptr().setID(to_std_string(id))

    def getModel(self, ):
        """
        Returns a reference to the parent steps.model.Model container object.
        
        Syntax::
        
            getModel()
        
        Arguments:
        None
        
        Return:
        steps.model.Model
        
        """
        return _py_Model.from_ptr(self.ptr().getModel())

    def getSReac(self, str id):
        """
        Returns a reference to the steps.model.SReac surface-reaction object
        with identifier sreac_id (if defined in the surface system).
        
        Syntax::
        
            getSReac(id)
        
        Arguments:
        string id
        
        Return:
        steps.model.SReac
        
        """
        return _py_SReac.from_ptr(self.ptr().getSReac(to_std_string(id)))

    def delSReac(self, str id):
        """
        Remove the steps.model.SReac surface-reaction object with identifier
        id from the surface system.
        
        Syntax::
        
            delSReac(id)
        
        Arguments:
        string id
        
        Return:
        None
        
        """
        self.ptr().delSReac(to_std_string(id))

    def getAllSReacs(self, ):
        """
        Returns a list of references to all steps.model.SReac surface-reaction
        objects defined in the surface system.
        
        Syntax::
        
            getAllSReacs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.SReac>
        
        """
        return _py_SReac.vector2list(self.ptr().getAllSReacs())

    def getDiff(self, str id):
        """
        Returns a reference to the steps.model.Diff diffusion-rule object with
        identifier diff_id (if defined in the surface system).
        
        Syntax::
        
            getDiff(diff_id)
        
        Arguments:
        string diff_id
        
        Return:
        steps.model.Diff
        
        """
        return _py_Diff.from_ptr(self.ptr().getDiff(to_std_string(id)))

    def delDiff(self, str id):
        """
        Remove the steps.model.Diff diffusion-rule object with identifier diff_id
        from the surface system.
        
        Syntax::
        
            delDiff(diff_id)
        
        Arguments:
        string diff_id
        
        Return:
        None
        
        """
        self.ptr().delDiff(to_std_string(id))

    def getAllDiffs(self, ):
        """
        Returns a list of references to all steps.model.Diff diffusion-rule objects
        defined in the surface system.
        
        Syntax::
        
            getAllDiffs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Diff>
        
        """
        return _py_Diff.vector2list(self.ptr().getAllDiffs())

    def getVDepTrans(self, str id):
        """
        Returns a reference to the steps.model.VDepTrans voltage-dependent transition object
        with identifier id (if defined in the surface system).
        
        Syntax::
        
            getVDepTrans(id)
        
        Arguments:
        string id
        
        Return:
        steps.model.VDepTrans
        
        """
        return _py_VDepTrans.from_ptr(self.ptr().getVDepTrans(to_std_string(id)))

    def delVDepTrans(self, str id):
        """
        Remove the steps.model.VDepTrans voltage-dependent transition object with identifier
        id from the surface system.
        
        Syntax::
        
            delVDepTrans(id)
        
        Arguments:
        string id
        
        Return:
        None
        
        """
        self.ptr().delVDepTrans(to_std_string(id))

    def getAllVDepTrans(self, ):
        """
        Returns a list of references to all steps.model.VDepTrans voltage-dependent transition
        objects defined in the surface system.
        
        Syntax::
        
            getAllVDepTrans()
        
        Arguments:
        None
        
        Return:
        list<steps.model.VDepTrans>
        
        """
        return _py_VDepTrans.vector2list(self.ptr().getAllVDepTrans())

    def getVDepSReac(self, str id):
        """
        Returns a reference to the steps.model.VDepSReac voltage-dependent surface reaction object
        with identifier id (if defined in the surface system).
        
        Syntax::
        
            getVDepSReac(id)
        
        Arguments:
        string id
        
        Return:
        steps.model.VDepSReac
        
        """
        return _py_VDepSReac.from_ptr(self.ptr().getVDepSReac(to_std_string(id)))

    def delVDepSReac(self, str id):
        """
        Remove the steps.model.VDepSReac voltage-dependent surface reaction object with identifier
        id from the surface system.
        
        Syntax::
        
            delVDepSReac(id)
        
        Arguments:
        string id
        
        Return:
        None
        
        """
        self.ptr().delVDepSReac(to_std_string(id))

    def getAllVDepSReacs(self, ):
        """
        Returns a list of references to all steps.model.VDepSReac voltage-dependent surface reaction 
        objects defined in the surface system.
        
        Syntax::
        
            getAllVDepSReacs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.VDepSReac>
        
        """
        return _py_VDepSReac.vector2list(self.ptr().getAllVDepSReacs())

    def getOhmicCurr(self, str id):
        """
        Returns a reference to the steps.model.OhmicCurr ohmic current object
        with identifier id (if defined in the surface system).
        
        Syntax::
        
            getOhmicCurr(id)
        
        Arguments:
        string id
        
        Return:
        steps.model.OhmicCurr
        
        """
        return _py_OhmicCurr.from_ptr(self.ptr().getOhmicCurr(to_std_string(id)))

    def delOhmicCurr(self, str id):
        """
        Remove the steps.model.OhmicCurr ohmic current object with identifier
        id from the surface system.
        
        Syntax::
        
            delOhmicCurr(id)
        
        Arguments:
        string id
        
        Return:
        None
        
        """
        self.ptr().delOhmicCurr(to_std_string(id))

    def getAllOhmicCurrs(self, ):
        """
        Returns a list of references to all steps.model.OhmicCurr ohmic current
        objects defined in the surface system.
        
        Syntax::
        
            getAllOhmicCurrs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.OhmicCurr>
        
        """
        return _py_OhmicCurr.vector2list(self.ptr().getAllOhmicCurrs())

    def getGHKcurr(self, str id):
        """
        Returns a reference to the steps.model.GHKcurr ghk current object
        with identifier id (if defined in the surface system).
        
        Syntax::
        
            getGHKcurr(id)
        
        Arguments:
        string id
        
        Return:
        steps.model.GHKcurr
        
        """
        return _py_GHKcurr.from_ptr(self.ptr().getGHKcurr(to_std_string(id)))

    def delGHKcurr(self, str id):
        """
        Remove the steps.model.GHKcurr ghk current object with identifier
        id from the surface system.
        
        Syntax::
        
            delGHKcurr(id)
        
        Arguments:
        string id
        
        Return:
        None
        
        """
        self.ptr().delGHKcurr(to_std_string(id))

    def getAllGHKcurrs(self, ):
        """
        Returns a list of references to all steps.model.GHKcurr ghk current
        objects defined in the surface system.
        
        Syntax::
        
            getAllGHKcurrs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.GHKcurrs>
        
        """
        return _py_GHKcurr.vector2list(self.ptr().getAllGHKcurrs())

    def getAllSpecs(self, ):
        """
        Returns a list of references to all steps.model.Spec species objects included
        in the surface system; that is all reactants and products in the surface
        reactions belonging to this surface system. No duplicate member is included.
        
        Syntax::
        
            getAllSpecs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
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
    id = property(getID, setID, doc="Identifier string of the surface system.")
    model = property(getModel, doc="Reference to parent model.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Volsys(_py__base):
    "Python wrapper class for Volsys"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef shared_ptr[Volsys] _autodealoc
    cdef Volsys *ptr(self):
        return <Volsys*> self._ptr

    def __init__(self, str id, _py_Model model):
        """
        Construction::
        
            v = steps.model.Volsys(id, mdl)
        
        Construct a volume system object with identifier string id and assign
        the object mdl as its parent model.
        
        Arguments:
        string id
        steps.model.Model mdl
        """
        self._ptr = new Volsys(to_std_string(id), model.ptr())
        #self._autodealoc.reset(self.ptr())

    def getID(self, ):
        """
        Get the identifier string of the volume system.
        
        Syntax::
        
            getID()
        
        Arguments:
        None
        
        Return:
        string
        
        """
        return from_std_string(deref(self.ptr()).getID())

    def setID(self, str newid):
        """
        Set the identifier string of the volume system.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        self.ptr().setID(to_std_string(newid))

    def getModel(self, ):
        """
        Returns a reference to the parent steps.model.Model container object.
        
        Syntax::
        
            getModel()
        
        Arguments:
        None
        
        Return:
        steps.model.Model
        
        """
        return _py_Model.from_ptr(deref(self.ptr()).getModel())

    def getReac(self, str id):
        """
        Returns a reference to the steps.model.Reac reaction-rule object with
        identifier string reac_id (if defined in the volume system).
        
        Syntax::
        
            getReac(reac_id)
        
        Arguments:
        string reac_id
        
        Return:
        steps.model.Reac
        
        """
        return _py_Reac.from_ptr(self.ptr().getReac(to_std_string(id)))

    def delReac(self, str id):
        """
        Remove the steps.model.Reac reaction-rule object with identifier reac_id
        (if defined) from the volume system.
        
        Syntax::
        
            delReac(reac_id)
        
        Arguments:
        string reac_id
        
        Return:
        None
        
        """
        return self.ptr().delReac(to_std_string(id))

    def getAllReacs(self, ):
        """
        Returns a list of references to all steps.model.Reac objects in this volume
        system; that is all reaction rules belonging to this volume system. No duplicate
        member is included.
        
        Syntax::
        
            getAllReacs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Reac>
        
        """
        return _py_Reac.vector2list(self.ptr().getAllReacs())

    def getDiff(self, str id):
        """
        Returns a reference to the steps.model.Diff diffusion-rule object with
        identifier diff_id (if defined in the volume system).
        
        Syntax::
        
            getDiff(diff_id)
        
        Arguments:
        string diff_id
        
        Return:
        steps.model.Diff
        
        """
        return _py_Diff.from_ptr(self.ptr().getDiff(to_std_string(id)))

    def delDiff(self, str id):
        """
        Remove the steps.model.Diff diffusion-rule object with identifier diff_id
        from the volume system.
        
        Syntax::
        
            delDiff(diff_id)
        
        Arguments:
        string diff_id
        
        Return:
        None
        
        """
        return self.ptr().delDiff(to_std_string(id))

    def getAllDiffs(self, ):
        """
        Returns a list of references to all steps.model.Diff diffusion-rule objects
        defined in the volume system.
        
        Syntax::
        
            getAllDiffs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Diff>
        
        """
        return _py_Diff.vector2list(self.ptr().getAllDiffs())

    def getAllSpecs(self, ):
        """
        Returns a list of references to all steps.model.Spec objects in this volume system;
        that is all reactants, products or diffusing species in the reaction and diffusion
        rules belonging to this volume system. No duplicate member is included.
        
        Syntax::
        
            getAllSpecs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
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
    id = property(getID, setID, doc="Identifier string of the volume system.")
    model = property(getModel, doc="Reference to parent model.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Diff(_py__base):
    "Python wrapper class for Diff"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Diff] _autodealoc
    cdef Diff *ptr(self):
        return <Diff*> self._ptr

    def __init__(self, str id, _py__base volsys_or_surfsys, _py_Spec lig, double dcst=0 ):
        """
        Construction::
        
            diff = steps.model.Diff(id, volsys, lig, dcst = 0.0)
        
        Construct a diffusion rule object with identifier string id applied to species
        lig and assign volsys as the parent volume system. Diffusion constant is
        set by dcst.
        
        Arguments:
        string id
        steps.model.Volsys volsys
        steps.model.Spec lig
        float dcst (default = 0.0)
        
        
        Construction::
        
            diff = steps.model.Diff(id, surfsys, lig, dcst = 0.0)
        
        Construct a diffusion rule object with identifier string id applied to species
        lig and assign surfsys as the parent surface system. Diffusion constant is
        set by dcst.
        
        Arguments:
        string id
        steps.model.Surfsys surfsys
        steps.model.Spec lig
        float dcst (default = 0.0)
        """

        if isinstance(volsys_or_surfsys, _py_Volsys):
            self._ptr = new Diff(to_std_string(id), <Volsys*>(volsys_or_surfsys._ptr), lig.ptr(), dcst)
        elif isinstance(volsys_or_surfsys, _py_Surfsys):
            self._ptr = new Diff(to_std_string(id), <Surfsys*>(volsys_or_surfsys._ptr), lig.ptr(), dcst)
        else:
            raise Exception("Wrong argument: model.Diff(std.string id, Volsys_or_Surfsys, Spec lig [, float dcst])")

    def getID(self, ):
        """
        Get the identifier string of the diffusion rule.
        
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
        Set the identifier string of the diffusion rule.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        self.ptr().setID(to_std_string(id))

    def getVolsys(self, ):
        """
        Returns a reference to the parent steps.model.Volsys volume system object if a volume diffusion object.
        
        Syntax::
        
            getVolsys()
        
        Arguments:
        None
        
        Return:
        steps.model.Volsys
        
        """
        return _py_Volsys.from_ptr(self.ptr().getVolsys())

    def getSurfsys(self, ):
        """
        Returns a reference to the parent steps.model.Surfsys surface system object if surface diffusion object.
        
        Syntax::
        
            getSurfys()
        
        Arguments:
        None
        
        Return:
        steps.model.Surfsys
        
        """
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        """
        Returns a reference to the parent steps.model.Model container object.
        
        Syntax::
        
            getModel()
        
        Arguments:
        None
        
        Return:
        steps.model.Model
        
        """
        return _py_Model.from_ptr(self.ptr().getModel())

    def getLig(self, ):
        """
        get a reference to the steps.model.Spec species object to which this
        diffusion rule is applied.
        
        Syntax::
        
            getLig()
        
        Arguments:
        None
        
        Return:
        steps.model.Spec
        
        """
        return _py_Spec.from_ptr(self.ptr().getLig())

    def setLig(self, _py_Spec lig):
        """
        Set a reference to the steps.model.Spec species object to which this
        diffusion rule is applied.
        
        Syntax::
        
            setLig(lig)
        
        Arguments:
        steps.model.Spec lig
        
        Return:
        None
        
        """
        self.ptr().setLig(lig.ptr())

    def getDcst(self, ):
        """
        Get the diffusion constant for the diffusion rule, in s.i. units.
        
        Syntax::
        
            getDcst()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().getDcst()

    def setDcst(self, double dcst):
        """
        Set the diffusion constant for the diffusion rule, in s.i. units.
        
        Syntax::
        
            setDcst(dcst)
        
        Arguments:
        float dcst
        
        Return:
        None
        
        """
        self.ptr().setDcst(dcst)

    def getAllSpecs(self, ):
        """
        Return a reference to the steps.model.Spec species object to which this
        diffusion rule is applied as a list of length 1.
        
        Syntax::
        
            getAllSpecs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec, length=1>
        
        """
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
    id      = property(getID, setID, doc="Identifier string of the diffusion.")
    model   = property(getModel, doc="Reference to parent model.")
    volsys  = property(getVolsys, doc="Reference to parent volume system.")
    surfsys = property(getSurfsys, doc="Reference to parent surface system.")
    dcst    = property(getDcst, setDcst, doc="Diffusion constant.")
    lig     = property(getLig, doc="Reference to diffusion species.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Reac(_py__base):
    "Python wrapper class for Reac"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Reac] _autodealoc
    cdef Reac *ptr(self):
        return <Reac*> self._ptr

    def __init__(self, str id, _py_Volsys volsys, list lhs=[], list rhs=[], double kcst=0):
        """
        Construction::
        
            reac = steps.model.Reac(id, volsys, lhs = [ ], rhs = [ ], kcst = 0.0)
        
        Construct a reaction rule object with identifier string id and assign
        volsys as the parent volume system. A list of left hand side reactants
        may be assigned with lhs, whilst a list of right hand side products may
        be assigned with rhs, the kinetic reaction rate constant is set by kcst.
        
        Arguments:
        string id
        steps.model.Volsys volsys
        list(steps.model.Spec) lhs (default = [ ])
        list(steps.model.Spec) rhs (default = [ ])
        float kcst (default = 0.0)
        """
        if id=="" or not volsys:
            raise Exception("React requires at least two arguments: std.string id and Volsys volsys")
        #Convert vectors
        cdef std.vector[Spec*] cpp_lhs
        cdef std.vector[Spec*] cpp_rhs
        for elem in lhs:
            cpp_lhs.push_back( (<_py_Spec>elem).ptr())
        for elem in rhs:
            cpp_rhs.push_back( (<_py_Spec>elem).ptr())

        self._ptr = new Reac(to_std_string(id), volsys.ptr(), cpp_lhs, cpp_rhs, kcst)

    def getID(self):
        """
        Get the identifier string of the reaction rule.
        
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
        Set the identifier string of the reaction rule.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        self.ptr().setID(to_std_string(id))

    def getModel(self):
        """
        Returns a reference to the parent steps.model.Model container object.
        
        Syntax::
        
            getModel()
        
        Arguments:
        None
        
        Return:
        steps.model.Model
        
        """
        return _py_Model.from_ptr(self.ptr().getModel())

    def getVolsys(self):
        """
        Returns a reference to the parent steps.model.Volsys volume system object.
        
        Syntax::
        
            getVolsys()
        
        Arguments:
        None
        
        Return:
        steps.model.Volsys
        
        """
        return _py_Volsys.from_ptr(self.ptr().getVolsys())

    def getLHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects on the
        left hand side of the reaction: the reactants.
        
        Syntax::
        
            getLHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getLHS())

    def setLHS(self, list lhs):
        """
        Set a list of references to steps.model.Spec species objects on the
        left hand side of the reaction: the reactants.
        
        Syntax::
        
            setLHS(lhs)
        
        Arguments:
        list<steps.model.Spec> lhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(lhs, &vec)
        self.ptr().setLHS(vec)

    def getRHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects on the
        right hand side of the reaction: the reactants.
        
        Syntax::
        
            getRHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getRHS())

    def setRHS(self, list rhs):
        """
        Set a list of references to steps.model.Spec species objects on the
        right hand side of the reaction: the reactants.
        
        Syntax::
        
            setRHS(rhs)
        
        Arguments:
        list<steps.model.Spec> rhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(rhs, &vec)
        self.ptr().setRHS(vec)

    def getAllSpecs(self, ):
        """
        Returns a list of references to all steps.model.Spec species objects in
        the reaction; that is all reactants and products. No duplicate member
        is included.
        
        Syntax::
        
            getAllSpecs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getAllSpecs())

    def getKcst(self):
        """
        Get the kinetic reaction rate constant, in s.i. units,
        where the actual units depend on the order of the reaction.
        
        Syntax::
        
            getKcst()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().getKcst()

    def setKcst(self, double kcst):
        """
        Set the kinetic reaction rate constant, in s.i. units,
        where the actual units depend on the order of the reaction.
        
        Syntax::
        
            setKcst(kcst)
        
        Arguments:
        float kcst
        
        Return:
        None
        
        """
        self.ptr().setKcst(kcst)

    def getOrder(self):
        """
        Returns the order of this reaction.
        
        Syntax::
        
            getOrder()
        
        Arguments:
        None
        
        Return:
        int
        
        """
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
    id = property(getID, setID, doc="Identifier string of the reaction.")
    model = property(getModel, doc="Reference to parent model.")
    volsys = property(getVolsys, doc="Reference to parent volume system.")
    kcst = property(getKcst, setKcst, doc="Reaction constant.")
    lhs = property(getLHS, setLHS, doc="Left hand side reactants.")
    rhs = property(getRHS, setRHS, doc="Right hand side reactants.")
    order = property(getOrder, doc="Order of the reaction.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_SReac(_py__base):
    "Python wrapper class for SReac"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[SReac] _autodealoc
    cdef SReac *ptr(self):
        return <SReac*> self._ptr

    def __init__(self, str id, _py_Surfsys surfsys, list olhs=[], list ilhs=[], list slhs=[], list irhs=[], list srhs=[], list orhs=[], double kcst=0):
        """
        Construction::
        
            sreac = steps.model.SReac(id, surfsys, ilhs = [ ], olhs = [ ], slhs = [ ], irhs = [ ], orhs = [ ], srhs = [ ], kcst = 0.0)
        
        Construct a surface reaction rule object with identifier string
        id and assign surfsys as the parent surface system. A list of
        left hand reactants are assigned with ilhs, olhs and slhs
        (default for each is an empty list). A list of right hand side
        products are assigned with irhs, orhs and srhs (default for each
        is an empty list). The kinetic reaction rate constant is set with kcst.
        
        
        Arguments:
        string id
        steps.model.Surfsys surfsys
        list(steps.model.Spec) ilhs (default = [ ])
        list(steps.model.Spec) olhs (default = [ ])
        list(steps.model.Spec) slhs (default = [ ])
        list(steps.model.Spec) irhs (default = [ ])
        list(steps.model.Spec) orhs (default = [ ])
        list(steps.model.Spec) srhs (default = [ ])
        float kcst (default = 0.0)
        """
        cdef std.vector[Spec*] _olhs, _ilhs, _slhs, _irhs, _srhs, _orhs
        _py_Spec.list2vector(olhs, &_olhs)
        _py_Spec.list2vector(ilhs, &_ilhs)
        _py_Spec.list2vector(slhs, &_slhs)
        _py_Spec.list2vector(irhs, &_irhs)
        _py_Spec.list2vector(srhs, &_srhs)
        _py_Spec.list2vector(orhs, &_orhs)
        # Instantiate
        self._ptr = new SReac( to_std_string(id), surfsys.ptr(), _olhs, _ilhs, _slhs, _irhs, _srhs, _orhs, kcst )

    def getID(self, ):
        """
        Get the identifier string of the surface reaction rule.
        
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
        Set the identifier string of the surface reaction rule.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        return self.ptr().setID(to_std_string(id))

    def getSurfsys(self, ):
        """
        Returns a reference to the parent steps.model.Surfsys surface system object.
        
        Syntax::
        
            getSurfsys()
        
        Arguments:
        None
        
        Return:
        steps.model.Surfsys
        
        """
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        """
        Returns a reference to the parent steps.model.Model container object of parent surface system object.
        
        Syntax::
        
            getModel()
        
        Arguments:
        None
        
        Return:
        steps.model.Model
        
        """
        return _py_Model.from_ptr(self.ptr().getModel())

    #obsolete def getInner(self, ):
    #    return self.ptr().getInner()
        
    # obsolete def getOuter(self, ):
    #    return self.ptr().getOuter()

    def getOLHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the left hand side outer volume reactants.
        
        Syntax::
        
            getOLHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getOLHS())

    def setOLHS(self, list olhs):
        """
        Set a list of references to steps.model.Spec species objects;
        the left hand side outer volume reactants.
        
        Syntax::
        
            setOLHS(olhs)
        
        Arguments:
        list<steps.model.Spec) olhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(olhs, &vec)
        self.ptr().setOLHS(vec)

    def getILHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the left hand side inner volume reactants.
        
        Syntax::
        
            getILHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getILHS())

    def setILHS(self, list ilhs):
        """
        Set a list of references to steps.model.Spec species objects;
        the left hand side inner volume reactants.
        
        Syntax::
        
            setILHS(ilhs)
        
        Arguments:
        list<steps.model.Spec> ilhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(ilhs, &vec)
        self.ptr().setILHS(vec)

    def getSLHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the left hand side surface reactants.
        
        Syntax::
        
            getSLHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getSLHS())

    def setSLHS(self, list slhs):
        """
        Set a list of references to steps.model.Spec species objects;
        the left hand side surface reactants.
        
        Syntax::
        
            setSLHS(slhs)
        
        Arguments:
        list<steps.model.Spec> slhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(slhs, &vec)
        self.ptr().setSLHS(vec)

    def getIRHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the right hand side inner volume reactants.
        
        Syntax::
        
            getIRHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getIRHS())

    def setIRHS(self, list irhs):
        """
        Set a list of references to steps.model.Spec species objects;
        the right hand side inner volume reactants.
        
        Syntax::
        
            setIRHS(irhs)
        
        Arguments:
        list<steps.model.Spec> irhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(irhs, &vec)
        self.ptr().setIRHS(vec)

    def getSRHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the right hand side surface reactants.
        
        Syntax::
        
            getSRHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getSRHS())

    def setSRHS(self, list srhs):
        """
        Set a list of references to steps.model.Spec species objects;
        the right hand side surface reactants.
        
        Syntax::
        
            setSRHS(srhs)
        
        Arguments:
        list<steps.model.Spec> srhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(srhs, &vec)
        self.ptr().setSRHS(vec)

    def getORHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the right hand side outer volume reactants.
        
        Syntax::
        
            getORHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getORHS())

    def setORHS(self, list orhs):
        """
        Get a list of references to steps.model.Spec species objects;
        the right hand side outer volume reactants.
        
        Syntax::
        
            setORHS(orhs)
        
        Arguments:
        list<steps.model.Spec> orhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(orhs, &vec)
        self.ptr().setORHS(vec)

    def getOrder(self, ):
        """
        Returns the order of this surface reaction.
        
        Syntax::
        
            getOrder()
        
        Arguments:
        None
        
        Return:
        int
        
        """
        return self.ptr().getOrder()

    def getKcst(self, ):
        """
        Get the kinetic reaction rate constant, in s.i. units,
        where the actual units depend on the order of the surface reaction.
        
        Syntax::
        
            getKcst()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().getKcst()

    def setKcst(self, double kcst):
        """
        Set the kinetic reaction rate constant, in s.i. units,
        where the actual units depend on the order of the surface reaction.
        
        Syntax::
        
            setKcst(kcst)
        
        Arguments:
        float kcst
        
        Return:
        None
        
        """
        self.ptr().setKcst(kcst)

    def getAllSpecs(self, ):
        """
        Returns a list of references to all steps.model.Spec species objects in
        the surface reaction; that is all reactants and products. No duplicate member
        is included.
        
        Syntax::
        
            getAllSpecs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
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
    id      = property(getID, setID, doc="Identifier string of the surface reaction.")
    model   = property(getModel, doc="Reference to parent model.")
    surfsys = property(getSurfsys, doc="Reference to parent surface system.")
    kcst    = property(getKcst, setKcst, doc="Reaction constant.")
    olhs    = property(getOLHS, setOLHS, doc="Left hand side reactants in outer compartment.")
    ilhs    = property(getILHS, setILHS, doc="Left hand side reactants in inner compartment.")
    slhs    = property(getSLHS, setSLHS, doc="Left hand side reactants on surface.")
    orhs    = property(getORHS, setORHS, doc="Right hand side reactants in outer compartment.")
    irhs    = property(getIRHS, setIRHS, doc="Right hand side reactants in inner compartment.")
    order   = property(getOrder, doc="Order of the reaction.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_VDepTrans(_py__base):
    "Python wrapper class for VDepTrans"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[VDepTrans] _autodealoc
    cdef VDepTrans *ptr(self):
        return <VDepTrans*> self._ptr

    def __init__(self, str id, _py_Surfsys surfsys, _py_ChanState src, _py_ChanState dst, std.vector[double] ratetab, double vmin, double vmax, double dv, uint tablesize):
        """
        Construction::
        
            vdeptrans = steps.model.VDepTrans(id, surfsys, src, dst, rate=<function>)
        
        Construct a voltage-dependent transition object with identifier string id
        and assign surfsys as the parent surface system. The 'source'
        channel state is assigned with src and the 'destination' channel state
        is assigned with dst. A function that returns the transition rate in /s
        at any voltage (in volts) is supplied with rate.
        
        Arguments:
        string id
        steps.model.Surfsys surfsys
        steps.model.ChanState src
        steps.model.ChanState dst
        function rate
        
        """
        self._ptr = new VDepTrans(to_std_string(id), surfsys.ptr(), src.ptrx(), dst.ptrx(), ratetab, vmin, vmax, dv, tablesize)

    def getID(self, ):
        """
        Get the identifier string of the voltage-dependent transition.
        
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
        Set the identifier string of the voltage-dependent transition.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        self.ptr().setID(to_std_string(id))

    def getSurfsys(self, ):
        """
        Returns a reference to the parent steps.model.Surfsys surface system object.
        
        Syntax::
        
            getSurfsys()
        
        Arguments:
        None
        
        Return:
        steps.model.Surfsys
        
        """
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        """
        Returns a reference to the parent steps.model.Model container object of parent surface system object.
        
        Syntax::
        
            getModel()
        
        Arguments:
        None
        
        Return:
        steps.model.Model
        
        """
        return _py_Model.from_ptr(self.ptr().getModel())

    def getChan(self, ):
        """
        Returns a reference to the steps.model.Chan container object of source and destination channel states.
        
        Syntax::
        
            getChan()
        
        Arguments:
        None
        
        Return:
        steps.model.Chan
        
        """
        return _py_Chan.from_ptr(self.ptr().getChan())

    def getSrc(self, ):
        """
        Returns a reference to the 'source' (left-hand side) steps.model.ChanState object.
        
        Syntax::
        
            getSrc()
        
        Arguments:
        None
        
        Return:
        steps.model.ChanState
        
        """
        return _py_ChanState.from_ptr(self.ptr().getSrc())

    def setSrc(self, _py_ChanState src):
        """
        Set the 'source' (left-hand side) channel state.
        
        Syntax::
        
            setSrc(src)
        
        Arguments:
        steps.model.ChanState src
        
        Return:
        None
        
        """
        self.ptr().setSrc(src.ptrx())

    def getDst(self, ):
        """
        Returns a reference to the 'destination' (right-hand side) steps.model.ChanState object.
        
        Syntax::
        
            getDst()
        
        Arguments:
        None
        
        Return:
        steps.model.ChanState
        
        """
        return _py_ChanState.from_ptr(self.ptr().getDst())

    def setDst(self, _py_ChanState dst):
        """
        Set the 'destination' (right-hand side) channel state.
        
        Syntax::
        
            setDst(dst)
        
        Arguments:
        steps.model.ChanState dst
        
        Return:
        None
        
        """
        self.ptr().setDst(dst.ptrx())

    def getRate(self, ):
        """
        Return a list of transition rates in the default voltage range.
        
        Syntax::
        
            getRate()
        
        Arguments:
        None
        
        Return:
        list<float>
        
        """
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
    id      = property(getID, setID, doc="Identifier string of the voltage-dependent transition.")
    model   = property(getModel, doc="Reference to parent model.")
    surfsys = property(getSurfsys, doc="Reference to parent surface system.")
    src     = property(getSrc, setSrc, doc=" Reference to the channel state object that describes the 'source' channel. ") 
    dst     = property(getDst, setDst, doc=" Reference to the channel state object that describes the 'destination' channel. ") 


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_VDepSReac(_py__base):
    "Python wrapper class for VDepSReac"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[VDepSReac] _autodealoc
    cdef VDepSReac *ptr(self):
        return <VDepSReac*> self._ptr

    def __init__(self, str id, _py_Surfsys surfsys, list olhs=[], list ilhs=[], list slhs=[], list irhs=[], list srhs=[], list orhs=[], std.vector[double] ktab=[], double vmin=0, double vmax=0, double dv=0, uint tablesize=0):
        """
        Construction::
        
            vdepsreac = steps.model.VDepSReac(id, surfsys, ilhs = [ ], olhs = [ ], slhs = [ ], irhs = [ ], orhs = [ ], srhs = [ ], k = <function>, vrange = [-150.0e-3, 100.0e-3, 1.0e-4] )
            
        Construct a voltage-dependent reaction object with identifier string id
        and assign surfsys as the parent surface system. A list of
        left hand reactants are assigned with ilhs, olhs and slhs
        (default for each is an empty list). A list of right hand side
        products are assigned with irhs, orhs and srhs (default for each
        is an empty list).
        A function that returns the kinetic reaction 'constant' in ordinary, Molar
        units at any voltage (in volts) is supplied with argument k.
        A 'voltage range' over which to calculate the reaction rate is
        optionally provided by the argument vrange (in Volts) as a list:
        [minimum voltage, maximum voltage, voltage step]
        
        Arguments:
        string id
        steps.model.Surfsys surfsys
        list(steps.model.Spec) ilhs (default = [ ])
        list(steps.model.Spec) olhs (default = [ ])
        list(steps.model.Spec) slhs (default = [ ])
        list(steps.model.Spec) irhs (default = [ ])
        list(steps.model.Spec) orhs (default = [ ])
        list(steps.model.Spec) srhs (default = [ ])
		function k
        list vrange (default = [-150.0e-3, 100.0e-3, 1.0e-4])
        """
        cdef std.vector[Spec*] _olhs, _ilhs, _slhs, _irhs, _srhs, _orhs
        _py_Spec.list2vector(olhs, &_olhs)
        _py_Spec.list2vector(ilhs, &_ilhs)
        _py_Spec.list2vector(slhs, &_slhs)
        _py_Spec.list2vector(irhs, &_irhs)
        _py_Spec.list2vector(srhs, &_srhs)
        _py_Spec.list2vector(orhs, &_orhs)
        # Instantiate
        self._ptr = new VDepSReac(to_std_string(id), surfsys.ptr(), _olhs, _ilhs, _slhs, _irhs, _srhs, _orhs, ktab, vmin, vmax, dv, tablesize )      # We create an object
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        """
        Get the identifier string of the voltage-dependent reaction.
        
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
        Set the identifier string of the voltage-dependent reaction.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        self.ptr().setID(to_std_string(id))

    def getSurfsys(self, ):
        """
        Returns a reference to the parent steps.model.Surfsys surface system object.
        
        Syntax::
        
            getSurfsys()
        
        Arguments:
        None
        
        Return:
        steps.model.Surfsys
        
        """
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        """
        Returns a reference to the parent steps.model.Model container object of parent surface system object.
        
        Syntax::
        
            getModel()
        
        Arguments:
        None
        
        Return:
        steps.model.Model
        
        """
        return _py_Model.from_ptr(self.ptr().getModel())

    # Obsolete def getInner(self, ):
    #    return self.ptr().getInner()

    # obsolete def getOuter(self, ):
    #    return self.ptr().getOuter()

    def getOLHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the left hand side outer volume reactants.
        
        Syntax::
        
            getOLHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getOLHS())

    def setOLHS(self, list olhs):
        """
        Set a list of references to steps.model.Spec species objects;
        the left hand side outer volume reactants.
        
        Syntax::
        
            setOLHS(olhs)
        
        Arguments:
        list<steps.model.Spec> olhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(olhs, &vec)
        self.ptr().setOLHS(vec)

    def getILHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the left hand side inner volume reactants.
        
        Syntax::
        
            getILHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getILHS())

    def setILHS(self, list ilhs):
        """
        Set a list of references to steps.model.Spec species objects;
        the left hand side inner volume reactants.
        
        Syntax::
        
            setILHS(ilhs)
        
        Arguments:
        list<steps.model.Spec> ilhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(ilhs, &vec)
        self.ptr().setILHS(vec)

    def getSLHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the left hand side surface reactants.
        
        Syntax::
        
            getSLHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getSLHS())

    def setSLHS(self, list slhs):
        """
        Set a list of references to steps.model.Spec species objects;
        the left hand side surface reactants.
        
        Syntax::
        
            setSLHS(slhs)
        
        Arguments:
        list<steps.model.Spec> slhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(slhs, &vec)
        self.ptr().setSLHS(vec)

    def getIRHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the right hand side inner volume reactants.
        
        Syntax::
        
            getIRHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getIRHS())

    def setIRHS(self, list irhs):
        """
        Set a list of references to steps.model.Spec species objects;
        the right hand side inner volume reactants.
        
        Syntax::
        
            setIRHS(irhs)
        
        Arguments:
        list<steps.model.Spec> irhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(irhs, &vec)
        self.ptr().setIRHS(vec)

    def getSRHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the right hand side surface reactants.
        
        Syntax::
        
            getSRHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getSRHS())

    def setSRHS(self, list srhs):
        """
        Set a list of references to steps.model.Spec species objects;
        the right hand side surface reactants.
        
        Syntax::
        
            setSRHS(srhs)
        
        Arguments:
        list<steps.model.Spec> srhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(srhs, &vec)
        self.ptr().setSRHS(vec)

    def getORHS(self, ):
        """
        Get a list of references to steps.model.Spec species objects;
        the right hand side outer volume reactants.
        
        Syntax::
        
            getORHS()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
        return _py_Spec.vector2list(self.ptr().getORHS())

    def setORHS(self, list orhs):
        """
        Get a list of references to steps.model.Spec species objects;
        the right hand side outer volume reactants.
        
        Syntax::
        
            setORHS(orhs)
        
        Arguments:
        list<steps.model.Spec> orhs
        
        Return:
        None
        
        """
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(orhs, &vec)
        self.ptr().setORHS(vec)

    def getOrder(self, ):
        """
        Returns the order of this voltage-dependent reaction.
        
        Syntax::
        
            getOrder()
        
        Arguments:
        None
        
        Return:
        int
        
        """
        return self.ptr().getOrder()

    def getK(self, ):
        """
        Return a list of reaction 'constants' in the default voltage range.
        
        Syntax::
        
            getK()
        
        Arguments:
        None
        
        Return:
        list<float>
        
        """
        return self.ptr().getK()

    def getAllSpecs(self, ):
        """
        Returns a list of references to all steps.model.Spec species objects in
        the voltage-dependent reaction; that is all reactants and products. No duplicate member
        is included.
        
        Syntax::
        
            getAllSpecs()
        
        Arguments:
        None
        
        Return:
        list<steps.model.Spec>
        
        """
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
    id      = property(getID, setID, doc="Identifier string of the voltage-dependent reaction.")
    model   = property(getModel, doc="Reference to parent model.")
    surfsys = property(getSurfsys, doc="Reference to parent surface system.")
    olhs    = property(getOLHS, setOLHS, doc="left hand side reactants in outer compartment.")
    ilhs    = property(getILHS, setILHS, doc="Left hand side reactants in inner compartment.")
    slhs    = property(getSLHS, setSLHS, doc="Left hand side reactants on surface.")
    orhs    = property(getORHS, setORHS, doc="Right hand side reactants in outer compartment.")
    irhs    = property(getIRHS, setIRHS, doc="Right hand side reactants in inner compartment.")
    srhs    = property(getSRHS, setSRHS, doc="Right hand side reactants on surface.")
    order   = property(getOrder, doc="Order of the voltage-dependent reaction.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_OhmicCurr(_py__base):
    "Python wrapper class for OhmicCurr"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[OhmicCurr] _autodealoc
    cdef OhmicCurr *ptr(self):
        return <OhmicCurr*> self._ptr

    def __init__(self, str id, _py_Surfsys surfsys, _py_ChanState chanstate, double erev, double g):
        """
        Construction::
        
            ohmiccurr = steps.model.OhmicCurr(id, surfsys, chanstate, erev, g)
        
        Construct an ohmic curernt object with identifier string id
        and assign surfsys as the parent surface system. Assign to channel state
        chanstate, set the reversal potential to erev (in volts) and the single-channel
        conductance to g (in Siemens).
        
        Arguments:
        string id
        steps.model.Surfsys surfsys
        steps.model.ChanState chanstate
        float erev
        float g
        
        """
        self._ptr = new OhmicCurr(to_std_string(id), surfsys.ptr(), chanstate.ptrx(), erev, g)
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        """
        Get the identifier string of the ohmic current.
        
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
        Set the identifier string of the ohmic current.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        self.ptr().setID(to_std_string(id))

    def getSurfsys(self, ):
        """
        Returns a reference to the parent steps.model.Surfsys surface system object.
        
        Syntax::
        
            getSurfsys()
        
        Arguments:
        None
        
        Return:
        steps.model.Surfsys
        
        """
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        """
        Returns a reference to the parent steps.model.Model container object of parent surface system object.
        
        Syntax::
        
            getModel()
        
        Arguments:
        None
        
        Return:
        steps.model.Model
        
        """
        return _py_Model.from_ptr(self.ptr().getModel())

    def getChanState(self, ):
        """
        Returns a reference to the steps.model.ChanState channel state object.
        
        Syntax::
        
            getChanState()
        
        Arguments:
        None
        
        Return:
        steps.model.ChanState
        
        """
        return _py_ChanState.from_ptr(self.ptr().getChanState())

    def setChanState(self, _py_ChanState chanstate):
        """
        Set the channel state for this ohmic current.
        
        Syntax::
        
            setChanState(chanstate)
        
        Arguments:
        steps.model.ChanState chanstate
        
        Return:
        None
        
        """
        self.ptr().setChanState(chanstate.ptrx())

    def getERev(self, ):
        """
        Returns the reversal potential of this ohmic current in volts.
        
        Syntax::
        
            getERev()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().getERev()

    def setERev(self, double erev):
        """
        Set the reveral potential for this ohmic current in volts.
        
        Syntax::
        
            setERev(erev)
        
        Arguments:
        float erev
        
        Return:
        None
        
        """
        self.ptr().setERev(erev)

    def getG(self, ):
        """
        Returns the single-channel conductance for this ohmic current in Siemens.
        
        Syntax::
        
            getG()
        
        Arguments:
        None
        
        Return:
        float
        
        """
        return self.ptr().getG()

    def setG(self, double g):
        """
        Set the single-channel conductance for this ohmic current in Siemens.
        
        Syntax::
        
            setG(g)
        
        Arguments:
        float g
        
        Return:
        None
        
        """
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
    id      = property(getID, setID, doc="Identifier string of the ohmic current.")
    model   = property(getModel, doc="Reference to parent model.")
    surfsys = property(getSurfsys, doc="Reference to parent surface system.")
    erev    = property(getERev, setERev, doc=" The reversal potential (in volts). ")
    g       = property(getG, setG, doc=" The single-channel conductance (in Siemens). ")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_GHKcurr(_py__base):
    "Python wrapper class for GHKcurr"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[GHKcurr] _autodealoc
    cdef GHKcurr *ptr(self):
        return <GHKcurr*> self._ptr

    def __init__(self, str id, _py_Surfsys surfsys, _py_ChanState chanstate, _py_Spec ion, bool computeflux=True, double virtual_oconc=-1.0, double vshift=0.0):
        """
        Construction::
        
            ghkcurr = steps.model.GHKcurr(id, surfsys, chanstate, ion, computeflux=True)
        
        Construct a ghk current object with identifier string id
        and assign surfsys as the parent surface system. Assign to channel state
        chanstate, set the species that describes the current with ion- this
        species object must have a valence specified. If computeflux flag is
        set to True then the current will result in movement of ions between compartments,
        if False the current will be calculated but will not correspond to a real ion flux.
        A 'virtual outer concentration' can be specified so that the outer compartment does 
        not have to be explicitly simulated (if it retains default negative value then the outer
        concentration of ion must be simulated).
        
        Arguments:
        string id
        steps.model.Surfsys surfsys
        steps.model.ChanState chanstate
        steps.model.Spec ion
        bool computeflux
        float virtual_oconc (default = -1.0)
        
        NOTE: function setP or setPInfo must be called on the object before creating simulation object.
        """
        self._ptr = new GHKcurr(to_std_string(id), surfsys.ptr(), chanstate.ptrx(), ion.ptr(), computeflux, virtual_oconc, vshift)
        #self._autodealoc.reset(self.ptr()) # So we need to ensure its destroyed too

    def getID(self, ):
        """
        Get the identifier string of the ghk current.
        
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
        Set the identifier string of the ghk current.
        
        Syntax::
        
            setID(name)
        
        Arguments:
        string name
        
        Return:
        None
        
        """
        self.ptr().setID(to_std_string(id))

    def getSurfsys(self, ):
        """
        Returns a reference to the parent steps.model.Surfsys surface system object.
        
        Syntax::
        
            getSurfsys()
        
        Arguments:
        None
        
        Return:
        steps.model.Surfsys
        
        """
        return _py_Surfsys.from_ptr(self.ptr().getSurfsys())

    def getModel(self, ):
        """
        Returns a reference to the parent steps.model.Model container object of parent surface system object.
        
        Syntax::
        
            getModel()
        
        Arguments:
        None
        
        Return:
        steps.model.Model
        
        """
        return _py_Model.from_ptr(self.ptr().getModel())

    def getChanState(self, ):
        """
        Returns a reference to the steps.model.ChanState channel state object.
        
        Syntax::
        
            getChanState()
        
        Arguments:
        None
        
        Return:
        steps.model.ChanState
        
        """
        return _py_ChanState.from_ptr(self.ptr().getChanState())

    def setChanState(self, _py_ChanState chanstate):
        """
        Set the channel state for this ghk current.
        
        Syntax::
        
            setChanState(chanstate)
        
        Arguments:
        steps.model.ChanState chanstate
        
        Return:
        None
        
        """
        self.ptr().setChanState(chanstate.ptrx())

    def getIon(self, ):
        """
        Returns a reference to steps.model.Spec object- the ion of this ghk current.
        
        Syntax::
        
            getIon()
        
        Arguments:
        None
        
        Return:
        steps.model.Spec
        
        """
        return _py_Spec.from_ptr(self.ptr().getIon())

    def setIon(self, _py_Spec ion):
        """
        Set the ion for this ghk current.
        
        Syntax::
        
            setIon(ion)
        
        Arguments:
        steps.model.Spec ion
        
        Return:
        None
        
        """
        self.ptr().setIon(ion.ptr())

    def setPInfo(self, double g, double V, double T, double oconc, double iconc):
        """
        Supply information from a channel measurement in order to find the permeability.
        A measured single-channel coonductance (in Siemens) should be supplied, along with the potential (in volts),
        temperature (in Kelvins), the 'outer' concentration and 'inner' concentration of ion (in molar units).
        
        Syntax::
        
            setPInfo(g, V, T, oconc, iconc)
        
        Arguments:
        float g (the conductance)
        float V (the voltage)
        float T (the temperatur
        float oconc (the 'outer' concentration)
        float iconc (the 'inner' concentration)
        
        Return:
        None
        
        """
        self.ptr().setPInfo(g, V, T, oconc, iconc)

    def setP(self, double p):
        """
        Set the single-channel permeability directly (units: cubic meters / second).
        
        Syntax::
        
            setP(p)
        
        Arguments:
        float p (the single-channel permeability)
        
        Return:
        None
        
        """
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
    id      = property(getID, setID, doc="Identifier string of the ghk current.")
    model   = property(getModel, doc="Reference to parent model.")
    surfsys = property(getSurfsys, doc="Reference to parent surface system.")
    ion     = property(getIon, setIon, doc=" The current ion. ")
    chanState = property(getChanState, setChanState)
