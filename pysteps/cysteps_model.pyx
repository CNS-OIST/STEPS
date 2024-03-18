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

import cython

from steps_model cimport *
from steps_common cimport *

import warnings

from enum import Enum
class _py_Immobilization(Enum):
    IMMOBILIZING = steps_model.IMMOBILIZING
    MOBILIZING   = steps_model.MOBILIZING
    NO_EFFECT    = steps_model.NO_EFFECT

cdef inline list string_flat_set_to_list(const flat_set[std.string]& fs):
    """
    Return a python list of strings from a flat_set of strings.
    """
    ret = list()
    cdef flat_set[std.string].const_iterator it = fs.const_begin()
    while it != fs.const_end():
        ret.append(from_std_string(deref(it)))
        cython.operator.preincrement(it)
    return ret

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
        return _py_Spec.from_ptr(&self.ptr().getSpec(to_std_string(id)))

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
        return _py_Chan.from_ptr(&self.ptr().getChan(to_std_string(id)))

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
        return _py_Volsys.from_ptr(&self.ptr().getVolsys(to_std_string(id)))

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
        return _py_Surfsys.from_ptr(&self.ptr().getSurfsys(to_std_string(id)))

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

    def getVesicle(self, str id):
        return _py_Vesicle.from_ptr(&self.ptr().getVesicle(to_std_string(id)))

    def getAllVesicles(self):
        return _py_Vesicle.vector2list(self.ptr().getAllVesicles())

    def getRaft(self, str id):
        return _py_Raft.from_ptr(&self.ptr().getRaft(to_std_string(id)))

    def getAllRafts(self):
        return _py_Raft.vector2list(self.ptr().getAllRafts())

    def getLinkSpec(self, str id):
        return _py_LinkSpec.from_ptr(&self.ptr().getLinkSpec(to_std_string(id)))

    def getAllLinkSpecs(self):
        return _py_LinkSpec.vector2list(self.ptr().getAllLinkSpecs())

    def getVesSurfsys(self, str id):
        return _py_VesSurfsys.from_ptr(&self.ptr().getVesSurfsys(to_std_string(id)))

    def getAllVesSurfsyss(self):
        return _py_VesSurfsys.vector2list(self.ptr().getAllVesSurfsyss())

    def getRaftsys(self, str id):
        return _py_Raftsys.from_ptr(&self.ptr().getRaftsys(to_std_string(id)))

    def getAllRaftsyss(self):
        return _py_Raftsys.vector2list(self.ptr().getAllRaftsyss())

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
        self._ptr = new Spec(to_std_string(id), deref(model.ptr()), valence)      # We create an object
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
        return _py_Model.from_ptr(&self.ptr().getModel())

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

    @staticmethod
    cdef list flat_set2list(flat_set[Spec*] specVec):
        return [ _py_Spec.from_ptr(elem) for elem in specVec ]

    ## Converting dicts ##
    @staticmethod
    cdef std.map[SpecP, SpecP] *dict2map(dict specDict, std.map[SpecP, SpecP] *dstMap):
        for key in specDict:
            val = specDict[key]
            assert isinstance(key, _py_Spec), "Wrong type of spec: " + str(type(key))
            assert isinstance(val, _py_Spec), "Wrong type of spec: " + str(type(val))
            dstMap.insert(std.pair[SpecP, SpecP]((<_py_Spec>key).ptr(), (<_py_Spec>val).ptr()))
        return dstMap

    @staticmethod
    cdef dict map2dict(std.map[SpecP, SpecP] specMap):
        specDict = {}
        for elems in specMap:
            specDict[_py_Spec.from_ptr(elems.first)] = _py_Spec.from_ptr(elems.second)
        return specDict

    ## properties ##
    id = property(getID, setID, doc="Identifier string of the species.")
    model = property(getModel, doc="Reference to parent model.")


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Complex(_py__base):
    "Python wrapper class for Complex"
# ----------------------------------------------------------------------------------------------------------------------
    #cdef unique_ptr[Spec] _autodealoc
    cdef Complex *ptr(self):
        return <Complex*> self._ptr

    def __init__(self, str id, _py_Model model, uint nbSubunits, uint nbStates):
        self._ptr = new Complex(to_std_string(id), deref(model.ptr()), nbSubunits, nbStates)

    def getID(self, ):
        return from_std_string(self.ptr().getID())

    def getModel(self, ):
        return _py_Model.from_ptr(&self.ptr().getModel())

    @staticmethod
    cdef _py_Complex from_ptr(Complex *ptr):
        cdef _py_Complex obj = _py_Complex.__new__(_py_Complex)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_Complex from_ref(const Complex &ref):
        return _py_Complex.from_ptr(<Complex*>&ref)


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
        self._ptr = new Chan(to_std_string(id), deref(model.ptr()))      # We create an object
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
        return _py_Model.from_ptr(&self.ptr().getModel())

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
        return _py_ChanState.from_ptr(&self.ptr().getChanState(to_std_string(id)))

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
        self._ptr = new ChanState(to_std_string(id), deref(model.ptr()), deref(chan.ptr()))      # We create an object

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
        return _py_Chan.from_ptr(&self.ptrx().getChan())

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
        self._ptr = new Surfsys(to_std_string(id), deref(model.ptr()))

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
        return _py_Model.from_ptr(&self.ptr().getModel())

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
        return _py_SReac.from_ptr(&self.ptr().getSReac(to_std_string(id)))

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

    def getAllComplexSReacs(self, ):
        """
        Returns a list of references to all steps.model.ComplexSReac complex surface-reaction
        objects defined in the surface system.

        Syntax::

            getAllComplexSReacs()

        Arguments:
        None

        Return:
        list<steps.model.ComplexSReac>

        """
        return _py_ComplexSReac.vector2list(self.ptr().getAllComplexSReacs())

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
        return _py_Diff.from_ptr(&self.ptr().getDiff(to_std_string(id)))

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
        return _py_VDepSReac.from_ptr(&self.ptr().getVDepSReac(to_std_string(id)))

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
        return _py_OhmicCurr.from_ptr(&self.ptr().getOhmicCurr(to_std_string(id)))

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
        return _py_GHKcurr.from_ptr(&self.ptr().getGHKcurr(to_std_string(id)))

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

    def getRaftGen(self, id):
        """
        Returns a reference to the steps.model.RaftGen raft generation object
        with identifier id (if defined in the surface system).

        Syntax::

            getRaftGen(id)

        Arguments:
        string id

        Return:
        steps.model.RaftGen

        """
        return _py_RaftGen.from_ptr(&self.ptr().getRaftGen(to_std_string(id)))

    def getAllRaftGens(self):
        """
        Returns a list of references to all steps.model.RaftGen raft generation
        objects defined in the surface system.

        Syntax::

            getAllRaftGens()

        Arguments:
        None

        Return:
        list<steps.model.RaftGen>

        """
        return _py_RaftGen.vector2list(self.ptr().getAllRaftGens())

    def getEndocytosis(self, id):
        """
        Returns a reference to the steps.model.Endocytosis endocytosis object
        with identifier id (if defined in the surface system).

        Syntax::

            getEndocytosis(id)

        Arguments:
        string id

        Return:
        steps.model.Endocytosis

        """        
        return _py_Endocytosis.from_ptr(&self.ptr().getEndocytosis(to_std_string(id)))

    def getAllEndocytosis(self):
        """
        Returns a list of references to all steps.model.Endocytosis endocytosis
        objects defined in the surface system.

        Syntax::

            getAllEndocytosis()

        Arguments:
        None

        Return:
        list<steps.model.Endocytosis>

        """
        return _py_Endocytosis.vector2list(self.ptr().getAllEndocytosis())


    def getAllSpecs(self, ):
        """
        Returns a list of references to all steps.model.Spec species objects included
        in the surface system; that is all reactants, products and channel states
        that appear in the interactions belonging to this surface system, regardless
        of where they appear in the geometry (i.e. includes all species in referenced
        patches and compartments). No duplicate member is included.

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
        self._ptr = new Volsys(to_std_string(id), deref(model.ptr()))
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
        return _py_Model.from_ptr(&self.ptr().getModel())

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
        return _py_Reac.from_ptr(&self.ptr().getReac(to_std_string(id)))

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
        return _py_Reac.flat_set2list(self.ptr().getAllReacs())

    def getAllComplexReacs(self):
        """
        Returns a list of references to all steps.model.ComplexReac complex reaction objects
        defined in the volume system.

        Syntax::

            getAllComplexReacs()

        Arguments:
        None

        Return:
        list<steps.model.ComplexReac>

        """
        return _py_ComplexReac.flat_set2list(self.ptr().getAllComplexReacs())

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
        return _py_Diff.from_ptr(&self.ptr().getDiff(to_std_string(id)))

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
        return _py_Diff.flat_set2list(self.ptr().getAllDiffs())

    def getVesBind(self, id):
        """
        Returns a reference to the steps.model.VesBind vesicle binding object with
        identifier vesbind_id (if defined in the volume system).

        Syntax::

            getVesBind(vesbind_id)

        Arguments:
        string vesbind_id

        Return:
        steps.model.VesBind

        """
        return _py_VesBind.from_ptr(&self.ptr().getVesBind(to_std_string(id)))

    def getAllVesBinds(self):
        """
        Returns a list of references to all steps.model.VesBind vesicle binding objects
        defined in the volume system.

        Syntax::

            getAllVesBinds()

        Arguments:
        None

        Return:
        list<steps.model.VesBind>

        """
        return _py_VesBind.flat_set2list(self.ptr().getAllVesBinds())

    def getVesUnbind(self, id):
        """
        Returns a reference to the steps.model.VesUnbind vesicle unbinding object with
        identifier vesunbind_id (if defined in the volume system).

        Syntax::

            getVesUnbind(vesunbind_id)

        Arguments:
        string vesunbind_id

        Return:
        steps.model.VesUnbind

        """
        return _py_VesUnbind.from_ptr(&self.ptr().getVesUnbind(to_std_string(id)))

    def getAllVesUnbinds(self):
        """
        Returns a list of references to all steps.model.VesUnbind vesicle unbinding objects
        defined in the volume system.

        Syntax::

            getAllVesUnbinds()

        Arguments:
        None

        Return:
        list<steps.model.VesUnbind>

        """
        return _py_VesUnbind.flat_set2list(self.ptr().getAllVesUnbinds())


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
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

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
            self._ptr = new Diff(to_std_string(id), deref(<Volsys*>(volsys_or_surfsys._ptr)), deref(lig.ptr()), dcst)
        elif isinstance(volsys_or_surfsys, _py_Surfsys):
            self._ptr = new Diff(to_std_string(id), deref(<Surfsys*>(volsys_or_surfsys._ptr)), deref(lig.ptr()), dcst)
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
        return _py_Model.from_ptr(&self.ptr().getModel())

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
        return _py_Spec.from_ptr(&self.ptr().getLig())

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
        self.ptr().setLig(deref(lig.ptr()))

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

    @staticmethod
    cdef list flat_set2list(flat_set[Diff*] vec):
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

        self._ptr = new Reac(to_std_string(id), deref(volsys.ptr()), cpp_lhs, cpp_rhs, kcst)

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
        return _py_Model.from_ptr(&self.ptr().getModel())

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
        return _py_Volsys.from_ptr(&self.ptr().getVolsys())

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
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

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

    @staticmethod
    cdef list flat_set2list(flat_set[Reac*] vec):
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
        self._ptr = new SReac( to_std_string(id), deref(surfsys.ptr()), _olhs, _ilhs, _slhs, _irhs, _srhs, _orhs, kcst )

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
        return _py_Surfsys.from_ptr(&self.ptr().getSurfsys())

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
        return _py_Model.from_ptr(&self.ptr().getModel())

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
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

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

    @staticmethod
    cdef list flat_set2list(flat_set[SReac*] vec):
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
        self._ptr = new VDepSReac(to_std_string(id), deref(surfsys.ptr()), _olhs, _ilhs, _slhs, _irhs, _srhs, _orhs, ktab, vmin, vmax, dv, tablesize )      # We create an object
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
        return _py_Surfsys.from_ptr(&self.ptr().getSurfsys())

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
        return _py_Model.from_ptr(&self.ptr().getModel())

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
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

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


# ======================================================================================================================
# Python bindings to namespace steps::model
# ======================================================================================================================

## Enums ##
cimport steps_model

_py_COMPLEX_FILTER_MAX_VALUE = steps_model.COMPLEX_FILTER_MAX_VALUE

# ----------------------------------------------------------------------------------------------------------------------

cdef std.vector[std.vector[SubunitStateFilter]] _get_filters(filts):
    cdef std.vector[std.vector[SubunitStateFilter]] filters
    filters.resize(len(filts))
    for i, filt in enumerate(filts):
        filters[i].resize(len(filt))
        for j, (_min, _max) in enumerate(filt):
            filters[i][j].min = _min
            filters[i][j].max = _max
    return filters

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_ComplexEvent(_py__base):
    "Python wrapper class for ComplexUpdateEvent"

    cdef ComplexEvent *ptr(self):
        return <ComplexEvent*> self._ptr
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_ComplexUpdateEvent(_py_ComplexEvent):
    "Python wrapper class for ComplexUpdateEvent"
    cdef ComplexEvent *ptr(self):
        return <ComplexEvent*> self._ptr

    def __init__(self, str comp, filts, std.vector[uint] reac, upds, destLoc=None):
        cdef steps_model.ComplexLocation cpp_destloc
        if destLoc is None:
            cpp_destloc = steps_model.COMP
        else:
            cpp_destloc = destLoc.value
        cdef std.vector[ComplexUpdate] updates
        updates.reserve(len(upds))
        for req, upd in upds:
            updates.push_back(ComplexUpdate(req, upd))
        self._ptr = new ComplexUpdateEvent(to_std_string(comp), _get_filters(filts), reac, updates, cpp_destloc)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_ComplexCreateEvent(_py_ComplexEvent):
    "Python wrapper class for ComplexCreateEvent"
    cdef ComplexEvent *ptr(self):
        return <ComplexEvent*> self._ptr

    def __init__(self, comp, std.vector[uint] init):
        self._ptr = new ComplexCreateEvent(to_std_string(comp), init)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_ComplexDeleteEvent(_py_ComplexEvent):
    "Python wrapper class for ComplexDeleteEvent"
    cdef ComplexEvent *ptr(self):
        return <ComplexEvent*> self._ptr

    def __init__(self, comp, filts):
        self._ptr = new ComplexDeleteEvent(to_std_string(comp), _get_filters(filts))
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_ComplexReac(_py__base):
    "Python wrapper class for ComplexReac"
# ----------------------------------------------------------------------------------------------------------------------
    cdef ComplexReac *ptr(self):
        return <ComplexReac*> self._ptr

    def __init__(self, str id, _py_Volsys volsys, list lhs=[], list rhs=[], list compEvs=[], double kcst=0):
        #Convert vectors
        cdef std.vector[Spec*] cpp_lhs
        cdef std.vector[Spec*] cpp_rhs
        for elem in lhs:
            cpp_lhs.push_back( (<_py_Spec>elem).ptr())
        for elem in rhs:
            cpp_rhs.push_back( (<_py_Spec>elem).ptr())

        cdef std.vector[ComplexEvent*] cpp_compEvs
        for ev in compEvs:
            cpp_compEvs.push_back((<_py_ComplexEvent>ev).ptr())

        self._ptr = new ComplexReac(to_std_string(id), deref(volsys.ptr()), cpp_lhs, cpp_rhs, cpp_compEvs, kcst)

    def getID(self):
        return from_std_string(self.ptr().getID())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getVolsys(self):
        return _py_Volsys.from_ptr(&self.ptr().getVolsys())

    def getLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getLHS())

    def getRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getRHS())

    def getAllSpecs(self, ):
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, double kcst):
        self.ptr().setKcst(kcst)

    def getOrder(self):
        return self.ptr().getOrder()

    @staticmethod
    cdef _py_ComplexReac from_ptr(ComplexReac *ptr):
        cdef _py_ComplexReac obj = _py_ComplexReac.__new__(_py_ComplexReac)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_ComplexReac from_ref(const ComplexReac &ref):
        return _py_ComplexReac.from_ptr(<ComplexReac*>&ref)

    @staticmethod
    cdef list flat_set2list(flat_set[ComplexReac*] vec):
        return [ _py_ComplexReac.from_ptr(elem) for elem in vec ]

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_ComplexSReac(_py__base):
    "Python wrapper class for ComplexReac"
# ----------------------------------------------------------------------------------------------------------------------
    cdef ComplexSReac *ptr(self):
        return <ComplexSReac*> self._ptr

    def __init__(self, str id, _py_Surfsys surfsys, list ilhs=[], list slhs=[], list olhs=[], list irhs=[], list srhs=[], list orhs=[], list icompEvs=[], list scompEvs=[], list ocompEvs=[], double kcst=0):
        cdef std.vector[Spec*] _ilhs, _slhs, _olhs, _irhs, _srhs, _orhs
        _py_Spec.list2vector(ilhs, &_ilhs)
        _py_Spec.list2vector(slhs, &_slhs)
        _py_Spec.list2vector(olhs, &_olhs)
        _py_Spec.list2vector(irhs, &_irhs)
        _py_Spec.list2vector(srhs, &_srhs)
        _py_Spec.list2vector(orhs, &_orhs)

        cdef std.vector[ComplexEvent*] _icompEvs, _scompEvs, _ocompEvs
        for ev in icompEvs:
            _icompEvs.push_back((<_py_ComplexEvent>ev).ptr())
        for ev in scompEvs:
            _scompEvs.push_back((<_py_ComplexEvent>ev).ptr())
        for ev in ocompEvs:
            _ocompEvs.push_back((<_py_ComplexEvent>ev).ptr())

        self._ptr = new ComplexSReac(to_std_string(id), deref(surfsys.ptr()), _ilhs, _slhs, _olhs, _irhs, _srhs, _orhs, _icompEvs, _scompEvs, _ocompEvs, kcst )

    def getID(self, ):
        return from_std_string(self.ptr().getID())

    def getSurfsys(self, ):
        return _py_Surfsys.from_ptr(&self.ptr().getSurfsys())

    def getModel(self, ):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getOLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getOLHS())

    def getILHS(self, ):
        return _py_Spec.vector2list(self.ptr().getILHS())

    def getSLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getSLHS())

    def getIRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getIRHS())

    def getSRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getSRHS())

    def getORHS(self, ):
        return _py_Spec.vector2list(self.ptr().getORHS())

    def getAllSpecs(self, ):
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, double kcst):
        self.ptr().setKcst(kcst)

    def getOrder(self):
        return self.ptr().getOrder()

    @staticmethod
    cdef _py_ComplexSReac from_ptr(ComplexSReac *ptr):
        cdef _py_ComplexSReac obj = _py_ComplexSReac.__new__(_py_ComplexSReac)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef _py_ComplexSReac from_ref(const ComplexSReac &ref):
        return _py_ComplexSReac.from_ptr(<ComplexSReac*>&ref)

    @staticmethod
    cdef list vector2list(std.vector[ComplexSReac*] vec):
        return [ _py_ComplexSReac.from_ptr(elem) for elem in vec ]


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

        Construct an ohmic current object with identifier string id
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
        self._ptr = new OhmicCurr(to_std_string(id), deref(surfsys.ptr()), deref(chanstate.ptrx()), erev, g)

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
        return _py_Surfsys.from_ptr(&self.ptr().getSurfsys())

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
        return _py_Model.from_ptr(&self.ptr().getModel())

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
        return _py_ChanState.from_ptr(&self.ptr().getChanState())

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
        self.ptr().setChanState(deref(chanstate.ptrx()))

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
        self._ptr = new GHKcurr(to_std_string(id), deref(surfsys.ptr()), deref(chanstate.ptrx()), deref(ion.ptr()), computeflux, virtual_oconc, vshift)
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
        return _py_Surfsys.from_ptr(&self.ptr().getSurfsys())

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
        return _py_Model.from_ptr(&self.ptr().getModel())

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
        return _py_ChanState.from_ptr(&self.ptr().getChanState())

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
        self.ptr().setChanState(deref(chanstate.ptrx()))

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
        return _py_Spec.from_ptr(&self.ptr().getIon())

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
        self.ptr().setIon(deref(ion.ptr()))

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

    def getP(self):
        """
        Get the single-channel permeability (units: cubic meters / second).

        Syntax::

            getP()

        Arguments:
        None

        Return:
        float
        """
        return self.ptr().getP()

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


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Vesicle(_py__base):
    "Python wrapper class for Vesicle"
# ----------------------------------------------------------------------------------------------------------------------
    cdef Vesicle *ptr(self):
        return <Vesicle*> self._ptr

    def __init__(self, str id, _py_Model model, double diameter, double dcst=0):
        self._ptr = new Vesicle(to_std_string(id), deref(model.ptr()), diameter, dcst)

    def getID(self):
        return from_std_string(self.ptr().getID())

    def addVesSurfsys(self, name):
        self.ptr().addVesSurfsys(to_std_string(name))

    def getVesSurfsys(self):
        return string_flat_set_to_list(self.ptr().getVesSurfsys())

    def getDiameter(self):
        return self.ptr().getDiameter()

    def getDcst(self):
        return self.ptr().getDcst()

    def setDcst(self, double dcst):
        return self.ptr().setDcst(dcst)

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    @staticmethod
    cdef _py_Vesicle from_ptr(Vesicle *ptr):
        if ptr == NULL:
            return None
        cdef _py_Vesicle obj = _py_Vesicle.__new__(_py_Vesicle)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[Vesicle*] vec):
        return [ _py_Vesicle.from_ptr(elem) for elem in vec ]


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_VesSurfsys(_py__base):
    "Python wrapper class for VesSurfsys"
# ----------------------------------------------------------------------------------------------------------------------
    cdef VesSurfsys *ptr(self):
        return <VesSurfsys*> self._ptr

    def __init__(self, str id, _py_Model model):
        self._ptr = new VesSurfsys(to_std_string(id), deref(model.ptr()))

    def getID(self):
        return from_std_string(self.ptr().getID())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getVesSReac(self, id):
        return _py_VesSReac.from_ptr(&self.ptr().getVesSReac(to_std_string(id)))

    def getAllVesSReacs(self):
        return _py_VesSReac.vector2list(self.ptr().getAllVesSReacs())

    def getVesSDiff(self, id):
        return _py_VesSDiff.from_ptr(&self.ptr().getVesSDiff(to_std_string(id)))

    def getAllVesSDiffs(self):
        return _py_VesSDiff.vector2list(self.ptr().getAllVesSDiffs())

    def getExocytosis(self, id):
        return _py_Exocytosis.from_ptr(&self.ptr().getExocytosis(to_std_string(id)))

    def getAllExocytosis(self):
        return _py_Exocytosis.vector2list(self.ptr().getAllExocytosis())

    def getAllSpecs(self):
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_VesSurfsys from_ptr(VesSurfsys *ptr):
        if ptr == NULL:
            return None
        cdef _py_VesSurfsys obj = _py_VesSurfsys.__new__(_py_VesSurfsys)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[VesSurfsys*] vec):
        return [ _py_VesSurfsys.from_ptr(elem) for elem in vec ]


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Raft(_py__base):
    "Python wrapper class for Raft"
# ----------------------------------------------------------------------------------------------------------------------
    cdef Raft *ptr(self):
        return <Raft*> self._ptr

    def __init__(self, str id, _py_Model model, double diameter, double dcst=0):
        self._ptr = new Raft(to_std_string(id), deref(model.ptr()), diameter, dcst)

    def getID(self):
        return from_std_string(self.ptr().getID())

    def addRaftsys(self, name):
        self.ptr().addRaftsys(to_std_string(name))

    def getRaftsys(self):
        return string_flat_set_to_list(self.ptr().getRaftsys())

    def getDiameter(self):
        return self.ptr().getDiameter()

    def getDcst(self):
        return self.ptr().getDcst()

    def setDcst(self, double dcst):
        return self.ptr().setDcst(dcst)

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    @staticmethod
    cdef _py_Raft from_ptr(Raft *ptr):
        if ptr == NULL:
            return None
        cdef _py_Raft obj = _py_Raft.__new__(_py_Raft)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[Raft*] vec):
        return [ _py_Raft.from_ptr(elem) for elem in vec ]


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Raftsys(_py__base):
    "Python wrapper class for Raftsys"
# ----------------------------------------------------------------------------------------------------------------------
    cdef Raftsys *ptr(self):
        return <Raftsys*> self._ptr

    def __init__(self, str id, _py_Model model):
        self._ptr = new Raftsys(to_std_string(id), deref(model.ptr()))

    def getID(self):
        return from_std_string(self.ptr().getID())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getRaftSReac(self, id):
        return _py_RaftSReac.from_ptr(&self.ptr().getRaftSReac(to_std_string(id)))

    def getAllRaftSReacs(self):
        return _py_RaftSReac.vector2list(self.ptr().getAllRaftSReacs())

    def getRaftEndocytosis(self, id):
        return _py_RaftEndocytosis.from_ptr(&self.ptr().getRaftEndocytosis(to_std_string(id)))

    def getAllRaftEndocytosiss(self):
        return _py_RaftEndocytosis.vector2list(self.ptr().getAllRaftEndocytosiss())

    def getRaftDis(self, id):
        return _py_RaftDis.from_ptr(&self.ptr().getRaftDis(to_std_string(id)))

    def getAllRaftDiss(self):
        return _py_RaftDis.vector2list(self.ptr().getAllRaftDiss())

    def getAllSpecs(self):
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_Raftsys from_ptr(Raftsys *ptr):
        cdef _py_Raftsys obj = _py_Raftsys.__new__(_py_Raftsys)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[Raftsys*] vec):
        return [ _py_Raftsys.from_ptr(elem) for elem in vec ]


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Endocytosis(_py__base):
    "Python wrapper class for Endocytosis"
# ----------------------------------------------------------------------------------------------------------------------
    cdef Endocytosis *ptr(self):
        return <Endocytosis*> self._ptr

    def __init__(self, str id, _py_Surfsys ssys, _py_Vesicle irhs=None, _py_Vesicle orhs=None, list deps=[], double kcst=0):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(deps, &vec)

        self._ptr = new Endocytosis(to_std_string(id), deref(ssys.ptr()), irhs.ptr() if irhs is not None else <Vesicle *>NULL,
            orhs.ptr() if orhs is not None else <Vesicle *>NULL, vec, kcst
        )

    def getID(self):
        return from_std_string(self.ptr().getID())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getInner(self):
        return self.ptr().getInner()

    def getSpecDeps(self):
        return _py_Spec.vector2list(self.ptr().getSpecDeps())

    def getIRHS(self):
        return _py_Vesicle.from_ptr(&self.ptr().getIRHS())

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, kcst):
        return self.ptr().setKcst(kcst)

    def getAllSpecs(self):
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

    def getSurfsys(self):
        return _py_Surfsys.from_ptr(&self.ptr().getSurfsys())

    @staticmethod
    cdef _py_Endocytosis from_ptr(Endocytosis *ptr):
        cdef _py_Endocytosis obj = _py_Endocytosis.__new__(_py_Endocytosis)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[Endocytosis*] vec):
        return [ _py_Endocytosis.from_ptr(elem) for elem in vec ]


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_Exocytosis(_py__base):
    "Python wrapper class for Exocytosis"
# ----------------------------------------------------------------------------------------------------------------------
    cdef Exocytosis *ptr(self):
        return <Exocytosis*> self._ptr

    def __init__(self, str id, _py_VesSurfsys ssys, list deps=[], _py_Raft raft=None, double kcst=0, bool kiss_and_run=False, dict kiss_and_run_spec_changes={}, double kiss_and_run_partial_release=1.0):
        cdef std.vector[Spec*] dep_vec
        _py_Spec.list2vector(deps, &dep_vec)
        cdef std.map[SpecP, SpecP] knr_map
        _py_Spec.dict2map(kiss_and_run_spec_changes, &knr_map)

        self._ptr = new Exocytosis(
            to_std_string(id), deref(ssys.ptr()), dep_vec, raft.ptr() if raft is not None else <Raft *>NULL, kcst, kiss_and_run, knr_map, kiss_and_run_partial_release
        )

    def getID(self):
        return from_std_string(self.ptr().getID())

    def getSpecDeps(self):
        return _py_Spec.vector2list(self.ptr().getSpecDeps())

    def getVesSurfsys(self):
        return _py_VesSurfsys.from_ptr(&self.ptr().getVesSurfsys())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getRaft(self):
        return _py_Raft.from_ptr(self.ptr().getRaft())

    def getKissAndRun(self):
        return self.ptr().getKissAndRun()

    def getKissAndRunSpecChanges(self):
        return _py_Spec.map2dict(self.ptr().getKissAndRunSpecChanges())

    def getKissAndRunPartRelease(self):
        return self.ptr().getKissAndRunPartRelease()

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, kcst):
        return self.ptr().setKcst(kcst)

    @staticmethod
    cdef _py_Exocytosis from_ptr(Exocytosis *ptr):
        cdef _py_Exocytosis obj = _py_Exocytosis.__new__(_py_Exocytosis)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[Exocytosis*] vec):
        return [ _py_Exocytosis.from_ptr(elem) for elem in vec ]


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_RaftEndocytosis(_py__base):
    "Python wrapper class for RaftEndocytosis"
# ----------------------------------------------------------------------------------------------------------------------
    cdef RaftEndocytosis *ptr(self):
        return <RaftEndocytosis*> self._ptr

    def __init__(self, str id, _py_Raftsys ssys, _py_Vesicle irhs, _py_Vesicle orhs=None, list deps=[], double kcst=0):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(deps, &vec)

        self._ptr = new RaftEndocytosis(to_std_string(id), deref(ssys.ptr()), irhs.ptr() if irhs is not None else <Vesicle *>NULL,
            orhs.ptr() if orhs is not None else <Vesicle *>NULL, vec, kcst
        )

    def getID(self):
        return from_std_string(self.ptr().getID())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getRaftsys(self):
        return _py_Raftsys.from_ptr(&self.ptr().getRaftsys())

    def getInner(self):
        return self.ptr().getInner()

    def getSpecDeps(self):
        return _py_Spec.vector2list(self.ptr().getSpecDeps())

    def getRHS(self):
        return _py_Vesicle.from_ptr(&self.ptr().getRHS())

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, kcst):
        return self.ptr().setKcst(kcst)

    def getAllSpecs(self):
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_RaftEndocytosis from_ptr(RaftEndocytosis *ptr):
        cdef _py_RaftEndocytosis obj = _py_RaftEndocytosis.__new__(_py_RaftEndocytosis)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[RaftEndocytosis*] vec):
        return [ _py_RaftEndocytosis.from_ptr(elem) for elem in vec ]


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_RaftGen(_py__base):
    "Python wrapper class for RaftGen"
# ----------------------------------------------------------------------------------------------------------------------
    cdef RaftGen *ptr(self):
        return <RaftGen*> self._ptr

    def __init__(self, str id, _py_Surfsys ssys, list deps, _py_Raft raft, double kcst=0):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(deps, &vec)

        self._ptr = new RaftGen(
            to_std_string(id), deref(ssys.ptr()), vec, deref(raft.ptr()), kcst
        )

    def getID(self):
        return from_std_string(self.ptr().getID())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getSurfsys(self):
        return _py_Surfsys.from_ptr(&self.ptr().getSurfsys())

    def getSpecDeps(self):
        return _py_Spec.vector2list(self.ptr().getSpecSignature())

    def getRaft(self):
        return _py_Raft.from_ptr(&self.ptr().getRaft())

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, kcst):
        return self.ptr().setKcst(kcst)

    def getAllSpecs(self):
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_RaftGen from_ptr(RaftGen *ptr):
        cdef _py_RaftGen obj = _py_RaftGen.__new__(_py_RaftGen)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[RaftGen*] vec):
        return [ _py_RaftGen.from_ptr(elem) for elem in vec ]


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_RaftDis(_py__base):
    "Python wrapper class for RaftDis"
# ----------------------------------------------------------------------------------------------------------------------
    cdef RaftDis *ptr(self):
        return <RaftDis*> self._ptr

    def __init__(self, str id, _py_Raftsys ssys, list antideps, double kcst=0):
        cdef std.vector[Spec*] vec
        _py_Spec.list2vector(antideps, &vec)

        self._ptr = new RaftDis(
            to_std_string(id), deref(ssys.ptr()), vec, kcst
        )

    def getID(self):
        return from_std_string(self.ptr().getID())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getRaftsys(self):
        return _py_Raftsys.from_ptr(&self.ptr().getRaftsys())

    def getRaft(self):
        return _py_Raft.from_ptr(self.ptr().getRaft())

    def getSpecAntiDeps(self):
        return _py_Spec.vector2list(self.ptr().getSpecSignature())

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, kcst):
        return self.ptr().setKcst(kcst)

    def getAllSpecs(self):
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_RaftDis from_ptr(RaftDis *ptr):
        cdef _py_RaftDis obj = _py_RaftDis.__new__(_py_RaftDis)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[RaftDis*] vec):
        return [ _py_RaftDis.from_ptr(elem) for elem in vec ]


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_LinkSpec(_py__base):
    "Python wrapper class for LinkSpec"
# ----------------------------------------------------------------------------------------------------------------------
    cdef LinkSpec *ptr(self):
        return <LinkSpec*> self._ptr

    def __init__(self, str id, _py_Model model, double dcst=0):
        self._ptr = new LinkSpec(to_std_string(id), deref(model.ptr()), dcst)

    def getID(self, ):
        return from_std_string(self.ptr().getID())

    def getDcst(self):
        return self.ptr().getDcst()

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    @staticmethod
    cdef _py_LinkSpec from_ptr(LinkSpec *ptr):
        cdef _py_LinkSpec obj = _py_LinkSpec.__new__(_py_LinkSpec)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef std.vector[LinkSpec*] *list2vector(list specList, std.vector[LinkSpec*] *dstVec):
        for item in specList:
            assert isinstance(item, _py_LinkSpec), "Wrong type of spec: " + str(type(item))
            dstVec.push_back( (<_py_LinkSpec>item).ptr())
        return dstVec

    @staticmethod
    cdef list vector2list(std.vector[LinkSpec*] specVec):
        return [ _py_LinkSpec.from_ptr(elem) for elem in specVec ]

    @staticmethod
    cdef list flat_set2list(flat_set[LinkSpec*] specVec):
        return [ _py_LinkSpec.from_ptr(elem) for elem in specVec ]

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_VesSReac(_py__base):
    "Python wrapper class for VesSReac"
# ----------------------------------------------------------------------------------------------------------------------
    cdef VesSReac *ptr(self):
        return <VesSReac*> self._ptr

    def __init__(self, str id, _py_VesSurfsys surfsys,
            list olhs=[], list slhs=[], list vlhs=[], list llhs=[],
            list lrhs=[], list vrhs=[], list srhs=[], list orhs=[], list irhs=[],
            list vdeps=[], double kcst=0, immobilization=_py_Immobilization.NO_EFFECT,
            double max_dist=-1, immobility=None
        ):
        if immobility is not None:
            ShowDeprecationWarning(
                item='The `immobility` keyword argument',
                replacement='`immobilization`',
                version='6.0'
            )
            immobilization = _py_Immobilization(immobility)
        if not isinstance(immobilization, _py_Immobilization):
            raise TypeError(f'Expected an Immobilization value, got {immobilization} instead.')
        cdef std.vector[Spec*] _olhs, _slhs, _vlhs, _vrhs, _srhs, _orhs, _irhs, _vdeps
        cdef std.vector[LinkSpec*] _llhs, _lrhs
        _py_Spec.list2vector(olhs, &_olhs)
        _py_Spec.list2vector(slhs, &_slhs)
        _py_Spec.list2vector(vlhs, &_vlhs)
        _py_LinkSpec.list2vector(llhs, &_llhs)
        _py_LinkSpec.list2vector(lrhs, &_lrhs)
        _py_Spec.list2vector(vrhs, &_vrhs)
        _py_Spec.list2vector(srhs, &_srhs)
        _py_Spec.list2vector(orhs, &_orhs)
        _py_Spec.list2vector(irhs, &_irhs)
        _py_Spec.list2vector(vdeps, &_vdeps)
        # Instantiate
        self._ptr = new VesSReac(
            to_std_string(id), deref(surfsys.ptr()),
            _olhs, _slhs, _vlhs, _llhs,
            _lrhs, _vrhs, _srhs, _orhs, _irhs,
            _vdeps, kcst, immobilization.value, max_dist
        )

    def getID(self, ):
        return from_std_string(self.ptr().getID())

    def getVesSurfsys(self):
        return _py_VesSurfsys.from_ptr(&self.ptr().getVesSurfsys())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getOLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getOLHS())

    def getSLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getSLHS())

    def getVLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getVLHS())

    def getLLHS(self, ):
        return _py_LinkSpec.vector2list(self.ptr().getLLHS())

    def getLRHS(self, ):
        return _py_LinkSpec.vector2list(self.ptr().getLRHS())

    def getVRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getVRHS())

    def getSRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getSRHS())

    def getORHS(self, ):
        return _py_Spec.vector2list(self.ptr().getORHS())

    def getIRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getIRHS())

    def getVDeps(self, ):
        return _py_Spec.vector2list(self.ptr().getVDeps())

    def getOrder(self):
        return self.ptr().getOrder()

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, double kcst):
        self.ptr().setKcst(kcst)

    def getImmobilization(self):
        return _py_Immobilization(self.ptr().getImmobilization())

    def getMaxDistance(self):
        return self.ptr().getMaxDistance()

    def getAllSpecs(self):
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

    def getAllLinkSpecs(self):
        return _py_LinkSpec.flat_set2list(self.ptr().getAllLinkSpecs())

    @staticmethod
    cdef _py_VesSReac from_ptr(VesSReac *ptr):
        cdef _py_VesSReac obj = _py_VesSReac.__new__(_py_VesSReac)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[VesSReac*] vec):
        return [ _py_VesSReac.from_ptr(elem) for elem in vec ]


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_RaftSReac(_py__base):
    "Python wrapper class for RaftSReac"
# ----------------------------------------------------------------------------------------------------------------------
    cdef RaftSReac *ptr(self):
        return <RaftSReac*> self._ptr

    def __init__(self, str id, _py_Raftsys surfsys,
            list ilhs=[], list olhs=[], list slhs=[], list rlhs=[],
            list rrhs=[], list srhs=[], list orhs=[], list irhs=[],
            list rdeps=[], list anti_rdeps=[], double kcst=0,
            immobilization=_py_Immobilization.NO_EFFECT, immobility=None
        ):
        if immobility is not None:
            ShowDeprecationWarning(
                item='The `immobility` keyword argument',
                replacement='`immobilization`',
                version='6.0'
            )
            immobilization = _py_Immobilization(immobility)
        if not isinstance(immobilization, _py_Immobilization):
            raise TypeError(f'Expected an Immobilization value, got {immobilization} instead.')
        cdef std.vector[Spec*] _ilhs, _olhs, _slhs, _rlhs, _rrhs, _srhs, _orhs, _irhs, _rdeps, _anti_rdeps
        _py_Spec.list2vector(ilhs, &_ilhs)
        _py_Spec.list2vector(olhs, &_olhs)
        _py_Spec.list2vector(slhs, &_slhs)
        _py_Spec.list2vector(rlhs, &_rlhs)
        _py_Spec.list2vector(rrhs, &_rrhs)
        _py_Spec.list2vector(srhs, &_srhs)
        _py_Spec.list2vector(orhs, &_orhs)
        _py_Spec.list2vector(irhs, &_irhs)
        _py_Spec.list2vector(rdeps, &_rdeps)
        _py_Spec.list2vector(anti_rdeps, &_anti_rdeps)
        # Instantiate
        self._ptr = new RaftSReac(
            to_std_string(id), deref(surfsys.ptr()),
            _ilhs, _olhs, _slhs, _rlhs,
            _rrhs, _srhs, _orhs, _irhs,
            _rdeps, _anti_rdeps, kcst, immobilization.value
        )

    def getID(self, ):
        return from_std_string(self.ptr().getID())

    def getRaftsys(self):
        return _py_Raftsys.from_ptr(&self.ptr().getRaftsys())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getInner(self):
        return self.ptr().getInner()

    def getOuter(self):
        return self.ptr().getOuter()

    def getILHS(self, ):
        return _py_Spec.vector2list(self.ptr().getILHS())

    def getOLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getOLHS())

    def getSLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getSLHS())

    def getRLHS(self, ):
        return _py_Spec.vector2list(self.ptr().getRsLHS())

    def getRRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getRsRHS())

    def getSRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getSRHS())

    def getORHS(self, ):
        return _py_Spec.vector2list(self.ptr().getORHS())

    def getIRHS(self, ):
        return _py_Spec.vector2list(self.ptr().getIRHS())

    def getRDeps(self, ):
        return _py_Spec.vector2list(self.ptr().getRsDeps())

    def getAntiRDeps(self, ):
        return _py_Spec.vector2list(self.ptr().getAntiRsDeps())

    def getOrder(self):
        return self.ptr().getOrder()

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, double kcst):
        self.ptr().setKcst(kcst)

    def getImmobilization(self):
        return _py_Immobilization(self.ptr().getImmobilization())

    def getAllSpecs(self):
        return _py_Spec.flat_set2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_RaftSReac from_ptr(RaftSReac *ptr):
        cdef _py_RaftSReac obj = _py_RaftSReac.__new__(_py_RaftSReac)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[RaftSReac*] vec):
        return [ _py_RaftSReac.from_ptr(elem) for elem in vec ]

    @staticmethod
    cdef list flat_set2list(flat_set[RaftSReac*] vec):
        return [ _py_RaftSReac.from_ptr(elem) for elem in vec ]

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_VesSDiff(_py__base):
    "Python wrapper class for VesSDiff"
# ----------------------------------------------------------------------------------------------------------------------
    cdef VesSDiff *ptr(self):
        return <VesSDiff*> self._ptr

    def __init__(self, str id, _py_VesSurfsys vssys, _py_Spec lig, double dcst=0):
        self._ptr = new VesSDiff(to_std_string(id), deref(vssys.ptr()), deref(lig.ptr()), dcst)

    def getID(self):
        return from_std_string(self.ptr().getID())

    def getVesSurfsys(self):
        return _py_VesSurfsys.from_ptr(&self.ptr().getVesSurfsys())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getLig(self):
        return _py_Spec.from_ptr(&self.ptr().getLig())

    def getDcst(self):
        return self.ptr().getDcst()

    def setDcst(self, double dcst):
        return self.ptr().setDcst(dcst)

    def getAllSpecs(self):
        return _py_Spec.vector2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_VesSDiff from_ptr(VesSDiff *ptr):
        cdef _py_VesSDiff obj = _py_VesSDiff.__new__(_py_VesSDiff)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[VesSDiff*] vec):
        return [ _py_VesSDiff.from_ptr(elem) for elem in vec ]


# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_VesBind(_py__base):
    "Python wrapper class for VesBind"
# ----------------------------------------------------------------------------------------------------------------------
    cdef VesBind *ptr(self):
        return <VesBind*> self._ptr

    def __init__(self, str id, _py_Volsys volsys,
            _py_Vesicle ves1, _py_Spec r1, _py_Vesicle ves2, _py_Spec r2,
            _py_LinkSpec l1, _py_LinkSpec l2, double length_max, double length_min,
            list vdeps1=[], list vdeps2=[], list ldeps1=[], list ldeps2=[],
            double kcst=0, immobilization=_py_Immobilization.NO_EFFECT, immobility=None
        ):
        if immobility is not None:
            ShowDeprecationWarning(
                item='The `immobility` keyword argument',
                replacement='`immobilization`',
                version='6.0'
            )
            immobilization = _py_Immobilization(immobility)
        if not isinstance(immobilization, _py_Immobilization):
            raise TypeError(f'Expected an Immobilization value, got {immobilization} instead.')
        cdef std.vector[Spec*] _vdeps1, _vdeps2
        cdef std.vector[LinkSpec*] _ldeps1, _ldeps2
        _py_Spec.list2vector(vdeps1, &_vdeps1)
        _py_Spec.list2vector(vdeps2, &_vdeps2)
        _py_LinkSpec.list2vector(ldeps1, &_ldeps1)
        _py_LinkSpec.list2vector(ldeps2, &_ldeps2)
        # Instantiate
        self._ptr = new VesBind(
            to_std_string(id), deref(volsys.ptr()),
            deref(ves1.ptr()), deref(r1.ptr()), deref(ves2.ptr()), deref(r2.ptr()), deref(l1.ptr()), deref(l2.ptr()),
            length_max, length_min, _vdeps1, _vdeps2, _ldeps1, _ldeps2,
            kcst, immobilization.value
        )

    def getID(self):
        return from_std_string(self.ptr().getID())

    def getVolsys(self):
        return _py_Volsys.from_ptr(&self.ptr().getVolsys())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getVesicle1(self):
        return _py_Vesicle.from_ptr(self.ptr().getReactants1().first)

    def getVesicle2(self):
        return _py_Vesicle.from_ptr(self.ptr().getReactants2().first)

    def getReactant1(self):
        return _py_Spec.from_ptr(self.ptr().getReactants1().second)

    def getReactant2(self):
        return _py_Spec.from_ptr(self.ptr().getReactants2().second)

    def getProduct1(self):
        return _py_LinkSpec.from_ptr(self.ptr().getProducts1().second)

    def getProduct2(self):
        return _py_LinkSpec.from_ptr(self.ptr().getProducts2().second)

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, double kcst):
        self.ptr().setKcst(kcst)

    def getLengthMin(self):
        return self.ptr().getLengthMin()

    def getLengthMax(self):
        return self.ptr().getLengthMax()

    def getImmobilization(self):
        return _py_Immobilization(self.ptr().getImmobilization())

    def getVDeps1(self):
        return _py_Spec.vector2list(self.ptr().getVDeps1())

    def getVDeps2(self):
        return _py_Spec.vector2list(self.ptr().getVDeps2())

    def getLDeps1(self):
        return _py_LinkSpec.vector2list(self.ptr().getLDeps1())

    def getLDeps2(self):
        return _py_LinkSpec.vector2list(self.ptr().getLDeps2())

    def getAllSpecs(self):
        return _py_Spec.vector2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_VesBind from_ptr(VesBind *ptr):
        cdef _py_VesBind obj = _py_VesBind.__new__(_py_VesBind)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[VesBind*] vec):
        return [ _py_VesBind.from_ptr(elem) for elem in vec ]

    @staticmethod
    cdef list flat_set2list(flat_set[VesBind*] vec):
        return [ _py_VesBind.from_ptr(elem) for elem in vec ]

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_VesUnbind(_py__base):
    "Python wrapper class for VesUnbind"
# ----------------------------------------------------------------------------------------------------------------------
    cdef VesUnbind *ptr(self):
        return <VesUnbind*> self._ptr

    def __init__(self, str id, _py_Volsys volsys,
            _py_LinkSpec l1, _py_LinkSpec l2,
            _py_Vesicle ves1, _py_Spec p1, _py_Vesicle ves2, _py_Spec p2,
            double kcst=0, immobilization=_py_Immobilization.NO_EFFECT, immobility=None
        ):
        if immobility is not None:
            ShowDeprecationWarning(
                item='The `immobility` keyword argument',
                replacement='`immobilization`',
                version='6.0'
            )
            immobilization = _py_Immobilization(immobility)
        if not isinstance(immobilization, _py_Immobilization):
            raise TypeError(f'Expected an Immobilization value, got {immobilization} instead.')
        # Instantiate
        self._ptr = new VesUnbind(
            to_std_string(id), deref(volsys.ptr()),
            deref(l1.ptr()), deref(l2.ptr()), deref(ves1.ptr()), deref(p1.ptr()), deref(ves2.ptr()), deref(p2.ptr()),
            kcst, immobilization.value
        )

    def getID(self):
        return from_std_string(self.ptr().getID())

    def getVolsys(self):
        return _py_Volsys.from_ptr(&self.ptr().getVolsys())

    def getModel(self):
        return _py_Model.from_ptr(&self.ptr().getModel())

    def getVesicle1(self):
        return _py_Vesicle.from_ptr(self.ptr().getProducts1().first)

    def getVesicle2(self):
        return _py_Vesicle.from_ptr(self.ptr().getProducts2().first)

    def getLink1(self):
        return _py_LinkSpec.from_ptr(self.ptr().getLinks1().second)

    def getLink2(self):
        return _py_LinkSpec.from_ptr(self.ptr().getLinks2().second)

    def getProduct1(self):
        return _py_Spec.from_ptr(self.ptr().getProducts1().second)

    def getProduct2(self):
        return _py_Spec.from_ptr(self.ptr().getProducts2().second)

    def getKcst(self):
        return self.ptr().getKcst()

    def setKcst(self, double kcst):
        self.ptr().setKcst(kcst)

    def getImmobilization(self):
        return _py_Immobilization(self.ptr().getImmobilization())

    def getAllSpecs(self):
        return _py_Spec.vector2list(self.ptr().getAllSpecs())

    @staticmethod
    cdef _py_VesUnbind from_ptr(VesUnbind *ptr):
        cdef _py_VesUnbind obj = _py_VesUnbind.__new__(_py_VesUnbind)
        obj._ptr = ptr
        return obj

    @staticmethod
    cdef list vector2list(std.vector[VesUnbind*] vec):
        return [ _py_VesUnbind.from_ptr(elem) for elem in vec ]

    @staticmethod
    cdef list flat_set2list(flat_set[VesUnbind*] vec):
        return [ _py_VesUnbind.from_ptr(elem) for elem in vec ]
