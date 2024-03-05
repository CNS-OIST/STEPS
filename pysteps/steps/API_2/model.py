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

import collections
import copy
from enum import Enum
import inspect
import itertools
import numbers
import warnings

from steps import stepslib

from . import geom as ngeom
from . import utils as nutils


try:
    COMPLEX_FILTER_MAX_VALUE = stepslib._py_COMPLEX_FILTER_MAX_VALUE
except AttributeError as ex:
    # If the installed STEPS version does not support complexes
    COMPLEX_FILTER_MAX_VALUE = 0
    stepslib._py_Complex = None
    stepslib._py_ComplexReac = None

__all__ = [
    'Model',
    'VolumeSystem',
    'SurfaceSystem',
    'VesicleSurfaceSystem',
    'RaftSurfaceSystem',
    'Vesicle',
    'Raft',
    'Species',
    'LinkSpecies',
    'Complex',
    'Channel',
    'Reaction',
    'Diffusion',
    'Endocytosis',
    'Exocytosis',
    'RaftEndocytosis',
    'RaftGen',
    'RaftDis',
    'VesicleBind',
    'VesicleUnbind',
    'Current',
    'OhmicCurr',
    'GHKCurr',
    'Location',
    'In',
    'Out',
    'Surf',
    'VesSurf',
    'RaftSurf',
    'ReactionManager',
    'ReactionElement',
    'ReactingElement',
    'ReactionSide',
    'SubUnit',
    'SubUnitState',
    'SubUnitSelector',
    'ComplexState',
    'ComplexSelector',
    'NoOrdering',
    'StrongOrdering',
    'RotationalSymmetryOrdering',
    'XDepFunc',
    'VDepFunc',
    'CompDepFunc',
    'VDepRate',
    'CompDepRate',
    'CompDepCond',
    'CompDepP',
    'CompDepDcst',
    'IMMOBILIZING',
    'MOBILIZING',
    'NO_EFFECT',
]


###################################################################################################
# Enums

try:
    Immobilization = stepslib._py_Immobilization
except AttributeError:
    # Older versions of cython do not export this enum
    class Immobilization(Enum):
        IMMOBILIZING = 0
        MOBILIZING   = 1
        NO_EFFECT    = 2

# Aliases
IMMOBILIZING = Immobilization.IMMOBILIZING
MOBILIZING = Immobilization.MOBILIZING
NO_EFFECT = Immobilization.NO_EFFECT


###################################################################################################
# Model


@nutils.FreezeAfterInit
class Model(nutils.UsableObject, nutils.StepsWrapperObject):
    """Top-level container for the objects in a kinetic model

    Should be used as a context manager for the declaration of Volume systems, reactions,
    diffusion rules, etc::

        mdl = Model()
        with mdl:
            vsys = VolumeSystem.Create()
            ... # Declare other objects in mdl

    After having declared children objects in the ``with mdl:`` block, they can be accessed
    as attributes of ``mdl`` with their name (see :py:class:`steps.API_2.utils.NamedObject` and
    :py:func:`steps.API_2.utils.NamedObject.Create`)::

        mdl.vsys
    """

    def __init__(self, *args, _createObj=True, **kwargs):
        super().__init__(*args, **kwargs)
        self.stepsModel = self._createStepsObj() if _createObj else None
        self.volSysConstraints = []

    def _SetUpMdlDeps(self, geom):
        """Set up structures that depend on objects declared in the model."""
        # Start with vesicles because they can update rafts
        for pl in self._getChildrenOfType(Vesicle):
            pl._SetUpMdlDeps(self, geom)
        for pl in self._getChildrenOfType(Raft):
            pl._SetUpMdlDeps(self, geom)
        for pl in self._getChildrenOfType(VesicleSurfaceSystem):
            pl._SetUpMdlDeps(geom)
        for pl in self._getChildrenOfType(RaftSurfaceSystem):
            pl._SetUpMdlDeps(geom)

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsModel]

    @classmethod
    def _FromStepsObject(cls, obj):
        """Create the interface object from a STEPS object."""
        mdl = cls(_createObj=False)
        mdl.stepsModel = obj
        with mdl:
            for spec in obj.getAllSpecs():
                Species._FromStepsObject(spec, mdl)
            for chan in obj.getAllChans():
                Channel._FromStepsObject(chan, mdl)
            for lspec in obj.getAllLinkSpecs():
                LinkSpecies._FromStepsObject(lspec, mdl)
            for ves in obj.getAllVesicles():
                Vesicle._FromStepsObject(ves, mdl)
            for raft in obj.getAllRafts():
                Raft._FromStepsObject(raft, mdl)
            for vssys in obj.getAllVesSurfsyss():
                VesicleSurfaceSystem._FromStepsObject(vssys, mdl)
            for rssys in obj.getAllRaftsyss():
                RaftSurfaceSystem._FromStepsObject(rssys, mdl)
            for vsys in obj.getAllVolsyss():
                VolumeSystem._FromStepsObject(vsys, mdl)
            for ssys in obj.getAllSurfsyss():
                SurfaceSystem._FromStepsObject(ssys, mdl)
        return mdl

    def _addVolSysConstraint(self, reac, volSys, surfSys, loc):
        self.volSysConstraints.append((reac, volSys, surfSys, loc))

    def _createStepsObj(self):
        return stepslib._py_Model()


###################################################################################################
# Space systems


@nutils.FreezeAfterInit
class SpaceSystem(nutils.UsableObject, nutils.UsingObjects(Model), nutils.StepsWrapperObject):
    """Base class for volume and surface system classes

    :meta private:
    """

    def __init__(self, *args, _createObj=True, **kwargs):
        super().__init__(*args, **kwargs)
        (mdl,) = self._getUsedObjects()
        self.stepsSys = self._createStepsObj(mdl) if _createObj else None
        self.locations = []

    def _addLocation(self, l):
        self.locations.append(l)

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsSys]


@nutils.FreezeAfterInit
class VolumeSystem(SpaceSystem):
    """A container that groups reactions and diffusion rules located in a volume

    Should be used as a context manager for the declaration of reactions, diffusion rules, etc::

        ...
        with mdl:
            ...
            vsys = VolumeSystem.Create()
            with vsys:
                diff1 = Diffusion.Create(...)
                ... # Declare other objects in vsys

    After having declared children objects in the ``with vsys:`` block, they can be accessed
    as attributes of ``vsys`` with their name (see :py:class:`steps.API_2.utils.NamedObject` and
    :py:func:`steps.API_2.utils.NamedObject.Create`)::

        vsys.diff1
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _createStepsObj(self, mdl):
        return stepslib._py_Volsys(self.name, mdl.stepsModel)

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        vsys = cls(_createObj=False, name=obj.getID())
        vsys.stepsSys = obj
        with vsys:
            for diff in obj.getAllDiffs():
                Diffusion._FromStepsObject(diff, mdl)
            for reac in obj.getAllReacs():
                # need to convert to actual steps classes
                Reaction._FromStepsObject(reac, mdl)
            for bind in obj.getAllVesBinds():
                VesicleBind._FromStepsObject(bind, mdl)
            for unbind in obj.getAllVesUnbinds():
                VesicleUnbind._FromStepsObject(unbind, mdl)
        return vsys


@nutils.FreezeAfterInit
class SurfaceSystem(SpaceSystem):
    """A container that groups reactions and diffusion rules located in a surface

    Should be used as a context manager for the declaration of surface reactions, surface
    diffusion rules, etc::

        ...
        with mdl:
            ...
            ssys = SurfaceSystem.Create()
            with ssys:
                diff2 = Diffusion.Create(...)
                ... # Declare other objects in ssys

    After having declared children objects in the ``with ssys:`` block, they can be accessed
    as attributes of ``ssys`` with their name (see :py:class:`steps.API_2.utils.NamedObject` and
    :py:func:`steps.API_2.utils.NamedObject.Create`)::

        ssys.diff2
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _createStepsObj(self, mdl):
        return stepslib._py_Surfsys(self.name, mdl.stepsModel)

    def _addLocation(self, l):
        super()._addLocation(l)
        if isinstance(l, ngeom.Patch):
            if l.innerComp is not None:
                l.innerComp.addSystem(self, Location.IN)
            if l.outerComp is not None:
                l.outerComp.addSystem(self, Location.OUT)

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        ssys = cls(_createObj=False, name=obj.getID())
        ssys.stepsSys = obj
        with ssys:
            for diff in obj.getAllDiffs():
                Diffusion._FromStepsObject(diff, mdl)
            for sreac in obj.getAllSReacs():
                Reaction._FromStepsObject(sreac, mdl)
            for vdsreac in obj.getAllVDepSReacs():
                Reaction._FromStepsObject(vdsreac, mdl)
            # TODO not urgent: implement _FromStepsObject for GHK and Ohmic currents
            # for ghkc in obj.getAllGHKcurrs():
                # GHKCurr._FromStepsObject(ghkc, mdl)
            # for ohmc in obj.getAllOhmicCurrs():
                # OhmicCurr._FromStepsObject(ohmc, mdl)

        return ssys


@nutils.FreezeAfterInit
class VesicleSurfaceSystem(SurfaceSystem):
    """A container that groups reactions and diffusion rules located in a vesicle surface

    Should be used as a context manager for the declaration of vesicle surface reactions, vesicle surface
    diffusion rules, etc::

        ...
        with mdl:
            ...
            vssys = VesicleSurfaceSystem.Create()
            with vssys:
                diff3 = Diffusion.Create(...)
                ... # Declare other objects in vssys

    After having declared children objects in the ``with vssys:`` block, they can be accessed
    as attributes of ``vssys`` with their name (see :py:class:`steps.API_2.utils.NamedObject` and
    :py:func:`steps.API_2.utils.NamedObject.Create`)::

        vssys.diff3
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _createStepsObj(self, mdl):
        return stepslib._py_VesSurfsys(self.name, mdl.stepsModel)

    def _addLocation(self, l):
        super()._addLocation(l)

    def _SetUpMdlDeps(self, geom):
        """Set up structures that depend on objects declared in the model."""
        # Potentially all compartments can contain any vesicle
        for comp in geom.ALL(ngeom.Compartment):
            comp.addSystem(self, Location.OUT)

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        vssys = cls(_createObj=False, name=obj.getID())
        vssys.stepsSys = obj
        with vssys:
            for diff in obj.getAllVesSDiffs():
                Diffusion._FromStepsObject(diff, mdl)
            for vsreac in obj.getAllVesSReacs():
                Reaction._FromStepsObject(vsreac, mdl)
            for exc in obj.getAllExocytosis():
                Exocytosis._FromStepsObject(exc, mdl)

        return vssys


@nutils.FreezeAfterInit
class RaftSurfaceSystem(SurfaceSystem):
    """A container that groups reactions and diffusion rules located in a raft surface

    Should be used as a context manager for the declaration of raft surface reactions, raft endocytosis
    rules, etc::

        ...
        with mdl:
            ...
            rssys = RaftSurfaceSystem.Create()
            with rssys:
                rend1 = RaftEndocytosis.Create(...)
                ... # Declare other objects in rssys

    After having declared children objects in the ``with rssys:`` block, they can be accessed
    as attributes of ``rssys`` with their name (see :py:class:`steps.API_2.utils.NamedObject` and
    :py:func:`steps.API_2.utils.NamedObject.Create`)::

        rssys.rend1
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _createStepsObj(self, mdl):
        return stepslib._py_Raftsys(self.name, mdl.stepsModel)

    def _addLocation(self, l):
        super()._addLocation(l)

    def _SetUpMdlDeps(self, geom):
        """Set up structures that depend on objects declared in the model."""
        # Potentially all patches can contain any raft
        for patch in geom.ALL(ngeom.Patch):
            patch.addSystem(self)
            if patch.innerComp is not None:
                patch.innerComp.addSystem(self, Location.IN)
            if patch.outerComp is not None:
                patch.outerComp.addSystem(self, Location.OUT)

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        rssys = cls(_createObj=False, name=obj.getID())
        rssys.stepsSys = obj
        with rssys:
            for rsreac in obj.getAllRaftSReacs():
                Reaction._FromStepsObject(rsreac, mdl)
            for end in obj.getAllRaftEndocytosiss():
                RaftEndocytosis._FromStepsObject(end, mdl)
            for dis in obj.getAllRaftDiss():
                RaftDissolution._FromStepsObject(dis, mdl)

        return rssys


###################################################################################################
# Complexes and channels

####################################
# State Ordering Functions


def StrongOrdering(state):
    """Strong ordering function for complex states

    :param state: A complex state, i.e. a sequence of :py:class:`SubUnitState`
    :type state: Tuple[:py:class:`SubUnitState`]

    :returns: The equivalent unique ordering
    :rtype: Tuple[:py:class:`SubUnitState`]

    Take a tuple of :py:class:`SubUnitState` as input and return the unique equivalent based on
    the type of ordering chosen.

    StrongOrdering corresponds to no symmetry in the complex.
    """
    return tuple(state)


def NoOrdering(state):
    """No ordering function for complex states

    :param state: A complex state, i.e. a sequence of :py:class:`SubUnitState`
    :type state: Tuple[:py:class:`SubUnitState`]

    :returns: The equivalent unique ordering
    :rtype: Tuple[:py:class:`SubUnitState`]

    Take a tuple of :py:class:`SubUnitState` as input and return the unique equivalent based on
    the type of ordering chosen.

    NoOrdering considers that the subunits have no physical position, a state is just a set of
    SubUnitStates.
    """
    return tuple(sorted(state, key=lambda x: x.name))


def RotationalSymmetryOrdering(state):
    """Rotational symmetry ordering function for complex states

    :param state: A complex state, i.e. a sequence of :py:class:`SubUnitState`
    :type state: Tuple[:py:class:`SubUnitState`]

    :returns: The equivalent unique ordering
    :rtype: Tuple[:py:class:`SubUnitState`]

    Take a tuple of :py:class:`SubUnitState` as input and return the unique equivalent based on
    the type of ordering chosen.

    RotationalSymmetryOrdering represents central symmetry: state [a,b,c] is equivalent to
    [b,c,a] and [c,a,b].
    """
    # Number the states based on sorting their names
    srt = sorted(enumerate(state), key=lambda x: x[1].name)
    num = 0
    s2n = {state[srt[0][0]]: num}
    for i, s in enumerate(srt[1:]):
        if s[1] != srt[i][1]:
            num += 1
        s2n[state[s[0]]] = num
    # Compute all rotations
    n = len(state)
    rots = [tuple(state[(i + shift) % n] for i in range(n)) for shift in range(n)]
    # select the one that maximizes the number formed by the sequence of digits
    return max(rots, key=lambda state: sum(s2n[s] * (num + 1) ** p for p, s in enumerate(state)))


####################################
# Complexes


@nutils.FreezeAfterInit
class Complex(nutils.UsingObjects(Model)):
    """Multi-state complex class

    A complex is composed of one or more subunits, each subunit has a finite (preferably low)
    number of states. Subunit states can be used inside reactions so as to simplify the declaration
    of reactions for multi-subunits complexes.
    This class can also be used as a way to automatically declare complex states as species.

    :param subUnits: The list of :py:class:`SubUnit` of the complex. Can also be a list of
        :py:class:`SubUnitState` in which case a :py:class:`SubUnit` will automatically be created
        and the complex will be composed of only this :py:class:`SubUnit`.
    :type subUnits: List[:py:class:`SubUnit`] or List[:py:class:`SubUnitState`]
    :param statesAsSpecies: If True, all the possible states of the complex are declared as
        :py:class:`Species`. If not, and if your STEPS implementation supports it, the complex will be
        declared as a native STEPS multi-state complex.
    :type statesAsSpecies: bool
    :param order: A callable that defines a custom state ordering. It should map a tuple of
        :py:class:`SubUnitState` into a unique ordering. It defaults to and :py:func:`NoOrdering`,
        which considers that subunits have no spatial position. See :py:func:`StrongOrdering`,
        :py:func:`NoOrdering`, and :py:func:`RotationalSymmetryOrdering` for details.
    :type order: Callable[[Tuple[:py:class:`SubUnitState`]], Tuple[:py:class:`SubUnitState`]]

    Declaration::

        S1A, S1B, S1C, S2A, S2B = SubUnitState.Create()
        SU1, SU2 = SubUnit.Create(
            [S1A, S1B, S1C],
            [S2A, S2B]
        )
        Comp = Complex.Create([SU1, SU1, SU2])

    Usage in reactions:
        - Centered on complexes::

            # Get a single 'Comp' complex, when several complexes of the same type are used in
            # the same reaction, it allows to determine which complex undergoes which changes.
            C = Comp.get()

            # Ca binding to the first subunit of complex Comp in state S1A and changing its
            # state to S1B
            C[S1A, :, :] + Ca <r[1]> C[S1B, :, :]

        - Centered on subunit states::

            with Comp[...]:
                # Ca Binding to any subunit in state S1A and changing it to state S1B
                S1A + Ca <r[1]> S1B

                # Ca Binding to first subunit in state S1A and changing it to state S1B
                S1A[0] + Ca <r[1]> S1B[0]

            with Comp[...] as C1, Comp[..., S2B] as C2:
                S1A[C1] + S1B[C2] <r[1]> S1B[C1] + S1C[C2]

    """

    _elemStr = 'Complex'

    _SIMPATH_ONLY_CHILDREN = True

    def __init__(self, subUnits, *args, statesAsSpecies=False, order=NoOrdering, _createObj=True, **kwargs):
        """
        Declare a complex from a list of subunits.

        If statesAsSpecies is True, automatically declare each complex state as a steps Species.

        'order' allows the user to define a custom ordering, it should be callable and take
        a state (tuple of SubUnitStates) as argument and return a unique version of this state
        according to the custom ordering. It defaults to NoOrdering, which considers that subunits
        have no spatial position.
        """
        # Creating from a list of states instead of a list of subunits
        if all(isinstance(su, SubUnitState) for su in subUnits):
            subUnits = [SubUnit(subUnits)]

        self._subUnits = subUnits

        super().__init__(*args, **kwargs)

        for i, sub in enumerate(self._subUnits):
            sub._addComplex(self)
            # Add children
            self._addChildren(sub)
            for sus in sub._states:
                self._addChildren(sus)

        if order is not NoOrdering and not statesAsSpecies:
            raise Exception(
                f'Steps complexes need to be declared with "order=NoOrdering", no '
                f'other orders are accepted.'
            )
        if not statesAsSpecies and stepslib._py_Complex is None:
            warnings.warn(
                f'Your version of STEPS does not support multi-state complexes, '
                f'Complex {self.name} will instead be declared with '
                f'statesAsSpecies = True.'
            )
            statesAsSpecies = True

        self._order = order

        (mdl,) = self._getUsedObjects()
        self._statesAsSpecies = statesAsSpecies
        self.stepsComplex = None
        self._compStates = {}
        if _createObj:
            if statesAsSpecies:
                # Create all complex states as steps species
                self.stepsComplex = self._createStepsStates(mdl)
            else:
                self._initializePoolStructs()
                self.stepsComplex = self._createStepsObj(mdl)

    def _areStatesAsSpecies(self):
        return self._statesAsSpecies

    def _getAdditionalChildren(self):
        """
        Return a list of NamedObjects that should be added as children of the used object
        along with the current object.
        """
        res = set()
        for su in self._subUnits:
            res |= su._states
        return res

    def _initializePoolStructs(self):
        """Initialize structures related to real complexes."""
        allSus = set()
        for su in set(self._subUnits):
            allSus |= su._states

        self._poolElems = sorted(allSus, key=lambda x: x.name)
        self._sus2PoolInd = {sus: self._poolElems.index(sus) for sus in allSus}

    def _getSUSInd(self, sus):
        """Return the index os the sus SubUnitState."""
        if sus not in self._sus2PoolInd:
            raise Exception(f'SubUnitState {sus} is not associated with Complex {self}.')
        return self._sus2PoolInd[sus]

    def _createStepsObj(self, mdl):
        # Not implemented in steps yet
        return stepslib._py_Complex(self.name, mdl.stepsModel, len(self._subUnits), len(self._poolElems))

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        if self._statesAsSpecies:
            return [obj for state, obj in self._compStates.items()]
        else:
            return [self.stepsComplex]

    def _getAllElems(self, loc=None):
        """
        Return all elements that can directly be accessed by name from the parent physical
        location.
        """
        return [state for state, obj in self._compStates.items()]

    def _createStepsStates(self, mdl):
        """Create steps species for each complex state."""
        for comb in itertools.product(*[sorted(su._states, key=lambda x: x.name) for su in self._subUnits]):
            cs = ComplexState(self, comb)
            if cs not in self._compStates:
                spec = stepslib._py_Spec(cs.name, mdl.stepsModel)
                cs._setStepsObjects(spec)
                self._compStates[cs] = spec
        return None

    def __iter__(self):
        """Iterator over all the possible states of the complex.

        :returns: An iterator that iterates over all unique (as defined by the custom ordering)
            :py:class:`ComplexState` of the complex.
        :rtype: Iterator[:py:class:`ComplexState`]

        Usage::

            Comp1 = Complex.Create(...)

            for state in Comp1:
                # Do something with the state
                ...

        :meta public:
        """
        return iter(self.__getitem__(...))

    def _simPathWalkExpand(self):
        """Return an iterable of the elements that should be part of ResultSelector paths."""
        if self._statesAsSpecies:
            return self.__iter__()
        else:
            return [self.__getitem__(...)]

    def _simPathCombinerClass(self):
        """Return the class that needs to be used to combine expanded elements."""
        if self._statesAsSpecies:
            return nutils.SumSimPathComb
        else:
            return None

    def __getitem__(self, subStates, _selName=None):
        """Square-bracket operator for defining states or sets of states.

        :param subStates: A sequence of one or more :py:class:`SubUnitState` or
            :py:class:`SubUnitSelector` that corresponds to one or several :py:class:`ComplexState`.
        :type subStates: Tuple[Union[:py:class:`SubUnitState`, :py:class:`SubUnitSelector`], ...]

        :returns: A single complex state or a set of complex states matching the arguments.
        :rtype: Union[:py:class:`ComplexState`, :py:class:`ComplexSelector`]

        Return the :py:class:`ComplexState` corresponding to the combination of
        :py:class:`SubUnitState` given as argument. If the state is fully determined (the number and
        type of SubUnitStates must match the :py:class:`SubUnit` that have been declared when the
        complex was created), return the corresponding :py:class:`ComplexState`. Otherwise, return
        a :py:class:`ComplexSelector` representing all the possible states.

        The slice operator ('``:``') as well as ellipsis ('``...``') can be used in the same way as
        when slicing numpy arrays. See https://numpy.org/doc/1.18/reference/arrays.indexing.html
        for details.

        Valid examples::

            Comp[S1A, S1B, S2A] # Valid, returns the corresponding ComplexState
            Comp[S1A, :, :]     # Valid, returns a ComplexSelector in which the first SubUnit is
                                # in state S1A and the other subunits are in any state
            Comp[S1A, ...]      # Valid, the ellipsis operator ('...') is equivalent to filling
                                # the remaining subUnits with ':'
            Comp[..., S2A]      # Valid, returns a ComplexSelector in which the last SubUnit is in
                                # state S2A and the other subunits are in any state
            Comp[S1A|S1B, :, :] # Valid, first subunit is either in state S1A or S1B and the other
                                # subunits are in any state
            Comp[~S1C, :, :]    # Valid, first subunit is not in state S1C and the other subunits
                                # are in any state

        See :py:class:`SubUnitState` and :py:class:`SubUnitSelector` for the syntax relative to
        the last two examples.

        Invalid examples::

            Comp[S2A, S1B, S1A] # Invalid, the order of states does not match the one used for
                                # complex creation
            Comp[S1A]           # Invalid, only one subunit is specified while the complex has 3
            Comp[..., S1A, ...] # Invalid, cannot use more than one ellipsis operator


        :meta public:
        """
        if isinstance(subStates, ComplexState):
            subStates = tuple(subStates)
        if not isinstance(subStates, tuple):
            subStates = (subStates,)
        # If ellipsis is used, expand it first
        if len(subStates) > len(self._subUnits):
            raise KeyError(
                'The number of specified states is higher than the number of subunits ' 'in the complex.'
            )
        elif len(subStates) <= len(self._subUnits):
            if any(s is Ellipsis for s in subStates):
                if subStates.count(Ellipsis) > 1:
                    raise KeyError('Cannot use the ellipsis operator ("...") more than once.')
                # If the ellipsis operator was used correctly, fill the missing slots with ':'
                i = subStates.index(Ellipsis)
                m = len(self._subUnits) - len(subStates) + 1
                subStates = subStates[:i] + (slice(None),) * m + subStates[i + 1 :]
            elif len(subStates) < len(self._subUnits):
                raise KeyError(
                    'The number of states does not match the number of subunits. '
                    'A complex selector needs to either specify all subunits states '
                    'or use the Ellipsis operator ("...").'
                )
        # If the state is fully declared
        if all(isinstance(s, SubUnitState) for s in subStates):
            cs = ComplexState(self, subStates, _parentSelName=_selName)
            cs._stepsObj = self._compStates[cs] if self._statesAsSpecies else None
            return cs
        else:
            suSelectors = []
            for i, su in enumerate(subStates):
                if isinstance(su, SubUnitState):
                    if su._subUnit is not self._subUnits[i]:
                        raise KeyError(
                            f'SubUnitState {su} is not associated to SubUnit ' f'{self._subUnits[i]}.'
                        )
                    suSelectors.append(SubUnitSelector(self._subUnits[i], [su], _pos=i))
                elif isinstance(su, SubUnitSelector):
                    if su._isEmpty():
                        raise KeyError(
                            f'SubUnitSelector {su} is empty. Cannot build a complex '
                            f'selector from an empty SubUnitSelector.'
                        )
                    su2 = copy.copy(su)
                    su2._pos = i
                    suSelectors.append(su2)
                elif su == slice(None):
                    suSelectors.append(SubUnitSelector(self._subUnits[i], self._subUnits[i]._states, _pos=i))
                else:
                    raise KeyError(f'Unexpected entry in complex selector expression: {su}.')
            return ComplexSelector(self, suSelectors, name=_selName)

    def get(self):
        """Get a specific instance of the complex (i.e. a single molecule)

        :returns: An object that represents a specific instance of a complex, to be used in
            reactions. This object behaves like a :py:class:`Complex`.
        :rtype: :py:class:`_ComplexSelectorSpawner`
        """
        return _ComplexSelectorSpawner(ComplexSelector(self, []))

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return Complex._elemStr


@nutils.FreezeAfterInit
class SubUnit(nutils.NamedObject):
    """Class for representing complex subunits

    :param subStates: List of states in which the subunit can be.
    :type subStates: List[:py:class:`SubUnitState`]

    Does not necessarily represent a biological subunit, it can be used to group substates of
    a complex. For example, if a complex can bind to species, it could have a binding
    :py:class:`SubUnit` with :py:class:`SubUnitState` representing the species being bound.

    A :py:class:`SubUnit` can only be asscoiated to a single :py:class:`Complex`.
    """

    _elemStr = 'SU'

    def __init__(self, subStates, *args, **kwargs):
        """Initialize the SubUnit from a list of SubUnitState objects."""
        super().__init__(*args, **kwargs)
        self._orderedStates = subStates
        self._states = set(subStates)
        for state in self._states:
            state._setSubUnit(self)
        self._complex = None

    def __iter__(self):
        """Iterator over all the possible states of the subunit.

        :returns: An iterator that iterates over all :py:class:`SubUnitState` in the subunit
            state. The order is the same as the one defined during the creation of the subunit.
        :rtype: Iterator[:py:class:`SubUnitState`]

        Usage::

            >>> su = SubUnit.Create([h0, h1])
            >>> for subState in su:
            ...     print(subState.name)
            ...
            h0
            h1

        :meta public:
        """
        return iter(self._orderedStates)

    def _addComplex(self, comp):
        if self._complex is not None and self._complex != comp:
            raise Exception(f'A SubUnit can only be associated to one complex.')
        self._complex = comp

    def _areStatesAsSpecies(self):
        return self._complex._areStatesAsSpecies()

    def __hash__(self):
        return id(self)

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return SubUnit._elemStr


class _ComplexSelectorSpawner:
    """
    Utility class for spawning complex selectors with identical names.
    Having identical names means that they will represent the same instance of a complex in
    reactions.

    Usage:
        C = Comp.get() # Get the _ComplexSelectorSpawner

        C[S1A, ...] <r[1]> C[S2A, ...] # Signifies that the complexSelectors on both sides are
                                       # related to the same complex.
    """

    def __init__(self, cs):
        self._cs = cs

    def __getitem__(self, key):
        return self._cs._complex.__getitem__(key, _selName=self._cs.name)

    def __lshift__(self, other):
        return self._cs._complex.__getitem__(..., _selName=self._cs.name) << other


####################################
# Channels


@nutils.FreezeAfterInit
class Channel(Complex):
    """Ion channel class

    :parameters: See :py:class:`Complex` for a description of the parameters.

    Channels behave like complexes. This simplifies the declaration of reactions between channel
    states in cases like the Hodgkin-Huxley model of Na+ channel, in which gating variable m and h
    are independent and can thus be declared as separate SubUnits. All combinations are then
    automatically computed for creating the channel states.
    """

    def __init__(self, subUnits, *args, statesAsSpecies=True, **kwargs):
        if not statesAsSpecies:
            raise NotImplementedError()
        super().__init__(subUnits, *args, statesAsSpecies=True, **kwargs)

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        chanStates = obj.getAllChanStates()
        allSus = [SubUnitState() for state in chanStates]
        channel = cls(allSus, statesAsSpecies=True, _createObj=False, name=obj.getID())
        channel.stepsComplex = obj
        for state, sus in zip(chanStates, allSus):
            cs = ComplexState(channel, (sus,))
            cs._setStepsObjects(state)
            channel._compStates[cs] = state
        return channel

    def _createStepsStates(self, mdl):
        """Create steps ChanState objects for each channel state."""
        chan = stepslib._py_Chan(self.name, mdl.stepsModel)
        # Compute the channel state from the subunit states
        for comb in itertools.product(*[sorted(su._states, key=lambda x: x.name) for su in self._subUnits]):
            cs = ComplexState(self, comb)
            if cs not in self._compStates:
                chanState = stepslib._py_ChanState(cs.name, mdl.stepsModel, chan)
                cs._setStepsObjects(chanState)
                self._compStates[cs] = chanState
        return chan


###################################################################################################
# Reaction elements


@nutils.FreezeAfterInit
class ReactionElement(nutils.NamedObject):
    """Base class for objects that can be included in a reaction

    Can be combined to other :py:class:`ReactionElement` with the ``+`` operator to form a
    :py:class:`ReactionSide`.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _toReactionSide(self):
        return ReactionSide([ReactingElement(self, 1)])

    def __add__(self, other):
        """Add two elements and form a reaction side

        :param other: The other element
        :type other: :py:class:`ReactionElement`

        :returns: A reaction side containing itself and the other element it was added to.
        :rtype: :py:class:`ReactionSide`

        Usage::

            ...
            S1 + S2        <r[1]> ...
            S1 + 2*S2      <r[2]> ...
            S1 + (S2 + S3) <r[3]> ...
            ...

        :meta public:
        """
        if not isinstance(other, ReactionElement):
            raise TypeError(f'Cannot add {other} to a reaction.')
        reacSelf = self if isinstance(self, ReactingElement) else ReactingElement(self, 1)
        reacOther = other if isinstance(other, ReactingElement) else ReactingElement(other, 1)
        return ReactionSide([reacSelf, reacOther])

    def __rmul__(self, stoich):
        """Specify a stoichiometric coefficient

        :param stoich: The stoichimetric coefficient
        :type stoich: int

        :returns: A reacting element object that uses stoich as a stoichiometric coefficient.
        :rtype: :py:class:`ReactingElement`

        Usage::

            ...
            2*S1 <r[1]> ...
            # Equivalent to:
            S1 + S1 <r[1]> ...
            ...

        :meta public:
        """
        self._checkStoech(stoich)
        return ReactingElement(self, stoich)

    def _toReactionElem(self):
        return self

    def _checkStoech(self, stoich):
        if not isinstance(stoich, numbers.Integral) or stoich <= 0:
            raise TypeError(
                f'{stoich} is not an acceptable stoichiometric coefficient, it needs '
                f'to be a strictly positive integer.'
            )

    def _makeCopy(self):
        """
        Make a copy of the reaction element.
        For most elements, like species, self should be returned because there should not be more
        than one instance of this element.
        """
        return self

    # Location specifiers
    @property
    def i(self):
        """Get a version of the element located in the inside compartment

        Shorthand for :py:func:`In`.

        :type: :py:class:`ReactingElement`, read-only
        """
        return In(self)

    @property
    def o(self):
        """Get a version of the element located in the outside compartment

        Shorthand for :py:func:`Out`.

        :type: :py:class:`ReactingElement`, read-only
        """
        return Out(self)

    @property
    def s(self):
        """Get a version of the element located on the surface patch

        Shorthand for :py:func:`Surf`.

        :type: :py:class:`ReactingElement`, read-only
        """
        return Surf(self)

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        return {'obj_type': self._solverStr(), 'obj_id': self.name}

    @property
    def v(self):
        """Get a version of the element located on the surface of a vesicle

        Shorthand for :py:func:`VesSurf`.

        :type: :py:class:`ReactingElement`, read-only
        """
        return VesSurf(self)

    @property
    def r(self):
        """Get a version of the element located on the surface of a raft

        Shorthand for :py:func:`RaftSurf`.

        :type: :py:class:`ReactingElement`, read-only
        """
        return RaftSurf(self)


@nutils.FreezeAfterInit
class ReactingElement(ReactionElement):
    """A reaction element with a stoichiometric coefficient and a location

    All elements in a :py:class:`ReactionSide` are :py:class:`ReactingElement`. While
    :py:class:`ReactionElement` is the base class for all objects that could be included in a
    reaction, :py:class:`ReactingElement` is a wrapper that is actually included in reactions.
    It contains additional information about stoichiometry and location of the element.

    .. note::
        This class should not be instantiated by the user, it is only documented for clarity.
    """

    def __init__(self, elem, stoich, loc=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._elem = elem._toReactionElem()
        self.name = self._elem.name
        self.stoich = stoich
        self.loc = loc

    def _toReactionSide(self):
        return ReactionSide([self])

    def _SetLoc(self, loc):
        self.loc = loc

    def __rmul__(self, stoich):
        """Return a ReactingElement object that uses stoich as a stoichiometric coefficient."""
        self._checkStoech(stoich)
        return ReactingElement(self._elem, stoich * self.stoich, self.loc)

    def __repr__(self):
        if self.loc is None:
            return '{}{}'.format(self.stoich if self.stoich > 1 else '', str(self._elem))
        else:
            s = f'{self.stoich} ' if self.stoich > 1 else ''
            pos = {
                Location.IN: 'i', Location.OUT: 'o', Location.SURF: 's', Location.VESSURF: 'v',
                Location.RAFTSURF: 'r'
            }[self.loc]
            return f'{s}{self._elem}.{pos}'

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return self._elem._getStepsObjects()

    def __eq__(self, other):
        return isinstance(other, ReactingElement) and (self._elem, self.stoich, self.loc) == (
            other._elem,
            other.stoich,
            other.loc,
        )

    def __hash__(self):
        return hash((self._elem, self.stoich, self.loc))


####################################
# Location specifiers


class Location(Enum):
    """The different locations that reaction elements can be in"""

    IN = 1
    SURF = 2
    OUT = 3
    VESSURF = 4
    RAFTSURF = 5


ALL_LOCATIONS = [None, Location.IN, Location.SURF, Location.OUT, Location.VESSURF, Location.RAFTSURF]
ALL_SURF_LOCATIONS = [Location.SURF, Location.VESSURF, Location.RAFTSURF]


def _locSpecifier(elem, loc):
    """Return a localized ReactingElement with location 'loc'."""
    if isinstance(elem, ReactionElement) and not isinstance(elem, ReactingElement):
        return ReactingElement(elem, 1, loc)
    elif isinstance(elem, ReactionSide):
        elem._SetLoc(loc)
        return elem
    else:
        raise TypeError(f'Expected a ReactionElement or a ReactionSide, got {elem} instead.')


def In(elem):
    """Get a version of an element located in the inside compartment

    :param elem: The element
    :type elem: Union[:py:class:`ReactionElement`, :py:class:`ReactionSide`]

    :returns: A version of the element located in the inside compartment
    :rtype: Union[:py:class:`ReactingElement`, :py:class:`ReactionSide`]
    """
    return _locSpecifier(elem, Location.IN)


def Out(elem):
    """Get a version of an element located in the outside compartment

    :param elem: The element
    :type elem: Union[:py:class:`ReactionElement`, :py:class:`ReactionSide`]

    :returns: A version of the element located in the outside compartment
    :rtype: Union[:py:class:`ReactingElement`, :py:class:`ReactionSide`]
    """
    return _locSpecifier(elem, Location.OUT)


def Surf(elem):
    """Get a version of an element located on the surface patch

    :param elem: The element
    :type elem: Union[:py:class:`ReactionElement`, :py:class:`ReactionSide`]

    :returns: A version of the element located on the surface patch
    :rtype: Union[:py:class:`ReactingElement`, :py:class:`ReactionSide`]
    """
    return _locSpecifier(elem, Location.SURF)

def VesSurf(elem):
    """Get a version of an element located on the surface of a vesicle

    :param elem: The element
    :type elem: Union[:py:class:`ReactionElement`, :py:class:`ReactionSide`]

    :returns: A version of the element located on the surface of a vesicle
    :rtype: Union[:py:class:`ReactingElement`, :py:class:`ReactionSide`]
    """
    return _locSpecifier(elem, Location.VESSURF)

def RaftSurf(elem):
    """Get a version of an element located on the surface of a raft

    :param elem: The element
    :type elem: Union[:py:class:`ReactionElement`, :py:class:`ReactionSide`]

    :returns: A version of the element located on the surface of a raft
    :rtype: Union[:py:class:`ReactingElement`, :py:class:`ReactionSide`]
    """
    return _locSpecifier(elem, Location.RAFTSURF)


####################################
# Vesicles

class _RaftVesLocation(nutils.UsingObjects(Model),
        nutils.StepsWrapperObject,
        nutils.ParameterizedObject
    ):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.systems = []
        self.sysNames = []

    def _SetUpMdlDeps(self, mdl, geom):
        """Set up the structures that will allow species or reaction access from locations."""
        for sname, loc in self.sysNames:
            s = getattr(mdl, sname)
            self.systems.append((s, loc))
            # Add the reactions and the corresponding reactants to children
            for name, c in s.children.items():
                if _RaftVesLocation._canBeChild(self, c):
                    if self._canBeChild(c):
                        self._addChildren(c)
                    for re in c._getAllElems(loc):
                        if re.name not in self.children:
                            self._addChildren(re)

            s._addLocation(self)
        # Add all species as children since users can set any species on vesicles and rafts
        for spec in mdl.ALL(Species):
            if spec.name not in self.children:
                self._addChildren(spec)

    def _canBeChild(self, c):
        """Return whether c can be a child of self."""
        return isinstance(c, (Reaction, Diffusion, Exocytosis, RaftEndocytosis, RaftDis))

    def addSystem(self, sys, _loc=None):
        """Add a vesicle or raft surface system to the location

        :param sys: The volume or surface system to be added, or its name.
        :type sys: :py:class:`VesicleSurfaceSystem`, :py:class:`RaftSurfaceSystem` or `str`

        :returns: None
        """
        if isinstance(sys, SpaceSystem):
            self.sysNames.append((sys.name, _loc))
        elif isinstance(sys, str):
            self.sysNames.append((sys, _loc))
        else:
            raise TypeError(
                f'Expected a VesicleSurfaceSystem or a RaftSurfaceSystem, got {sys} instead.'
            )


class _VesicleSelection(nutils.SolverPathObject):
    """Class that represents a subselection of vesicles with additional information"""

    def __init__(self, ves, inside):
        self._ves = ves
        self._inside = inside

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return self._ves._solverStr() + ('Inner' if self._inside else 'Surface')

    def _solverId(self):
        """Return the id that is used to identify an object in the solver."""
        return self._ves._solverId()

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        return {'ves_loc': self.loc}

    def __getattr__(self, name):
        return getattr(self._ves, name)

    @property
    def children(self):
        return self._ves.children

    @property
    def ves(self):
        return self._ves

    @property
    def loc(self):
        return 'in' if self._inside else 'surf'

    def __repr__(self):
        return f"{self.ves}('{self.loc}')"


@nutils.FreezeAfterInit
class Vesicle(_RaftVesLocation):
    """Represents a vesicle object

    Vesicles are described by a spherical diameter and diffusion rate.
    It models the many behaviours and interactions of these complex biological
    entities. A vesicle may be formed by endocytosis or by direct user input,
    may exit by exocytosis or user input, may transport species on its surface or
    luminally, and may undergo many interactions with its environment. These
    interactions are individually described in each corresponding class.

    :param diameter: Diameter of the vesicle (in m)
    :type diameter: float
    :param dcst: Default diffusion coefficient of the vesicle (in m^2 s^-1)
    :param vesSurfSys: Optional vesicle surface system to be associated with the vesicle
    :type vesSurfSys: :py:class:`VesicleSurfaceSystem`
    """

    _locStr = 'Vesicle'

    def __init__(self, diameter, dcst=0, vesSurfSys=None, _createObj=True, **kwargs):
        super().__init__(**kwargs)
        (mdl,) = self._getUsedObjects()

        self._setParameter('Diameter', diameter, nutils.Units('m'))
        self._setParameter('Dcst', dcst, nutils.Units('m^2 s^-1'))

        if _createObj:
            self.stepsVesicle = stepslib._py_Vesicle(self.name, mdl.stepsModel, self.Diameter, self.Dcst)
        else:
            self.stepsVesicle = None

        if vesSurfSys is not None:
            self.addSystem(vesSurfSys)

    def _SetUpMdlDeps(self, mdl, geom):
        """Set up the structures that will allow species or reaction access from locations."""
        super()._SetUpMdlDeps(mdl, geom)
        # Add all relevant link species as children
        for vsys in mdl.ALL(VolumeSystem):
            for vbind in vsys.ALL(VesicleBind, VesicleUnbind):
                for ves, ls in vbind._getVesLinkSpecPairs():
                    if ves == ves and ls.name not in self.children:
                        self._addChildren(ls)

    def _canBeChild(self, c):
        """Return whether c can be a child of self."""
        return super()._canBeChild(c)

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsVesicle]

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return Vesicle._locStr

    def addSystem(self, sys, _loc=None):
        """Add a vesicle surface system to the vesicle

        :param sys: The vesicle surface system to be added, or its name.
        :type sys: Union[:py:class:`steps.API_2.model.VesicleSurfaceSystem`, str]

        :returns: None
        """
        super().addSystem(sys, Location.VESSURF)
        if _loc is None:
            if isinstance(sys, SpaceSystem):
                sys = sys.name
            self.stepsVesicle.addVesSurfsys(sys)

    def __call__(self, loc):
        """Get a version of the vesicle with added information

        Parentheses notation (function call) can be used on a vesicle to specify additional
        information. This should never be needed during model declaration but is frequently needed
        for simulation control and data saving (see :py:class:`steps.API_2.sim.SimPath`).

        :param loc: If `'surf'`, the simulation path represents the surface of the selected vesicles,
            if `'in'`, it represents the inside of the selected vesicles
        :type inside: str

        :returns: An object that represents the vesicle with the added information

        :meta public:
        """
        if loc not in {'surf', 'in'}:
            raise ValueError(f"Location string can only be 'surf' or 'in', got {loc} instead.")
        return _VesicleSelection(self, loc == 'in')

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        ves = cls(obj.getDiameter(), obj.getDcst(), _createObj=False, name=obj.getID())
        ves.stepsVesicle = obj
        for sysname in obj.getVesSurfsys():
            ves.addSystem(sysname)
        return ves

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m'))
    def Diameter(self):
        """Vesicle diameter (in m)

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`], read-only
        """
        pass

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m^2 s^-1'))
    def Dcst(self):
        """Diffusion constant (in m^2 s^-1)

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        pass

    @Dcst.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('m^2 s^-1'))
    def Dcst(self, val):
        if isinstance(val, numbers.Number):
            self.stepsVesicle.setDcst(val)
        else:
            raise TypeError(f'{val} is not a valid diffusion constant.')


@nutils.FreezeAfterInit
class Raft(_RaftVesLocation):
    """Represents a membrane raft

    The membrane analogy of vesicles, rafts effectively group surface species
    within a defined radius and exist in and may diffuse within patches. Rafts
    and the species they contain undergo special interactions such as raft
    endocytosis, raft generation and raft dissolution.

    Rafts can represent lipid rafts or other clusters. They occupy space
    defined by a diameter projected onto the surface, and can be mobile.

    :param diameter: Diameter of the raft (in m)
    :type diameter: float
    :param dcst: Default diffusion coefficient of the raft (in m^2 s^-1)
    :param raftSurfSys: Optional raft surface system to be associated with the raft
    :type raftSurfSys: :py:class:`RaftSurfaceSystem`
    """

    _locStr = 'Raft'

    def __init__(self, diameter, dcst=0, raftSurfSys=None, _createObj=True, **kwargs):
        super().__init__(**kwargs)
        (mdl,) = self._getUsedObjects()

        self._setParameter('Diameter', diameter, nutils.Units('m'))
        self._setParameter('Dcst', dcst, nutils.Units('m^2 s^-1'))

        if _createObj:
            self.stepsRaft = stepslib._py_Raft(self.name, mdl.stepsModel, self.Diameter, self.Dcst)
        else:
            self.stepsRaft = None

        if raftSurfSys is not None:
            self.addSystem(raftSurfSys)

    def _canBeChild(self, c):
        """Return whether c can be a child of self."""
        if not super()._canBeChild(c):
            return False
        if ((isinstance(c, Reaction) and not c._isRaftSurfReac()) or
            isinstance(c, (Diffusion, Exocytosis))
        ):
            return False
        return True

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsRaft]

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return Raft._locStr

    def addSystem(self, sys, _loc=None):
        """Add a raft surface system to the raft

        :param sys: The raft surface system to be added, or its name.
        :type sys: Union[:py:class:`steps.API_2.model.RaftSurfaceSystem`, str]

        :returns: None
        """
        super().addSystem(sys, Location.RAFTSURF)
        if _loc is None:
            if isinstance(sys, SpaceSystem):
                sys = sys.name
            self.stepsRaft.addRaftsys(sys)

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        raft = cls(obj.getDiameter(), obj.getDcst(), _createObj=False, name=obj.getID())
        raft.stepsRaft = obj
        for sysname in obj.getRaftsys():
            raft.addSystem(sysname)
        return raft

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m'))
    def Diameter(self):
        """Raft diameter (in m)

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`], read-only
        """
        pass

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m^2 s^-1'))
    def Dcst(self):
        """Raft diffusion constant (in m^2 s^-1)

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        pass

    @Dcst.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('m^2 s^-1'))
    def Dcst(self, val):
        if isinstance(val, numbers.Number):
            self.stepsRaft.setDcst(val)
        else:
            raise TypeError(f'{val} is not a valid diffusion constant.')


####################################


@nutils.FreezeAfterInit
class Species(
        nutils.UsingObjects(Model),
        ReactionElement,
        nutils.StepsWrapperObject,
        nutils.ParameterizedObject
    ):
    """A chemical species

    Species can be involved in :py:class:`Reaction`, :py:class:`Diffusion`, or other transport
    mechanisms.

    :param valence: Optional, the valence of the species (defaults to 0)
    :type valence: int
    """

    _elemStr = 'Spec'

    def __init__(self, *args, valence=0, _createObj=True, **kwargs):
        super().__init__(*args, **kwargs)
        (mdl,) = self._getUsedObjects()
        self.stepsSpecies = stepslib._py_Spec(self.name, mdl.stepsModel) if _createObj else None
        if valence != 0:
            self.valence = valence

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsSpecies]

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        spec = cls(_createObj=False, name=obj.getID())
        spec.stepsSpecies = obj
        val = obj.getValence()
        if val != 0:
            spec.valence = val
        return spec

    @property
    @nutils.ParameterizedObject.RegisterGetter()
    def valence(self):
        """The valence of the species (defaults to 0)

        :type: Union[int, :py:class:`steps.API_2.utils.Parameter`]
        """
        return self.stepsSpecies.getValence()

    @valence.setter
    @nutils.ParameterizedObject.RegisterSetter()
    def valence(self, val):
        self.stepsSpecies.setValence(val)

    def __repr__(self):
        return self.name

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return Species._elemStr


@nutils.FreezeAfterInit
class LinkSpecies(
        nutils.UsingObjects(Model),
        ReactionElement,
        nutils.StepsWrapperObject,
        nutils.ParameterizedObject
    ):
    """A chemical species that can form links between vesicles

    A link species is a special kind of chemical species that links two
    vesicles (by forming a bond with another link species). A link species is
    formed by a vesicle binding event and exists within a specified upper and
    lower bound of length for the duration of its existence. A link species
    may diffuse on a vesicle surface, but only within this length bound. Link
    species are destroyed by a vesicle unbinding event.

    :param dcst: Optional, diffusion coefficient of the vesicle (in m^2 s^-1, defaults to 0)
    :type dcst: float
    """

    _elemStr = 'LinkSpec'

    def __init__(self, dcst=0, _createObj=True, **kwargs):
        super().__init__(**kwargs)
        (mdl,) = self._getUsedObjects()

        self._setParameter('Dcst', dcst, nutils.Units('m^2 s^-1'))

        self.stepsSpecies = stepslib._py_LinkSpec(self.name, mdl.stepsModel, self.Dcst) if _createObj else None

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsSpecies]

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        spec = cls(obj.getDcst(), _createObj=False, name=obj.getID())
        spec.stepsSpecies = obj
        return spec

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m^2 s^-1'))
    def Dcst(self):
        """Diffusion constant (in m^2 s^-1)

        :type: float, read-only
        """
        pass

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return LinkSpecies._elemStr


@nutils.FreezeAfterInit
class SubUnitState(ReactionElement):
    """State of a :py:class:`SubUnit`

    Can take part in reactions (see :py:class:`Complex`) and be combined with ``|`` and ``~`` to
    form a :py:class:`SubUnitSelector`.

    In reactions, further information can be provided with the square bracket operator, see
    :py:func:`SubUnitState.__getitem__`.

    A :py:class:`SubUnitState` can only be associated to a single :py:class:`SubUnit`.

    .. warning::
        When used in reactions in conjunction with the `|` operator, the expression must be wrapped
        in parentheses to avoid incorrect operator precedence.
    """

    _elemStr = 'SUS'

    def __init__(self, _orig=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._subUnit = None
        # Each subunit state creates its own modified version, for use in reaction computation
        if _orig is None:
            self._modVers = SubUnitState(name=f'special%mod', _orig=self)
            self._orig = None
        else:
            self._modVers = None
            self._orig = _orig

    def _setSubUnit(self, su):
        if self._subUnit is None:
            self._subUnit = su
            if self._modVers is not None:
                self._modVers._setSubUnit(su)
        else:
            raise Exception('A SubUnitState can only be associated to one subunit.')

    def _getModif(self):
        return self._modVers if self._modVers is not None else self

    def _getOrig(self):
        return self if self._orig is None else self._orig

    def _toReactionElem(self):
        return self._toSelector()

    def _areStatesAsSpecies(self):
        return self._subUnit._areStatesAsSpecies()

    def __or__(self, other):
        """Get the union between two :py:class:`SubUnitState`

        :param other: The other state(s)
        :type other: Union[:py:class:`SubUnitState`, :py:class:`SubUnitSelector`]

        :returns: An object representing the union between the current state and the state(s)
            from 'other'.
        :rtype: :py:class:`SubUnitSelector`

        Usage::

            (S1A | S1B)
            # It can also be chained
            (S1A | S1B | S1C)

        If the associated :py:class:`SubUnit` was defined with ``S1A``, ``S1B``, and ``S1C`` as
        states, ``(S1A | S1B | S1C)`` is equivalent to ``:`` (all SubUnitStates).

        .. warning::
            Do not forget to wrap the overall expression in parentheses in order to avoid issues
            with operator precedence when used in reaction declaration.

        :meta public:
        """
        return self._toSelector() | other

    def __invert__(self):
        """Get all :py:class:`SubUnitState` except itself

        :returns: An object representing all states of the associated :py:class:`SubUnit`
            except itself
        :rtype: :py:class:`SubUnitSelector`

        Usage::

            ~S1A # All subunits except S1A

        If the associated :py:class:`SubUnit` was defined with ``S1A``, ``S1B``, and ``S1C`` as
        states, ``~S1A`` is equivalent to ``(S1B | S1C)``

        :meta public:
        """
        return SubUnitSelector(self._subUnit, [s for s in self._subUnit._states if s is not self])

    def _toSelector(self):
        return SubUnitSelector(self._subUnit, [self])

    def __getitem__(self, key):
        """Specify additional information

        Transform the subunit state into a :py:class:`SubUnitSelector` that only selects the
        current state, and call ``__getitem__`` on it.
        See :py:func:`SubUnitSelector.__getitem__` for details.

        :meta public:
        """
        return SubUnitSelector(self._subUnit, [self]).__getitem__(key)

    def __iter__(self):
        raise NotImplementedError()

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        return isinstance(other, SubUnitState) and self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        if self._areStatesAsSpecies():
            return ''
        else:
            return SubUnitState._elemStr

    def _solverId(self):
        """Return the id that is used to identify an object in the solver."""
        if self._areStatesAsSpecies():
            return tuple()
        else:
            return (self._subUnit._complex._getSUSInd(self),)

    def _solverModifier(self):
        """Return None or a function that will modify the output from the solver."""
        if self._areStatesAsSpecies():
            return lambda x: (x, self)
        else:
            return None


@nutils.FreezeAfterInit
class SubUnitSelector(ReactionElement):
    """Set of :py:class:`SubUnitState` associated to a :py:class:`SubUnit`

    Represents a selection among the states of a :py:class:`SubUnit`.

    .. note::
        Should not be instantiated by a user, the class is only documented for clarity.
    """

    _OR = 0
    _AND = 1

    def __init__(self, su, states, *args, _pos=None, _compSel=None, _ident=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._subUnit = su
        if any(sus._subUnit != su for sus in states):
            raise Exception(f'Some subunit states in {states} are not associated with subunit {su}.')
        self._states = frozenset(states)
        self._pos = _pos
        self._compSel = _compSel
        self._id = _ident

    def _setOp(self, op, other):
        if isinstance(other, SubUnitState):
            other = other._toSelector()
        if isinstance(other, SubUnitSelector):
            if other._subUnit is not self._subUnit:
                raise TypeError('Cannot combine subunit states from different subunits.')
            if other._compSel != self._compSel:
                raise Exception('Cannot combine subunit states from different complex selectors.')
            if other._pos != self._pos:
                raise Exception('Cannot combine subunit states from different positions.')
            tmpId = self._id if self._id is not None else other._id

            if op == SubUnitSelector._OR:
                return SubUnitSelector(
                    self._subUnit,
                    self._states | other._states,
                    _pos=self._pos,
                    _compSel=self._compSel,
                    _ident=tmpId,
                )
            elif op == SubUnitSelector._AND:
                return SubUnitSelector(
                    self._subUnit,
                    self._states & other._states,
                    _pos=self._pos,
                    _compSel=self._compSel,
                    _ident=tmpId,
                )
            else:
                raise NotImplementedError()
        elif isinstance(other, (ReactionSide, ReactionElement)):
            # Common syntax error
            raise SyntaxError(
                'SubUnitSelectors with "|" or "&" operator need to be wrapped by '
                'parentheses when used in reactions.'
            )
        else:
            raise TypeError(f'Expected a SubUnitState or a SubUnitSelector, got {other} instead.')

    def __getitem__(self, key):
        """Specify additional information

        :param key: The information(s) (it is possible to provide several within the same square
            brackets)
        :type key: Union[int, str, :py:class:`ComplexSelector`]

        :returns: The subunit selector, annotated with the information
        :rtype: :py:class:`SubUnitSelector`

        This method allows the usage of square brackets notation for providing the following
        additional informations:

        - Position of the associated :py:class:`SubUnit` in the complex::

            with Comp[...]:
                # Ca binds to the first subunit (index 0) when in state S1A
                S1A[0] + Ca <r[1]> S1B[0]

        - Association with a specific complex within a reaction::

            with Comp[...] as C1, Comp[...] as C2:
                # Two instances of Comp react, S1A in the first one is modified into S1B
                # while S1B in the second is modified to S1C.
                S1A[C1] + S1B[C2] <r[1]> S1B[C1] + S1C[C2]

        - Custom identifier (useful when a reaction is ambiguous)::

            with Comp[...]:
                # A subunit in state S1A reacts with another subunit is state S1B. Since
                # both states are associated to the same subunit (and this subunit is
                # present twice in the complex), the following reaction could be read in
                # two different ways:
                #     S1A + S1B >r[1]> S1B + S1C
                # We could have (S1A, S1B) -> (S1B, S1C) or (S1A, S1B) -> (S1C, S1B). In
                # the first case both subunit changed while in the second, only the first
                # one changed. Custom identifiers allow us to specify this reaction
                # without ambiguities:
                S1A['a'] + S1B['b'] >r[1]> S1B['a'] + S1C['b']

        Note that it is possible to provide several information at the same time. For example::

            with Comp[...] as C1, Comp[...] as C2:
                # Two instances of Comp react, in the first one (S1A, S1B) -> (S1B, S1C) while
                # in the second one S1B -> S1A
                S1A[C1, 'a'] + S1B[C1, 'b'] + S1B[C2] >r[1]> S1B[C1, 'a'] + S1C[C1, 'b'] + S1A[C2]


        :meta public:
        """
        if not isinstance(key, tuple):
            key = (key,)
        if len(key) == 0 or len(key) > 3:
            raise KeyError(
                'SubUnitStates can only be indexed with one complex selector, one '
                'integer indicating their position, or one identifier string.'
            )
        res = copy.copy(self)
        for v in key:
            if isinstance(v, _ComplexSelectorSpawner):
                v = v[...]
            if isinstance(v, ComplexSelector):
                if res._compSel is not None:
                    raise Exception('The SubUnitSelector is already associated with a complex selector.')
                res._compSel = v
            elif isinstance(v, numbers.Integral):
                if res._pos is not None:
                    raise Exception('The SubUnitSelector is already associated with a position.')
                res._pos = v
            elif isinstance(v, str):
                if res._id is not None:
                    raise Exception('The SubUnitSelector is already associated with an identifier.')
                res._id = v
            else:
                raise TypeError(
                    f'Only complex selectors, strings, and integers can be used as a '
                    f'key, got {v} instead.'
                )
        return res

    def _isEmpty(self):
        return len(self._states) == 0

    def _isFull(self):
        """Return True if the SubUnitSelector covers all states of its corresponding subunit."""
        return len(self._states) == len(self._subUnit._states)

    def _checkPos(self, comp):
        """Raise an exception is the position is not compatible with complex 'comp'."""
        n = len(comp._subUnits)
        if self._pos is not None and not (0 <= self._pos < n):
            raise Exception(f'{self._pos} is not a valid position, it needs to be >= 0 and < {n}')

    def __or__(self, other):
        """Get the union between the states of two :py:class:`SubUnitSelector`

        :param other: The other subunit selector, has to be associated to the same
            :py:class:`SubUnit`. Can also be a single :py:class:`SubUnitState`.
        :type other: Union[:py:class:`SubUnitState`, :py:class:`SubUnitSelector`]

        :returns: A subunit selector representing the union between the current selector states and
            the state(s) from 'other'.
        :rtype: :py:class:`SubUnitSelector`

        See :py:class:`SubUnitState.__or__` for usage examples.

        .. warning::
            Do not forget to wrap the overall expression in parentheses in order to avoid issues
            with operator precedence when used in reaction declaration.

        :meta public:
        """
        return self._setOp(SubUnitSelector._OR, other)

    def __and__(self, other):
        """Get the intersection between the states of two :py:class:`SubUnitSelector`

        :param other: The other subunit selector, has to be associated to the same
            :py:class:`SubUnit`. Can also be a single :py:class:`SubUnitState`.
        :type other: Union[:py:class:`SubUnitState`, :py:class:`SubUnitSelector`]

        :returns: A subunit selector representing the intersection between the current selector
            states and the state(s) from 'other'.
        :rtype: :py:class:`SubUnitSelector`

        Usage::

            AorB = S1A | S1B
            BorC = S1B | S1C
            (AorB & BorC)            # Equivalent to S1B
            (AorB & (~S1A) & (~S1C)) # Equivalent to S1B

        If the intersection is empty, it returns an empty :py:class:`SubUnitSelector` which
        might raise exceptions if later used in a reaction declaration.

        .. warning::
            Do not forget to wrap the overall expression in parentheses in order to avoid issues
            with operator precedence when used in reaction declaration.

        :meta public:
        """
        return self._setOp(SubUnitSelector._AND, other)

    def __invert__(self):
        """Get all :py:class:`SubUnitState` except the ones in self

        :returns: An object representing all states of the associated :py:class:`SubUnit`
            minus the ones present in self
        :rtype: :py:class:`SubUnitSelector`

        Usage::

            ~(S1A | S1C) # All subunits except S1A and S1C

        If the associated :py:class:`SubUnit` was defined with ``S1A``, ``S1B``, and ``S1C`` as
        states, ``~(S1A | S1C)`` is equivalent to ``S1B``

        :meta public:
        """
        return SubUnitSelector(self._subUnit, self._subUnit._states - self._states)

    def __iter__(self):
        return iter(sorted(self._states, key=lambda x: x.name))

    def __contains__(self, sus):
        """Test whether a subunit state is present in the subunit selector

        :param sus: The subunit state
        :type sus: :py:class:`SubUnitState`

        :returns: True if the given subunit state is in the current selector, False otherwise.
        :rtype: bool

        Usage::

            >>> AorB = S1A | S1B
            >>> S1A in AorB
            True
            >>> S1C in AorB
            False

        :meta public:
        """
        return isinstance(sus, SubUnitState) and sus._subUnit is self._subUnit and sus in self._states

    def __repr__(self):
        if self._states == self._subUnit._states:
            return ':'
        elif len(self._states) <= len(self._subUnit._states) / 2:
            return '|'.join(map(str, self._states))
        else:
            # Express it as a negation if it is more concise
            if len(self._states) == len(self._subUnit._states) - 1:
                return f'~{~self}'
            else:
                return f'~({~self})'

    def __eq__(self, other):
        return isinstance(other, SubUnitSelector) and (
            self._states,
            self._compSel._complex if self._compSel is not None else None,
        ) == (other._states, other._compSel._complex if other._compSel is not None else None)

    def __hash__(self):
        return hash((self._states, self._compSel._complex if self._compSel is not None else None))

    def _makeCopy(self):
        """Make a copy of the reaction element."""
        return copy.copy(self)


@nutils.FreezeAfterInit
class ComplexState(ReactionElement):
    """Fully specified state of a complex

    Can be used in the same way as a :py:class:`ReactionElement`.

    .. note::
        A :py:class:`ComplexState` object should not be directly instantiated by the user, instead
        it is obtained by using the square bracket notation on a :py:class:`Complex` object::

            Comp[S1B, S1B, S2A] # Returns the complex state for which the two S1 subunits are
                                # in state S1B and the S2 subunit is in state S2A.

        See :py:func:`Complex.__getitem__` for details.
    """

    _SIMPATH_ONLY_CHILDREN = True

    def __init__(
        self, comp, state, *args, _stepsObj=None, _parentCompSel=None, _parentSelName=None, **kwargs
    ):
        self._comp = comp
        # check that the subunit states are compatible with the subunits
        for su, st in zip(self._comp._subUnits, state):
            if st not in su._states:
                raise TypeError(f'{st} is not a state of subunit {su}.')

        self._state = self._reOrderState(self._comp._order(state))
        self._stepsObj = _stepsObj

        self._parentCompSel = _parentCompSel
        self._parentSelName = _parentSelName

        super().__init__(*args, name=self._getStateName(), **kwargs)

    def _getReferenceObject(self):
        """
        Return the object this object was derived from. Useful for getting the complex associated
        with a complex selector, etc.
        """
        return self._comp

    def _getComplexName(self):
        return self._comp.name

    def _getStateName(self):
        """Return the name of the state."""
        return self._comp.name + '_' + '_'.join(map(str, self._state))

    def _setStepsObjects(self, obj):
        self._stepsObj = obj

    def _areStatesAsSpecies(self):
        return self._comp._areStatesAsSpecies()

    def _getStepsObjects(self):
        if self._stepsObj is None and self._comp._areStatesAsSpecies():
            return [self._comp._compStates[self]]
        else:
            return [self._stepsObj]

    def _reOrderState(self, state):
        """
        Transforms the output of the ordering function in such a way that subunit states are 
        associated to the correct subunits, i.e. reorder subunit states to match subunit order.
        Conserve relative ordering of subunit states linked to identical subunits.
        """
        subunits = copy.copy(self._comp._subUnits)
        res = [None] * len(subunits)
        for sus in state:
            i = subunits.index(sus._subUnit)
            res[i] = sus
            subunits[i] = None

        if any(sus is None for sus in res) or any(su is not None for su in subunits):
            raise Exception(
                f'An invalid complex state was used. There is probably an issue with the ordering function.'
            )

        return tuple(res)

    def _toComplexSelector(self):
        """
        Return a complex selector describing this complex state, keep the link to the original
        complexSelector that was used to instantiante the state, if it exists.
        """
        stateVals = zip(range(len(self._state)), self._comp._subUnits, self._state)
        subSels = [SubUnitSelector(su, [s], _pos=i) for i, su, s in stateVals]
        if self._parentCompSel is not None:
            res = copy.copy(self._parentCompSel)
            res._initSubSels(subSels)
            return res
        elif self._parentSelName is not None:
            return ComplexSelector(self._comp, subSels, name=self._parentSelName)
        else:
            return ComplexSelector(self._comp, subSels)

    def _toUnorderedVect(self):
        """
        Return a subunit state pool vector representing the complex state
        """
        vec = [0] * len(self._comp._poolElems)
        for sus in self:
            vec[self._comp._sus2PoolInd[sus]] += 1
        return tuple(vec)

    def _toUnorderedFilter(self):
        """
        Return a filter representing the complex state
        """
        vec = [(0, 0) for i in range(len(self._comp._poolElems))]
        for sus in self:
            _min, _max = vec[self._comp._sus2PoolInd[sus]]
            vec[self._comp._sus2PoolInd[sus]] = (_min + 1, _max + 1)
        return (tuple(vec),)

    def Count(self, sus):
        """Count the number of subunits that match a given subunit selector

        :param sus: A subunit state or a subunit selector
        :type sus: Union[:py:class:`SubUnitSelector`, :py:class:`SubUnitState`]

        :returns: The number of subunits in the complex state that match ``sus``.
        :rtype: int

        Usage::

            >>> state = Comp[S1A, S1B, S2A]
            >>> state.Count(S1B)
            1
            >>> state.Count(S1C)
            0
            >>> state.Count(S1A | S1B)
            2
        """
        if isinstance(sus, SubUnitState):
            sus = sus._toSelector()
        if not isinstance(sus, SubUnitSelector):
            raise TypeError(f'Expected a SubUnitState or a SubUnitSelector, got {sus} instead.')
        return sum(1 for s in self._state if s in sus)

    def __iter__(self):
        """Iterator over all the subunit states of the complex state.

        :returns: An iterator that iterates over all :py:class:`SubUnitState` in the complex
            state. The order is the same as the one defined during the creation of the associated
            :py:class:`Complex`.
        :rtype: Iterator[:py:class:`SubUnitState`]

        Usage::

            >>> state = Comp[S1A, S1B, S2A]
            >>> for subState in state:
            ...     print(subState.name)
            ...
            S1A
            S1B
            S2A

        :meta public:
        """
        return iter(self._state)

    def __eq__(self, other):
        return isinstance(other, ComplexState) and other._comp == self._comp and other._state == self._state

    def __hash__(self):
        return hash((self._comp, self._state))

    def __or__(self, other):
        """Get a complex selector representing the union between the states of self and other

        :param other: The other complex state(s), has to be associated to the same
            :py:class:`Complex`.
        :type other: Union[:py:class:`ComplexState`, :py:class:`ComplexSelector`]

        :returns: An object representing the union between the current state and the complex
            state(s) from 'other'.
        :rtype: :py:class:`ComplexSelector`

        Usage::

            Comp[S1A, S1A, S2A] | Comp[S1B, S1B, S2A] # Represents only these two states

            # Note that Comp[S1A | S1B, S1A | S1B, S2A] is equivalent to more than two states
            Comp[S1A, S1A, S2A] | Comp[S1A, S1B, S2A] | Comp[S1B, S1B, S2A]

            # If the custom order of the complex (See :py:class:`Complex`) is set to
            # StrongOrdering, it is equivalent to 4 states:
            Comp[S1A, S1A, S2A] | Comp[S1A, S1B, S2A] | Comp[S1B, S1A, S2A] | Comp[S1B, S1B, S2A]

        .. warning::
            Do not forget to wrap the overall expression in parentheses in order to avoid issues
            with operator precedence when used in reaction declaration.

        :meta public:
        """
        if isinstance(other, ComplexState):
            return self._toComplexSelector().__or__(other._toComplexSelector())
        elif isinstance(other, ComplexSelector):
            return self._toComplexSelector().__or__(other)
        else:
            raise TypeError('ComplexStates can only be combined with ComplexSelectors or ' 'ComplexStates.')

    def __repr__(self):
        return f'{self._comp.name}[' + ', '.join(map(str, self)) + ']'

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        if self._areStatesAsSpecies():
            return Species._elemStr
        else:
            return Complex._elemStr

    def _solverId(self):
        """Return the id that is used to identify an object in the solver."""
        if self._areStatesAsSpecies():
            return (self.name,)
        else:
            return (
                self._comp.name,
                self._toUnorderedFilter(),
            )

    def _solverModifier(self):
        """Return None or a function that will modify the output from the solver."""
        if self._areStatesAsSpecies():

            def modif(x):
                if isinstance(x, tuple):
                    return x[0] * self.Count(x[1]._toSelector())
                else:
                    return x

            return modif
        else:
            return None

    def _simPathCheckParent(self):
        """
        Determines whether the object needs to be in the parent's children in order to be added to
        a simulation path
        """
        return self._areStatesAsSpecies()


class _ComplexReactionElement(ReactionElement):
    def _areStatesAsSpecies(self):
        pass


@nutils.FreezeAfterInit
class ComplexSelector(_ComplexReactionElement, nutils.MultiUsable):
    """Partially specified state of a complex

    Can be used in the same way as a :py:class:`ReactionElement`.

    .. note::
        A :py:class:`ComplexSelector` object should not be directly instantiated by the user,
        instead it is obtained by using the square bracket notation on a :py:class:`Complex`
        object in a way that does not fully specify the whole state::

            Comp[S1B, S1B, :] # Returns the complex state for which the two S1 subunits are
                              # in state S1B and the S2 subunit is in any state.

        See :py:func:`Complex.__getitem__` for details.
    """

    _OR = 0
    _AND = 1

    def __init__(self, comp, subSels, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._complex = comp
        self._initSubSels(subSels)

    def _getReferenceObject(self):
        """
        Return the object this object was derived from. Useful for getting the complex associated
        with a complex selector, etc.
        """
        return self._complex

    def _areStatesAsSpecies(self):
        return self._complex._areStatesAsSpecies()

    def _getComplexName(self):
        return self._complex.name

    def _initSubSels(self, subSels):
        self._subSels = set(
            [tuple(SubUnitSelector(su, su._states, _pos=i) for i, su in enumerate(self._complex._subUnits))]
        )
        for ss in subSels:
            self._addSubSel(ss)

    def _addSubSel(self, ss):
        newSubSels = set()
        for row in self._subSels:
            newRow = row[: ss._pos] + (row[ss._pos] & ss,) + row[ss._pos + 1 :]
            if not newRow[ss._pos]._isEmpty():
                newSubSels.add(newRow)
        self._subSels = newSubSels

    def _isEmpty(self):
        """Return True if no state can match the ComplexSelector."""
        return len(self._subSels) == 0

    def _toUnorderedFilter(self):
        """
        Return a filter vector representing the complex selector
        If necessary, the returned vector contains several chained filter vectors.
        """
        totvec = []
        for row in self._subSels:
            filt = []
            for sus in self._complex._poolElems:
                minv, maxv, realMax = 0, 0, 0
                for ss in row:
                    if sus in ss:
                        maxv += 1
                        if len(ss._states) == 1:
                            minv += 1
                    if sus in ss._subUnit._states:
                        realMax += 1
                filt.append((minv, maxv if maxv < realMax else COMPLEX_FILTER_MAX_VALUE))
            totvec.append(tuple(filt))
        return tuple(totvec)

    def _setOp(self, op, other):
        """
        Return the complex selector resulting from applying the operation op to self and other.
        """
        if isinstance(other, ComplexState):
            other = other._toComplexSelector()
        if not isinstance(other, ComplexSelector):
            raise TypeError('ComplexSelector objects can only be combined with other ' '_Complexselectors.')
        if other._complex is not self._complex:
            raise Exception('Cannot combine _Complexselectors associated to different complexes.')

        cs = copy.copy(self)

        if op == ComplexSelector._OR:
            cs._subSels = cs._subSels | other._subSels
        elif op == ComplexSelector._AND:
            cs._subSels = set()
            for ss1 in self._subSels:
                for ss2 in other._subSels:
                    inter = tuple(su1 & su2 for su1, su2 in zip(ss1, ss2))
                    if not any(su._isEmpty() for su in inter):
                        cs._subSels.add(inter)
        else:
            raise NotImplementedError()

        return cs

    def _getAllStates(self):
        """Return the set of all states that match the complex selector."""
        states = set()
        for ss in self._subSels:
            for comb in itertools.product(*[sus._states for sus in ss]):
                states.add(ComplexState(self._complex, comb, _parentCompSel=self))
        return sorted(states, key=lambda x: x.name)

    def _getMatchingStates(self, rcs):
        """
        Return the pairs of states needed to declare a reaction in which self is on the lhs and
        rcs on the rhs.
        """
        l2rstates = {}
        for lss in self._subSels:
            for comb in itertools.product(*[sorted(sus._states, key=lambda x: x.name) for sus in lss]):
                # For each state in the lhs complex, try to find the corresponding state
                # in the rhs complex
                destComb = None
                moreGeneralError = 0
                for rss in rcs._subSels:
                    tmpDestComb = []
                    ok = True
                    # For each subunit
                    for i, s in enumerate(comb):
                        # If there are identifiers and they are not matching, just skip
                        if rss[i]._id != lss[i]._id:
                            ok = False
                            break
                        rstates = rss[i]._states
                        lstates = lss[i]._states
                        if len(lstates) < len(rstates):
                            # If the destination state is more general than the source state
                            moreGeneralError += 1
                            ok = False
                            break
                        if len(rstates) == 1:
                            # If the destination state modifies the subunit
                            tmpDestComb.append(next(iter(rstates)))
                        elif s in rstates:
                            # If the destination state can contain the source state
                            tmpDestComb.append(s)
                        else:
                            # Otherwise, stop considering
                            ok = False
                            break
                    if ok:
                        if destComb is None:
                            destComb = tmpDestComb
                        elif ComplexState(rcs._complex, destComb) != ComplexState(rcs._complex, tmpDestComb):
                            # Only raise an exception if the states corresponding to the
                            # combinations are different
                            raise Exception(
                                f'ComplexSelector {self} can be matched with ComplexSelector '
                                f'{rcs} in more than one way. The reaction is ambiguous. State '
                                f'{comb} could at least be matched to {destComb} or to '
                                f'{tmpDestComb}. Try adding position or identification to '
                                f'SubUnitStates.'
                            )
                if destComb is None:
                    # Try to give a specific error message
                    if moreGeneralError == len(rcs._subSels):
                        raise Exception(
                            f'ComplexSelector {self} cannot be matched with ComplexSelector {rcs} '
                            f'because the latter is more general than the former. This can happen '
                            f'when the sides of a reaction using SubUnitStates are not balanced.'
                        )
                    else:
                        # If not possible, give a generic one
                        raise Exception(
                            f'ComplexSelector {self} cannot be matched with ComplexSelector {rcs} '
                            f'for state {comb}. The reaction is undefined.'
                        )
                lstate = ComplexState(self._complex, comb, _parentCompSel=self)
                rstate = ComplexState(rcs._complex, destComb, _parentCompSel=rcs)
                l2rstates.setdefault(lstate, set()).add(rstate)
        allPairs = []
        for ls, rss in l2rstates.items():
            for rs in sorted(rss, key=lambda s: s.name):
                allPairs.append((ls, rs, 1 / len(rss)))
        return allPairs

    def __or__(self, other):
        """Get a complex selector representing the union between the states of self and other

        :param other: The other complex state(s), has to be associated to the same
            :py:class:`Complex`.
        :type other: Union[:py:class:`ComplexSelector`, :py:class:`ComplexState`]

        :returns: A complex selector representing the union between the states in the current
            selector and the state(s) from 'other'.
        :rtype: :py:class:`ComplexSelector`

        Usage:
            States that have the first two subunits as ``S1A`` *or* the first as ``S1A`` and third
            as ``S2A``::

                Comp[S1A, S1A, :] | Comp[S1A, :, S2A]

            Beware, the above expression is **NOT** equivalent to ``Comp[S1A, :, :]``::

                >>> for state in (Comp[S1A, S1A, :] | Comp[S1A, :, S2A]):
                ...     print(state.name)
                ...
                Comp_S1A_S1A_S2A
                Comp_S1A_S1A_S2B
                Comp_S1A_S1B_S2A
                Comp_S1A_S1C_S2A
                >>> for state in Comp[S1A, :, :]:
                ...     print(state.name)
                ...
                Comp_S1A_S1A_S2A
                Comp_S1A_S1A_S2B
                Comp_S1A_S1B_S2A
                Comp_S1A_S1B_S2B
                Comp_S1A_S1C_S2A
                Comp_S1A_S1C_S2B

            In general, an OR operation between two :py:class:`ComplexSelector` is **NOT**
            equivalent to doing OR operations between their :py:class:`SubUnitSelector`.

        .. warning::
            Do not forget to wrap the overall expression in parentheses in order to avoid issues
            with operator precedence when used in reaction declaration.

        :meta public:
        """
        return self._setOp(ComplexSelector._OR, other)

    def __and__(self, other):
        """Get a complex selector representing the intersection between the states of self and other

        :param other: The other complex state(s), has to be associated to the same
            :py:class:`Complex`.
        :type other: Union[:py:class:`ComplexSelector`, :py:class:`ComplexState`]

        :returns: A complex selector representing the intersection between the states in the
            current selector and the state(s) from 'other'.
        :rtype: :py:class:`ComplexSelector`

        Usage:
            States that have the first two subunits as ``S1A`` *and* the first as ``S1A`` and third
            as ``S2A``::

                Comp[S1A, S1A, :] & Comp[S1A, :, S2A]

            The above expression corresponds to only one state: ``Comp[S1A, S1A, S2A]``. Note that,
            in constrast to :py:func:`__or__`, this does correspond to a series of AND operations
            between the :py:class:`SubUnitSelector`.

        If the intersection is empty, it returns an empty :py:class:`ComplexSelector` which
        might raise exceptions if later used in a reaction declaration.

        .. warning::
            Do not forget to wrap the overall expression in parentheses in order to avoid issues
            with operator precedence when used in reaction declaration.

        :meta public:
        """
        return self._setOp(ComplexSelector._AND, other)

    def __lshift__(self, subSels):
        """Get a complex selector representing the insertion of one or several subunit selector(s).

        This operator is mainly used with complexes that have a :py:func:`StrongOrdering`, meaning
        that subunits have specific positions. With this ordering, ``Comp[S1A, S1B, S2A]`` is
        not the same as ``Comp[S1B, S1A, S2A]`` (the order of the first two subunits was changed).

        Since, with this order function, the specific positions of the subunits matter, it would be
        cumbersome to write ``Comp[S1A, :, :] | Comp[:, S1A, :]`` to express "a SU1 subunit is in
        state S1A, in any position". And it would become harder and harder as the number of
        identical subunits grows. The ``<<`` operator helps expressing this type of selectors.

        :param subSels: One or several subunit states or selectors. If several, they should be
            provided in the form of a :py:class:`ReactionSide` (``S1A + S1B`` for example).
        :type subSels: Union[:py:class:`SubUnitState`, :py:class:`SubUnitSelector`,
            :py:class:`ReactionSide`]

        :returns: A complex selector to which subunit selectors present in subSels were
            added in all possible positions.
        :rtype: :py:class:`ComplexSelector`

        Assuming we have a complex that has 3 identical subunits, declared like so::

            with mdl:
                S3A, S3B = SubUnitState.Create()
                S3U = SubUnit.Create([S3A, S3B])
                Comp2 = Complex([S3U, S3U, S3U], order=StrongOrdering, statesAsSpecies=True)

        Usage of the ``<<`` operator::

            (Comp2[...] << S3A)       # Equivalent to:
                                      # (Comp2[S3A,  : ,  : ] |
                                      #  Comp2[ : , S3A,  : ] |
                                      #  Comp2[ : ,  : , S3A])

            (Comp2[...] << S3A + S3B) # Equivalent to:
                                      # (Comp2[S3A, S3B,  : ] |
                                      #  Comp2[S3A,  : , S3B] |
                                      #  Comp2[S3B, S3A,  : ] |
                                      #  Comp2[ : , S3A, S3B] |
                                      #  Comp2[S3B,  : , S3A] |
                                      #  Comp2[ : , S3B, S3A])

        .. warning::
            Do not forget to wrap the overall expression in parentheses in order to avoid issues
            with operator precedence when used in reaction declaration.

        :meta public:
        """
        if isinstance(subSels, (ReactionElement, ReactingElement)):
            subSels = subSels._toReactionSide()
        if not isinstance(subSels, ReactionSide):
            raise TypeError(f'Expected reaction elements, got {subSels} instead.')
        if not all(isinstance(el._elem, SubUnitSelector) for el in subSels.reacElems):
            raise TypeError(f'Some elements in "{subSels}" are not subunit states or selectors.')

        return self._mergeSubSels(subSels, enforceNoLocation=True)

    def __iter__(self):
        """Iterator over all the possible states of the complex selector.

        :returns: An iterator that iterates over all unique (as defined by the custom ordering)
            :py:class:`ComplexState` of the complex selector.
        :rtype: Iterator[:py:class:`ComplexState`]

        Usage::

            >>> selector = Comp[S1A, :, S2A]
            >>> for state in selector:
            ...     print(state.name)
            ...
            Comp_S1A_S1A_S2A
            Comp_S1A_S1B_S2A
            Comp_S1A_S1C_S2A

        :meta public:
        """
        return iter(self._getAllStates())

    def __len__(self):
        """Return the number of states covered by the complex selector.

        :returns: The number of unique states covered by the complex selector.
        :rtype: int

        Usage::

            >>> print(len(Comp[S1A, :, S2A]))
            3

        :meta public:
        """
        return len(self._getAllStates())

    def _simPathWalkExpand(self):
        """Return an iterable of the elements that should be part of ResultSelector paths."""
        if self._areStatesAsSpecies():
            return self.__iter__()
        else:
            return [self]

    def _simPathCombinerClass(self):
        """Return the class that needs to be used to combine expanded elements."""
        if self._areStatesAsSpecies():
            return nutils.SumSimPathComb
        else:
            return None

    def _simPathCheckParent(self):
        """
        Determines whether the object needs to be in the parent's children in order to be added to
        a simulation path
        """
        return False

    def _mergeSubSels(self, subSels, enforceNoLocation=False):
        """
        Return a complex selector to which SubUnitSelectors present in subSels were added in all
        possible positions if they were unordered or in their position if specified.
        """
        if enforceNoLocation and subSels._hasLocation():
            raise Exception('Location cues are not accepted when building a ComplexSelector.')
        cs = self
        for re in subSels.reacElems:
            for s in range(re.stoich):
                newcs = None
                ss = re._elem
                if ss._pos is None:
                    # If the subunitselector has no specified positions, generate the different
                    # positions available
                    for p in cs._getFreeSlots(ss):
                        ss2 = copy.copy(ss)
                        ss2._pos = p
                        if newcs is None:
                            newcs = cs._applysSubSel(ss2)
                        else:
                            newcs |= cs._applysSubSel(ss2)
                else:
                    newcs = cs._applysSubSel(ss)
                if newcs is None or newcs._isEmpty():
                    raise Exception(f'Complex selector {cs} is incompatible with SubUnitSelector ' f'{ss}.')
                cs = newcs
        return cs

    def _applysSubSel(self, ss):
        res = copy.copy(self)
        res._subSels = set()
        for row in self._subSels:
            if len(row[ss._pos]._states) > len(ss._states) and not (row[ss._pos] & ss)._isEmpty():
                newSus = row[ss._pos] & ss
                res._subSels.add(row[: ss._pos] + (newSus,) + row[ss._pos + 1 :])
        return res

    def _getFreeSlots(self, ss):
        """
        Return the subunit positions where ss could be applied. Only considers the positions that
        are completely undetermined.
        """
        for i in range(len(self._complex._subUnits)):
            for st in self._subSels:
                if st[i]._subUnit is ss._subUnit and st[i]._isFull():
                    yield i
                    break

    def _getFullyDeterminedState(self):
        """
        If the ComplexSelector only represents a single fully determined state, return this state,
        otherwise, return None.
        """
        fullState = None
        for row in self._subSels:
            if any(len(ss._states) != 1 for ss in row):
                return None
            state = tuple(next(iter(ss._states)) for ss in row)
            if fullState is not None and state != fullState:
                return None
            fullState = state
        return ComplexState(self._complex, fullState, _parentCompSel=self)

    def _getUpdateVects(self, rcs):
        """
        Return a list of pool vectors tuples corresponding to the complex reaction.
        """
        if len(rcs._subSels) > 1:
            raise Exception(
                f'Cannot declare a reaction with composite complex selectors on the right hand '
                f'side: {rcs}'
            )

        rrow = next(iter(rcs._subSels))
        comp = self._complex
        n = len(comp._poolElems)

        filt = self._toUnorderedFilter()
        updates = set()
        for lrow in self._subSels:
            # For each pair of SubUnitSelector rows
            lds, rds = [], []
            for lss, rss in zip(lrow, rrow):
                if lss != rss and not lss._isFull():
                    if len(rss._states) == 1:
                        rsus = next(iter(rss))
                        lds.append(list(lss._states))
                        rds.append(rsus)
                    elif any(lsus not in rss for lsus in lss):
                        # Should never happen
                        raise Exception()

            # For each possible state of the considered SubUnits
            for comb in itertools.product(*lds):
                req = [0] * n
                upd = [0] * n
                for lsus, rsus in zip(comb, rds):
                    req[comp._sus2PoolInd[lsus]] += 1
                    upd[comp._sus2PoolInd[lsus]] += -1
                    upd[comp._sus2PoolInd[rsus]] += 1
                updates.add((tuple(req), tuple(upd)))

        react = [0] * n
        return [(filt, react, sorted(updates))]

    @staticmethod
    def _removeMoreSpecific(lst, func):
        """
        Takes a list of identical size iterables and returns a simplified list by removing the
        elements that are strictly more specific than other elements. The function 'func' should
        take two elements and return True if the first one is more specific than the other.
        """
        finallst = []
        while len(lst) > 0:
            delInd = 0
            for i, vect1 in enumerate(lst):
                totLst = finallst + lst[0:i] + lst[i + 1 :]
                if any(func(vect1, vect2) for vect2 in totLst):
                    # If vect1 is always more specific than at least one other vector, remove vect1
                    break
                delInd += 1
            finallst += lst[0:delInd]
            lst = lst[delInd + 1 :]
        return finallst

    def __repr__(self):
        return ' | '.join(f'{self._complex.name}[' + ', '.join(map(str, row)) + ']' for row in self._subSels)

    def __eq__(self, other):
        return isinstance(other, ComplexSelector) and self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        if self._areStatesAsSpecies():
            return Species._elemStr
        else:
            return Complex._elemStr

    def _solverId(self):
        """Return the id that is used to identify an object in the solver."""
        if self._areStatesAsSpecies():
            raise Exception('_solverId should never be called when statesAsSpecied is True.')
        return (
            self._complex.name,
            self._toUnorderedFilter(),
        )


@nutils.FreezeAfterInit
class _ComplexReactants(_ComplexReactionElement):
    """
    Represents the combination of a complex selectors and a list of reactants involved in a
    reaction. The reactants are SubUnitSelectors and can be associated with any compatible position
    in the complex selector.
    """

    def __init__(self, filt, reactants, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._filter = filt
        self._reactants = []
        for re in reactants.reacElems:
            self._reactants += [re._elem] * re.stoich

    def _areStatesAsSpecies(self):
        return self._filter._areStatesAsSpecies()

    def _hasSameReactants(self, cr):
        """Return True if cr has the same reactants in the same stoichiometry."""
        return set((r, self._reactants.count(r)) for r in self._reactants) == set(
            (r, self._reactants.count(r)) for r in cr._reactants
        )

    def _getMatchingStates(self, rcs):
        """
        Return the pairs of states needed to declare a reaction in which self is on the lhs and
        rcs on the rhs.
        """
        # Compute the pairs of reactants
        reactantPairs = self._getMatchingReactants(rcs)
        # Order them from the smallest to the biggest to ensure correct positions computation
        reactantPairs = sorted(reactantPairs, key=lambda x: (len(x[0]._states), hash(x[0])))

        groupedStatePairs = {}
        for filtSubSels in self._filter._subSels:
            # Find all possible positions for the pairs of reactants
            for posReacPairs in _ComplexReactants._getAllReactantsPositions(filtSubSels, reactantPairs):
                # Get resulting src and dst rows
                src = filtSubSels
                dst = filtSubSels
                reactantPos = set()
                for lrp, rrp, p in posReacPairs:
                    lrp = copy.copy(lrp)
                    lrp._pos = p
                    rrp = copy.copy(rrp)
                    rrp._pos = p
                    src = src[0:p] + (lrp & src[p],) + src[p + 1 :]
                    if (rrp & dst[p])._isEmpty():
                        dst = dst[0:p] + (rrp,) + dst[p + 1 :]
                    else:
                        dst = dst[0:p] + (rrp & dst[p],) + dst[p + 1 :]
                    reactantPos.add(p)
                # Generate states from src
                for comb in itertools.product(*[sorted(sus._states, key=lambda x: x.name) for sus in src]):
                    # Compute dst state
                    dstState, dstModState = [], []
                    for srcState, srcSubSel, dstSubSel, p in zip(comb, src, dst, range(len(comb))):
                        if len(dstSubSel._states) == 1:
                            sus = next(iter(dstSubSel._states))
                            dstState.append(sus)
                            # Tag the destination state as modified, useful if destination is the same as src
                            dstModState.append(sus._getModif() if p in reactantPos else sus)
                        elif srcSubSel == dstSubSel:
                            dstState.append(srcState)
                            # Tag the destination state as modified, useful if destination is the same as src
                            dstModState.append(srcState._getModif() if p in reactantPos else srcState)
                        else:
                            raise Exception(
                                f'Cannot tranform state {srcState} associated with '
                                f'SubUnitSelector {srcSubSel} to SubUnitSelector '
                                f'{dstSubSel}, reaction is ambiguous.'
                            )
                    # Get the corresponding generic state (identical if using ordered Complexes)
                    srcos = tuple(comb)
                    dstos = tuple(dstState)
                    dstoms = tuple(dstModState)
                    srcgs = ComplexState(self._filter._complex, srcos, _parentCompSel=self)
                    dstgs = ComplexState(rcs._filter._complex, dstos, _parentCompSel=rcs)
                    # Add the state pairs to the pair counting structure
                    key = (srcgs, dstgs)
                    groupedStatePairs.setdefault(key, {}).setdefault(srcos, set()).add(dstoms)

        if len(groupedStatePairs) == 0:
            raise Exception(
                f'It was not possible to find a combination of subunit states position compatible with the '
                f'complex selector.'
            )

        allPairs = []
        # Count the pairs to determine the required rate multiplier
        for groupedPair, pairs in groupedStatePairs.items():
            # Count total number of subpairs
            totSubPairs = sum(len(dstStates) for srcos, dstStates in pairs.items())
            rateMult = totSubPairs / len(pairs.keys())
            allPairs.append(groupedPair + (rateMult,))

        return allPairs

    @staticmethod
    def _getAllReactantsPositions(filtSubSels, reactantPairs):
        """
        Compute all the possible positions for reactant pairs given a filter.
        """
        if len(reactantPairs) == 0:
            yield []
        else:
            lrp, rrp = reactantPairs[0]
            # Get the possible positions for the first pair
            if lrp._pos is None:
                positions = list(_ComplexReactants._getFreeSlotsForReactantsPair(filtSubSels, lrp, rrp))
            else:
                if (filtSubSels[lrp._pos] & lrp)._isEmpty():
                    positions = []
                else:
                    positions = [lrp._pos]

            # Recursively compute the other positions
            for p in positions:
                empty = lrp & ~lrp
                newFilt = filtSubSels[0:p] + (empty,) + filtSubSels[p + 1 :]
                for subRPpos in _ComplexReactants._getAllReactantsPositions(newFilt, reactantPairs[1:]):
                    if len(subRPpos) == len(reactantPairs[1:]) :
                        yield [(lrp, rrp, p)] + subRPpos

    @staticmethod
    def _getFreeSlotsForReactantsPair(filtSubSels, lrp, rrp):
        """
        Return the positions in filtSubSels for which both lrp and rrp could fit. Only considers
        the positions that are completely undetermined.
        """
        for i, ss in enumerate(filtSubSels):
            unpos = copy.copy(ss)
            unpos._pos = None
            if lrp._subUnit is ss._subUnit and not (lrp & unpos)._isEmpty():
                yield i

    def _getMatchingReactants(self, rcs):
        """
        Check that the reaction is balanced and return the pairs of reactants needed to
        Generate the lhs and rhs states.
        """
        # Group by subunit
        lhsSU2React = {}
        for r in self._reactants:
            lhsSU2React.setdefault(r._subUnit, []).append(r)
        rhsSU2React = {}
        for r in rcs._reactants:
            rhsSU2React.setdefault(r._subUnit, []).append(r)

        # Check that the same subUnits are covered by both sides
        if set(lhsSU2React.keys()) != set(rhsSU2React.keys()):
            raise Exception(
                f'The Subunits involved in the LHS of the reaction ({lhsSU2React.keys()}) '
                f'are different from those involved in the RHS ({rhsSU2React.keys()})'
            )

        reactantPairs = []
        # For each subunit
        for su, lhsReacts in lhsSU2React.items():
            # First perform basic checks
            if su not in rhsSU2React:
                raise Exception(f'SubUnitSelectors {lhsReacts} have no match in the RHS of the reaction.')
            rhsReacts = rhsSU2React[su]
            if len(lhsReacts) != len(rhsReacts):
                raise Exception(
                    f'The number of SubUnitSelectors associated to SubUnit {su} should '
                    f'be the same on both sides of the reaction.'
                )
            # Try to find an unambiguous matching
            for lhsr in lhsReacts:
                tmpDest = []
                for rhsr in rhsReacts:
                    # If there are identifiers or positions and they are not matching, just skip
                    if (
                        lhsr._id == rhsr._id
                        and lhsr._pos == rhsr._pos
                        and len(lhsr._states) >= len(rhsr._states)
                    ):
                        tmpDest.append(rhsr)
                # Check that the matching is possible
                if len(tmpDest) == 0:
                    raise Exception(
                        f'SubunitSelector {lhsr} from the LHS of the reaction could '
                        f'not be matched with any SubUnitSelector from the RHS.'
                    )
                if len(set(tmpDest)) > 1:
                    raise Exception(
                        f'The reaction is ambiguous: SubUnitSelector {lhsr} from the LHS '
                        f'could be matched with any of {tmpDest} from the RHS.'
                    )
                # Add the pair
                reactantPairs.append((lhsr, tmpDest[0]))
                # Remove therhs one that just got added
                rhsReacts.remove(tmpDest[0])
            # If there are SubUnitSelector from the rhs that were not added
            if len(rhsReacts) > 0:
                raise Exception(
                    f'SubUnitSelectors from the RHS {rhsReacts} could not be matched '
                    f'with any SubUnitSelectors from the LHS.'
                )

        # TODO Not urgent: Throw exception if reactantPairs is empty?
        return reactantPairs

    def _getUpdateVects(self, rcs):
        """
        Return a list of pool vectors tuples corresponding to the complex reaction.
        """
        comp = self._filter._complex
        n = len(comp._poolElems)
        allFilts = self._filter._toUnorderedFilter()

        allVects = []

        for filt in allFilts:
            allReactants = [[0] * n]
            allUpdates = [([0] * n, [0] * n)]
            # Get matching reactants
            for lhsr, rhsr in self._getMatchingReactants(rcs):
                newReactants = []
                newUpdates = []
                for lsus in sorted(lhsr._states, key=lambda x: x.name):
                    for reac, (req, upd) in zip(allReactants, allUpdates):
                        nr = copy.copy(reac)
                        nreq, nupd = copy.copy(req), copy.copy(upd)
                        ind = comp._sus2PoolInd[lsus]

                        nr[ind] += 1
                        nreq[ind] += 1
                        nupd[ind] += -1
                        if len(rhsr._states) == 1:
                            rsus = next(iter(rhsr._states))
                        elif lsus in rhsr._states:
                            rsus = lsus
                        else:
                            # This should never happen
                            raise Exception()
                        nupd[comp._sus2PoolInd[rsus]] += 1

                        newReactants.append(nr)
                        newUpdates.append((nreq, nupd))

                allReactants = newReactants
                allUpdates = newUpdates

            # Update the filters with the reactants
            for reac, (req, upd) in zip(allReactants, allUpdates):
                newFilt = []
                for i, r in enumerate(reac):
                    # SubUnitStates from the filter can also be reactants
                    newFilt.append((max(filt[i][0], r), filt[i][1]))
                allVects.append(([newFilt], reac, [(req, upd)]))

        # Remove duplicates or more specific
        def func(a, b):
            return all(
                ra == rb and a[0][0][i][0] >= b[0][0][i][0] for i, ra, rb in zip(range(len(a[2])), a[2], b[2])
            )

        return ComplexSelector._removeMoreSpecific(allVects, func)

    def _getReferenceObject(self):
        """
        Return the object this object was derived from. Useful for getting the complex associated
        with a complex selector, etc.
        """
        return self._filter._getReferenceObject()

    @property
    def _complex(self):
        return self._filter._complex

    def __eq__(self, other):
        return isinstance(other, _ComplexReactants) and self._filter == other._filter

    def __hash__(self):
        return hash(self._filter)


###################################################################################################
# Reactions


@nutils.FreezeAfterInit
class ReactionSide:
    """Unordered collection of reactants

    Represents one of the sides of a reaction. Reactants can be added with `+`.
    Reactants are unordered.

    .. note::
        This class should not be instantiated by the user, it is only documented for clarity.
    """

    def __init__(self, elems, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reacElems = list(elems)

    def _add(self, other):
        res = copy.copy(self)
        if other is None or self.reacElems is None:
            raise Exception('Cannot add reactant to None.')
        if isinstance(other, ReactingElement):
            res.reacElems = self.reacElems + [other]
        elif isinstance(other, ReactionElement):
            res.reacElems = self.reacElems + [ReactingElement(other, 1)]
        elif isinstance(other, ReactionSide):
            res.reacElems = self.reacElems + other.reacElems
        else:
            raise TypeError(f'Cannot add {other} to a reaction side.')
        return res

    def __add__(self, other):
        """Get the reaction side resulting from the addition of other to self

        :param other: The object to be added to the reaction side
        :type other: Union[:py:class:`ReactionElement`, :py:class:`ReactionSide`]

        :returns: The reaction side resulting from the addition of other to self
        :rtype: :py:class:`ReactionSide`

        :meta public:
        """
        return self._add(other)

    def __radd__(self, other):
        return self._add(other)

    def __or__(self, other):
        """
        Raise an exception for the common case in which a SubUnitSelector made from '|' operator
        was not wrapped in parentheses.
        """
        if isinstance(other, SubUnitState):
            raise SyntaxError(
                'SubUnitSelectors with "|" operator need to be wrapped by '
                'parentheses when used in reactions.'
            )
        else:
            raise SyntaxError('Cannot use "|" operator in reactions except for combining SubUnitStates.')

    def _SetLoc(self, loc):
        for s in self.reacElems:
            s.loc = loc

    def _GetStepsElems(self, loc=None):
        res = []
        for re in self.reacElems:
            if loc is None or re.loc == loc:
                res += re._getStepsObjects() * re.stoich
        return res

    def _getElemsOfType(self, *cls):
        for e in self.reacElems:
            for c in cls:
                if isinstance(e._elem, c):
                    yield e
                    break

    def _hasElemsOfType(self, *cls):
        try:
            next(self._getElemsOfType(*cls))
            return True
        except StopIteration:
            return False

    def _isEmpty(self):
        return len(self.reacElems) == 0

    def _getAllElems(self, loc=None):
        res = set()
        for re in self.reacElems:
            if loc is None or re.loc == loc:
                res.add(re._elem)
                if isinstance(re._elem, ComplexState):
                    res.add(re._elem._comp)
                if isinstance(re._elem, ComplexSelector):
                    res.add(re._elem._complex)
        return res

    def _isSurfaceSide(self):
        if len(self.reacElems) == 0:
            return None
        elif all(s.loc is None for s in self.reacElems):
            return False
        elif all(s.loc is not None for s in self.reacElems):
            return True
        else:
            raise Exception(
                'Cannot mix localized and non localized reacElems in the same reaction side.'
            )

    def _hasSingleLocation(self):
        return len(set(s.loc for s in self.reacElems)) < 2

    def _hasVolumeLocation(self):
        return all(s.loc not in ALL_SURF_LOCATIONS for s in self.reacElems)

    def _hasLocation(self):
        return any(s.loc is not None for s in self.reacElems)

    def _getLocation(self):
        if len(self.reacElems) > 0 and self._hasSingleLocation():
            return self.reacElems[0].loc
        else:
            return None

    def __repr__(self):
        if len(self.reacElems) > 0:
            return ' + '.join(map(str, self.reacElems))
        else:
            return 'Void'

    def __iter__(self):
        return iter(self.reacElems)

    def __len__(self):
        return len(self.reacElems)

    def _toReactionSide(self):
        return self

    def _getFlatList(self):
        """Return an equivalent list of reacting elements that all have stoich = 1"""
        res = []
        for re in self:
            res += [ReactingElement(re._elem, 1, re.loc)] * re.stoich
        return res

    def _makeCopy(self):
        """Make a copy of the reaction side."""
        res = []
        for e in self.reacElems:
            e2 = copy.copy(e)
            e2._elem = e._elem._makeCopy()
            res.append(e2)
        return ReactionSide(res)

    def _getFinal(self, globCompSels):
        """
        Return a ReactionSide that only contains ComplexSelectors and Species. Transform all
        SubUnitSelectors into ComplexSelectors.
        """
        # Split the ReactionSide into ComplexSelectors, SubUnitSelectors, and Species
        compSelElems = list(self._getElemsOfType(ComplexSelector))
        specElems = list(self._getElemsOfType(Species, LinkSpecies))
        compStateElems = list(self._getElemsOfType(ComplexState))
        subSelElems = list(self._getElemsOfType(SubUnitSelector))

        # Transform complex states in complex selectors if they are real steps complexes
        tmpCS = []
        for cs in compStateElems:
            if not cs._elem._comp._areStatesAsSpecies():
                cs._elem = cs._elem._toComplexSelector()
                compSelElems.append(cs)
            else:
                tmpCS.append(cs)
        compStateElems = tmpCS

        # Automatically attribute ComplexSelectors to SubUnitSelectors
        compSels2subSels = {}
        for re in subSelElems:
            ss = re._elem
            if ss._compSel is None:
                if len(globCompSels) == 0:
                    raise Exception(
                        f'The reaction was declared outside of a ComplexSelector and '
                        f'SubUnitSelector {ss} was not associated with any '
                        f'ComplexSelector.'
                    )
                potCompSels = [c for c in globCompSels if ss._subUnit in c._complex._subUnits]
                if len(potCompSels) == 0:
                    raise Exception(
                        f'SubUnitSelector {ss} does not match any of the complex '
                        f'selectors: {globCompSels}.'
                    )
                if len(potCompSels) > 1:
                    raise Exception(
                        f'SubUnitSelector {ss} does not have a clear complex selector.'
                        f' It could be matched with any of these: {potCompSels}.'
                    )
                compSel = potCompSels[0]
            else:
                compSel = ss._compSel
            compSels2subSels.setdefault(compSel, ReactionSide([]))
            # Check that the positions cues are valid for this complex
            ss._checkPos(compSel._complex)
            # Clear compSel, it should not be used from this point on
            ss._compSel = None
            compSels2subSels[compSel] += re
        # Generate ComplexReactants objects
        resCs = ReactionSide([])
        for cs, subSels in compSels2subSels.items():
            if not subSels._hasSingleLocation():
                raise Exception(f'SubUnits from the same complex have different locations: {subSels}.')
            if cs in [e._elem for e in compSelElems]:
                raise Exception(
                    f'Cannot mix ComplexSelector and SubunitSelectors for the same '
                    f'ComplexSelector in a reaction side. Either use one complex '
                    f'selector or only SubUnitSelectors.'
                )
            compReacts = _ComplexReactants(cs, subSels)
            resCs += ReactionSide([ReactingElement(compReacts, 1, subSels._getLocation())])

        # Reassemble the ReactionSide
        finalSide = (
            ReactionSide(compSelElems) + resCs + ReactionSide(specElems) + ReactionSide(compStateElems)
        )

        # Check that each ComplexSelector is only present once
        allCompSels = list(finalSide._getElemsOfType(ComplexSelector))
        if len(allCompSels) > 0:
            compSels, stoichs = zip(*[(re._elem, re.stoich) for re in allCompSels])
            for cs, st in zip(compSels, stoichs):
                if st > 1 or compSels.count(cs) > 1:
                    raise Exception(
                        f'Complex selector {cs} was used more than once in the reaction side.'
                    )

        return finalSide


@nutils.FreezeAfterInit
class _SubReaction(nutils.NamedObject):
    """
    Wrapper for a single steps reaction object.
    """

    _elemStrDict = {
        stepslib._py_Reac: 'Reac',
        stepslib._py_SReac: 'SReac',
        stepslib._py_VDepSReac: 'VDepSReac',
        stepslib._py_VesSReac: 'VesSReac',
        stepslib._py_RaftSReac: 'RaftSReac',
        stepslib._py_ComplexReac: 'ComplexReac',
        stepslib._py_ComplexSReac: 'ComplexSReac',
    }

    def __init__(self, parent, stepsReac, lhs=None, rateMult=1, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent = parent
        self._stepsReac = stepsReac
        cls = stepsReac.__class__
        if cls not in _SubReaction._elemStrDict:
            raise TypeError(f'Unknown steps reaction type {cls}.')
        self._elemStr = _SubReaction._elemStrDict[cls]
        self._lhs = lhs
        self._rateMult = rateMult

    def _getRateUnits(self):
        """Return rate units of the reaction"""
        return self._parent._getRateUnits()

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self._stepsReac]

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return self._elemStr
    
    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        return {'obj_type': 'Reac', 'obj_id': self.name}

    def _solverSetValue(self, valName, v):
        """
        Return the value that should actually be set in the solver when value 'v' is given by
        the user.
        """
        if isinstance(v, CompDepFunc):
            if self._lhs is None:
                raise NotImplementedError()
            v = v(self._lhs)

        if valName == 'K':
            return self._rateMult * v
        else:
            return v


@nutils.FreezeAfterInit
class _SubReactionList(nutils.SolverPathObject, nutils.ParameterizedObject, list):
    """
    Represents a unidirectional reaction consisting in a list of single steps reaction objects.
    """

    def __init__(self, parent, reacName, LRP, compEvs, isFwd=True, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent = parent
        self._reacName = reacName
        self._LRP = LRP
        self._compEvs = compEvs
        self._isFwd = isFwd

        # Add the reaction name as parameter if the reaction is not anonymous
        if not isinstance(self._parent, AnonymousReaction):
            self._setParameter('Reaction ID', self._parent.name + ("['fwd']" if self._isFwd else "['bkw']"))
        # Add the complex selector as parameter, if any
        _, compSels, *_ = self._parent._usedObjects
        if len(compSels) > 0:
            self._setParameter('Complex', nutils._ParamList(compSels))

        self.K = 0

        self._reacElems = {side: {loc: set() for loc in ALL_LOCATIONS} for side in 'LR'}

    def append(self, item):
        if not isinstance(item, _SubReaction):
            raise TypeError('Only _SubReaction objects can be added to _SubReactionList.')
        super().append(item)

    def _getRateUnits(self):
        """Return rate units of the reaction"""
        units = set()

        # some reactions are volume-volume regardless of their content
        isVolVolReac = self._parent._isVesicleSurfReac()

        for l, r, rm in self._LRP:
            l = l._getFlatList()
            molExp = 1 - len(l)
            unit = ''
            if molExp != 0:
                if isVolVolReac or len(l) == 0 or any(s.loc != Location.SURF for s in l):
                    unit = f'M'
                else:
                    unit = f'(mol m^-2)'
                unit += f'^{molExp} ' if molExp != 1 else ' '
            units.add(unit + 's^-1')

        if len(units) > 1:
            raise Exception(f'Inconsistent rate units in reaction {self}: {units}')

        return nutils.Units(list(units)[0])

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=_getRateUnits)
    def K(self):
        """Reaction rate(s)

        If the reaction is reversible, this value corresponds to a tuple composed of the forward
        and backward rate. The rate can be a special function that depends e.g. on membrane potential
        (:py:class:`VDepRate`) or on complex states (:py:class:`CompDepRate`).

        Usage examples::

            # Bidirectional reactions
            SA + SB <r[1]> 2 * SC
            r[1].K = forward_rate, backward_rate

            # Unidirectional reactions
            SC >r[2]> SD
            r[1].K = single_rate

        The reaction rate(s) default(s) to 0.

        :type: Union[float, :py:class:`steps.API_2.utils.Parameter`, :py:class:`XDepFunc`]
        """
        pass

    @K.setter
    @nutils.ParameterizedObject.RegisterSetter(units=_getRateUnits)
    def K(self, rate):
        if not isinstance(rate, (numbers.Number, XDepFunc)):
            raise TypeError(f'{rate} cannot be used as a rate.')

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def Dependencies(self):
        """Reaction dependencies

        A set of species that need to be present in order for the reaction to happen but their amount
        will not have any effect on the reaction rate.

        If the reaction is reversible, this value corresponds to a tuple composed of the forward
        and backward dependencies. The dependencies can be a special function that depends e.g. on
        complex states (:py:class:`CompDepFunc`).

        Note that for now, only species on the surfaces of rafts and vesicles can be used as dependencies.

        The dependencies should be given in the shape of a reaction side, for example: `SA.v + SC.v`

        Usage examples::

            # Bidirectional reactions
            SA.v + SB.o <r[1]> SA.v + SB.v
            r[1].K = fwd_rate, bkw_rate

            # SC is a dependency for the forward reaction, SD for the backward:
            r[1].Dependencies = SC.v, SD.v

            # Unidirectional reactions
            SC.v >r[2]> SD.v
            r[1].K = single_rate
            r[1].Dependencies = SA.v

        The dependencies are empty by default.

        :type: Union[:py:class:`steps.API_2.ReactionElement`, :py:class:`steps.API_2.ReactionSide`,
            :py:class:`CompDepFunc`]
        """
        # Return a default value if the property was not set
        return None

    @Dependencies.setter
    @nutils.ParameterizedObject.RegisterSetter(units=None)
    def Dependencies(self, deps):
        if not isinstance(deps, (ReactionElement, ReactionSide, CompDepFunc, type(None))):
            raise TypeError(f'{deps} cannot be used as a dependency.')

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def AntiDependencies(self):
        """Reaction anti-dependencies

        A set of species that, if present, will prevent the reaction from happenning.

        If the reaction is reversible, this value corresponds to a tuple composed of the forward
        and backward anti-dependencies. The anti-dependencies can be a special function that depends e.g. on
        complex states (:py:class:`CompDepFunc`).

        Note that for now, only species on the surfaces of rafts can be used as anti-dependencies.

        The anti-dependencies should be given in the shape of a reaction side, for example: `SA.r + SC.r`
        In that example, the reaction will be prevented if either `SA` or `SC` are present in the raft.
        If e.g. `2*SA.R + SC.r` is used, the reaction will be prevented is at least two `SA` are present
        or at least one `SC`.

        Usage examples::

            # Bidirectional reactions
            SA.r + SB.s <r[1]> SA.r + SB.r
            r[1].K = fwd_rate, bkw_rate

            # SC will prevent the forward reaction, SD the backward:
            r[1].AntiDependencies = SC.r, SD.r

            # Unidirectional reactions
            SC.r >r[2]> SD.r
            r[1].K = single_rate
            r[1].AntiDependencies = SA.r

        The anti-dependencies are empty by default.

        :type: Union[:py:class:`steps.API_2.ReactionElement`, :py:class:`steps.API_2.ReactionSide`,
            :py:class:`CompDepFunc`]
        """
        # Return a default value if the property was not set
        return None

    @AntiDependencies.setter
    @nutils.ParameterizedObject.RegisterSetter(units=None)
    def AntiDependencies(self, deps):
        if not isinstance(deps, (ReactionElement, ReactionSide, CompDepFunc, type(None))):
            raise TypeError(f'{deps} cannot be used as an anti-dependency.')

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def Immobilization(self):
        """Reaction immobilization flag

        The immobilization flag is only applicable to reactions involving rafts or vesicles, it can
        either be :py:data:`IMMOBILIZING`, :py:data:`NO_EFFECT` or :py:data:`MOBILIZING`.

        If the reaction is reversible, this value corresponds to a tuple composed of the forward
        and backward immobilization flags. Immobilization flags can be a special function that depends e.g. on
        complex states (:py:class:`CompDepFunc`).

        Usage examples::

            # Bidirectional reactions
            SA.v + SB.o <r[1]> SA.v + SB.v
            r[1].K = fwd_rate, bkw_rate

            # Forward reaction immobilizes, backward reaction does not:
            r[1].Immobilization = IMMOBILIZING, NO_EFFECT

            # Unidirectional reactions
            SC.v >r[2]> SD.v
            r[1].K = single_rate
            r[1].Immobilization = IMMOBILIZING

        The immobilization flag defaults to :py:data:`NO_EFFECT`.

        :type: Union[:py:data:`steps.API_2.model.IMMOBILIZING`, :py:data:`steps.API_2.model.NO_EFFECT`,
            :py:data:`steps.API_2.model.MOBILIZING`, :py:class:`CompDepFunc`]
        """
        # Return a default value if the property was not set
        return NO_EFFECT

    @Immobilization.setter
    @nutils.ParameterizedObject.RegisterSetter(units=None)
    def Immobilization(self, val):
        if not isinstance(val, Immobilization):
            raise TypeError(f'Expected an immobilization flag, got {val} instead.')

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m'))
    def MaxDistance(self):
        """Maximum distance (in m) between vesicle and patch for the reaction to occur

        Note that only reactions involving vesicles can have a maximum distance parameter.

        If the reaction is reversible, this value corresponds to a tuple composed of the forward
        and backward maximum distances. Maximum distances can be a special function that depends e.g. on
        complex states (:py:class:`CompDepFunc`).

        Usage examples::

            # Bidirectional reactions
            SA.v + SB.s <r[1]> SA.v + SB.v
            r[1].K = fwd_rate, bkw_rate

            # Forward reaction has a maximum distance, backward reaction does not:
            r[1].MaxDistance = 0.1e-6, None

            # Unidirectional reactions
            SC.v + SD.s >r[2]> 2 * SD.v
            r[1].K = single_rate
            r[1].MaxDistance = 5e-7

        The maximum distance defaults to None, i.e. no maximum distance.

        :type: Union[float, :py:class:`CompDepFunc`]
        """
        # Return a default value if the property was not set
        return None

    @MaxDistance.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('m'))
    def MaxDistance(self, md):
        if md is not None and (isinstance(md, bool) or not isinstance(md, (numbers.Number, CompDepFunc))):
            raise TypeError(f'{md} cannot be used as maximum distance.')

    def __repr__(self):
        if self._isFwd:
            return f'{self._parent.lhs} → {self._parent.rhs}'
        else:
            return f'{self._parent.rhs} → {self._parent.lhs}'

    @classmethod
    def _getDisplayName(cls):
        """Name that will be used as column title during parameter export
        Defaults to class name.
        """
        return Reaction._getDisplayName()

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        res = []
        for sr in self:
            res += sr._getStepsObjects()
        return res

    def _simPathWalkExpand(self):
        """Return an iterable of the elements that should be part of ResultSelector paths."""
        return self

    def _getAllElems(self, loc, sides='LR'):
        elems = set()
        for side in sides:
            elems |= self._reacElems[side][loc]
        return elems

    def _decl(self):
        """
        Return a string that describes the declaration of the object. It returns its name by default
        but it can contain e.g. the file and line at which it was declared.
        """
        return self._parent._decl()

    def _addToSteps(self):
        """
        Declare and add all steps reactions corresponding to the (lhs, rhs) pairs in LRP.
        """
        mdl, globCompSels, volSys, surfSys = self._parent._usedObjects

        for i, reac in enumerate(self._LRP):
            l, r, rm = reac
            # Add reac elems
            for loc in ALL_LOCATIONS:
                self._reacElems['L'][loc] |= l._getAllElems(loc)
                self._reacElems['R'][loc] |= r._getAllElems(loc)
            # Try to declare complex reac, if applicable
            for j, compEvents in enumerate(self._compEvs):
                # Add complex reac elems
                for ce in compEvents:
                    for loc in ce.getLocations():
                        for side in ce.affectedSides():
                            self._reacElems[side][loc].add(ce.comp)

                # Declare reaction
                name = self._reacName + (f'_{i}' if len(self._LRP) > 1 else '')
                if len(compEvents) > 0:
                    name_c = name + (f'_{j}' if len(self._compEvs) > 1 else '')
                    res = self._declareComplexReac(
                        name_c, l, r, compEvents, mdl, volSys, surfSys, rm * self.K, rm,
                        self.Dependencies, self.AntiDependencies, self.Immobilization, self.MaxDistance
                    )
                else:
                    res = self._declareSimpleReac(
                        name, l, r, mdl, volSys, surfSys, rm * self.K, rm, self.Dependencies,
                        self.AntiDependencies, self.Immobilization, self.MaxDistance
                    )

                # Add the subreactions to self
                self.append(res)

    def _declareComplexReac(self, name, simpLHS, simpRHS, compEvents, mdl, volSys, surfSys, rate, rateMult,
            deps, antideps, immobilization, max_dist):
        """Declare complex reaction involving real complexes."""

        if deps is not None:
            raise Exception(f'Cannot use dependencies with complex reactions.')

        nutils._print(f'\tAdding STEPS complex reaction {name}: {compEvents}, rate = {rate}', 3)

        if self._parent._isSurfaceReac():
            ilhs = simpLHS._GetStepsElems(Location.IN)
            slhs = simpLHS._GetStepsElems(Location.SURF)
            olhs = simpLHS._GetStepsElems(Location.OUT)
            irhs = simpRHS._GetStepsElems(Location.IN)
            srhs = simpRHS._GetStepsElems(Location.SURF)
            orhs = simpRHS._GetStepsElems(Location.OUT)

            if isinstance(rate, VDepRate):
                raise NotImplementedError('Cannot use voltage dependent rates with STEPS complex reactions')

            icompEvs = [ce._stepsObj for ce in compEvents if ce.loc == Location.IN]
            scompEvs = [ce._stepsObj for ce in compEvents if ce.loc == Location.SURF]
            ocompEvs = [ce._stepsObj for ce in compEvents if ce.loc == Location.OUT]

            stepsReac = stepslib._py_ComplexSReac(
                name,
                surfSys.stepsSys,
                ilhs=ilhs,
                slhs=slhs,
                olhs=olhs,
                irhs=irhs,
                srhs=srhs,
                orhs=orhs,
                icompEvs=icompEvs,
                scompEvs=scompEvs,
                ocompEvs=ocompEvs,
                kcst=rate,
            )
        else:
            lhs = simpLHS._GetStepsElems()
            rhs = simpRHS._GetStepsElems()

            if isinstance(rate, CompDepRate):
                raise NotImplementedError(
                    f'{self._decl()}: Cannot use CompDepRate with STEPS complex '
                    f'reactions.'
                )


            compEvs = [ev._stepsObj for ev in compEvents]

            stepsReac = stepslib._py_ComplexReac(name, volSys.stepsSys, lhs=lhs, rhs=rhs, compEvs=compEvs, kcst=rate)

        return _SubReaction(self, stepsReac, rateMult=rateMult, name=name)

    def _declareSimpleReac(self, name, lhs, rhs, mdl, volSys, surfSys, rate, rateMult, deps, antideps,
            immobilization, max_dist):
        """Declare and return the steps reactions from the final sides lhs and rhs."""

        # Compute rates and other parameters in case of complex dependencies
        if isinstance(rate, CompDepRate):
            rate = rate(lhs)
        if isinstance(deps, CompDepFunc):
            deps = deps(lhs)
        if isinstance(immobilization, CompDepFunc):
            immobilization = immobilization(lhs)
        if isinstance(max_dist, CompDepFunc):
            max_dist = max_dist(lhs)

        nutils._print(f'\tAdding STEPS reaction {name}: {lhs} -> {rhs}, rate = {rate}, rateMult = {rateMult}', 3)

        if self._parent._isSurfaceReac():
            # Surface system reaction
            ilhs = lhs._GetStepsElems(Location.IN)
            slhs = lhs._GetStepsElems(Location.SURF)
            olhs = lhs._GetStepsElems(Location.OUT)
            irhs = rhs._GetStepsElems(Location.IN)
            srhs = rhs._GetStepsElems(Location.SURF)
            orhs = rhs._GetStepsElems(Location.OUT)

            if isinstance(surfSys, VesicleSurfaceSystem):
                # Vesicle surface reaction

                if len(ilhs) > 0:
                    raise Exception(
                        f'{self._decl()}: Vesicle reactions cannot involve species inside the vesicle on'
                        f' the LHS.'
                    )

                # On vesicle surface
                _vlhs = lhs._GetStepsElems(Location.VESSURF)
                _vrhs = rhs._GetStepsElems(Location.VESSURF)
                vlhs = [se for se in _vlhs if isinstance(se, stepslib._py_Spec)]
                llhs = [se for se in _vlhs if isinstance(se, stepslib._py_LinkSpec)]
                vrhs = [se for se in _vrhs if isinstance(se, stepslib._py_Spec)]
                lrhs = [se for se in _vrhs if isinstance(se, stepslib._py_LinkSpec)]

                # Check that no forbidden locations are present
                for elem in lhs:
                    if elem.loc not in [Location.IN, Location.OUT, Location.SURF, Location.VESSURF]:
                        raise Exception(
                            f'{self._decl()}: Element {elem} cannot appear in a vesicle surface reaction.'
                        )

                # Dependencies
                vdeps = []
                if deps is not None:
                    deps = deps._toReactionSide()
                    if any(elem.loc != Location.VESSURF for elem in deps):
                        raise Exception(
                            f'{self._decl()}: Dependencies of reactions declared in a '
                            f'VesicleSurfaceSystem need to be explicitely declared on the surface of the '
                            f'vesicle (e.g. "spec.v"). The dependencies given do not satisfy this '
                            f'constraint: {deps}.'
                        )

                    vdeps = deps._GetStepsElems(Location.VESSURF)

                # Immobilization
                if immobilization not in [IMMOBILIZING, NO_EFFECT, MOBILIZING]:
                    raise ValueError(
                        f'{self._decl()}: Unsupported immobilization flag. IMMOBILIZING, NO_EFFECT'
                        f', or MOBILIZING.')

                # Max distance
                if max_dist is None:
                    max_dist = -1
                elif max_dist < 0:
                    raise ValueError(
                        f'{self._decl()}: Maximum distance parameter should be a positive number, got '
                        f'{max_dist} instead.'
                    )

                stepsReac = stepslib._py_VesSReac(
                    name,
                    surfSys._getStepsObjects()[0],
                    olhs=olhs,
                    slhs=slhs,
                    vlhs=vlhs,
                    llhs=llhs,
                    lrhs=lrhs,
                    vrhs=vrhs,
                    srhs=srhs,
                    orhs=orhs,
                    irhs=irhs,
                    vdeps=vdeps,
                    kcst=rate,
                    immobilization=immobilization,
                    max_dist=max_dist,
                )
            elif isinstance(surfSys, RaftSurfaceSystem):
                # Raft surface reaction

                # On vesicle surface
                rlhs = lhs._GetStepsElems(Location.RAFTSURF)
                rrhs = rhs._GetStepsElems(Location.RAFTSURF)

                # Check that no forbidden locations are present
                for elem in lhs:
                    if elem.loc not in [Location.IN, Location.OUT, Location.SURF, Location.RAFTSURF]:
                        raise Exception(
                            f'{self._decl()}: Element {elem} cannot appear in a raft surface reaction.'
                        )

                # Dependencies
                rdeps = []
                if deps is not None:
                    deps = deps._toReactionSide()
                    if any(elem.loc != Location.RAFTSURF for elem in deps):
                        raise Exception(
                            f'{self._decl()}: Dependencies of reactions declared in a '
                            f'RaftSurfaceSystem need to be explicitely declared on the surface of the '
                            f'raft (e.g. "spec.r"). The dependencies given do not satisfy this '
                            f'constraint: {deps}.'
                        )

                    rdeps = deps._GetStepsElems(Location.RAFTSURF)

                # Anti Dependencies
                anti_rdeps = []
                if antideps is not None:
                    antideps = antideps._toReactionSide()
                    if any(elem.loc != Location.RAFTSURF for elem in antideps):
                        raise Exception(
                            f'{self._decl()}: Anti-dependencies of reactions declared in a '
                            f'RaftSurfaceSystem need to be explicitely declared on the surface of the '
                            f'raft (e.g. "spec.r"). The anti-dependencies given do not satisfy this '
                            f'constraint: {antideps}.'
                        )

                    anti_rdeps = antideps._GetStepsElems(Location.RAFTSURF)

                # Immobilization
                if immobilization not in [IMMOBILIZING, NO_EFFECT, MOBILIZING]:
                    raise ValueError(
                        f'{self._decl()}: Unsupported immobilization flag. IMMOBILIZING, NO_EFFECT'
                        f', or MOBILIZING.')

                stepsReac = stepslib._py_RaftSReac(
                    name,
                    surfSys._getStepsObjects()[0],
                    ilhs=ilhs,
                    olhs=olhs,
                    slhs=slhs,
                    rlhs=rlhs,
                    rrhs=rrhs,
                    srhs=srhs,
                    orhs=orhs,
                    irhs=irhs,
                    rdeps=rdeps,
                    anti_rdeps=anti_rdeps,
                    kcst=rate,
                    immobilization=immobilization,
                )
            else:
                # Normal surface reaction

                if isinstance(rate, VDepRate):
                    # Voltage dependent surface reaction
                    stepsReac = stepslib._py_VDepSReac(
                        name,
                        surfSys.stepsSys,
                        ilhs=ilhs,
                        slhs=slhs,
                        olhs=olhs,
                        irhs=irhs,
                        srhs=srhs,
                        orhs=orhs,
                        **rate._getVDepSReacParams(self._getRateUnits())
                    )
                else:
                    stepsReac = stepslib._py_SReac(
                        name,
                        surfSys.stepsSys,
                        ilhs=ilhs,
                        slhs=slhs,
                        olhs=olhs,
                        irhs=irhs,
                        srhs=srhs,
                        orhs=orhs,
                        kcst=rate,
                    )
        else:
            # Volume system reaction
            elhs = lhs._GetStepsElems()
            erhs = rhs._GetStepsElems()
            if isinstance(rate, VDepRate):
                # Undeclare the parent reaction to avoid further exception raising when exiting context
                # managers
                # self._parent._declared = False
                raise NotImplementedError(
                    f'{self._decl()}: Cannot declare a voltage-dependent volume '
                    f'reaction.'
                )
            else:
                stepsReac = stepslib._py_Reac(name, volSys.stepsSys, lhs=elhs, rhs=erhs, kcst=rate)

        # Check that there were no silently unused Parameters
        if not isinstance(stepsReac, (stepslib._py_VesSReac, stepslib._py_RaftSReac)):
            if deps is not None:
                raise Exception(
                    f'{self._decl()}: Dependencies can only be used with vesicle and raft surface reactions.'
                )
            if immobilization != NO_EFFECT:
                raise Exception(
                    f'{self._decl()}: Immobilization flags can only be used with vesicle and raft surface '
                    f'reactions.'
                )
        if not isinstance(stepsReac, stepslib._py_VesSReac):
            if max_dist is not None:
                raise Exception(
                    f'{self._decl()}: Maximum distances can only be used with vesicle surface reactions.'
                )

        return _SubReaction(self, stepsReac, lhs=lhs, rateMult=rateMult, name=name)


class _ComplexEvent:
    """Utility class for storing complex events"""
    def __init__(self, comp, loc):
        self.comp = comp
        self.loc = loc
        self._stepsObj = self._createStepsObj()

    def getLocations(self):
        return [self.loc]


class _ComplexCreateEvent(_ComplexEvent):
    def __init__(self, comp, loc, init):
        self.init = init
        super().__init__(comp, loc)

    def affectedSides(self):
        return 'R'

    def _createStepsObj(self):
        return stepslib._py_ComplexCreateEvent(self.comp.name, self.init)

    def __repr__(self):
        return f'DeleteEv(comp={self.comp}, loc={self.loc}, init={self.init})'


class _ComplexLHSEvent(_ComplexEvent):
    def __init__(self, comp, loc, filt):
        self.filt = filt
        super().__init__(comp, loc)


class _ComplexDeleteEvent(_ComplexLHSEvent):
    def affectedSides(self):
        return 'L'

    def _createStepsObj(self):
        return stepslib._py_ComplexDeleteEvent(self.comp.name, self.filt)

    def __repr__(self):
        return f'DeleteEv(comp={self.comp}, loc={self.loc}, filt={self.filt})'


class _ComplexUpdateEvent(_ComplexLHSEvent):
    def __init__(self, comp, loc, filt, reac, upd, destLoc=None):
        self.reac = reac
        self.upd = upd
        self.destLoc = destLoc
        super().__init__(comp, loc, filt)

    def affectedSides(self):
        return 'LR'

    def getLocations(self):
        return [self.loc, self.destLoc]

    def _createStepsObj(self):
        return stepslib._py_ComplexUpdateEvent(self.comp.name, self.filt, self.reac, self.upd, self.destLoc)

    def __repr__(self):
        return f'UpdateEv(comp={self.comp}, loc={self.loc}, filt={self.filt}, reac={self.reac}, '\
               f'upd={self.upd}, dest={self.destLoc})'


@nutils.FreezeAfterInit
class Reaction(
    nutils.UsingObjects((Model, ComplexSelector), (VolumeSystem, SurfaceSystem)),
    nutils.StepsWrapperObject,
    nutils.ParameterizedObject,
):
    """A (possibly reversible) reaction.

    A reaction needs to be declared within a :py:class:`Model` and within a
    :py:class:`VolumeSystem` or a :py:class:`SurfaceSystem` (i.e. inside a ``with`` block,
    see e.g. :py:class:`VolumeSystem`)

    A :py:class:`Reaction` is defined by a left hand side (lhs) and a right hand side (rhs).
    Both sides are instances of :py:class:`ReactionSide` and the overall reaction is declared
    with a custom syntax designed to look like standard chemical reactions notation.
    When :py:class:`Complex` are involved, a :py:class:`Reaction` object can represent several
    elementary reactions, see user guide for details.

    .. note::
        Reactions should always be instantiated from a :py:class:`ReactionManager`.

    Usage::

        r = ReactionManager()
        with vsys:
            # Simple forward reaction
            SA >r[1]> SB
            r[1].K = 1      # Set forward rate

            # Bidirectional reaction
            SA + SB <r[1]> SC + SD
            r[1].K = 2, 1   # Set forward and backward rates

            # Chained reactions
            SA <r[1]> SB >r[2]> SC
            r[1].K = 2, 1
            r[2].K = 3

    .. warning::
        The corresponding reactions are only truly added to steps after the `with vsys:` block
        is exited.

    """

    _FwdSpecifier = 'fwd'
    _BkwSpecifier = 'bkw'

    def __init__(self, _managerInfo=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.lhs = None
        self.rhs = None

        self._subReactions = collections.OrderedDict()
        self._subReactions[Reaction._FwdSpecifier] = None
        self._subReactions[Reaction._BkwSpecifier] = None

        self._declared = False
        self._added = False
        self._bidir = None

        self._declarationInfo = None

        self._managerInfo = _managerInfo
        self._usedObjects = list(self._getUsedObjects())

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        reac = cls(name=obj.getID())

        # Reaction sides
        lhs, rhs = ReactionSide([]), ReactionSide([])
        if isinstance(obj, stepslib._py_Reac):
            for lhsSpec in obj.getLHS():
                lhs += getattr(mdl, lhsSpec.getID())
            for rhsSpec in obj.getRHS():
                rhs += getattr(mdl, rhsSpec.getID())
        elif isinstance(obj, (
                stepslib._py_SReac, stepslib._py_VDepSReac, 
                stepslib._py_VesSReac, stepslib._py_RaftSReac)):
            if not isinstance(obj, stepslib._py_VesSReac):
                for lhsSpec in obj.getILHS():
                    lhs += getattr(mdl, lhsSpec.getID()).i
            for lhsSpec in obj.getSLHS():
                lhs += getattr(mdl, lhsSpec.getID()).s
            for lhsSpec in obj.getOLHS():
                lhs += getattr(mdl, lhsSpec.getID()).o
            for rhsSpec in obj.getIRHS():
                rhs += getattr(mdl, rhsSpec.getID()).i
            for rhsSpec in obj.getSRHS():
                rhs += getattr(mdl, rhsSpec.getID()).s
            for rhsSpec in obj.getORHS():
                rhs += getattr(mdl, rhsSpec.getID()).o
            if isinstance(obj, stepslib._py_VesSReac):
                for lhsSpec in obj.getVLHS() + obj.getLLHS():
                    lhs += getattr(mdl, lhsSpec.getID()).v
                for rhsSpec in obj.getVRHS() + obj.getLRHS():
                    rhs += getattr(mdl, rhsSpec.getID()).v
            if isinstance(obj, stepslib._py_RaftSReac):
                for lhsSpec in obj.getRLHS():
                    lhs += getattr(mdl, lhsSpec.getID()).r
                for rhsSpec in obj.getRRHS():
                    rhs += getattr(mdl, rhsSpec.getID()).r
        else:
            raise NotImplementedError(f'Cannot import from STEPS reaction {obj}')

        lhs > reac > rhs

        subreacLst = _SubReactionList(reac, obj.getID(), [(lhs, rhs, 1)], None)
        subreacLst.append(_SubReaction(subreacLst, obj, lhs=lhs, rateMult=1.0, name=obj.getID()))
        reac._subReactions[Reaction._FwdSpecifier] = subreacLst

        # Rate
        if not isinstance(obj, stepslib._py_VDepSReac):
            # TODO not urgent: set K for VDepSReacs
            reac.K = obj.getKcst()

        # Immobilization
        if isinstance(obj, (stepslib._py_VesSReac, stepslib._py_RaftSReac)):
            reac.Immobilization = obj.getImmobilization()

        # MaxDistance
        if isinstance(obj, stepslib._py_VesSReac) and obj.getMaxDistance() != -1.0:
            reac.MaxDistance = obj.getMaxDistance()

        # Dependencies
        if isinstance(obj, (stepslib._py_VesSReac, stepslib._py_RaftSReac)):
            deps = ReactionSide([])
            if isinstance(obj, stepslib._py_VesSReac):
                for depSpec in obj.getVDeps():
                    deps += getattr(mdl, depSpec.getID()).v
            else:
                for depSpec in obj.getRDeps():
                    deps += getattr(mdl, depSpec.getID()).r
            reac.Dependencies = deps

        # Anti Dependencies
        if isinstance(obj, (stepslib._py_RaftSReac, )):
            antideps = ReactionSide([])
            for adepSpec in obj.getAntiRDeps():
                deps += getattr(mdl, adepSpec.getID()).r
            reac.AntiDependencies = deps

        reac._added = True
        for loc in ALL_LOCATIONS:
            subreacLst._reacElems['L'][loc] |= lhs._getAllElems(loc)
            subreacLst._reacElems['R'][loc] |= rhs._getAllElems(loc)

        return reac

    @staticmethod
    def _findMatchingComplexInd(csl, rhs):
        """Return the index of the complexSelector corresponding to csl in rhs."""
        if csl.loc is None:
            return rhs.index(csl) if csl in rhs else None
        else:
            # Handle the case in which complex selectors have different location
            csrInds = [i for i, csr in enumerate(rhs) if (csl._elem, csl.stoich) == (csr._elem, csr.stoich)]
            return csrInds[0] if len(csrInds) > 0 else None

    def _expandSimpleReactions(self, lhs, rhs):
        """
        Expand ComplexSelectors whose complex has '_statesAsSpecies = True' to a list of
        lhs and rhs pairs.
        """
        cre = _ComplexReactionElement
        simpleCompLhs = [cs for cs in lhs._getElemsOfType(cre) if cs._elem._areStatesAsSpecies()]
        simpleCompRhs = [cs for cs in rhs._getElemsOfType(cre) if cs._elem._areStatesAsSpecies()]
        specLhs = lhs._getElemsOfType(Species, LinkSpecies, ComplexState)
        specRhs = rhs._getElemsOfType(
            Species, LinkSpecies, ComplexState
        )
        # TODO Later release: Using complexState in reaction when using C++ complex reacs: is it working?

        # Initialize the list of reaction pairs with species
        lhsRhsPairs = [(ReactionSide(specLhs), ReactionSide(specRhs), 1)]

        # Group complexReactionElements by complex
        comp2csls = {}
        for csl in simpleCompLhs:
            cmplx = csl._elem._getReferenceObject()
            comp2csls.setdefault(cmplx, []).append(csl)

        # Compute all reaction pairs for lhs ComplexSelectors
        for comp, csls in comp2csls.items():
            if len(csls) > 2:
                raise Exception(
                    f'Cannot use more than 2 instances of a complex in a reaction, '
                    f'the current reaction uses {comp} {len(csls)} times.'
                )

            allStatePairs = []
            for csl in csls:
                csrInd = Reaction._findMatchingComplexInd(csl, simpleCompRhs)
                if csrInd is not None:
                    # Compute state pairs if there is a matching ComplexSelector in the rhs
                    csr = simpleCompRhs[csrInd]
                    statePairs = []
                    for nl, nr, rm in csl._elem._getMatchingStates(csr._elem):
                        statePairs.append(
                            (
                                ReactingElement(nl, csl.stoich, csl.loc),
                                ReactingElement(nr, csr.stoich, csr.loc),
                                rm,
                            )
                        )
                    del simpleCompRhs[csrInd]
                else:
                    if isinstance(csl, _ComplexReactants):
                        raise Exception(
                            f'ComplexReactants {csl} has no matching ComplexReactants on '
                            f'the RHS of the reaction.'
                        )
                    # Otherwise, just add all the states of csl
                    allStates = csl._elem._getAllStates()
                    statePairs = [
                        (ReactingElement(nl, csl.stoich, csl.loc), ReactionSide([]), 1) for nl in allStates
                    ]
                allStatePairs.append(statePairs)

            # If 2 instances of the same complex are used, update the rate multipliers accordingly
            newASP = []
            if len(csls) == 2:
                sameReacs = (
                    all(isinstance(csl._elem, _ComplexReactants) for csl in csls)
                    and csls[0]._elem._hasSameReactants(csls[1]._elem)
                )
                allL1 = {l1: rm1 for l1, r1, rm1 in allStatePairs[0]}
                allL2 = {l2: rm2 for l2, r2, rm2 in allStatePairs[1]}
                for l1, r1, rm1 in allStatePairs[0]:
                    for l2, r2, rm2 in allStatePairs[1]:
                        if sameReacs and l1 in allL2 and l2 in allL1:
                            newASP.append((l1 + l2, r1 + r2, rm1 * rm2 / 2))
                        else:
                            newASP.append((l1 + l2, r1 + r2, rm1 * rm2))
            else:
                newASP = allStatePairs[0]

            # Update lhsRhsPairs with the matching states
            newLRP = []
            for l, r, rm in lhsRhsPairs:
                for nl, nr, nrm in newASP:
                    newLRP.append((l + nl, r + nr, rm * nrm))
            lhsRhsPairs = newLRP

        # Add reaction pairs from remaining rhs elements if they are completely defined
        # Otherwise, raise an exception.
        if len(simpleCompRhs) > 0:
            newLRP = []
            for l, r, rm in lhsRhsPairs:
                for csr in simpleCompRhs:
                    states = csr._elem._getAllStates()
                    if len(states) > 1:
                        raise Exception(
                            f'Complex selector {csr._elem} is used in the right hand side of a '
                            f'reaction but is not matching anything in the left hand side and is '
                            f'not fully defined. The reaction is ambiguous.'
                        )
                    newLRP.append((l, r + ReactingElement(states.pop(), csr.stoich, csr.loc), rm))
            lhsRhsPairs = newLRP

        return lhsRhsPairs

    def _expandComplexReactions(self, lhs, rhs):
        """
        Expand ComplexSelectors whose complex has '_statesAsSpecies = False' to a list of
        complex events.
        """
        cre = _ComplexReactionElement
        clhs = [cs for cs in lhs._getElemsOfType(cre) if not cs._elem._areStatesAsSpecies()]
        crhs = [cs for cs in rhs._getElemsOfType(cre) if not cs._elem._areStatesAsSpecies()]

        # Initialize the list of complex reaction events
        compEvents = [[]]

        for cl in clhs:
            comp = cl._elem._complex
            crInd = Reaction._findMatchingComplexInd(cl, crhs)
            if crInd is not None:
                cr = crhs[crInd]._elem
                events = [
                    _ComplexUpdateEvent(comp, cl.loc, filt, reac, upd, destLoc=crhs[crInd].loc)
                    for filt, reac, upd in cl._elem._getUpdateVects(cr)
                ]
                del crhs[crInd]
            else:
                if isinstance(cl, _ComplexReactants):
                    raise Exception(
                        f'ComplexReactants {cl} has no matching ComplexReactants on '
                        f'the RHS of the reaction.'
                    )

                events = [_ComplexDeleteEvent(comp, cl.loc, cl._elem._toUnorderedFilter())]

            # Update compEvents with the new events
            newCE = []
            for evs in compEvents:
                for e in events:
                    newCE.append(evs + [e])
            compEvents = newCE

        # Add creation events
        for cr in crhs:
            state = cr._elem._getFullyDeterminedState()
            if state is None:
                raise Exception(
                    f'Complex selector {cr._elem} is used on the RHS of the reaction and '
                    f'should lead to a new complex creation but is not fully defined. '
                    f'The reaction is ambiguous.'
                )
            for evs in compEvents:
                evs.append(_ComplexCreateEvent(state._comp, cr.loc, state._toUnorderedVect()))

        return compEvents

    def _parseReac(self):
        """
        Parse the reaction, expand complex reactions and setup data structures for declaring the
        reactions in steps.
        """
        mdl, globCompSels, volSys, surfSys = self._usedObjects

        # Check that no SubUnitStates are wrongly used
        if (self.lhs._isEmpty() and self.rhs._hasElemsOfType(SubUnitSelector)) or (
            self.rhs._isEmpty() and self.lhs._hasElemsOfType(SubUnitSelector)
        ):
            raise Exception('Cannot mix SubUnitState or SubUnitSelectors with empty reaction side.')

        # Check that the correct spatial systems are being used
        if self._isSurfaceReac() and surfSys is None:
            raise Exception(
                f'Reaction {self.name} is a surface reaction but is being declared '
                f'outside of a surface system.'
            )
        if not self._isSurfaceReac() and volSys is None:
            raise Exception(
                f'Reaction {self.name} is a volume reaction but is being declared '
                f'outside of a volume system.'
            )

        if not self._isSurfaceReac():
            # Check that if reacElems have a location, it's the same for all reacElems and both sides
            if self.lhs._hasLocation() or self.rhs._hasLocation():
                if self.lhs._hasSingleLocation() and self.rhs._hasSingleLocation():
                    if self.lhs._hasVolumeLocation() and self.rhs._hasVolumeLocation():
                        # Add a constraint to the model that will be checked later,
                        # to make sure that the volume system matches the location of the surface system
                        loc = self.lhs._getLocation() if self.lhs._hasLocation() else self.rhs._getLocation()
                        mdl._addVolSysConstraint(self, volSys, surfSys, loc)
                    else:
                        raise Exception(
                            f'Reaction {self.name} is a volume reaction but has a side with '
                            f'surface located component.'
                        )
                else:
                    raise Exception(
                        f'Reaction {self.name} is a volume reaction but has component in different locations.'
                    )


        # Transform SubUnitSelectors to ComplexSelectors or ComplexReactants
        lhs = self.lhs._getFinal(globCompSels)
        rhs = self.rhs._getFinal(globCompSels)

        # Transform ComplexSelectors whose complex has "statesAsSpecies = True" to species
        fwdLRP = self._expandSimpleReactions(lhs, rhs)
        fwdCompEvs = self._expandComplexReactions(lhs, rhs)
        nameFwd = self.name + Reaction._FwdSpecifier if self._bidir else self.name

        self._subReactions[Reaction._FwdSpecifier] = _SubReactionList(self, nameFwd, fwdLRP, fwdCompEvs, True)

        if self._bidir:
            bkwLRP = self._expandSimpleReactions(rhs, lhs)
            bkwCompEvs = self._expandComplexReactions(rhs, lhs)
            nameBkw = self.name + Reaction._BkwSpecifier

            self._subReactions[Reaction._BkwSpecifier] = _SubReactionList(self, nameBkw, bkwLRP, bkwCompEvs, False)

        # Retrieve filename and line number for easier debugging
        for i, fi in enumerate(inspect.stack()):
            if fi.function in ['__lt__', '__gt__']:
                i += 1
                break
        frameInfo =  inspect.stack()[i]
        self._declarationInfo = (frameInfo.filename, frameInfo.lineno)

        self._declared = True

    def _comparison(self, other, gt):
        if self._declared:
            if self._managerInfo is None:
                raise Exception('This reaction has already been declared.')
            # Return the new reaction
            # TODO not urgent: Manager info being a tuple hurts readability, make it a namedtuple
            p = self._managerInfo[0]._popNewReacArgs()
            if p is None:
                self._managerInfo[0]._setNewReacArgs(other, gt)
                # python exhibits special behavior for chained comparison operator
                # see https://docs.python.org/3/reference/expressions.html#comparisons
                # The first comparison needs to return True
                return True
            else:
                r = self._managerInfo[0]._getNew(self._managerInfo[1])
                r._comparison(*p)
                return r._comparison(other, gt)

        if isinstance(other, ReactionElement):
            other = other._toReactionSide()
        elif other is None:
            other = ReactionSide([])
        elif not isinstance(other, ReactionSide):
            raise Exception(f'{other} cannot be one of the sides of reaction {self}.')

        if self.lhs is None:
            self.lhs = other._makeCopy()
            self._bidir = gt
            return self
        else:
            if not gt:
                raise Exception('Used "<" instead of ">" after the reaction.')
            self.rhs = other._makeCopy()
            self._parseReac()
            return other

    def __gt__(self, other):
        return self._comparison(other, True)

    def __lt__(self, other):
        return self._comparison(other, False)

    def _getAllElems(self, loc, sides='LR'):
        if self._subReactions[Reaction._FwdSpecifier] is not None:
            res = self._subReactions[Reaction._FwdSpecifier]._getAllElems(loc, sides)
            if self._subReactions[Reaction._BkwSpecifier] is not None:
                res |= self._subReactions[Reaction._BkwSpecifier]._getAllElems(loc, sides)
            return res
        return set()

    def _isSurfaceReac(self):
        if self.lhs is None or self.rhs is None:
            return None
        lhsS, rhsS = self.lhs._isSurfaceSide(), self.rhs._isSurfaceSide()
        return (lhsS and rhsS) or (None in [lhsS, rhsS] and True in [lhsS, rhsS])

    def _isVesicleSurfReac(self):
        return isinstance(self._usedObjects[3], VesicleSurfaceSystem)

    def _isRaftSurfReac(self):
        return isinstance(self._usedObjects[3], RaftSurfaceSystem)

    def __getitem__(self, key):
        """Get sub-reactions from a reaction

        Square bracket notation can be used on a reaction to retrieve either the forward
        reaction(s), or the backward one(s). This should never be needed during model declaration
        but it is very important for simulation control and data saving (see
        :py:class:`steps.API_2.sim.SimPath`).

        :param key: Use ``'{{model.Reaction._FwdSpecifier}}'`` to access the forward
            reaction(s) or ``'{{model.Reaction._BkwSpecifier}}'`` to access the backward
            reaction(s).
        :type key: str

        :returns: The querried reaction(s).

        See :ref:`/API_2/STEPS_Tutorial_wm.ipynb` for usage example.

        :meta public:
        """
        if key not in self._subReactions:
            raise KeyError(f'No reaction can be accessed with key {key}.')
        return self._subReactions[key]

    def __repr__(self):
        return '{} {}-> {}'.format(self.lhs, '-' if not self._bidir else '<', self.rhs)

    def _decl(self):
        """
        Return a string that describes the declaration of the object. It returns its name by default
        but it can contain e.g. the file and line at which it was declared.
        """
        return f'{self.name} ({self._declarationInfo[0]}: {self._declarationInfo[1]})'

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        res = []
        for specif, srl in self._subReactions.items():
            if srl is not None:
                res += srl._getStepsObjects()
        return res

    def _simPathWalkExpand(self):
        """Return an iterable of the elements that should be part of ResultSelector paths."""
        return itertools.chain(
            *[srl._simPathWalkExpand() for specif, srl in self._subReactions.items() if srl is not None]
        )

    def _getSubParameterizedObjects(self):
        """Return all subobjects that can hold parameters."""
        return [srl for specif, srl in self._subReactions.items() if srl is not None]

    def _includeinParamTables(self):
        """Whether the object should be included in parameter tables"""
        # Only directed reactions should be included in parameter tables
        return False

    def _exitCallback(self, parent):
        """
        Method to be called the first time we get out of a context manager in which the object as been
        declared.
        """
        # Only declare the reaction if we are exiting the volume / surface system
        if isinstance(parent, SpaceSystem):
            if self._declared and not self._added:
                nutils._print(f'Adding pysteps reaction {self.name}: {self}', 2)

                self._added = True
                self._subReactions[Reaction._FwdSpecifier]._addToSteps()
                if self._bidir:
                    self._subReactions[Reaction._BkwSpecifier]._addToSteps()

    @staticmethod
    def _getProperty(prop):
        """Build properties for the Reaction class based on the property prop."""

        def getter(self):
            if not self._declared:
                raise Exception(f'Cannot get properties of a reaction that was not yet declared.')
            if self._bidir:
                return (
                    prop.fget(self._subReactions[Reaction._FwdSpecifier]),
                    prop.fget(self._subReactions[Reaction._BkwSpecifier]),
                )
            else:
                return prop.fget(self._subReactions[Reaction._FwdSpecifier])

        def setter(self, val):
            if not self._declared:
                raise Exception('Cannot set properties of a reaction that was not yet declared.')
            if self._added:
                raise Exception('Cannot set properties of a reaction that was already added to STEPS.')

            if self._bidir:
                if not (isinstance(val, tuple) and len(val) == 2):
                    raise Exception(
                        f'The reaction is reversible, two properties should be set, got {val} instead.'
                    )
                prop.fset(self._subReactions[Reaction._FwdSpecifier], val[0])
                prop.fset(self._subReactions[Reaction._BkwSpecifier], val[1])
            else:
                if not isinstance(val, tuple):
                    val = (val, )
                if len(val) > 1:
                    raise Exception(f'The reaction is unidirectional but more than one property was set')
                prop.fset(self._subReactions[Reaction._FwdSpecifier], val[0])

        return property(getter, setter)


# Map _SubReactionList properties to Reaction
for propName in dir(_SubReactionList):
    prop = getattr(_SubReactionList, propName)
    if isinstance(prop, property):
        newProp = Reaction._getProperty(prop)
        newProp.__doc__ = prop.__doc__
        setattr(Reaction, propName, newProp)


@nutils.FreezeAfterInit
class AnonymousReaction(Reaction):
    """
    :meta private:
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


@nutils.FreezeAfterInit
class ReactionManager:
    """Class used to instantiate reactions

    A :py:class:`ReactionManager` object allows to easily declare reactions by associating them
    to an identifier using the square bracket syntax::

        r = ReactionManager()  # Only one ReactionManager is needed

        S1 + S2 >r['R01']> S3  # String identifier for permanently naming the reaction

        S1 + S2 >r[1]> S3      # Integer identifier for reactions that do not need explicit
                               # naming

    String identifiers are used when the reaction will need to be accessed during simulation (see
    :py:class:`steps.API_2.sim.SimPath` for usage), they cannot be reused for different reactions.
    Integer identifiers are temporary identifiers and can be reused for declaring several un-named
    reactions::

        S1 + S2 >r[1]> S3    # First use of temporary identifier 1. r[1] now refers to this
        r[1].K = ...         # reaction until a new reaction is defined with r[1]. 

        S1 + S3 <r[1]> S4    # r[1] now refers to this new reaction.
        r[1].K = ...

    In general, only one :py:class:`ReactionManager` should be instantiated and given a short
    variable name (``r`` in all examples) for convenience. In some specific cases however, one might
    need several independent reaction managers so this class is not implemented as a singleton.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reactions = {}
        self._newReacArgs = None

    def __getitem__(self, i):
        """Square brackets syntax for instantiating or accessing a reaction

        :param i: Identifier of the reaction. Integer for a temporary id and string for a
            permanent one.
        :type i: Union[int, str]

        :returns: The reaction corresponding to the given identifier
        :rtype: :py:class:`Reaction`

        :meta public:
        """
        if i not in self.reactions:
            if isinstance(i, str):
                self.reactions[i] = Reaction(name=i, _managerInfo=(self, i))
            else:
                self.reactions[i] = AnonymousReaction(_managerInfo=(self, i))
        return self.reactions[i]

    def _getNew(self, i):
        del self.reactions[i]
        return self.__getitem__(i)

    def _setNewReacArgs(self, *args):
        self._newReacArgs = args

    def _popNewReacArgs(self):
        p = self._newReacArgs
        self._newReacArgs = None
        return p


###################################################################################################
# State-dependent rates


@nutils.FreezeAfterInit
class XDepFunc(nutils.NamedObject, nutils.ParameterizedObject):
    r"""Base class for functions depending on some simulation or object state

    This is a wrapper around a function that adds some functionalities.

    :param func: lambda function or standard function. There is not support for variable named
        arguments (\*\*kwargs).
    :type func: Callable[[Any, ...], Number]
    """

    _DEPENDS_ON_NAME = 'Depends on'

    def __init__(self, func, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not hasattr(func, '__call__'):
            raise TypeError(f'Expected a callable, got {func} instead.')
        self._func = func
        self._args = []
        self._mult = 1

        self._setParameter('Function', self._func)
        self._subParamObjs = set()

    def _getSubParameterizedObjects(self):
        """Return all subobjects that can hold parameters."""
        return sorted(self._subParamObjs, key=lambda x: x.name)

    def __mul__(self, m):
        """Multiplication with a constant or another function of the same class

        :param m: Constant number or another :py:class:`XDepFunc` of the same type
        :type m: Union[Number]

        :returns: A new function, which corresponds to the multiplication of both operands
        :rtype: :py:class:`XDepFunc`
        :meta public:
        """
        if isinstance(m, numbers.Number):
            ret = copy.copy(self)
            ret._mult *= m
            return ret
        else:
            raise TypeError(f'Cannot multiply {self} with {m}.')

    def __rmul__(self, m):
        """See :py:func:`__mul__`

        :meta public:
        """
        return self.__mul__(m)

    # TODO Not urgent: add more operations

    @property
    def func(self):
        """Access the function itself

        :type: Callable, read-only
        """
        return self._func

    def __call__(self, *lst, _recParams=True):
        """Call the wrapped function

        :meta public:
        """
        if _recParams:
            nutils.Parameter._startRecording()

        ret = self._mult * self._func(*lst)

        if _recParams:
            # Retrieve the parameters used when the function was called
            allParams = [param for param in nutils.Parameter._endRecordingAndGetUsage() if param._isNamed() and param._isUserDefined()]
            if len(allParams) > 0:
                depon = self._getParameter(XDepFunc._DEPENDS_ON_NAME)
                names = set([param.name for param in allParams])
                if depon is not None:
                    names |= set(depon.value.lst)
                    allParams = list(set(depon._dependencies) | set(allParams))
                nameList = nutils._ParamList(list(names))
                self._setParameter(XDepFunc._DEPENDS_ON_NAME, nutils.Parameter(nameList, _dependencies=allParams))

        if isinstance(ret, XDepFunc):
            self._subParamObjs.add(ret)

        return ret


@nutils.FreezeAfterInit
class VDepFunc(XDepFunc):
    """Hold a function and a voltage range for computing e.g. voltage-dependent reaction rates

    For efficiency reasons, the values of these functions are pre-computed for specific voltage
    values (by default between -150mV and 100mV with a 0.1mV step) and these values might be
    changed by providing the optional ``vrange`` parameter.

    :param func: lambda function or standard function that takes a membrane potential (in Volts)
        as parameter and returns a single numerical value.
    :type func: Callable[[float], Number]
    :param vrange: Voltage range at which the function should be computed. Should be a 3-tuple
        in the form ``(min_voltage, max_voltage, step_voltage)``, all expressed in Volts.
    :type vrange: Tuple[float, float, float]

    See :py:class:`VDepRate` for an example.
    """
    _DEFAULT_VMIN = -150e-3
    _DEFAULT_VMAX = 100e-3
    _DEFAULT_DV   = 1e-4

    def __init__(self, func, vrange=None, *args, **kwargs):
        super().__init__(func=func, *args, **kwargs)
        if vrange is None:
            vrange = [VDepFunc._DEFAULT_VMIN, VDepFunc._DEFAULT_VMAX, VDepFunc._DEFAULT_DV]
        if vrange[0] >= vrange[1]:
            raise Exception(f'The maximum voltage of the range is lower than the minimum.')
        if vrange[2] > vrange[1] - vrange[0]:
            raise Exception(f'The voltage step is greater than the voltage range.')

        self.vrange = vrange

        self._args = [self.vrange]

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('V'))
    def vrange(self):
        """The associated voltage range

        :type: Tuple[float, float, float]
        """
        pass

    @vrange.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('V'))
    def vrange(self, val):
        pass


@nutils.FreezeAfterInit
class CompDepFunc(XDepFunc):
    """Hold a function that depends on the fully specified state of a complex

    :param func: lambda function or standard function that takes one or several fully specified
        state(s) of complex(es) as parameter and returns a single numerical value.
    :type func: Callable[[:py:class:`ComplexState`, ...], Number]
    :param compsOrCompSels: A list of complexes or complex selectors in the same order as they
        should be provided to the wrapped function.
    :type compsOrCompSels: List[Union[:py:class:`Complex`, :py:class:`ComplexSelector`]]

    See :py:class:`CompDepRate` for an example.
    """

    def __init__(self, func, compsOrCompSels, *args, **kwargs):
        super().__init__(func=func, *args, **kwargs)
        if not hasattr(compsOrCompSels, '__iter__'):
            raise TypeError(f'Expected an iterable containing complexes, got {compsOrCompSels} instead.')
        self._complexes = []
        for cs in compsOrCompSels:
            if isinstance(cs, Complex):
                self._complexes.append(cs)
            elif isinstance(cs, ComplexSelector):
                self._complexes.append(cs._complex)
            elif isinstance(cs, ComplexState):
                self._complexes.append(cs._comp)
            else:
                raise TypeError(f'Expected a complex or a complex selector, got {cs} instead.')

        self._compsOrCompSels = compsOrCompSels

        self._args = [self._compsOrCompSels]

        self._setParameter('Complex', self._compsOrCompSels)

    def _getStatesFromReacSide(self, lhs):
        lhsCompStates = list(lhs._getElemsOfType(ComplexState))
        lhsCompSels = list(lhs._getElemsOfType(ComplexSelector))
        states = []
        for cs in self._compsOrCompSels:
            candidates = []
            if isinstance(cs, ComplexSelector):
                for re in lhsCompStates:
                    if (
                        isinstance(re._elem._parentCompSel, _ComplexReactants)
                        and re._elem._parentCompSel._filter == cs
                    ) or re._elem._parentCompSel == cs:
                        candidates.append(re._elem)
                for re in lhsCompSels:
                    if re._elem == cs:
                        for state in re._elem._getAllStates():
                            candidates.append(state)
            elif isinstance(cs, Complex):
                for re in lhsCompStates:
                    if re._elem._comp is cs:
                        candidates.append(re._elem)
                for re in lhsCompSels:
                    if re._elem._complex is cs:
                        for state in re._elem._getAllStates():
                            candidates.append(state)
            else:
                raise TypeError(f'Expected a Complex or a ComplexSelector, got {cs} instead.')
            if len(candidates) == 1:
                states.append(candidates[0])
            elif len(candidates) > 1:
                raise Exception(
                    f'More than one complex state in {lhs} corresponded to {cs} for '
                    f'determining a rate: {candidates}'
                )
            else:
                raise Exception(f'No complex states in {lhs} matched {cs}.')
        return states

    def __call__(self, *lst):
        if all(isinstance(e, ComplexState) for e in lst):
            return super().__call__(*lst)
        elif len(lst) == 1 and isinstance(lst[0], ReactionSide):
            return super().__call__(*self._getStatesFromReacSide(lst[0]))
        else:
            raise NotImplementedError()


@nutils.FreezeAfterInit
class VDepRate(VDepFunc):
    """Voltage dependent reaction rate

    See base class :py:class:`VDepFunc` for how to create.

    Usage::

        with CaPchan[...]:

            # VDepRate declaration
            # alpha_cap and beta_cap are normal functions taking the membrane potential in Volts
            # as parameter and returning the corresponding rate.

            # with custom voltage range:
            a_cap = VDepRate(alpha_cap, vrange = [-100.0e-3, 50e-3, 1e-4])

            # with default voltage range
            b_cap = VDepRate(beta_cap)

            m0.s <r[1]> m1.s <r[2]> m2.s <r[3]> m3.s
            r[1].K = 3 * a_cap, 1 * b_cap
            r[2].K = 2 * a_cap, 2 * b_cap
            r[3].K = 1 * a_cap, 3 * b_cap

            # The transitions from m0 to m3 have rates that are multiples of a_cap and b_cap
            # it is thus easier to first define these two VDepRate functions and reuse them
            # when declaring reactions.
    """

    def _getVDepSReacParams(self, units):
        """Return the parameters that are actually needed by the cython constructor"""
        vmin, vmax, dv = self.vrange
        tablesize = int((vmax - vmin) / dv) + 1

        # Retrieve the parameters used when the function is called with vmin
        # This can possibly miss some parameters if they are only involved in the computation
        # for some specific values of the membrane potential
        self(vmin)

        ktab = []
        for i in range(tablesize):
            ret = self(vmin + i * dv, _recParams=False)
            ktab.append(ret.valueIn(units) if isinstance(ret, nutils.Parameter) else ret)

        return dict(ktab=ktab, vmin=vmin, vmax=vmax, dv=dv, tablesize=tablesize)


@nutils.FreezeAfterInit
class CompDepRate(CompDepFunc):
    """Complex state dependent reaction rate

    See base class :py:class:`CompDepFunc` for how to create.

    Usage::

        with Comp[...]:
            S1A >r[1]> S1B

            rate = CompDepRate(lambda s: k1 + k2 * s.Count(S1B), [Comp])

            r[1].K = rate
            # The above reaction specifies that the subunit state S1A changes to S1B with
            # a basal rate of k1 and this rate is increased by k2 for each subunit in the complex
            # being in state S1B.

        with Comp1[...] as C1, Comp2[...] as C2:
            S1A[C1] + S3A[C2] >r[1]> S1B[C1] + S3B[C2]

            rate = CompDepRate(lambda s1, s2: k3 + k4 * s1.Count(S2A) * s2.Count(S4A), [C1, C2])

            r[1].K = rate
            # The above reaction takes place between a subunit of Comp1 in state S1A and a
            # subunit of Comp2 in state S3A. It takes place with a basal rate k3 and this rate
            # is increased by a value that depends on the number of subunits in state S2A in
            # Comp1 and the number of subunits in state S4A in Comp2.

    """

    pass


@nutils.FreezeAfterInit
class CompDepCond(CompDepFunc):
    """Complex state dependent conductance

    See base class :py:class:`CompDepFunc` for how to create.
    """

    pass


@nutils.FreezeAfterInit
class CompDepP(CompDepFunc):
    """Complex state dependent permeability

    See base class :py:class:`CompDepFunc` for how to create.
    """

    pass


@nutils.FreezeAfterInit
class CompDepDcst(CompDepFunc):
    """Complex state dependent diffusion constant

    See base class :py:class:`CompDepFunc` for how to create.
    """

    pass


###################################################################################################
# Diffusion


class _SubDiffusion(nutils.NamedObject):
    """Wrapper for a single steps diffusion object"""

    _DiffStr = 'Diff'
    _SDiffStr = 'SDiff'

    _DiffDirecStr = 'direction_tet'
    _SDiffDirecStr = 'direction_tri'

    def __init__(self, stepsDiff, parent, keyElem=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._stepsDiff = stepsDiff
        self._elemStr = _SubDiffusion._SDiffStr if parent._isSurfaceDiff() else _SubDiffusion._DiffStr
        self._keyElem = keyElem

        self._direcStr = (
            _SubDiffusion._SDiffDirecStr if parent._isSurfaceDiff() else _SubDiffusion._DiffDirecStr
        )
        self._direc = None

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self._stepsDiff]

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return self._elemStr

    def _solverKeywordParams(self):
        """Return the additional keyword parameters that should be passed to the solver"""
        if self._direc is not None:
            return {self._direcStr: self._direc._solverId()[0]}
        else:
            return {}

    def __call__(self, direc=None):
        """Return a specialized version of the subdiffusion"""
        sd = copy.copy(self)
        sd._direc = direc
        return sd

    def _solverSetValue(self, valName, v):
        """
        Return the value that should actually be set in the solver when value 'v' is given by
        the user.
        """
        if valName == 'D' and isinstance(v, CompDepDcst):
            return v(self._keyElem)
        else:
            return v


class _SubDiffusionList(nutils.SolverPathObject):
    """List of _SubDiffusion objects"""

    def __init__(self, lst, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._diffLst = lst

    def __iter__(self):
        return iter(self._diffLst)

    def _simPathWalkExpand(self):
        """Return an iterable of the elements that should be part of ResultSelector paths."""
        return self.__iter__()

    def __call__(self, direc=None):
        return _SubDiffusionList([d.__call__(direc) for d in self._diffLst])


@nutils.FreezeAfterInit
class Diffusion(
        nutils.UsingObjects(Model, (VolumeSystem, SurfaceSystem)),
        nutils.StepsWrapperObject, nutils.ParameterizedObject
    ):
    """A diffusion rule for a chemical species in a volume or on a surface

    :param elem: The element that subject to diffusion
    :type elem: Union[:py:class:`Species`, :py:class:`ComplexState`, :py:class:`Complex`,
        :py:class:`ComplexSelector`]
    :param Dcst: The diffusion constant in S.I. units
    :type Dcst: Union[float, :py:class:`CompDepDcst`]

    Usage::

        with vsys:

            # Volume diffusion of Species S1
            # The diffusion rule is named 'diffS1' for later access during simulation
            diffS1 = Diffusion.Create(S1, 0.2e-9)

        with ssys:

            # Surface diffusion of Species S2
            # The diffusion rule is not named and will thus not be accessible during simulation
            Diffusion(S2, 0.08e-12)
    """

    def __init__(self, elem, Dcst, *args, _createObj=True, **kwargs):
        super().__init__(*args, **kwargs)
        mdl, volSys, surfSys = self._getUsedObjects()
        if volSys is not None and surfSys is not None:
            raise Exception(
                'A diffusion rule should only apply to a volume or a surface system, but not both.'
            )

        if not isinstance(elem, nutils.NamedObject):
            raise TypeError(f'Expected a named ligand, got {elem} instead.')

        self.sys = volSys if surfSys is None else surfSys
        self._stepsDiffs = {}
        self._elem = elem._getReferenceObject()

        if _createObj:
            if isinstance(self.sys, VesicleSurfaceSystem):
                cls = stepslib._py_VesSDiff
            elif self.sys.__class__ in [VolumeSystem, SurfaceSystem]:
                cls = stepslib._py_Diff
            else:
                raise Exception(f'Cannot declare a diffusion rule in {self.sys}')

            if (
                isinstance(elem, (Complex, ComplexSelector))
            ) and elem._areStatesAsSpecies():
                for state in elem:
                    name = f'{self.name}_{state.name}'
                    self._stepsDiffs[state] = _SubDiffusion(
                        cls(name, self.sys._getStepsObjects()[0], state._getStepsObjects()[0]),
                        self,
                        keyElem=state,
                        name=name,
                    )
            elif isinstance(elem, (Species, ComplexState)):
                self._stepsDiffs[elem] = _SubDiffusion(
                    cls(self.name, self.sys._getStepsObjects()[0], elem._getStepsObjects()[0]),
                    self,
                    keyElem=elem,
                    name=self.name,
                )
            else:
                raise TypeError(f'Cannot declare a diffusion rule for {elem}.')

        # Setting diffusion constant
        self.Dcst = Dcst

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [sd._getStepsObjects()[0] for _, sd in self._stepsDiffs.items()]

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        lig = getattr(mdl, obj.getLig().getID())
        diff = cls(lig, obj.getDcst(), _createObj=False, name=obj.getID())
        diff._stepsDiffs[lig] = _SubDiffusion(obj, diff, keyElem=lig, name=diff.name)
        return diff

    def _isSurfaceDiff(self):
        return isinstance(self.sys, SurfaceSystem)

    def __getitem__(self, key):
        """Get sub diffusion rules from a diffusion object

        Square bracket notation can be used on a diffusion rule when it was defined for a complex
        or a complex selector. This should never be needed during model declaration but it is very
        important for simulation control and data saving (see :py:class:`steps.API_2.sim.SimPath`).

        :param key: A sequence of one or more :py:class:`SubUnitState` or
            :py:class:`SubUnitSelector` that corresponds to one or several :py:class:`ComplexState`.
        :type key: Tuple[Union[:py:class:`SubUnitState`, :py:class:`SubUnitSelector`], ...]

        :returns: The diffusion rule associated to ``key``

        :meta public:
        """
        if isinstance(self._elem, Complex):
            e = self._elem.__getitem__(key)
            if isinstance(e, ComplexState):
                return self._stepsDiffs[e]
            elif isinstance(e, ComplexSelector):
                return _SubDiffusionList([self._stepsDiffs[s] for s in e])
        raise NotImplementedError()

    def __call__(self, direc=None):
        """Get a version of the diffusion rule with added information

        Parentheses notation (function call) can be used on a diffusion rule to represent a
        diffusion rule further specified with additional information. This should never be needed
        during model declaration but can become handy for simulation control and data saving
        (see :py:class:`steps.API_2.sim.SimPath`).

        :param direc: A direction for the diffusion rule
        :type key: :py:class:`Reference`

        :returns: A version of the diffusion rule that is specific to direction ``direc``.

        :meta public:
        """
        return _SubDiffusionList([sd.__call__(direc=direc) for _, sd in self._stepsDiffs.items()])

    def _getAllElems(self, loc):
        if loc is None:
            return [self._elem]
        else:
            if self._isSurfaceDiff():
                return [self._elem] if loc == Location.SURF else []
            else:
                return []

    def _simPathWalkExpand(self):
        """Return an iterable of the elements that should be part of ResultSelector paths."""
        for _, sd in self._stepsDiffs.items():
            yield sd

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m^2 s^-1'))
    def Dcst(self):
        """Diffusion constant (in m^2 s^-1)

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        pass

    @Dcst.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('m^2 s^-1'))
    def Dcst(self, val):
        if isinstance(val, numbers.Number):
            for _, sd in self._stepsDiffs.items():
                sd._stepsDiff.setDcst(val)
        elif isinstance(val, CompDepFunc):
            for cs, sd in self._stepsDiffs.items():
                sd._stepsDiff.setDcst(val(cs))
        else:
            raise TypeError(f'{val} is not a valid diffusion constant.')


###################################################################################################
# Vesicle and Raft creation / deletion


@nutils.FreezeAfterInit
class Endocytosis(
        nutils.UsingObjects(Model, SurfaceSystem),
        nutils.StepsWrapperObject, nutils.ParameterizedObject
    ):
    """A specific type of interaction that models endocytosis of a vesicle

    An endocytosis event models the process of vesicle endocytosis by creating
    a new vesicle within a compartment. A vesicle will be created at a given
    rate when an optional species signature is met within an endocytic zone of
    a patch (see :py:class:`steps.API_2.geom.EndocyticZone`).

    :param vesicle: Vesicle that should appear
    :type vesicle: :py:class:`Vesicle`
    :param kcst: Rate constant of the reaction in s^-1
    :type kcst: float
    :param deps: Species dependencies
    :type deps: Union[None, :py:class:`ReactionElement`, :py:class:`ReactionSide`]
    :param inside: If true, the vesicle should appear in the inside compartment, otherwise it should appear
        in the outside compartment.
    :type inside: bool
    """

    def __init__(self, vesicle, kcst=0, deps=None, inside=True, _createObj=True, **kwargs):
        super().__init__(**kwargs)
        mdl, surfSys = self._getUsedObjects()

        if isinstance(deps, ReactionElement):
            deps = deps._toReactionSide()
        elif deps is None:
            deps = ReactionSide([])
        elif not isinstance(deps, ReactionSide):
            raise TypeError(f'Expected species dependency, got {deps} instead.')

        if any(elem.loc not in [None, Location.SURF] for elem in deps):
            raise Exception(f'Species dependencies cannot include species in volumes.')

        if not isinstance(vesicle, Vesicle):
            raise TypeError(f'Expected a vesicle, got {vesicle} instead.')

        self._vesicle = vesicle
        self._setParameter('Dependencies', deps)
        self._setParameter('K', kcst, nutils.Units('s^-1'))

        if _createObj:
            deps = self.Dependencies._GetStepsElems()
            if inside:
                self.stepsEndo = stepslib._py_Endocytosis(
                    self.name, surfSys._getStepsObjects()[0], vesicle._getStepsObjects()[0], None,
                    deps, self.K
                )
            else:
                self.stepsEndo = stepslib._py_Endocytosis(
                    self.name, surfSys._getStepsObjects()[0], None, vesicle._getStepsObjects()[0],
                    deps, self.K
                )
        else:
            self.stepsEndo = None

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsEndo]

    def _getAllElems(self, loc=None):
        """
        Return all elements that can directly be accessed by name from the parent physical
        location.
        """
        return self.Dependencies._getAllElems(loc)

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        return {'obj_type': self._solverStr(), 'obj_id': self.name}

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        ves = getattr(mdl, obj.getIRHS().getID())
        deps = ReactionSide([])
        for spec in obj.getSpecDeps():
            deps += getattr(mdl, spec.getID())

        endo = cls(ves, obj.getKcst(), deps=deps, inside=obj.getInner(), _createObj=False, name=obj.getID())
        endo.stepsEndo = obj

        return endo

    @property
    def Vesicle(self):
        """Vesicle associated with the endocytosis reaction

        :type: :py:class:`Vesicle`, read-only
        """
        return self._vesicle

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def Dependencies(self):
        """Species dependencies for the endocytosis reaction

        :type: :py:class:`ReactionSide`, read-only
        """
        pass

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('s^-1'))
    def K(self):
        """Endocytosis reaction constant (in s^-1)

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        pass

    @K.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('s^-1'))
    def K(self, val):
        if isinstance(val, numbers.Number):
            self.stepsEndo.setKcst(val)
        else:
            raise TypeError(f'{val} is not a valid reaction constant.')


@nutils.FreezeAfterInit
class Exocytosis(
        nutils.UsingObjects(Model, VesicleSurfaceSystem),
        nutils.StepsWrapperObject, nutils.ParameterizedObject
    ):
    """A specific type of interaction that models exocytosis of a vesicle

    An exocytosis event models the process of vesicle exocytosis. By default the vesicle
    is destroyed, species in the vesicle surface are deposited in the patch at
    the location at which exocytosis occurs, and species inside the vesicle lumen
    are deposited in the opposite compartment.


    :param kcst: Rate constant of the exocytosis reaction in s^-1
    :type kcst: float
    :param deps: Species dependencies on the vesicle surface
    :type deps: Union[None, :py:class:`ReactionElement`, :py:class:`ReactionSide`]
    :param raft: If this argument is given the exocytosis event will create a
        raft on the patch at the location that exocytosis occurs. Any species on
        the vesicle surface will be deposited in the raft instead of in the patch.
    :type raft: :py:class:`Raft`
    :param kissAndRun: This parameter can be set to True to model a kiss-and-run
        exocytosis event. This does not result in collapse of the vesicle,
        instead the vesicle is maintained in position after the release of
        vesicle lumen contents into the opposite compartment. The vesicle
        surface species are maintained on the vesicle.
    :type kissAndRun: bool
    :param knrSpecChanges: To aid kiss-and-run modeling, this argument can be a
        list of tuples of two :py:class:`Species`; when the kiss-and-run event
        happens, all vesicle surface species from the first element of the tuple
        get changed to species of the second element. This feature can be used to
        specify any species changes that take place on the vesicle surface upon a
        kiss-and-run exocytosis event. This can be useful for modeling maturation
        of a complex for which undocking may be dependent.
    :type knrSpeChanges: List[Tuple[:py:class:`Species`, :py:class:`Species`]]
    :param knrPartialRelease: To aid kiss-and-run modeling, this argument can
        be used to specify partial release of any inner species that are released
        upon a kiss-and-run exocytosis event. 1.0 means all species are released,
        0.0 means none of the species are released. If not specified, it defaults
        to 1.0.
    :type knrPartialRelease: float
    """

    def __init__(self, kcst=0, deps=None, raft=None, kissAndRun=False, knrSpecChanges=None,
            knrPartialRelease=1, _createObj=True, **kwargs):
        super().__init__(**kwargs)
        mdl, surfSys = self._getUsedObjects()

        if isinstance(deps, ReactionElement):
            deps = deps._toReactionSide()
        elif deps is None:
            deps = ReactionSide([])
        elif not isinstance(deps, ReactionSide):
            raise TypeError(f'Expected species dependency, got {deps} instead.')

        if any(elem.loc not in [None, Location.VESSURF] for elem in deps):
            raise Exception(f'Species dependencies cannot include species in volumes.')

        if raft is not None and not isinstance(raft, Raft):
            raise TypeError(f'Expected a raft, got {raft} instad.')

        self._raft = raft
        self._setParameter('Dependencies', deps)
        self._setParameter('KissAndRun', kissAndRun)
        self._setParameter('KnRSpecChanges', knrSpecChanges)
        self._setParameter('KnRPartialRelease', knrSpecChanges)
        self._setParameter('K', kcst, nutils.Units('s^-1'))

        if _createObj:
            if self.KnRSpecChanges is not None:
                knrSpecChanges = {
                    spec1._getStepsObjects()[0]: spec2._getStepsObjects()[0]
                        for spec1, spec2 in self.KnRSpecChanges
                }
            else:
                knrSpecChanges = {}
            deps = self.Dependencies._GetStepsElems()
            self.stepsExo = stepslib._py_Exocytosis(
                self.name, surfSys._getStepsObjects()[0], deps,
                raft._getStepsObjects()[0] if raft is not None else None,
                self.K, self.KissAndRun, knrSpecChanges, knrPartialRelease
            )
        else:
            self.stepsExo = None

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsExo]

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        return {'obj_type': self._solverStr(), 'obj_id': self.name}

    def _getAllElems(self, loc=None):
        """
        Return all elements that can directly be accessed by name from the parent physical
        location.
        """
        return self.Dependencies._getAllElems(loc)

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        if obj.getRaft() is not None:
            raft = getattr(mdl, obj.getRaft().getID())
        else:
            raft = None
        deps = ReactionSide([])
        for spec in obj.getSpecDeps():
            deps += getattr(mdl, spec.getID())
        knrSpecChanges = list(obj.getKissAndRunSpecChanges().items())

        exo = cls(
            obj.getKcst(), deps=deps, raft=raft, kissAndRun=obj.getKissAndRun(),
            knrSpecChanges=knrSpecChanges, knrPartialRelease=obj.getKissAndRunPartRelease(),
            _createObj=False, name=obj.getID()
        )
        exo.stepsExo = obj

        return exo

    @property
    def Raft(self):
        """Raft optionally associated with the exocytosis reaction

        :type: Union[:py:class:`Raft`, None], read-only
        """
        return self._raft

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def Dependencies(self):
        """Species dependencies for the exocytosis reaction

        :type: :py:class:`ReactionSide`, read-only
        """
        pass

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def KissAndRun(self):
        """Whether the exocytosis is a kiss-and-run exocytosis

        :type: bool, read-only
        """
        pass

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def KnRSpecChanges(self):
        """Species changes upon kiss-and-run ecoxytosis reaction

        :type: List[Tuple[:py:class:`Species`, :py:class:`Species`]], read-only
        """
        pass

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def KnRPartialRelease(self):
        """Fraction of inner species that are released upon a kiss-and-run exocytosis event

        :type: float, read-only
        """
        pass

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('s^-1'))
    def K(self):
        """Exocytosys reaction constant (in s^-1)

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        pass

    @K.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('s^-1'))
    def K(self, val):
        if isinstance(val, numbers.Number):
            self.stepsExo.setKcst(val)
        else:
            raise TypeError(f'{val} is not a valid reaction constant.')


@nutils.FreezeAfterInit
class RaftEndocytosis(
        nutils.UsingObjects(Model, RaftSurfaceSystem),
        nutils.StepsWrapperObject, nutils.ParameterizedObject
    ):
    """A specific type of interaction that models endocytosis of a raft to form a vesicle

    A raft endocytosis event models the process of vesicle endocytosis by
    creating a new vesicle within a compartment from a raft. A vesicle will be
    created at a given rate when an optional species signature is met within
    the raft. 

    :param vesicle: Vesicle that should appear
    :type vesicle: :py:class:`Vesicle`
    :param kcst: Rate constant of the reaction
    :type kcst: float
    :param deps: Species dependencies for endocytosis on the raft
    :type deps: Union[None, :py:class:`ReactionElement`, :py:class:`ReactionSide`]
    :param inside: If true, the vesicle should appear in the inside compartment, otherwise it should appear
        in the outside compartment.
    :type inside: bool
    """

    def __init__(self, vesicle, kcst=0, deps=None, inside=True, _createObj=True, **kwargs):
        super().__init__(**kwargs)
        mdl, surfSys = self._getUsedObjects()

        if isinstance(deps, ReactionElement):
            deps = deps._toReactionSide()
        elif deps is None:
            deps = ReactionSide([])
        elif not isinstance(deps, ReactionSide):
            raise TypeError(f'Expected species dependency, got {deps} instead.')

        if any(elem.loc not in [None, Location.RAFTSURF] for elem in deps):
            raise Exception(f'Species dependencies cannot include species in volumes.')

        if not isinstance(vesicle, Vesicle):
            raise TypeError(f'Expected a vesicle, got {vesicle} instad.')

        self._vesicle = vesicle
        self._setParameter('Dependencies', deps)
        self._setParameter('K', kcst, nutils.Units('s^-1'))

        if _createObj:
            deps = self.Dependencies._GetStepsElems()
            if inside:
                self.stepsRaftEndo = stepslib._py_RaftEndocytosis(
                    self.name, surfSys._getStepsObjects()[0], vesicle._getStepsObjects()[0], None,
                    deps, self.K
                )
            else:
                self.stepsRaftEndo = stepslib._py_RaftEndocytosis(
                    self.name, surfSys._getStepsObjects()[0], None, vesicle._getStepsObjects()[0],
                    deps, self.K
                )
        else:
            self.stepsRaftEndo = None

    def _getAllElems(self, loc=None):
        """
        Return all elements that can directly be accessed by name from the parent physical
        location.
        """
        return self.Dependencies._getAllElems(loc)

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsRaftEndo]

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        return {'obj_type': self._solverStr(), 'obj_id': self.name}

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        ves = getattr(mdl, obj.getRHS().getID())
        deps = ReactionSide([])
        for spec in obj.getSpecDeps():
            deps += getattr(mdl, spec.getID())

        endo = cls(ves, obj.getKcst(), deps=deps, inside=obj.getInner(), _createObj=False, name=obj.getID())
        endo.stepsRaftEndo = obj

        return endo

    @property
    def Vesicle(self):
        """Vesicle associated with the raft endocytosis reaction

        :type: :py:class:`Vesicle`, read-only
        """
        return self._vesicle

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def Dependencies(self):
        """Species dependencies for the raft endocytosis reaction

        :type: :py:class:`ReactionSide`, read-only
        """
        pass

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('s^-1'))
    def K(self):
        """RaftEndocytosys reaction constant (in s^-1)

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        pass

    @K.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('s^-1'))
    def K(self, val):
        if isinstance(val, numbers.Number):
            self.stepsRaftEndo.setKcst(val)
        else:
            raise TypeError(f'{val} is not a valid reaction constant.')


class RaftGen(
        nutils.UsingObjects(SurfaceSystem),
        nutils.StepsWrapperObject, nutils.ParameterizedObject
    ):
    """A specific type of interaction that models a raft genesis event

    Generate a raft at a given rate when a defined species signature is met
    within a patch triangle.

    :param raft: Raft that should be created
    :type raft: :py:class:`Raft`
    :param kcst: Rate constant
    :type kcst: float
    :param deps: Species dependencies on the patch surface
    :type deps: Union[None, :py:class:`ReactionElement`, :py:class:`ReactionSide`]
    """

    def __init__(self, raft, kcst=0, deps=None, _createObj=True, **kwargs):
        super().__init__(**kwargs)
        surfSys, = self._getUsedObjects()

        if isinstance(deps, ReactionElement):
            deps = deps._toReactionSide()
        elif deps is None:
            deps = ReactionSide([])
        elif not isinstance(deps, ReactionSide):
            raise TypeError(f'Expected species dependency, got {deps} instead.')

        if any(elem.loc not in [None, Location.SURF] for elem in deps):
            raise Exception(f'Species dependencies cannot include species in volumes.')

        if not isinstance(raft, Raft):
            raise TypeError(f'Expected a raft, got {raft} instad.')

        self._raft = raft
        self._setParameter('Dependencies', deps)
        self._setParameter('K', kcst, nutils.Units('s^-1'))

        if _createObj:
            deps = self.Dependencies._GetStepsElems()
            self.stepsRaftGen = stepslib._py_RaftGen(
                self.name, surfSys._getStepsObjects()[0], deps,
                raft._getStepsObjects()[0], self.K
            )
        else:
            self.stepsRaftGen = None

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsRaftGen]

    def _getAllElems(self, loc=None):
        """
        Return all elements that can directly be accessed by name from the parent physical
        location.
        """
        return [self.Raft] + list(self.Dependencies._getAllElems(loc))

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        if obj.getRaft() is not None:
            raft = getattr(mdl, obj.getRaft().getID())
        else:
            raft = None
        deps = ReactionSide([])
        for spec in obj.getSpecDeps():
            deps += getattr(mdl, spec.getID())

        exo = cls(raft, obj.getKcst(), deps=deps, _createObj=False, name=obj.getID())
        exo.stepsRaftGen = obj

        return exo

    @property
    def Raft(self):
        """Raft associated with the genesis reaction

        :type: :py:class:`Raft`, read-only
        """
        return self._raft

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def Dependencies(self):
        """Species dependencies for the raft genesis reaction

        :type: :py:class:`ReactionSide`, read-only
        """
        pass

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('s^-1'))
    def K(self):
        """Raft genesis reaction constant (in s^-1)

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        pass

    @K.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('s^-1'))
    def K(self, val):
        if isinstance(val, numbers.Number):
            self.stepsRaftGen.setKcst(val)
        else:
            raise TypeError(f'{val} is not a valid reaction constant.')


@nutils.FreezeAfterInit
class RaftDis(
        nutils.UsingObjects(RaftSurfaceSystem),
        nutils.StepsWrapperObject, nutils.ParameterizedObject
    ):
    """A specific type of interaction that models a raft dissolution event

    A raft dissolution event results in removal of a raft at a given rate when
    the population of species within a raft are **at or below** an
    anti-dependencies signature. Any remaining species in the raft are
    inserted into patch triangles.

    :param kcst: Rate constant in s^-1
    :type kcst: float
    :param antideps: Species anti-dependencies on the raft surface. The raft
        dissolution cannot happen if species count on the raft are **strictly
        bigger** than at least one of the anti-dependencies.
    :type deps: Union[None, :py:class:`ReactionElement`, :py:class:`ReactionSide`]
    """

    def __init__(self, kcst=0, antideps=None, _createObj=True, **kwargs):
        super().__init__(**kwargs)
        surfSys, = self._getUsedObjects()

        if isinstance(antideps, ReactionElement):
            antideps = antideps._toReactionSide()
        elif antideps is None:
            antideps = ReactionSide([])
        elif not isinstance(antideps, ReactionSide):
            raise TypeError(f'Expected species dependency, got {antideps} instead.')

        if any(elem.loc not in [None, Location.RAFTSURF] for elem in antideps):
            raise Exception(f'Species dependencies cannot include species in volumes.')

        self._setParameter('AntiDependencies', antideps)
        self._setParameter('K', kcst, nutils.Units('s^-1'))

        if _createObj:
            antideps = self.AntiDependencies._GetStepsElems()
            self.stepsRaftDis = stepslib._py_RaftDis(
                self.name, surfSys._getStepsObjects()[0], antideps, self.K
            )
        else:
            self.stepsRaftDis = None

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsRaftDis]

    def _getAllElems(self, loc=None):
        """
        Return all elements that can directly be accessed by name from the parent physical
        location.
        """
        return self.AntiDependencies._getAllElems(loc)

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""
        antideps = ReactionSide([])
        for spec in obj.getSpecAntiDeps():
            antideps += getattr(mdl, spec.getID())

        exo = cls(obj.getKcst(), antideps=antideps, _createObj=False, name=obj.getID())
        exo.stepsRaftDis = obj

        return exo

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def AntiDependencies(self):
        """Species antidependencies for the raft dissolution reaction

        :type: :py:class:`ReactionSide`, read-only
        """
        pass

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('s^-1'))
    def K(self):
        """Raft dissolution reaction constant (in s^-1)

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        pass

    @K.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('s^-1'))
    def K(self, val):
        if isinstance(val, numbers.Number):
            self.stepsRaftDis.setKcst(val)
        else:
            raise TypeError(f'{val} is not a valid reaction constant.')


@nutils.FreezeAfterInit
class VesicleBind(
        nutils.UsingObjects(VolumeSystem),
        nutils.StepsWrapperObject, nutils.ParameterizedObject
    ):
    """Vesicle binding reaction in a volume system

    A vesicle binding event binds two vesicles by creating link species between them.
    The types of vesicle are arbitrarily termed 1 and 2, and may be the same.
    A reactant on vesicle surface of vesicle 1 binds with a reactant on vesicle
    2, and link species are created, one on vesicle 1 and one on vesicle 2.
    Internally these two link species will be associated with each other and
    exist within a defined length bound.

    :param vesicles: The two vesicles involved in the binding reaction (in a tuple)
    :type vesicles: Tuple[:py:class:`Vesicle`, :py:class:`Vesicle`]
    :param reactants: The species involved in the binding reaction for each vesicle (in a tuple,
        the order should match the vesicles parameter)
    :type reactants: Tuple[:py:class:`Species`, :py:class:`Species`]
    :param linkProducts: The link species resulting from the binding for each vesicle (in a tuple,
        the order should match the vesicles parameter)
    :type linkProducts: Tuple[:py:class:`LinkSpecies`, :py:class:`LinkSpecies`]
    :param lenMin: Minimum length of the link between the two vesicles, in m
    :type lenMin: float
    :param lenMax: Minimum length of the link between the two vesicles, in m
    :type lenMax: float
    :param kcst: Rate constant in M^-1 s^-1
    :type kcst: float
    :param immobilization: immobilization flag for the binding reaction (see
        :py:attr:`VesicleBind.Immobilization`)
    :type immobilization: Union[:py:data:`steps.API_2.model.IMMOBILIZING`,
        :py:data:`steps.API_2.model.NO_EFFECT`, :py:data:`steps.API_2.model.MOBILIZING`]
    :param deps: Dependencies for the binding reaction (in a tuple, the order should match the vesicle
        parameter)
    :type deps: Tuple[Union[None, :py:class:`ReactionElement`, :py:class:`ReactionSide`],
        Union[None, :py:class:`ReactionElement`, :py:class:`ReactionSide`]]
    """

    def __init__(
            self, vesicles, reactants, linkProducts, lenMin, lenMax, 
            kcst=0, immobilization=NO_EFFECT, deps=(None, None), **kwargs
        ):
        super().__init__(**kwargs)
        self._volSys, = self._getUsedObjects()

        self._stepsVesBind = None
        self._declared = False
        self._added = False

        # Check that the arguments are of correct types
        for val, cls in [(vesicles, Vesicle), (reactants, Species), (linkProducts, LinkSpecies)]:
            if not (isinstance(val, tuple) and len(val) == 2 and all(isinstance(v, cls) for v in val)):
                raise TypeError(f'Expected a 2-tuple of {cls.__name__}, got {val} instead.')

        self._vesicles = vesicles
        self._reactants = reactants
        self._linkProducts = linkProducts

        self.LengthMin = lenMin
        self.LengthMax = lenMax
        self.K = kcst
        self.Immobilization = immobilization
        self.Dependencies = deps

        self._declared = True

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self._stepsVesBind] if self._stepsVesBind is not None else []

    def _getAllElems(self, loc=None):
        """
        Return all elements that can directly be accessed by name from the parent physical
        location.
        """
        return list(self._vesicles)

    def _getVesLinkSpecPairs(self):
        """Return pairs of (vesicle, link species) that are involved in the reaction"""
        return zip(self._vesicles, self._linkProducts)

    def _checkNotAdded(self):
        """Raise an exception if the vesbind reaction was already added to STEPS"""
        if self._added:
            raise Exception(
                'Cannot set properties of a vesicle bind reaction that was already added to STEPS.'
            )

    def _exitCallback(self, parent):
        """
        Method to be called the first time we get out of a context manager in which the object as been
        declared.
        """
        if not self._declared or self._added:
            return

        allDeps = []
        for deps in self.Dependencies:
            vdeps = []
            ldeps = []
            if deps is not None:
                deps = deps._toReactionSide()
                if any(elem.loc != Location.VESSURF for elem in deps):
                    raise Exception(
                        f'{self._decl()}: Dependencies should all be located on the vesicle surface '
                        f'(spec.v): {deps}'
                    )
                for elem in deps._GetStepsElems():
                    if isinstance(elem, stepslib._py_Spec):
                        vdeps.append(elem)
                    elif isinstance(elem, stepslib._py_LinkSpec):
                        ldeps.append(elem)
                    else:
                        raise TypeError(f'{self._decl()}: Unexpected element in dependencies: {elem}')
            allDeps.append((vdeps, ldeps))

        self._stepsVesBind = stepslib._py_VesBind(
            self.name, self._volSys._getStepsObjects()[0],
            self._vesicles[0]._getStepsObjects()[0], self._reactants[0]._getStepsObjects()[0],
            self._vesicles[1]._getStepsObjects()[0], self._reactants[1]._getStepsObjects()[0],
            self._linkProducts[0]._getStepsObjects()[0], self._linkProducts[1]._getStepsObjects()[0],
            self.LengthMax, self.LengthMin,
            allDeps[0][0], allDeps[1][0],
            allDeps[0][1], allDeps[1][1],
            self.K, self.Immobilization
        )

        self._added = True

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""

        ves1 = getattr(mdl, obj.getVesicle1().getID())
        ves2 = getattr(mdl, obj.getVesicle2().getID())
        r1 = getattr(mdl, obj.getReactant1().getID())
        r2 = getattr(mdl, obj.getReactant2().getID())
        l1 = getattr(mdl, obj.getProduct1().getID())
        l2 = getattr(mdl, obj.getProduct2().getID())
        dep1, dep2 = ReactionSide([]), ReactionSide([])
        for spec in obj.getVDeps1():
            dep1 += getattr(mdl, spec.getID()).v
        for spec in obj.getVDeps2():
            dep2 += getattr(mdl, spec.getID()).v
        for spec in obj.getLDeps1():
            dep1 += getattr(mdl, spec.getID()).v
        for spec in obj.getLDeps2():
            dep2 += getattr(mdl, spec.getID()).v

        vb = cls(
            (ves1, ves2), (r1, r2), (l1, l2),
            obj.getLengthMin(), obj.getLengthMax(), obj.getKcst(), obj.getImmobilization(),
            (dep1, dep2), name=obj.getID()
        )
        vb._stepsVesBind = obj
        vb._added = True

        return vb

    @property
    def Vesicles(self):
        """Vesicles associated with the binding reaction

        :type: Tuple[:py:class:`Vesicle`, :py:class:`Vesicle`], read-only
        """
        return self._vesicles

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m'))
    def LengthMin(self):
        """Minimum length of the link (in m)

        :type: float
        """
        pass

    @LengthMin.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('m'))
    def LengthMin(self, val):
        self._checkNotAdded()
        if not isinstance(val, numbers.Number):
            raise TypeError(f'{val} is not a valid minimum link length.')

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m'))
    def LengthMax(self):
        """Maximum length of the link (in m)

        :type: float
        """
        pass

    @LengthMax.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('m'))
    def LengthMax(self, val):
        self._checkNotAdded()
        if not isinstance(val, numbers.Number):
            raise TypeError(f'{val} is not a valid maximum link length.')

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('M^-1 s^-1'))
    def K(self):
        """Binding reaction constant (in M^-1 s^-1)

        :type: float
        """
        return 0

    @K.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('M^-1 s^-1'))
    def K(self, val):
        self._checkNotAdded()
        if not isinstance(val, numbers.Number):
            raise TypeError(f'{val} is not a valid reaction constant.')

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def Immobilization(self):
        """Immobilization flag (:py:data:`IMMOBILIZING`, :py:data:`NO_EFFECT` or :py:data:`MOBILIZING`).

        Defaults to :py:data:`NO_EFFECT`.

        :type: Union[:py:data:`steps.API_2.model.IMMOBILIZING`, :py:data:`steps.API_2.model.NO_EFFECT`,
            :py:data:`steps.API_2.model.MOBILIZING`]
        """
        return NO_EFFECT

    @Immobilization.setter
    @nutils.ParameterizedObject.RegisterSetter(units=None)
    def Immobilization(self, val):
        if not isinstance(val, Immobilization):
            raise TypeError(f'Expected an immobilization flag, got {val} instead.')
        self._checkNotAdded()

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def Dependencies(self):
        """Vesicle binding dependencies

        A set of species or link species that need to be present in order for the binding reaction to happen.

        The dependencies should be given in the shape of a reaction side, for example: `SA.v + SC.v`

        The dependencies are empty by default.

        :type: Tuple[Union[None, :py:class:`ReactionElement`, :py:class:`ReactionSide`],
            Union[None, :py:class:`ReactionElement`, :py:class:`ReactionSide`]]
        """
        # Return a default value if the property was not set
        return None, None

    @Dependencies.setter
    @nutils.ParameterizedObject.RegisterSetter(units=None)
    def Dependencies(self, deps):
        self._checkNotAdded()
        if (not isinstance(deps, tuple) or 
                not all(d is None or isinstance(d, (ReactionElement, ReactionSide)) for d in deps)):
            raise TypeError(f'{deps} cannot be used as a dependency.')


@nutils.FreezeAfterInit
class VesicleUnbind(
        nutils.UsingObjects(VolumeSystem),
        nutils.StepsWrapperObject, nutils.ParameterizedObject
    ):
    """Vesicle unbinding reaction in a volume system

    A vesicle unbinding event unbinds two vesicles bound by two link species.
    The types of vesicle are arbitrarily termed 1 and 2, and may be the same.
    Upon application of this reaction, the link species on vesicle 1 becomes a
    species on the vesicle 1 surface and the link species on vesicle 2 becomes a
    species on the vesicle 2 surface.

    :param vesicles: The two vesicles involved in the binding reaction (in a tuple)
    :type vesicles: Tuple[:py:class:`Vesicle`, :py:class:`Vesicle`]
    :param linkReactants: The link species involved in the link for each vesicle (in a tuple,
        the order should match the vesicles parameter)
    :type linkReactants: Tuple[:py:class:`LinkSpecies`, :py:class:`LinkSpecies`]
    :param products: The species resulting from the unbinding for each vesicle (in a tuple,
        the order should match the vesicles parameter)
    :type products: Tuple[:py:class:`Species`, :py:class:`Species`]
    :param kcst: Rate constant in s^-1
    :type kcst: float
    :param immobilization: immobilization flag for the binding reaction (see 
        :py:attr:`VesicleUnbind.Immobilization`)
    :type immobilization: Union[:py:data:`steps.API_2.model.IMMOBILIZING`,
        :py:data:`steps.API_2.model.NO_EFFECT`, :py:data:`steps.API_2.model.MOBILIZING`]
    """

    def __init__(self, vesicles, linkReactants, products, kcst=0, immobilization=NO_EFFECT, **kwargs):
        super().__init__(**kwargs)
        self._volSys, = self._getUsedObjects()

        self._stepsVesUnbind = None
        self._declared = False
        self._added = False

        # Check that the arguments are of correct types
        for val, cls in [(vesicles, Vesicle), (products, Species), (linkReactants, LinkSpecies)]:
            if not (isinstance(val, tuple) and len(val) == 2 and all(isinstance(v, cls) for v in val)):
                raise TypeError(f'Expected a 2-tuple of {cls.__name__}, got {val} instead.')

        self._vesicles = vesicles
        self._linkReactants = linkReactants
        self._products = products

        self.K = kcst
        self.Immobilization = immobilization

        self._declared = True

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self._stepsVesUnbind] if self._stepsVesUnbind is not None else []

    def _getAllElems(self, loc=None):
        """
        Return all elements that can directly be accessed by name from the parent physical
        location.
        """
        return list(self._vesicles)

    def _getVesLinkSpecPairs(self):
        """Return pairs of (vesicle, link species) that are involved in the reaction"""
        return zip(self._vesicles, self._linkReactants)

    def _checkNotAdded(self):
        """Raise an exception if the vesunbind reaction was already added to STEPS"""
        if self._added:
            raise Exception(
                'Cannot set properties of a vesicle unbind reaction that was already added to STEPS.'
            )

    def _exitCallback(self, parent):
        """
        Method to be called the first time we get out of a context manager in which the object as been
        declared.
        """
        if not self._declared or self._added:
            return

        self._stepsVesUnbind = stepslib._py_VesUnbind(
            self.name, self._volSys._getStepsObjects()[0],
            self._linkReactants[0]._getStepsObjects()[0], self._linkReactants[1]._getStepsObjects()[0],
            self._vesicles[0]._getStepsObjects()[0], self._products[0]._getStepsObjects()[0],
            self._vesicles[1]._getStepsObjects()[0], self._products[1]._getStepsObjects()[0],
            self.K, self.Immobilization
        )

        self._added = True

    @classmethod
    def _FromStepsObject(cls, obj, mdl):
        """Create the interface object from a STEPS object."""

        ves1 = getattr(mdl, obj.getVesicle1().getID())
        ves2 = getattr(mdl, obj.getVesicle2().getID())
        l1 = getattr(mdl, obj.getLink1().getID())
        l2 = getattr(mdl, obj.getLink2().getID())
        p1 = getattr(mdl, obj.getProduct1().getID())
        p2 = getattr(mdl, obj.getProduct2().getID())

        vb = cls((ves1, ves2), (l1, l2), (p1, p2), obj.getKcst(), obj.getImmobilization(), name=obj.getID())
        vb._stepsVesUnbind = obj
        vb._added = True

        return vb

    @property
    def Vesicles(self):
        """Vesicles associated with the binding reaction

        :type: Tuple[:py:class:`Vesicle`, :py:class:`Vesicle`], read-only
        """
        return self._vesicles

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('s^-1'))
    def K(self):
        """Binding reaction constant (in s^-1)

        :type: float
        """
        return 0

    @K.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('s^-1'))
    def K(self, val):
        self._checkNotAdded()
        if not isinstance(val, numbers.Number):
            raise TypeError(f'{val} is not a valid reaction constant.')

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=None)
    def Immobilization(self):
        """Immobilization flag (:py:data:`IMMOBILIZING`, :py:data:`NO_EFFECT` or :py:data:`MOBILIZING`)

        Defaults to :py:data:`NO_EFFECT`.

        :type: Union[:py:data:`steps.API_2.model.IMMOBILIZING`, :py:data:`steps.API_2.model.NO_EFFECT`,
            :py:data:`steps.API_2.model.MOBILIZING`]
        """
        return NO_EFFECT

    @Immobilization.setter
    @nutils.ParameterizedObject.RegisterSetter(units=None)
    def Immobilization(self, val):
        if not isinstance(val, Immobilization):
            raise TypeError(f'Expected an immobilization flag, got {val} instead.')
        self._checkNotAdded()


###################################################################################################
# Currents


class _SubCurrent(nutils.NamedObject):
    """Wrapper class for the steps OhmicCurr or GHKCurr object"""

    def __init__(self, parent, obj, *args, **kwargs):
        super().__init__(*args, name=obj.getID(), **kwargs)
        self._parent = parent
        self.stepsCurrent = obj

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsCurrent]

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return self._parent.__class__._currStr

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        mtdt = {'obj_type': self._parent.__class__._currStr, 'obj_id': self.name}
        # Add information about parent comp or patch
        for key, val in self._parent._simPathAutoMetaData().items():
            mtdt['parent_' + key] = val
        return mtdt


class _SubCurrentList:
    """List of _SubCurrent objects"""

    def __init__(self, lst, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._currentLst = lst

    def __iter__(self):
        return iter(self._currentLst)

    def _simPathWalkExpand(self):
        """Return an iterable of the elements that should be part of ResultSelector paths."""
        return self.__iter__()

    def _simPathCombinerClass(self):
        """Return the class that needs to be used to combine expanded elements."""
        return nutils.SumSimPathComb


@nutils.FreezeAfterInit
class Current(nutils.UsingObjects(SurfaceSystem), nutils.StepsWrapperObject, nutils.ParameterizedObject):
    """Base class representing a current associated with one or more channel states

    All :py:class:`Current` objects need to be declared inside a :py:class:`SurfaceSystem`.

    .. note::
        A :py:class:`Current` object should not be directly instantiated by the user,
        only its subclasses should be instantiated. It is only documented for clarity.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._currents = {}
        self._complex = None
        self._added = False
        self._declared = False

    def _getAllElems(self, loc):
        return (
            list(set(state._comp for state, curr in self._currents.items()))
            + [state for state, curr in self._currents.items()]
            + [curr for state, curr in self._currents.items()]
        )

    def _getComplexStates(self, states):
        self._setParameter('Opened state', states)
        if isinstance(states, Complex):
            return states, states[...]._getAllStates()
        if isinstance(states, ComplexSelector):
            return states._complex, states._getAllStates()
        elif isinstance(states, ComplexState):
            return states._comp, [states]
        else:
            raise TypeError(f'Expected a ComplexSelector or a ComplexState, got {states} instead.')

    def _getCurrentFromState(self, state):
        if state not in self._currents:
            raise Exception(f'ComplexState {state} is not associated to a _SubCurrent in {self}.')
        return self._currents[state]

    def __getitem__(self, key):
        """Get sub current from a current object

        Square bracket notation can be used on a current to retrieve the specific current object
        corresponding to a specific channel state. This should never be needed during model
        declaration but it is very important for simulation control and data saving (see
        :py:class:`steps.API_2.sim.SimPath`). If a :py:class:`ComplexSelector` covering several complex
        states is used, an iterable list of subcurrents is returned.

        :param key: A sequence of one or more :py:class:`SubUnitState` or
            :py:class:`SubUnitSelector` that corresponds to one or several :py:class:`ComplexState`.
        :type key: Tuple[Union[:py:class:`SubUnitState`, :py:class:`SubUnitSelector`], ...]

        :returns: The sub current(s) associated to ``key``

        :meta public:
        """
        state = self._complex.__getitem__(key)
        if isinstance(state, ComplexState):
            return self._getCurrentFromState(state)
        elif isinstance(state, ComplexSelector):
            return _SubCurrentList([self._getCurrentFromState(s) for s in state._getAllStates()])
        else:
            raise NotImplementedError()

    def __iter__(self):
        return iter([curr for state, curr in self._currents.items()])

    def _simPathWalkExpand(self):
        """Return an iterable of the elements that should be part of ResultSelector paths."""
        return self.__iter__()

    def _simPathCombinerClass(self):
        """Return the class that needs to be used to combine expanded elements."""
        if self._complex._statesAsSpecies:
            return nutils.SumSimPathComb
        else:
            return None

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [obj.stepsCurrent for state, obj in self._currents.items()]

    def _createStepsObj(self):
        """Create the actual STEPS objects and populate self._currents"""
        raise NotImplementedError()

    def _exitCallback(self, parent):
        """
        Method to be called the first time we get out of a context manager in which the object as been
        declared.
        """
        self._added = True
        self._createStepsObj()

    def _checkNotAdded(self):
        """Raise an exception if the currents were already added"""
        if self._added:
            raise Exception(f'Cannot modify the properties of a current that was already added to STEPS.')

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        return {'obj_type': self.__class__._currStr, 'obj_id': self.name}


@nutils.FreezeAfterInit
class OhmicCurr(Current):
    """An ohmic current through a channel in one or more conducting states

    :param states: The conducting state(s) of the channel
    :type states: Union[:py:class:`ComplexState`, :py:class:`ComplexSelector`]
    :param conduct: A constant conductance in Siemens or a :py:class:`CompDepCond` function that
        returns a conductance in Siemens as a function of a complex state. Defaults to 0 Siemens.
    :type conduct: Union[float, :py:class:`CompDepCond`]
    :param rev_pot: Reversal potential in Volts. Defaults to 0 Volts.
    :type rev_pot: float

    An ohmic current object is a simple current based on a (possibly state specific) fixed value
    for single-channel conductance and a constant reversal potential. An ohmic current
    does not result in movement of ions between compartments, but simply
    contributes a continuous current to the EField solver for every conducting channel
    on a membrane surface that contains the ohmic current.

    Usage::

        with ssys:

            # Constant conductance:

            OC_SK = OhmicCurr.Create(SKchan[sko1|sko2], SK_G, SK_rev)
            # The above line declares an ohmic current named OC_SK for channel SKchan when its
            # single subunit is in states sko1 or sko2.

            # Complex dependent conductance:

            cond = CompDepCond(lambda s: k1 + k2 * s.Count(S1A), [Chan1])

            OC_Chan1 = OhmicCurr.Create(Chan1[:, :, S2A], cond, e_rev)

            # The above 2 lines declare an ohmic current from channel Chan1 when its 3rd subunit
            # is in state S2A. Thanks to the CompDepCond function, the conductance changes 
            # depending on the specific states of the other subunits. In this case, there is a
            # basal conductance k1 and this conductance is increased by k2 for each subunit in
            # state S1A.

    .. warning::
        The corresponding currents are only truly added to steps after the `with ssys:` block
        is exited.

    """

    _currStr = 'Ohmic'

    def __init__(self, states, conduct=0, rev_pot=0, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._complex, self._states = self._getComplexStates(states)

        # Properties
        self.G = conduct
        self.ERev = rev_pot

        self._declared = True

    def _createStepsObj(self):
        """Create the actual STEPS objects and populate self._currents"""
        if not self._declared:
            return

        surfSys = self._getParentOfType(SurfaceSystem)
        conduct = self.G
        rev_pot = self.ERev

        if isinstance(conduct, CompDepFunc):
            if len(conduct._complexes) != 1 or conduct._complexes[0] is not self._complex:
                raise Exception('The CompDepFunc used as conductance is not compatible with the channel.')
            condStates = [(s, conduct(s)) for s in self._states]
        else:
            condStates = [(s, conduct) for s in self._states]

        for state, g in condStates:
            chanState = self._complex._compStates[state]
            self._currents[state] = _SubCurrent(
                self,
                stepslib._py_OhmicCurr(
                    f'{self.name}_{state.name}', surfSys.stepsSys, chanstate=chanState, g=g, erev=rev_pot
                ),
            )

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('V'))
    def ERev(self):
        """Reversal potential in Volts

        Can only be set in the `with ssys:` block in which the current has been declared.

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        pass

    @ERev.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('V'))
    def ERev(self, val):
        self._checkNotAdded()
        if not isinstance(val, numbers.Number):
            raise TypeError(f'Expected a float, got {val} instead')

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('S'))
    def G(self):
        """Conductance in Siemens

        Can only be set in the `with ssys:` block in which the current has been declared.

        :type: Union[float, :py:class:`CompDepCond`]
        """
        pass

    @G.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('S'))
    def G(self, val):
        self._checkNotAdded()
        if not isinstance(val, (numbers.Number, CompDepCond)):
            raise TypeError(f'Expected a float or a CompDepCond object, got {val} instead')


@nutils.FreezeAfterInit
class GHKCurr(Current):
    """A Goldman-Hodgkin-Katz current through a channel in one or more conducting states

    :param states: The conducting state(s) of the channel.
    :type states: Union[:py:class:`ComplexState`, :py:class:`ComplexSelector`]
    :param spec: The species involved in the current, must have a non zero
        :py:func:`Species.valence`.
    :type spec: :py:class:`Species`
    :param P: A single channel permeability in m^3/s, or a :py:class:`PInfo` object, or a
        :py:class:`CompDepP` function that returns a permeability as a function of a complex state.
    :type P: Union[float, :py:class:`PInfo`, :py:class:`CompDepP`]
    :param computeflux: If True, then the current will result in movement of ions between
        compartments, if False the current will be calculated but will not correspond to a real
        ion flux.
    :type computeflux: bool
    :param virtual_oconc: A 'virtual outer concentration' that can be specified so that the outer
        compartment does not have to be explicitly simulated (if it retains default `None`
        value then the outer concentration of ion must be simulated).
    :type virtual_oconc: float
    :param vshift:
    :type vshift: float

    Each GHK current in the simulation is solved within the SSA with a rate determined
    from the simulation state, i.e. membrane potential, 'outer' and 'inner'
    concentration of the ion and temperature (which is fixed), and constant parameters
    ion valence and permeability. The user can supply a call to :py:func:`PInfo` instead of a
    permeability, STEPS will compute the corresponding permeability.

    A GHK current involves, optionally, a real transfer of ions between comparments separated
    by a membrane in the STEPS simulation. It is important that these ions implement diffusion
    if compartments are not well-mixed.

    Usage::

        with mdl:

            Ca = Species.Create()
            Ca.valence = 2

            with ssys:

                # Constant permeability:

                GC_CaP = GHKCurr.Create(CaPchan[m3], Ca, CaP_P, computeflux = True)
                # The above line declares GHK current named OC_CaP for Ca2+ ions and for
                # channel CaPchan when its single subunit is in state m3.

                # Complex dependent permeability:

                P = CompDepP(lambda s: p1 + p2 * s.Count(S1A), [Chan2])

                GC_Chan2 = GHKCurr.Create(Chan2[:, :, S2A], Ca, P, computeflux = True)

                # The above 2 lines declare GHK current for Ca2+ ions and for channel Chan2
                # when its 3rd subunit is in state S2A. Thanks to the CompDepP function,
                # the permeability changes depending on the specific states of the other
                # subunits. In this case, there is a basal permeability p1 and this
                # permeability is increased by p2 for each subunit in state S1A.

    """

    _currStr = 'GHK'

    class _GHKCurrPInfo(nutils.ParameterizedObject):
        """Simple wrapper for providing Pinfo instead of permeability"""

        def __init__(self, g, V, T, oconc, iconc, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.g = g
            self.V = V
            self.T = T
            self.oconc = oconc
            self.iconc = iconc
            self._valence = None

        def _setValence(self, valence):
            self._valence = valence

        @property
        def value(self):
            if self._valence is None:
                raise Exception(
                    f'The PInfo object was not associated with any GHK current, cannot compute permeability.'
                )
            return stepslib._py_permeability(self.g, self.V, self._valence, self.T, self.iconc, self.oconc)

    _GHKCurrPInfo.RegisterParameter('g', nutils.Units('S'))
    _GHKCurrPInfo.RegisterParameter('V', nutils.Units('V'))
    _GHKCurrPInfo.RegisterParameter('T', nutils.Units('K'))
    _GHKCurrPInfo.RegisterParameter('oconc', nutils.Units('M'))
    _GHKCurrPInfo.RegisterParameter('iconc', nutils.Units('M'))

    def __init__(self, states, spec, P, computeflux=True, virtual_oconc=None, vshift=0, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if not isinstance(spec, Species):
            raise TypeError(f'Expected a Species, got {spec} instead.')
        self._spec = spec
        self._complex, self._states = self._getComplexStates(states)
        self._computeflux = computeflux
        self._vshift = vshift

        # Properties
        self.P = copy.copy(P)
        if isinstance(P, GHKCurr._GHKCurrPInfo):
            self.P._setValence(spec.valence)

        self.VOConc = virtual_oconc

        self._declared = True

    def _createStepsObj(self):
        """Create the actual STEPS objects and populate self._currents"""
        if not self._declared:
            return

        surfSys = self._getParentOfType(SurfaceSystem)
        P = self.P
        virtual_oconc = self.VOConc if self.VOConc is not None else -1

        if isinstance(P, CompDepFunc):
            if len(P._complexes) != 1 or P._complexes[0] is not self._complex:
                raise Exception('The CompDepFunc used as permeability is not compatible with the channel.')
            condStates = [(s, P(s)) for s in self._states]
        else:
            condStates = [(s, P) for s in self._states]

        for state, perm in condStates:
            chanState = self._complex._compStates[state]
            self._currents[state] = _SubCurrent(
                self,
                stepslib._py_GHKcurr(
                    f'{self.name}_{state.name}',
                    surfSys.stepsSys,
                    chanState,
                    self._spec.stepsSpecies,
                    computeflux=self._computeflux,
                    virtual_oconc=virtual_oconc,
                    vshift=self._vshift,
                ),
            )
            if isinstance(perm, numbers.Number):
                self._currents[state].stepsCurrent.setP(perm)
            elif isinstance(perm, GHKCurr._GHKCurrPInfo):
                self._currents[state].stepsCurrent.setPInfo(
                    perm.g, perm.V, perm.T, perm.oconc, perm.iconc
                )
            else:
                raise TypeError(f'Expected a permeability, got {perm} instead.')

    @classmethod
    def PInfo(cls, g, V, T, oconc, iconc):
        """Provide information from channel measurement instead of permeability

        :param g: A measured single-channel conductance in Siemens
        :type g: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        :param V: The potential in Volts
        :type V: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        :param T: The temperature in Kelvin (may be different from the STEPS simulation
            temperature)
        :type T: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        :param oconc: The 'outer' concentration of the ion in molar units
        :type oconc: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        :param iconc: The 'inner' concentration of the ion in molar units
        :type iconc: Union[float, :py:class:`steps.API_2.nutils.Parameter`]

        :returns: An aggregate object. Does not return the permeability.

        Usage::

            Pinfos = GHKCurr.PInfo(g=20e-12, V=-22e-3, T=293.15, oconc=4e-3, iconc=155e-3)

            GC_CaP = GHKCurr.Create(CaPchan[m3], Ca, Pinfos, computeflux = True)

        The actual permeability can be accessed with `GC_CaP.P.value`.
        """
        return GHKCurr._GHKCurrPInfo(g, V, T, oconc, iconc)

    def _getAllElems(self, loc):
        if loc != Location.SURF:
            return super()._getAllElems(loc) + [self._spec]
        else:
            return super()._getAllElems(loc)

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('M'))
    def VOConc(self):
        """Virtual outer concentration, in Molars

        Can only be set in the `with ssys:` block in which the current has been declared.

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        pass

    @VOConc.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('M'))
    def VOConc(self, val):
        self._checkNotAdded()

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m^3 s^-1'))
    def P(self):
        """Single channel permeability in m^3 s^-1

        Can only be set in the `with ssys:` block in which the current has been declared.

        :type: Union[float, :py:class:`steps.API_2.nutils.Parameter`, :py:class:`PInfo`, :py:class:`CompDepP`]
        """
        pass

    @P.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('m^3 s^-1'))
    def P(self, val):
        self._checkNotAdded()

