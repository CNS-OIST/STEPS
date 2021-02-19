####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

import collections
import copy
from enum import Enum
import inspect
import itertools
import numbers
import warnings

import steps.API_1.model as smodel

from . import utils as nutils
from . import geom as ngeom

try:
    from steps.API_1.model import (
        COMPLEX_REAC_UPDATE,
        COMPLEX_REAC_CREATE,
        COMPLEX_REAC_DELETE,
        COMPLEX_FILTER_MAX_VALUE,
    )
except ImportError:
    # If the installed STEPS version does not support complexes
    COMPLEX_REAC_UPDATE, COMPLEX_REAC_CREATE, COMPLEX_REAC_DELETE, COMPLEX_FILTER_MAX_VALUE = 0, 0, 0, 0
    smodel.Complex = None
    smodel.ComplexReac = None
    smodel.stepslib._py_Complex = None
    smodel.stepslib._py_ComplexReac = None

__all__ = [
    'Model',
    'VolumeSystem',
    'SurfaceSystem',
    'Species',
    'Complex',
    'Channel',
    'Reaction',
    'Diffusion',
    'Current',
    'OhmicCurr',
    'GHKCurr',
    'Location',
    'In',
    'Out',
    'Surf',
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
]

###################################################################################################
# Model


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
            for vsys in obj.getAllVolsyss():
                VolumeSystem._FromStepsObject(vsys, mdl)
            for ssys in obj.getAllSurfsyss():
                SurfaceSystem._FromStepsObject(ssys, mdl)
            for chan in obj.getAllChans():
                Channel._FromStepsObject(chan, mdl)
        return mdl

    def _addVolSysConstraint(self, reac, volSys, surfSys, loc):
        self.volSysConstraints.append((reac, volSys, surfSys, loc))

    def _createStepsObj(self):
        return smodel.Model()


###################################################################################################
# Space systems


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
        return smodel.Volsys(self.name, mdl.stepsModel)

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
        return vsys


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
        return smodel.Surfsys(self.name, mdl.stepsModel)

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
            for ghkc in obj.getAllGHKcurrs():
                GHKCurr._FromStepsObject(ghkc, mdl)
            for ohmc in obj.getAllOhmicCurrs():
                OhmicCurr._FromStepsObject(ohmc, mdl)

        return ssys


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

    def __init__(self, subUnits, *args, statesAsSpecies=False, order=NoOrdering, **kwargs):
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
        if not statesAsSpecies and smodel.Complex is None:
            warnings.warn(
                f'Your version of STEPS does not support multi-state complexes, '
                f'Complex {self.name} will instead be declared with '
                f'statesAsSpecies = True.'
            )
            statesAsSpecies = True

        self._order = order

        (mdl,) = self._getUsedObjects()
        self._statesAsSpecies = statesAsSpecies
        self._compStates = {}
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
        return smodel.Complex(self.name, mdl.stepsModel, len(self._subUnits), len(self._poolElems))

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
        for comb in itertools.product(*[su._states for su in self._subUnits]):
            cs = ComplexState(self, comb)
            if cs not in self._compStates:
                spec = smodel.Spec(cs.name, mdl.stepsModel)
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
            Comp[S1A, :, :]     # Valid, returns a ComplexSelector in which the first SubUnit is in
                                # state S1A and the other subunits are in any state
            Comp[S1A, ...]      # Valid, the ellipsis operator ('...') is equivalent to filling the
                                # remaining subUnits with ':'
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

    def _createStepsStates(self, mdl):
        """Create steps ChanState objects for each channel state."""
        chan = smodel.Chan(self.name, mdl.stepsModel)
        # Compute the channel state from the subunit states
        for comb in itertools.product(*[su._states for su in self._subUnits]):
            cs = ComplexState(self, comb)
            if cs not in self._compStates:
                chanState = smodel.ChanState(cs.name, mdl.stepsModel, chan)
                cs._setStepsObjects(chanState)
                self._compStates[cs] = chanState
        return chan


###################################################################################################
# Reaction elements


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
            s = self.stoich if self.stoich > 1 else ''
            return f'{self.loc}({s}{self._elem})'

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


ALL_LOCATIONS = [None, Location.IN, Location.SURF, Location.OUT]


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


####################################


class Species(nutils.UsingObjects(Model), ReactionElement, nutils.StepsWrapperObject):
    """A chemical species

    Species can be involved in :py:class:`Reaction`, :py:class:`Diffusion`, or other transport
    mechanisms.

    :param valence: Optional, the valence of the species (defaults to 0)
    :type valence: int
    """

    def __init__(self, *args, valence=0, _createObj=True, **kwargs):
        super().__init__(*args, **kwargs)
        (mdl,) = self._getUsedObjects()
        self.stepsSpecies = smodel.Spec(self.name, mdl.stepsModel) if _createObj else None
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
        return spec

    @property
    def valence(self):
        """The valence of the species (defaults to 0)

        :type: int
        """
        return self.stepsSpecies.getValence()

    @valence.setter
    def valence(self, val):
        return self.stepsSpecies.setValence(val)

    def __repr__(self):
        return self.name


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

    def _toUnorderedFilter(self):
        """
        Return a filter vector representing the complex state
        """
        vec = [0] * (2 * len(self._comp._poolElems))
        for sus in self:
            vec[2 * self._comp._sus2PoolInd[sus]] += 1
            vec[2 * self._comp._sus2PoolInd[sus] + 1] += 1
        return vec

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

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        if self._areStatesAsSpecies():
            return ''
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


class _ComplexReactionElement(ReactionElement):
    def _areStatesAsSpecies(self):
        pass


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
            for sus in self._complex._poolElems:
                minv, maxv, realMax = 0, 0, 0
                for ss in row:
                    if sus in ss:
                        maxv += 1
                        if len(ss._states) == 1:
                            minv += 1
                    if sus in ss._subUnit._states:
                        realMax += 1
                totvec += [minv, maxv if maxv < realMax else COMPLEX_FILTER_MAX_VALUE]
        return totvec

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
        allStates = set()
        for lss in self._subSels:
            for comb in itertools.product(*[list(sus._states) for sus in lss]):
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
                allStates.add(
                    (
                        ComplexState(self._complex, comb, _parentCompSel=self),
                        ComplexState(rcs._complex, destComb, _parentCompSel=rcs),
                        1,
                    )
                )
        return allStates

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
        comp = self._complex
        n = len(comp._poolElems)

        totFilt = []
        totUpd = []
        for lrow in self._subSels:
            for rrow in rcs._subSels:
                # For each pair of SubUnitSelector rows
                lds, rds = [], []
                # Do not consider subUnits that are unchanged and unconstrained
                for lss, rss in zip(lrow, rrow):
                    if not lss._isFull() or not rss._isFull():
                        lds.append(list(lss._states))
                        rds.append(list(rss._states))

                # For each possible state of the considered SubUnits
                for comb in itertools.product(*lds):
                    nf = []
                    for i in range(n):
                        nf += [0, COMPLEX_FILTER_MAX_VALUE]

                    nu = [0] * n
                    for lsus, rss in zip(comb, rds):
                        # Compute the filt and upd vects
                        nf[2 * comp._sus2PoolInd[lsus]] += 1
                        if len(rss) == 1:
                            rsus = rss[0]
                            nu[comp._sus2PoolInd[lsus]] += -1
                            nu[comp._sus2PoolInd[rsus]] += 1
                        elif lsus not in rss:
                            # Should never happen
                            raise Exception()
                    # Append to the final list
                    totFilt.append(nf)
                    totUpd.append(nu)

        # Remove duplicates or more specific
        def func(a, b):
            return all(
                ra == rb and a[0][2 * i] >= b[0][2 * i] for i, ra, rb in zip(range(len(a[2])), a[2], b[2])
            )

        totReacts = [[0] * n] * len(totFilt)
        return ComplexSelector._removeMoreSpecific(list(zip(totFilt, totReacts, totUpd)), func)

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
            return ''
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
        reactantPairs = sorted(reactantPairs, key=lambda x: len(x[0]._states))

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
                for comb in itertools.product(*[list(sus._states) for sus in src]):
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
            # allDstStates = set()
            # for srcos, dstStates in pairs.items():
            # allDstStates |= dstStates
            # totSubPairs = len(allDstStates)
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
        # TODO Not urgent: Might need to update this code to take into account the change in how 
        # complex selectors interact with subunit reactants, subunit reactants can now match complex
        # selectors part that are not fully undetermined.

        comp = self._filter._complex
        n = len(comp._poolElems)
        allFilts = self._filter._toUnorderedFilter()
        vects = []
        for i in range(0, len(allFilts), 2 * n):
            vects.append(allFilts[i : (i + 2 * n)])

        allVects = []

        for filt in vects:
            allReactants = [[0] * n]
            allUpdates = [[0] * n]
            # Get matching reactants
            for lhsr, rhsr in self._getMatchingReactants(rcs):
                newReactants = []
                newUpdates = []
                for lsus in lhsr._states:
                    for reac, upd in zip(allReactants, allUpdates):
                        nr = copy.copy(reac)
                        nu = copy.copy(upd)
                        ind = comp._sus2PoolInd[lsus]

                        nr[ind] += 1
                        nu[ind] += -1
                        if len(rhsr._states) == 1:
                            rsus = next(iter(rhsr._states))
                        elif lsus in rhsr._states:
                            rsus = lsus
                        else:
                            # This should never happen
                            raise Exception()
                        nu[comp._sus2PoolInd[rsus]] += 1

                        newReactants.append(nr)
                        newUpdates.append(nu)

                allReactants = newReactants
                allUpdates = newUpdates

            # Update the filters with the reactants
            for reac, upd in zip(allReactants, allUpdates):
                newFilt = []
                for i, r in enumerate(reac):
                    newFilt += [filt[2 * i] + r, filt[2 * i + 1]]
                allVects.append((newFilt, reac, upd))

        # Remove duplicates or more specific
        def func(a, b):
            return all(
                ra == rb and a[0][2 * i] >= b[0][2 * i] for i, ra, rb in zip(range(len(a[2])), a[2], b[2])
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
            raise SyntaxError('Cannot use "|" operator in reactions except for combining ' 'SubUnitStates.')

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
        return all(s.loc != Location.SURF for s in self.reacElems)

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
        specElems = list(self._getElemsOfType(Species))
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
            if compSel not in compSels2subSels:
                compSels2subSels[compSel] = ReactionSide([])
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


class _SubReaction(nutils.NamedObject):
    """
    Wrapper for a single steps reaction object.
    """

    _elemStrDict = {
        smodel.Reac: 'Reac',
        smodel.SReac: 'SReac',
        smodel.VDepSReac: 'VDepSReac',
        smodel.ComplexReac: 'ComplexReac',
        smodel.stepslib._py_Reac: 'Reac',
        smodel.stepslib._py_SReac: 'SReac',
        smodel.stepslib._py_VDepSReac: 'VDepSReac',
        smodel.stepslib._py_ComplexReac: 'ComplexReac',
    }

    def __init__(self, stepsReac, lhs=None, rateMult=1, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._stepsReac = stepsReac
        cls = stepsReac.__class__
        if cls not in _SubReaction._elemStrDict:
            raise TypeError(f'Unknown steps reaction type {cls}.')
        self._elemStr = _SubReaction._elemStrDict[cls]
        self._lhs = lhs
        self._rateMult = rateMult

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self._stepsReac]

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return self._elemStr

    def _solverSetValue(self, v):
        """
        Return the value that should actually be set in the solver when value 'v' is given by
        the user.
        """
        if isinstance(v, CompDepFunc):
            if self._lhs is None:
                raise NotImplementedError()
            else:
                return self._rateMult * v(self._lhs)
        else:
            return self._rateMult * v


class _SubReactionList(nutils.SolverPathObject, list):
    """
    Represents a unidirectional reaction consisting in a list of single steps reaction objects.
    """

    def __init__(self, parent, reacName, LRP, compEvs, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parent = parent
        self._reacName = reacName
        self._LRP = LRP
        self._compEvs = compEvs

        # Properties
        self._rate = 0

        self._reacElems = {loc: set() for loc in ALL_LOCATIONS}

    def append(self, item):
        if not isinstance(item, _SubReaction):
            raise TypeError('Only _SubReaction objects can be added to _SubReactionList.')
        super().append(item)

    @property
    def K(self):
        """Reaction rate(s)

        If the reaction is bidirectional, this value corresponds to a tuple composed of the forward
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

        :type: Union[float, :py:class:`XDepFunc`]
        """
        return self._rate

    @K.setter
    def K(self, rate):
        if not isinstance(rate, (numbers.Number, XDepFunc)):
            raise TypeError(f'{rate} cannot be used as a rate.')
        self._rate = rate

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        res = []
        for sr in self:
            res += sr._getStepsObjects()
        return res

    def _simPathWalkExpand(self):
        """Return an iterable of the elements that should be part of ResultSelector paths."""
        return self

    def _getAllElems(self, loc):
        return self._reacElems[loc]

    def _addToSteps(self):
        """Add all the reactions to steps."""
        """
        Declare and add all steps reactions corresponding to the (lhs, rhs) pairs in LRP.
        """
        mdl, globCompSels, volSys, surfSys = self._parent._usedObjects

        for i, reac in enumerate(self._LRP):
            l, r, rm = reac
            # Add reac elems
            for loc in ALL_LOCATIONS:
                self._reacElems[loc] |= l._getAllElems(loc) | r._getAllElems(loc)
            # Try to declare complex reac, if applicable
            for j, compEvents in enumerate(self._compEvs):
                # Add complex reac elems
                for ce in compEvents:
                    # TODO Later release: Modify once surface complex reacs are implemented
                    for loc in ALL_LOCATIONS:
                        self._reacElems[loc].add(ce[1])

                # Declare reaction
                name = self._reacName + (f'_{i}' if len(self._LRP) > 1 else '')
                if len(compEvents) > 0:
                    name_c = name + (f'_{j}' if len(self._compEvs) > 1 else '')
                    res = self._declareComplexReac(
                        name_c, l, r, compEvents, mdl, volSys, surfSys, rm * self._rate, rm
                    )
                else:
                    res = self._declareSimpleReac(name, l, r, mdl, volSys, surfSys, rm * self._rate, rm)

                # Add the subreactions to self
                self.append(res)

    def _declareComplexReac(self, name, simpLHS, simpRHS, compEvents, mdl, volSys, surfSys, rate, rateMult):
        """Declare complex reaction involving real complexes."""

        # TMP Only implement volume reactions for now
        lhs = simpLHS._GetStepsElems()
        rhs = simpRHS._GetStepsElems()
        compEvs = []
        for ev in compEvents:
            compEvs.append(ev[0:1] + (ev[1].name,) + ev[2:])

        if isinstance(rate, CompDepRate):
            raise NotImplementedError(
                f'{self._parent._getDeclarationDescr()}: Cannot use CompDepRate with STEPS complex '
                f'reactions.'
            )

        nutils._print(f'Adding STEPS complex reaction {name}.', 2)

        stepsReac = smodel.ComplexReac(name, volSys.stepsSys, lhs=lhs, rhs=rhs, compEvs=compEvs, kcst=rate)
        return _SubReaction(stepsReac, rateMult=rateMult, name=name)
        # END TMP

    def _declareSimpleReac(self, name, lhs, rhs, mdl, volSys, surfSys, rate, rateMult):
        """Declare and return the steps reactions from the final sides lhs and rhs."""

        nutils._print(f'Adding STEPS reaction {name}: {lhs} --> {rhs}, rate = {rate}', 3)

        # Compute rates in case of complex dependent rates
        if isinstance(rate, CompDepRate):
            rate = rate(lhs)

        if self._parent._isSurfaceReac():
            # Surface system reaction
            ilhs = lhs._GetStepsElems(Location.IN)
            slhs = lhs._GetStepsElems(Location.SURF)
            olhs = lhs._GetStepsElems(Location.OUT)
            irhs = rhs._GetStepsElems(Location.IN)
            srhs = rhs._GetStepsElems(Location.SURF)
            orhs = rhs._GetStepsElems(Location.OUT)
            if isinstance(rate, VDepRate):
                stepsReac = smodel.VDepSReac(
                    name,
                    surfSys.stepsSys,
                    ilhs=ilhs,
                    slhs=slhs,
                    olhs=olhs,
                    irhs=irhs,
                    srhs=srhs,
                    orhs=orhs,
                    k=rate._func,
                    vrange=rate._vrange,
                )
            else:
                stepsReac = smodel.SReac(
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
            return _SubReaction(stepsReac, lhs=lhs, rateMult=rateMult, name=name)
        else:
            # Volume system reaction
            elhs = lhs._GetStepsElems()
            erhs = rhs._GetStepsElems()
            if isinstance(rate, VDepRate):
                # Undeclare the parent reaction to avoid further exception raising when exiting context
                # managers
                self._parent._declared = False
                raise NotImplementedError(
                    f'{self._parent._getDeclarationDescr()}: Cannot declare a voltage-dependent volume '
                    f'reaction.'
                )
            else:
                stepsReac = smodel.Reac(name, volSys.stepsSys, lhs=elhs, rhs=erhs, kcst=rate)
                return _SubReaction(stepsReac, lhs=lhs, rateMult=rateMult, name=name)



class Reaction(
    nutils.UsingObjects((Model, ComplexSelector), (VolumeSystem, SurfaceSystem)), nutils.StepsWrapperObject
):
    """A (possibly bidirectional) reaction.

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

        lhs, rhs = ReactionSide([]), ReactionSide([])
        if isinstance(obj, smodel.stepslib._py_Reac):
            for lhsSpec in obj.getLHS():
                lhs += getattr(mdl, lhsSpec.getID())
            for rhsSpec in obj.getRHS():
                rhs += getattr(mdl, rhsSpec.getID())
        elif isinstance(obj, smodel.stepslib._py_SReac):
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
        else:
            raise NotImplementedError(f'Cannot import from STEPS reaction {obj}')

        lhs > reac > rhs
        reac.K = obj.getKcst()
        subreacLst = _SubReactionList(reac, obj.getID(), [(lhs, rhs, 1)], None)
        subreacLst.append(_SubReaction(obj, lhs=lhs, rateMult=1.0, name=obj.getID()))
        reac._subReactions[Reaction._FwdSpecifier] = subreacLst
        reac._added = True
        for loc in ALL_LOCATIONS:
            subreacLst._reacElems[loc] |= lhs._getAllElems(loc) | rhs._getAllElems(loc)
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
        specLhs = lhs._getElemsOfType(Species, ComplexState)
        specRhs = rhs._getElemsOfType(
            Species, ComplexState
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
                    (COMPLEX_REAC_UPDATE, comp, filt, reac, upd)
                    for filt, reac, upd in cl._elem._getUpdateVects(cr)
                ]
                del crhs[crInd]
            else:
                if isinstance(cl, _ComplexReactants):
                    raise Exception(
                        f'ComplexReactants {cl} has no matching ComplexReactants on '
                        f'the RHS of the reaction.'
                    )

                events = [(COMPLEX_REAC_DELETE, comp, cl._elem._toUnorderedFilter())]

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
                evs.append((COMPLEX_REAC_CREATE, state._comp, state._toUnorderedFilter()))

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

        self._subReactions[Reaction._FwdSpecifier] = _SubReactionList(self, nameFwd, fwdLRP, fwdCompEvs)

        if self._bidir:
            bkwLRP = self._expandSimpleReactions(rhs, lhs)
            bkwCompEvs = self._expandComplexReactions(rhs, lhs)
            nameBkw = self.name + Reaction._BkwSpecifier

            self._subReactions[Reaction._BkwSpecifier] = _SubReactionList(self, nameBkw, bkwLRP, bkwCompEvs)

        # Retrieve filename and line number for easier debugging
        frameInfo = inspect.stack()[3]
        self._declarationInfo = (frameInfo.filename, frameInfo.lineno)

        self._declared = True

    def _comparaison(self, other, gt):
        if self._declared:
            if self._managerInfo is None:
                raise Exception('This reaction has already been declared.')
            # Return the new reaction
            # TODO not urgent: Manager info being a tuple hurts readability, make it a namedtuple
            p = self._managerInfo[0]._popNewReacArgs()
            if p is None:
                self._managerInfo[0]._setNewReacArgs(other, gt)
                # python exhibits special behavior for chained comparaison operator
                # see https://docs.python.org/3/reference/expressions.html#comparisons
                # The first comparaison needs to return True
                return True
            else:
                r = self._managerInfo[0]._getNew(self._managerInfo[1])
                r._comparaison(*p)
                return r._comparaison(other, gt)

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
        return self._comparaison(other, True)

    def __lt__(self, other):
        return self._comparaison(other, False)

    def _getAllElems(self, loc):
        if self._subReactions[Reaction._FwdSpecifier] is not None:
            res = self._subReactions[Reaction._FwdSpecifier]._getAllElems(loc)
            if self._subReactions[Reaction._BkwSpecifier] is not None:
                res |= self._subReactions[Reaction._BkwSpecifier]._getAllElems(loc)
            return res
        return set()

    def _isSurfaceReac(self):
        if self.lhs is None or self.rhs is None:
            return None
        lhsS, rhsS = self.lhs._isSurfaceSide(), self.rhs._isSurfaceSide()
        return (lhsS and rhsS) or (None in [lhsS, rhsS] and True in [lhsS, rhsS])

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

        See :ref:`/API_2/Interface_Tutorial_1_wm.ipynb` for usage example.

        :meta public:
        """
        if key not in self._subReactions:
            raise KeyError(f'No reaction can be accessed with key {key}.')
        return self._subReactions[key]

    def __repr__(self):
        return '{} {}-> {}'.format(self.lhs, '-' if not self._bidir else '<', self.rhs)

    def _getDeclarationDescr(self):
        return f'{self._declarationInfo[0]}: {self._declarationInfo[1]}'

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

    def _exitCallback(self, parent):
        """Method to be called when we get out of a context manager."""
        if self._declared and not self._added:
            nutils._print(f'Adding pysteps reaction {self.name}: {self}', 2)

            self._subReactions[Reaction._FwdSpecifier]._addToSteps()
            if self._bidir:
                self._subReactions[Reaction._BkwSpecifier]._addToSteps()
            self._added = True

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
                        f'The reaction is bidirectional, two properties should be set, got {val} instead.'
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


class AnonymousReaction(Reaction):
    """
    :meta private:
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class ReactionManager:
    """Class used to instantiate reactions

    A :py:class:`ReactionManager` object allows to easily declare reactions by associating them
    to an identifier using the square bracket syntax::

        r = ReactionManager()  # Only one ReactionManager is needed

        S1 + S2 >r['R01']> S3  # String identifier for permanently naming the reaction

        S1 + S2 >r[1]> S3      # Integer identifier for reactions that do not need explicit naming

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


class XDepFunc:
    """Base class for functions depending on some simulation or object state

    This is a wrapper around a function that adds some functionalities.

    :param func: lambda function or standard function. There is not support for variable named
        arguments (\*\*kwargs).
    :type func: Callable[[Any, ...], Number]
    """

    def __init__(self, func, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not hasattr(func, '__call__'):
            raise TypeError(f'Expected a callable, got {func} instead.')
        self._func = func
        self._args = []

    def __mul__(self, m):
        """Multiplication with a constant or another function of the same class

        :param m: Constant number or another :py:class:`XDepFunc` of the same type
        :type m: Union[Number, :py:class:`XDepFunc`]

        :returns: A new function, which corresponds to the multiplication of both operands
        :rtype: :py:class:`XDepFunc`
        :meta public:
        """
        if isinstance(m, numbers.Number):
            return self.__class__(lambda *lst: m * self._func(*lst), *self._args)
        elif isinstance(m, self.__class__):
            return self.__class__(lambda *lst: m(*lst) * self._func(*lst), *self._args)
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

    def __call__(self, *lst):
        """Call the wrapped function

        :meta public:
        """
        return self._func(*lst)


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

    def __init__(self, func, vrange=None, *args, **kwargs):
        super().__init__(func=func, *args, **kwargs)
        if vrange is None:
            vrange = [smodel._VoltageTable.vmin, smodel._VoltageTable.vmax, smodel._VoltageTable.dv]
        self._vrange = vrange

        self._args = [self._vrange]

    @property
    def vrange(self):
        """Access the voltage range

        :type: Tuple[float, float, float], read-only
        """
        return self._vrange


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

    pass


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
            # The above reaction takes place between a subunit of Comp1 in state S1A and a subunit
            # of Comp2 in state S3A. It takes place with a basal rate k3 and this rate is increased
            # by a value that depends on the number of subunits in state S2A in Comp1 and the
            # number of subunits in state S4A in Comp2.

    """

    pass


class CompDepCond(CompDepFunc):
    """Complex state dependent conductance

    See base class :py:class:`CompDepFunc` for how to create.
    """

    pass


class CompDepP(CompDepFunc):
    """Complex state dependent permeability

    See base class :py:class:`CompDepFunc` for how to create.
    """

    pass


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

    def _solverSetValue(self, v):
        """
        Return the value that should actually be set in the solver when value 'v' is given by
        the user.
        """
        if isinstance(v, CompDepDcst):
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


class Diffusion(nutils.UsingObjects(Model, (VolumeSystem, SurfaceSystem)), nutils.StepsWrapperObject):
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
        self.sys = volSys if surfSys is None else surfSys
        self._stepsDiffs = {}
        self._elem = elem._getReferenceObject()
        self._Dcst = Dcst

        if _createObj:
            if (
                isinstance(elem, (Complex, ComplexSelector))
            ) and elem._areStatesAsSpecies():
                for state in elem:
                    name = f'{self.name}_{state.name}'
                    self._stepsDiffs[state] = _SubDiffusion(
                        smodel.Diff(name, self.sys.stepsSys, state._getStepsObjects()[0]),
                        self,
                        keyElem=state,
                        name=name,
                    )
            elif isinstance(elem, (Species, ComplexState)):
                self._stepsDiffs[elem] = _SubDiffusion(
                    smodel.Diff(self.name, self.sys.stepsSys, elem._getStepsObjects()[0]),
                    self,
                    keyElem=elem,
                    name=self.name,
                )
            else:
                raise TypeError(f'Cannot declare a diffusion rule for {elem}.')

            # Setting diffusion constant
            if isinstance(Dcst, numbers.Number):
                for _, sd in self._stepsDiffs.items():
                    sd._stepsDiff.setDcst(Dcst)
            elif isinstance(Dcst, CompDepFunc):
                for cs, sd in self._stepsDiffs.items():
                    sd._stepsDiff.setDcst(Dcst(cs))
            else:
                raise TypeError(f'{Dcst} is not a valid diffusion constant.')

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


class Current(nutils.UsingObjects(SurfaceSystem), nutils.StepsWrapperObject):
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

    def _getAllElems(self, loc):
        return (
            list(set(state._comp for state, curr in self._currents.items()))
            + [state for state, curr in self._currents.items()]
            + [curr for state, curr in self._currents.items()]
        )

    def _getComplexStates(self, states):
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
        return [obj for state, obj in self._currents.items()]


class OhmicCurr(Current):
    """An ohmic current through a channel in one or more conducting states

    :param states: The conducting state(s) of the channel
    :type states: Union[:py:class:`ComplexState`, :py:class:`ComplexSelector`]
    :param conduct: A constant conductance in Siemens or a :py:class:`CompDepCond` function that
        returns a conductance in Siemens as a function of a complex state.
    :type conduct: Union[float, :py:class:`CompDepCond`]
    :param rev_pot: Reversal potential in Volts
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

            # The above 2 lines declare an ohmic current from channel Chan1 when its 3rd subunit is
            # in state S2A. Thanks to the CompDepCond function, the conductance changes depending
            # on the specific states of the other subunits. In this case, there is a basal
            # conductance k1 and this conductance is increased by k2 for each subunit in state S1A.

    """

    _currStr = 'Ohmic'

    def __init__(self, states, conduct, rev_pot, *args, **kwargs):
        super().__init__(*args, **kwargs)
        (surfSys,) = self._getUsedObjects()

        self._complex, states = self._getComplexStates(states)
        if isinstance(conduct, CompDepFunc):
            if len(conduct._complexes) != 1 or conduct._complexes[0] is not self._complex:
                raise Exception('The CompDepFunc used as conductance is not compatible with the channel.')
            condStates = [(s, conduct(s)) for s in states]
        else:
            condStates = [(s, conduct) for s in states]

        for state, g in condStates:
            chanState = self._complex._compStates[state]
            self._currents[state] = _SubCurrent(
                self,
                smodel.OhmicCurr(
                    f'{self.name}_{state.name}', surfSys.stepsSys, chanstate=chanState, g=g, erev=rev_pot
                ),
            )


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
        compartment does not have to be explicitly simulated (if it retains default negative
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

    class _GHKCurrPInfo:
        """Simple wrapper for providing Pinfo instead of permeability"""

        def __init__(self, *args):
            self.args = args

    def __init__(self, states, spec, P, computeflux=True, virtual_oconc=-1, vshift=0, *args, **kwargs):
        super().__init__(*args, **kwargs)
        (surfSys,) = self._getUsedObjects()
        if not isinstance(spec, Species):
            raise TypeError(f'Expected a Species, got {spec} instead.')
        self._spec = spec

        self._complex, states = self._getComplexStates(states)
        if isinstance(P, CompDepFunc):
            if len(P._complexes) != 1 or P._complexes[0] is not self._complex:
                raise Exception('The CompDepFunc used as permeability is not compatible with the channel.')
            condStates = [(s, P(s)) for s in states]
        else:
            condStates = [(s, P) for s in states]

        for state, perm in condStates:
            chanState = self._complex._compStates[state]
            self._currents[state] = _SubCurrent(
                self,
                smodel.GHKcurr(
                    f'{self.name}_{state.name}',
                    surfSys.stepsSys,
                    chanState,
                    self._spec.stepsSpecies,
                    computeflux=computeflux,
                    virtual_oconc=virtual_oconc,
                    vshift=vshift,
                ),
            )
            if isinstance(perm, numbers.Number):
                self._currents[state].stepsCurrent.setP(perm)
            elif isinstance(perm, GHKCurr._GHKCurrPInfo):
                self._currents[state].stepsCurrent.setPInfo(*perm.args)
            else:
                raise TypeError(f'Expected a permeability, got {perm} instead.')

    @classmethod
    def PInfo(cls, g, V, T, oconc, iconc):
        """Provide information from channel measurement instead of permeability

        :param g: A measured single-channel conductance in Siemens
        :type g: float
        :param V: The potential in Volts
        :type V: float
        :param T: The temperature in Kelvin (may be different from the STEPS simulation
            temperature)
        :type T: float
        :param oconc: The 'outer' concentration of the ion in molar units
        :type oconc: float
        :param iconc: The 'inner' concentration of the ion in molar units
        :type iconc: float

        :returns: An aggregate object. Does not return the permeability.

        Usage::

            Pinfos = GHKCurr.PInfo(g=20e-12, V=-22e-3, T=293.15, oconc=4e-3, iconc=155e-3)

            GC_CaP = GHKCurr.Create(CaPchan[m3], Ca, Pinfos, computeflux = True)
        """
        return GHKCurr._GHKCurrPInfo(g, V, T, oconc, iconc)

    def _getAllElems(self, loc):
        if loc != Location.SURF:
            return super()._getAllElems(loc) + [self._spec]
        else:
            return super()._getAllElems(loc)
