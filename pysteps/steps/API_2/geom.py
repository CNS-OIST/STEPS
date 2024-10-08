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

import contextlib
import copy
import functools
import numbers
import operator
import pickle
import types
import warnings

import numpy

from steps import stepslib

import steps.API_1.utilities.meshio as smeshio
import steps.API_1.utilities.geom_decompose as sgdecomp
import steps.API_1.utilities.metis_support as smetis
import steps.API_1.utilities.morph_support as smorph

from . import utils as nutils
from . import model as nmodel

__all__ = [
    'Geometry',
    'Compartment',
    'Patch',
    'Membrane',
    'EndocyticZone',
    'ROI',
    'TetMesh',
    'DistMesh',
    'DiffBoundary',
    'SDiffBoundary',
    'MeshPartition',
    'LinearMeshPartition',
    'MetisPartition',
    'GmshPartition',
    'MorphPartition',
    'Reference',
    'TetReference',
    'TriReference',
    'BarReference',
    'VertReference',
    'RefList',
    'TetList',
    'TriList',
    'BarList',
    'VertList',
    'Point',
    'BoundingBox',
    'Morph',
]


###################################################################################################

ELEM_VERTEX     = stepslib._py_ElementType.ELEM_VERTEX
ELEM_TRI        = stepslib._py_ElementType.ELEM_TRI
ELEM_TET        = stepslib._py_ElementType.ELEM_TET
ELEM_UNDEFINED  = stepslib._py_ElementType.ELEM_UNDEFINED

UNKNOWN_TET     = stepslib.UNKNOWN_TET
UNKNOWN_TRI     = stepslib.UNKNOWN_TRI
UNKNOWN_VERT    = stepslib.UNKNOWN_VERT
INDEX_NUM_BYTES = stepslib.INDEX_NUM_BYTES
INDEX_DTYPE     = numpy.uint32 if INDEX_NUM_BYTES == 4 else numpy.uint64

###################################################################################################


@nutils.FreezeAfterInit
class Geometry(nutils.UsableObject, nutils.StepsWrapperObject):
    """Top-level geometry container

    A number of compartment objects and patches objects may be grouped.
    Should be used as a context manager for the declaration of compartments etc::

        geom = Geometry()
        with geom:
            comp1 = Compartment.Create(vsys, volume)

            ... # Declare other objects in geom

    After having declared children objects in the ``with mesh:`` block, they can be accessed
    as attributes of ``geom`` with their name (see :py:class:`steps.API_2.utils.NamedObject` and
    :py:func:`steps.API_2.utils.NamedObject.Create`)::

        geom.comp1
    """

    def __init__(self, *args, _createObj=True, **kwargs):
        super().__init__(*args, **kwargs)
        self.stepsGeom = stepslib._py_Geom() if _createObj else None

        self._local = False
        self._owned = None
        self._callKwargs = {}
        self._lstArgs = {}

    def _SetUpMdlDeps(self, mdl):
        """Set up structures that depend on objects declared in the model."""
        # Start with patches because they can update compartments
        for pl in self._getChildrenOfType(Patch):
            pl._SetUpMdlDeps(mdl)
        for pl in self._getChildrenOfType(Membrane):
            pl._SetUpMdlDeps(mdl)
        for pl in self._getChildrenOfType(Compartment):
            pl._SetUpMdlDeps(mdl)

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsGeom]

    @classmethod
    def _FromStepsObject(cls, obj):
        """Create the interface object from a STEPS object."""
        geom = cls(_createObj=False)
        geom.stepsGeom = obj
        with geom:
            for comp in obj.getAllComps():
                Compartment._FromStepsObject(comp, geom)
            for patch in obj.getAllPatches():
                Patch._FromStepsObject(patch, geom)
        return geom


class _PhysicalLocation(nutils.UsingObjects(Geometry), nutils.StepsWrapperObject):
    """Base class for physical locations

    Compartment, Patch, Membrane, etc. derive from this class
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.systems = []
        self.sysNames = []
        (self._geom,) = self._getUsedObjects()

    def addSystem(self, sys, _loc=None):
        """Add a surface or a volume system to the location

        :param sys: The volume or surface system to be added, or its name.
        :type sys: :py:class:`steps.API_2.model.VolumeSystem`, :py:class:`steps.API_2.model.SurfaceSystem`\
        or `str`

        :returns: None
        """
        if isinstance(sys, nmodel.SpaceSystem):
            self.sysNames.append((sys.name, _loc))
        elif isinstance(sys, str):
            self.sysNames.append((sys, _loc))
        else:
            raise TypeError(f'Expected a VolumeSystem or a SurfaceSystem, got a {type(sys)} instead.')

    def _SetUpMdlDeps(self, mdl):
        """Set up the structures that will allow species or reaction access from locations."""
        for sname, loc in self.sysNames:
            s = getattr(mdl, sname)
            self.systems.append((s, loc))
            # Add the reactions and the corresponding reactants to children
            for name, c in s.children.items():
                if _PhysicalLocation._canBeChild(self, c):
                    if self._canBeChild(c):
                        self._addChildren(c)
                    for re in c._getAllElems(loc):
                        if re.name not in self.children:
                            self._addChildren(re)

            s._addLocation(self)

    def _canBeChild(self, c):
        """Return wether c can be a child of self."""
        return isinstance(c, (
            nmodel.Reaction, nmodel.Diffusion, nmodel.Current, nmodel.Complex, nmodel.VesicleBind,
            nmodel.VesicleUnbind, nmodel.Endocytosis, nmodel.RaftGen
        ))

    def _createStepsObj(self, geom):
        pass

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        return {'loc_type': self.__class__._locStr, 'loc_id': self.name}

    @property
    def _callKwargs(self):
        return self._geom._callKwargs

    @property
    def _local(self):
        return self._geom._local

    @property
    def _owned(self):
        return self._geom._owned


@nutils.FreezeAfterInit
class Compartment(_PhysicalLocation, nutils.ParameterizedObject, nutils.Facade):
    """Base class for compartment objects

    The same class is used to declare compartments in well-mixed models and in 3D tetrahedral
    meshes.

    Compartments in well-mixed geometry can be declared in the following ways::

        with geom:
            wmComp1 = Compartment.Create()
            wmComp2 = Compartment.Create(vsys)
            wmComp3 = Compartment.Create(vol=volume)
            wmComp4 = Compartment.Create(vsys, volume)

    :param vsys: Optional, the volume system associated with this compartment.
    :type vsys: :py:class:`steps.API_2.model.VolumeSystem` or :py:class:`str`
    :param vol: Optional, volume of the compartment
    :type vol: float

    Compartments in a tetrahedral mesh can be declared in the following ways::

        with mesh:
            tetComp1 = Compartment.Create(tetLst)
            tetComp2 = Compartment.Create(tetLst, vsys)

    :param tetLst: List of tetrahedrons associated with the compartment.
    :type tetLst: :py:class:`TetList` or any argument that can be used to build one
    :param vsys: Optional, the volume system associated with this compartment.
    :type vsys: :py:class:`steps.API_2.model.VolumeSystem` or :py:class:`str`

    The volume of a compartment in a tetrahedral mesh is the total volume of the encapsulated
    tetrahedrons.

    Compartments in a distributed tetrahedral mesh can be declared in the following ways::

        with mesh:
            distComp1 = Compartment.Create()
            distComp2 = Compartment.Create(vsys, conductivity=g)
            distComp3 = Compartment.Create(physicalTag=10)
            distComp4 = Compartment.Create(tetLst, vsys)

    :param vsys: Optional, the volume system associated with this compartment.
    :type vsys: :py:class:`steps.API_2.model.VolumeSystem` or :py:class:`str`
    :param conductivity: Optional, the volume conductivity of the compartment.
    :type conductivity: float
    :param physicalTag: Optional, the index of the physical tag used to build the compartment
        (see :py:class:`DistMesh`).
    :type physicalTag: int
    :param tetLst: Optional, List of tetrahedrons associated with the compartment. If it is ommited
        or an empty list is given, the compartment name will be used to find tetrahedrons tagged with it.
    :type tetLst: :py:class:`TetList` or any argument that can be used to build one

    Some of the methods documented below are only available for tetrahedral compartments, they are
    grouped accordingly.
    """

    _locStr = 'Comp'

    def __new__(cls, *args, **kwargs):
        # Compartment is a facade to its subclasses, the user always uses the Compartment class
        # but the appropriate subclass is used for object instantiation.

        if cls is not Compartment:
            return super(Compartment, cls).__new__(cls)

        geom = nutils.UsableObject._getUsedObjectOfClass(Geometry)
        if geom is None:
            raise Exception(f'Compartments need to be declared in a geometry.')

        if geom.__class__ == Geometry:
            return _WmCompartment.__new__(_WmCompartment, *args, **kwargs)
        elif geom.__class__ == TetMesh:
            return _TetCompartment.__new__(_TetCompartment, *args, **kwargs)
        elif geom.__class__ == DistMesh:
            return _DistCompartment.__new__(_DistCompartment, *args, **kwargs)
        else:
            raise NotImplementedError(f'Compartments cannot be created in {type(geom)}')

    @classmethod
    def _getDisplayName(cls):
        return Compartment.__name__

    @classmethod
    def _FromStepsObject(cls, obj, geom):
        if stepslib._STEPS_USE_DIST_MESH and isinstance(obj, stepslib._py_DistComp):
            return _DistCompartment._FromStepsObject(obj, geom)
        elif isinstance(obj, stepslib._py_TmComp):
            return _TetCompartment._FromStepsObject(obj, geom)
        elif isinstance(obj, stepslib._py_Comp):
            return _WmCompartment._FromStepsObject(obj, geom)
        else:
            raise NotImplementedError(f'Cannot create API2 version of {type(obj)}')

    def __init__(self, vsys, _createObj=True, **kwargs):
        super().__init__(**kwargs)

        self.stepsComp = self._createStepsObj(self._geom) if _createObj else None
        if vsys is not None:
            self.addSystem(vsys)

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsComp]

    def _createStepsObj(self, geom):
        raise NotImplementedError()

    def _SetUpMdlDeps(self, mdl):
        """Set up the structures that will allow species or reaction access from locations."""
        super()._SetUpMdlDeps(mdl)
        for ves in mdl._getChildrenOfType(nmodel.Vesicle):
            self._addChildren(ves)

    def _canBeChild(self, c):
        """Return wether c can be a child of self."""
        if not super()._canBeChild(c):
            return False
        if ((isinstance(c, nmodel.Reaction) and c._isSurfaceReac()) or
            (isinstance(c, nmodel.Diffusion) and c._isSurfaceDiff()) or
            isinstance(c, (nmodel.Endocytosis, nmodel.RaftGen))
        ):
            return False
        return True

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return Compartment._locStr

    def addSystem(self, vsys, _loc=None):
        """Add a volume system to the compartment

        :param vsys: The volume system to be added, or its name.
        :type vsys: :py:class:`steps.API_2.model.VolumeSystem` or `str`

        :returns: None
        """
        super().addSystem(vsys, _loc)
        if _loc is None:
            if isinstance(vsys, str):
                self.stepsComp.addVolsys(vsys)
            elif vsys.__class__ is nmodel.VolumeSystem:
                self.stepsComp.addVolsys(vsys.name)

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m^3'))
    @nutils.IgnoresGhostElements
    def Vol(self):
        """Volume of the compartment

        :type: Union[float, :py:class:`steps.API_2.utils.Parameter`], read-only for tetrahedral compartments
        """
        return self.stepsComp.getVol(**self._callKwargs)


@nutils.FreezeAfterInit
class _WmCompartment(Compartment):

    def __init__(self, vsys=None, vol=None, **kwargs):
        super().__init__(vsys=vsys, **kwargs)
        if vol is not None:
            self.Vol = vol

    def _createStepsObj(self, geom):
        return stepslib._py_Comp(self.name, geom.stepsGeom)

    @classmethod
    def _FromStepsObject(cls, obj, geom):
        """Create the interface object from a STEPS object."""
        comp = cls(_createObj=False, name=obj.getID())
        comp.stepsComp = obj
        # setup volume systems
        for vsysname in comp.stepsComp.getVolsys():
            comp.addSystem(vsysname)

        return comp

    @Compartment.Vol.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('m^3'))
    def Vol(self, v):
        self.stepsComp.setVol(v)


@nutils.FreezeAfterInit
class _BaseTetCompartment(Compartment):
    _FACADE_TITLE_STR = 'Only available for compartments defined in tetrahedral meshes'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def surface(self):
        """All triangles in the surface of the compartment

        :type: :py:class:`TriList`, read-only
        """
        return self.tets.surface

    @property
    def tets(self):
        """All tetrahedrons in the compartment

        :type: :py:class:`TetList`, read-only
        """
        raise NotImplementedError()

    @property
    @nutils.IgnoresGhostElements
    def bbox(self):
        """The bounding box of the compartment

        :type: :py:class:`steps.API_2.geom.BoundingBox`, read-only
        """
        return BoundingBox(
            Point(*self.stepsComp.getBoundMin(**self._callKwargs)),
            Point(*self.stepsComp.getBoundMax(**self._callKwargs))
        )


@nutils.FreezeAfterInit
class _TetCompartment(_BaseTetCompartment):
    _FACADE_TITLE_STR = None

    def __init__(self, tetLst, vsys=None, **kwargs):
        self._tetLst = tetLst
        super().__init__(vsys=vsys, **kwargs)

    def _createStepsObj(self, geom):
        self._tetLst = TetList._toRefList(self._tetLst, geom)
        return stepslib._py_TmComp(self.name, geom.stepsMesh, self._tetLst.indices)

    @classmethod
    def _FromStepsObject(cls, obj, geom):
        """Create the interface object from a STEPS object."""
        tetLst = TetList(obj.getAllTetIndices(), geom)

        comp = cls(tetLst, _createObj=False, name=obj.getID())
        comp.stepsComp = obj
        # setup volume systems
        for vsysname in comp.stepsComp.getVolsys():
            comp.addSystem(vsysname)

        return comp

    @_BaseTetCompartment.tets.getter
    def tets(self):
        return self._tetLst


@nutils.FreezeAfterInit
class _DistCompartment(_BaseTetCompartment):
    _FACADE_TITLE_STR = 'Only available for compartments defined in distributed tetrahedral meshes'

    def __init__(self, *args, conductivity=0, physicalTag=None, **kwargs):
        self._physicalTag = physicalTag

        # Allow any ordering of vsys and tetLst
        vsys, tetLst, args, kwargs = nutils.extractArgs(args, kwargs, [
            ('vsys', (nmodel.VolumeSystem, str), None),
            ('tetLst', (TetList, list, tuple), None),
        ])

        self._tetLst = tetLst
        super().__init__(vsys=vsys, **kwargs)

        self.Conductivity = conductivity

    def _createStepsObj(self, geom):
        if self._tetLst is None:
            comp = stepslib._py_DistComp(self.name, geom.stepsMesh, None, self._physicalTag)
        else:
            tetLst = TetList._toRefList(self._tetLst, geom)
            if tetLst.isLocal() and tetLst._owned != False:
                raise Exception(
                    f'The compartment is being created with a local tetrahedron list but without '
                    f'including non-owned elements. It should instead be created inside a with '
                    f'mesh.asLocal(owned=False) block.'
                )
            comp = stepslib._py_DistComp(
                self.name, geom.stepsMesh, tetLst.indices, self._physicalTag, local=tetLst.isLocal()
            )
        del self._tetLst
        return comp

    @_BaseTetCompartment.tets.getter
    def tets(self):
        return TetList(
            self.stepsComp.getAllTetIndices(**self._geom._callKwargs),
            mesh=self._geom, _immutable=True, **self._geom._lstArgs
        )

    @_BaseTetCompartment.surface.getter
    @nutils.IgnoresGhostElements
    def surface(self):
        return TriList(
            self.stepsComp.getSurfTris(**self._geom._callKwargs),
            mesh=self._geom, _immutable=True, **self._geom._lstArgs
        )

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('S m^-1'))
    def Conductivity(self):
        """Conductivity of the compartment

        :type: Union[float, :py:class:`steps.API_2.utils.Parameter`]
        """
        return self.stepsComp.getConductivity()

    @Conductivity.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('S m^-1'))
    def Conductivity(self, v):
        self.stepsComp.setConductivity(v)


@nutils.FreezeAfterInit
class Patch(_PhysicalLocation, nutils.UsableObject, nutils.ParameterizedObject, nutils.Facade):
    """Base class for patch objects

    A patch is a piece of 2D surface surrounding (part of) a 3D compartment, which may be
    connected to another compartment. The same class is used to declare patches in well-mixed
    models and in 3D tetrahedral meshes.

    Patches in well-mixed geometry can be declared in the following ways::

        with geom:
            patch1 = Patch.Create(inner)
            patch2 = Patch.Create(inner, outer)
            patch3 = Patch.Create(inner, outer, ssys)
            patch4 = Patch.Create(inner, None , ssys)
            patch5 = Patch.Create(inner, None , None, area)
            patch6 = Patch.Create(inner, outer, None, area)
            patch7 = Patch.Create(inner, outer, ssys, area)

    :param inner: Inner compartment, see below
    :type inner: :py:class:`steps.API_2.geom.Compartment`
    :param outer: Optional, Outer compartment, see below
    :type outer: :py:class:`steps.API_2.geom.Compartment`
    :param ssys: Optional, the surface system associated with this patch.
    :type ssys: :py:class:`steps.API_2.model.SurfaceSystem` or :py:class:`str`
    :param area: Optional, surface area of the patch
    :type area: float

    Patches in a tetrahedral mesh can be declared in the following ways::

        with mesh:
            patch1 = Patch.Create(triLst, inner)
            patch2 = Patch.Create(triLst, inner, outer)
            patch3 = Patch.Create(triLst, inner, None , ssys)
            patch4 = Patch.Create(triLst, inner, outer, ssys)

    :param triLst: List of triangles associated with the patch.
    :type triLst: :py:class:`TriList` or any argument that can be used to build one
    :param inner: Inner compartment, see below
    :type inner: :py:class:`steps.API_2.geom.Compartment`
    :param outer: Optional, Outer compartment, see below
    :type outer: :py:class:`steps.API_2.geom.Compartment`
    :param ssys: Optional, the surface system associated with this patch.
    :type ssys: :py:class:`steps.API_2.model.SurfaceSystem` or :py:class:`str`

    Patches in a distributed tetrahedral mesh can be declared in the following ways::

        with mesh:
            comp1 = Patch.Create(inner)
            comp2 = Patch.Create(inner, ssys=ssys)
            comp3 = Patch.Create(inner, outer, ssys)
            comp4 = Patch.Create(inner, physicalTag=10)
            comp5 = Patch.Create(triList, inner, ssys=ssys)

    :param inner: Inner compartment, see below
    :type inner: :py:class:`steps.API_2.geom.Compartment`
    :param outer: Optional, Outer compartment, see below
    :type outer: :py:class:`steps.API_2.geom.Compartment`
    :param ssys: Optional, the surface system associated with this patch.
    :type ssys: :py:class:`steps.API_2.model.SurfaceSystem` or :py:class:`str`
    :param physicalTag: Optional, the index of the physical tag used to build the patch
        (see :py:class:`DistMesh`).
    :type physicalTag: int
    :param triLst: Optional, list of triangles associated with the patch. If an empy list is given
        the patch will be created from triangles that are tagged with the patch name.
    :type triLst: :py:class:`TriList` or any argument that can be used to build one

    **Relationship between Compartments and Patches**

    It is necessary to explain the inner/outer relationship between compartments
    and patches. When a patch object is created. it is necessary to arbitrarily label the
    compartment(s) "inner" and "outer" (if a patch is connected to only one compartment
    then the compartment must be labelled "inner" by convention). This is necessary
    in order to fully describe the surface reaction rules. Accordingly, compartments
    also store a list of connections, "inner" patches and "outer" patches. So if a
    patch1 is created with comp1 as it's "inner" compartment, comp1 knows patch1 as
    an "outer" patch. The labelling is purely defined when creating the Patch objects,
    bearing in mind the stoichiometry defined in the surface reaction objects. This may
    seem a little confusing at first, but will become clearer when experience is gained
    with these objects.

    For patches in tetrahedral meshes, the area is the total area of the encapsulated triangles.

    Some of the methods documented below are only available for tetrahedral patches, they are
    grouped accordingly.
    """

    _locStr = 'Patch'

    def __new__(cls, *args, **kwargs):
        # Patch is a facade to its subclasses, the user always uses the Patch class
        # but the appropriate subclass is used for object instantiation.

        if cls is not Patch:
            return super(Patch, cls).__new__(cls)

        geom = nutils.UsableObject._getUsedObjectOfClass(Geometry)
        if geom is None:
            raise Exception(f'Patches need to be declared in a geometry.')

        if geom.__class__ == Geometry:
            return _WmPatch.__new__(_WmPatch, *args, **kwargs)
        elif geom.__class__ == TetMesh:
            return _TetPatch.__new__(_TetPatch, *args, **kwargs)
        elif geom.__class__ == DistMesh:
            return _DistPatch.__new__(_DistPatch, *args, **kwargs)
        else:
            raise NotImplementedError(f'Patches cannot be created in {type(geom)}')

    @classmethod
    def _getDisplayName(cls):
        return Patch.__name__

    @classmethod
    def _FromStepsObject(cls, obj, geom):
        if stepslib._STEPS_USE_DIST_MESH and isinstance(obj, stepslib._py_DistPatch):
            return _DistPatch._FromStepsObject(obj, geom)
        elif isinstance(obj, stepslib._py_TmPatch):
            return _TetPatch._FromStepsObject(obj, geom)
        elif isinstance(obj, stepslib._py_Patch):
            return _WmPatch._FromStepsObject(obj, geom)
        else:
            raise NotImplementedError(f'Cannot create API2 version of {type(obj)}')

    def __init__(self, inner, outer=None, ssys=None, _createObj=True, **kwargs):
        super().__init__(**kwargs)

        if not isinstance(inner, Compartment):
            raise TypeError(f'Expected an inner Compartment, got {inner} instead.')
        if outer is not None and not isinstance(outer, Compartment):
            raise TypeError(f'Expected a Compartment, got {outer} instead.')

        self._innerComp = inner
        self._outerComp = outer

        icomp = self._innerComp.stepsComp if isinstance(self._innerComp, Compartment) else None
        ocomp = self._outerComp.stepsComp if isinstance(self._outerComp, Compartment) else None
        self.stepsPatch = self._createStepsObj(icomp, ocomp, self._geom) if _createObj else None

        if ssys is not None:
            self.addSystem(ssys)

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsPatch]

    def addSystem(self, ssys):
        """Add a surface system to the patch

        :param ssys: The surface system to be added, or its name.
        :type ssys: :py:class:`steps.API_2.model.SurfaceSystem` or `str`

        :returns: None
        """
        super().addSystem(ssys, nmodel.Location.SURF)
        if isinstance(ssys, str):
            self.stepsPatch.addSurfsys(ssys)
        elif ssys.__class__ is nmodel.SurfaceSystem:
            self.stepsPatch.addSurfsys(ssys.name)

    def _createStepsObj(self, icomp, ocomp, geom):
        raise NotImplementedError()

    def _canBeChild(self, c):
        """Return wether c can be a child of self."""
        if not super()._canBeChild(c):
            return False
        if ((isinstance(c, nmodel.Reaction) and not c._isSurfaceReac()) or
            (isinstance(c, nmodel.Diffusion) and not c._isSurfaceDiff()) or
            isinstance(c, (nmodel.VesicleBind, nmodel.VesicleUnbind))
        ):
            return False
        return True

    def _SetUpMdlDeps(self, mdl):
        """Set up the structures that will allow species or reaction access from locations."""
        super()._SetUpMdlDeps(mdl)
        for raft in mdl._getChildrenOfType(nmodel.Raft):
            self._addChildren(raft)
        for endoZone in self._getChildrenOfType(EndocyticZone):
            endoZone._SetUpMdlDeps(mdl)

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return Patch._locStr

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('m^2'))
    @nutils.IgnoresGhostElements
    def Area(self):
        """Surface area of the patch

        :type: Union[float, :py:class:`steps.API_2.utils.Parameter`], read-only for tetrahedral patches
        """
        return self.stepsPatch.getArea(**self._callKwargs)

    @property
    def innerComp(self):
        """Inner compartment

        :type: :py:class:`Compartment`
        """
        return self._innerComp

    @property
    def outerComp(self):
        """Outer compartment

        :type: :py:class:`Compartment`
        """
        return self._outerComp


@nutils.FreezeAfterInit
class _BasePatch(Patch):
    @classmethod
    def _FromStepsObject(cls, obj, geom):
        """Create the interface object from a STEPS object."""
        args = cls._ArgsFromStepsObject(obj, geom)
        patch = cls(*args, _createObj=False, name=obj.getID())
        patch.stepsPatch = obj

        # setup surface systems
        for ssysname in patch.stepsPatch.getSurfsys():
            patch.addSystem(ssysname)

        return patch

    @classmethod
    def _ArgsFromStepsObject(cls, obj, geom):
        icomp = getattr(geom, obj.icomp.getID())
        ocomp = getattr(geom, obj.ocomp.getID()) if obj.ocomp is not None else None
        return (icomp, ocomp)


@nutils.FreezeAfterInit
class _WmPatch(_BasePatch):
    def __init__(self, inner, outer=None, ssys=None, area=None, **kwargs):
        super().__init__(inner, outer, ssys=ssys, **kwargs)
        if area is not None:
            self.Area = area

    def _createStepsObj(self, icomp, ocomp, geom):
        return stepslib._py_Patch(self.name, geom.stepsGeom, icomp, ocomp)

    @Patch.Area.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('m^2'))
    def Area(self, v):
        self.stepsPatch.setArea(v)


@nutils.FreezeAfterInit
class _BaseTetPatch(_BasePatch):
    _FACADE_TITLE_STR = 'Only available for patches defined in tetrahedral meshes'

    def __init__(self, triLst, *args, **kwargs):
        self._triLst = triLst
        super().__init__(*args, **kwargs)

    @property
    @nutils.IgnoresGhostElements
    def bbox(self):
        """The bounding box of the patch

        :type: :py:class:`steps.API_2.geom.BoundingBox`, read-only
        """
        return BoundingBox(
            Point(*self.stepsPatch.getBoundMin(**self._callKwargs)),
            Point(*self.stepsPatch.getBoundMax(**self._callKwargs))
        )

    @property
    def tris(self):
        """All triangles in the patch

        :type: :py:class:`TriList`, read-only
        """
        raise NotImplementedError()

    @property
    def edges(self):
        """Perimeter of the triangles in the patch

        I.e. all bars that are not shared by two triangles in the patch.

        :type: :py:class:`BarList`, read-only
        """
        return self._triLst.edges


@nutils.FreezeAfterInit
class _TetPatch(_BaseTetPatch):
    _FACADE_TITLE_STR = None

    def __init__(self, triLst, inner, outer=None, ssys=None, **kwargs):
        super().__init__(triLst, inner, outer, ssys=ssys, **kwargs)

    def _createStepsObj(self, icomp, ocomp, geom):
        if icomp is None and ocomp is None:
            raise Exception(
                f'Cannot declare a Patch in a tetrahedral mesh with no neighboring compartment.'
            )
        self._triLst = TriList._toRefList(self._triLst, geom)
        return stepslib._py_TmPatch(self.name, geom.stepsMesh, self._triLst.indices, icomp, ocomp)

    @classmethod
    def _FromStepsObject(cls, obj, geom):
        """Create the interface object from a STEPS object."""
        patch = super()._FromStepsObject(obj, geom)

        with patch:
            for subObj in obj.getAllEndocyticZones():
                EndocyticZone._FromStepsObject(subObj, geom)

        return patch

    @classmethod
    def _ArgsFromStepsObject(cls, obj, geom):
        return (TriList(obj.getAllTriIndices(), geom), ) + _BasePatch._ArgsFromStepsObject(obj, geom)

    @_BaseTetPatch.tris.getter
    def tris(self):
        return self._triLst


@nutils.FreezeAfterInit
class _DistPatch(_BaseTetPatch):
    _FACADE_TITLE_STR = None

    def __init__(self, *args, physicalTag=None, **kwargs):
        self._physicalTag = physicalTag

        # Allow any ordering of inner, outer, ssys and triLst
        inner, outer, ssys, triLst, args, kwargs = nutils.extractArgs(args, kwargs, [
            ('inner', Compartment, None),
            ('outer', Compartment, None),
            ('ssys', (nmodel.SurfaceSystem, str), None),
            ('triLst', (TriList, list, tuple), None),
        ])

        self._triLst = triLst
        super().__init__(triLst, inner, outer, ssys=ssys, **kwargs)

    def _createStepsObj(self, icomp, ocomp, geom):
        if self._triLst is None:
            patch = stepslib._py_DistPatch(self.name, geom.stepsMesh, None, icomp, ocomp, self._physicalTag)
        else:
            triLst = TriList._toRefList(self._triLst, geom)
            if triLst.isLocal() and triLst._owned != False:
                raise Exception(
                    f'The patch is being created with a local triangle list but without '
                    f'including non-owned elements. It should instead be created inside a with '
                    f'mesh.asLocal(owned=False) block.'
                )
            patch = stepslib._py_DistPatch(
                self.name, geom.stepsMesh, triLst.indices, icomp, ocomp, self._physicalTag, local=triLst.isLocal()
            )
        del self._triLst
        return patch

    @_BaseTetPatch.tris.getter
    def tris(self):
        return TriList(
            self.stepsPatch.getAllTriIndices(**self._geom._callKwargs),
            mesh=self._geom, _immutable=True, **self._geom._lstArgs
        )

    @_BaseTetPatch.edges.getter
    def edges(self):
        raise NotImplementedError()


@nutils.FreezeAfterInit
class EndocyticZone(nutils.UsingObjects(Patch)):
    """A set of triangles that model a zone in which endocytosis reactions can happen

    Endocytosis reactions declared in surface systems do not happen at any point of the corresponding
    patches, they only happen in endocytic zones. By default all endocytosis events are active in this
    zone, but they can be deactivated during simulation.

    :param tris: The list of triangle
    :type tris: :py:class:`TriList`

    Endocytic zones are declared by using a patch as a context manager:

        with patch1:

            zone1 = EndocyticZone.Create(triLst)

    Note that the triangles in `triLst` all need to be part of `patch1`.
    """
    def __init__(self, tris, _createObj=True, *args, **kwargs):
        super().__init__(*args, **kwargs)
        (patch,) = self._getUsedObjects()

        if not isinstance(patch, _TetPatch):
            raise TypeError('Endocytic zones can only be declared in patches from tetrahedral meshes.')

        if not isinstance(tris, TriList):
            mesh = patch._getParentOfType(Geometry)
            if not isinstance(mesh, TetMesh):
                raise TypeError(f'{patch} was not declared in a tetrahedral mesh.')
            tris = TriList(tris, mesh=mesh)

        self._patch = patch
        self._tris = tris

        if _createObj:
            self.stepsEndoZone = stepslib._py_EndocyticZone(
                self.name, self._patch._getStepsObjects()[0], self._tris.indices
            )

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsEndoZone]

    @classmethod
    def _FromStepsObject(cls, obj, geom):
        """Create the interface object from a STEPS object."""
        endoZone = cls(obj.getAllTriIndices(), _createObj=False, name=obj.getID())
        endoZone.stepsEndoZone = obj

        return endoZone

    def _SetUpMdlDeps(self, mdl):
        """Set up the structures that will allow species or reaction access from locations."""
        for endo in self._patch._getChildrenOfType(nmodel.Endocytosis):
            self._addChildren(endo)

    @property
    def tris(self):
        """All triangles in the endocytic zone

        :type: :py:class:`TriList`, read-only
        """
        return self._tris


@nutils.FreezeAfterInit
class Membrane(_PhysicalLocation, nutils.Facade):
    """A set of patches on which membrane potential will be computed

    This class provides annotation for a group of triangles that comprise
    a surface describing a membrane in a Tetmesh. This may be the same
    description of one or several TmPatches in order for voltage-dependent
    transitions, currents and so on to be inserted in the membrane.
    A Membrane object must be available if membrane potential calculation is to be
    performed.

    By default STEPS only adds the inner compartments of the patches to the conduction volume,
    if other compartments should also be added to the conduction volume, they should be supplied
    through the supplementary_comps keyword parameter. 

    :param patches: List of patches whose triangles will be added to the membrane
        (only one patch for distributed meshes)
    :type patches: :py:class:`list`

    The following parameters are only available for non-distributed meshes:

    :param verify: Perform some checks on suitability of membrane (see below)
    :type verify: :py:class:`bool`
    :param opt_method: Optimization method (see below)
    :type opt_method: int
    :param search_percent: See below
    :type search_percent: float
    :param opt_file_name: Full path to a membrane optimization file
    :type opt_file_name: str
    :param supplementary_comps: Supplementary compartments to be considered part of the conduction volume
    :type supplementary_comps: List[:py:class:`Compartment`]

    Perform some checks on suitability of membrane if *verify* is True: these checks will print
    warnings if membrane forms an open surface or if any triangle is found to have
    more than 3 neighbours.
    Specify optimization method with *opt_method*:

    1. principal axis ordering (quick to set up but usually results in slower simulation than
       method 2).
    2. breadth first search (can be time-consuming to set up, but usually faster simulation.

    If breadth first search is chosen then argument *search_percent* can specify the number of
    starting points to search for the lowest bandwidth.
    If a filename (with full path) is given in optional argument *opt_file_name* the membrane
    optimization will be loaded from file, which was saved previously for this membrane with
    :py:func:`steps.API_2.sim.Simulation.saveMembOpt`
    """

    _locStr = 'Memb'

    def __new__(cls, *args, **kwargs):
        # Membrane is a facade to its subclasses, the user always uses the Membrane class
        # but the appropriate subclass is used for object instantiation.

        if cls is not Membrane:
            return super(Membrane, cls).__new__(cls)

        geom = nutils.UsableObject._getUsedObjectOfClass(Geometry)
        if geom is None:
            raise Exception(f'Membranes need to be declared in a geometry.')

        if geom.__class__ == TetMesh:
            return _TetMembrane.__new__(_TetMembrane, *args, **kwargs)
        elif geom.__class__ == DistMesh:
            return _DistMembrane.__new__(_DistMembrane, *args, **kwargs)
        else:
            raise NotImplementedError(f'Membranes cannot be created in {type(geom)}')

    @classmethod
    def _getDisplayName(cls):
        return Membrane.__name__

    @classmethod
    def _FromStepsObject(cls, obj, geom):
        raise NotImplementedError(f'Cannot create API2 version of {type(obj)}')

    def __init__(
        self, patches, *args, **kwargs
    ):
        super().__init__(*args, **kwargs)
        (geom,) = self._getUsedObjects()
        self.mesh = geom
        self.patches = patches
        self.stepsMemb = self._createStepsObj(geom)

    def addSystem(self, _):
        """
        Even though a Membrane is a PhysicalLocation, surface systems should be added to the
        patches, not through the membrane.
        :meta private:
        """
        raise NotImplementedError(f'Cannot add space systems to a Membrane.')

    def _SetUpMdlDeps(self, mdl):
        """Set up the structures that will allow species or reaction access from locations."""
        for p in self.patches:
            self.sysNames += p.sysNames
        self.sysNames = list(set(self.sysNames))

        super()._SetUpMdlDeps(mdl)

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsMemb]

    def _createStepsObj(self, geom):
        raise NotImplementedError()

    @property
    def open(self):
        """Whether the membrane is open (contains holes) or not

        :type: `bool`, read-only
        """
        return self.stepsMemb.open()

    @property
    def tris(self):
        """List of triangles associated with the membrane

        :type: :py:class:`steps.API_2.geom.TriList`, read-only
        """
        res = TriList([], mesh=self.mesh)
        for p in self.patches:
            res += p.tris
        return res

    @property
    def Area(self):
        """Surface area of the membrane

        :type: float, read-only
        """
        return sum(p.Area for p in self.patches)

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return Membrane._locStr


@nutils.FreezeAfterInit
class _TetMembrane(Membrane):

    def __init__(
        self, patches, verify=False, opt_method=1, search_percent=100.0, opt_file_name='',
        supplementary_comps=[], *args, **kwargs
    ):
        self._tmpParams = dict(
            verify=verify,
            opt_method=opt_method,
            search_percent=search_percent,
            opt_file_name=opt_file_name,
            supplementary_comps=[comp._getStepsObjects()[0] for comp in supplementary_comps],
        )
        super().__init__(patches, *args, **kwargs)

    def _createStepsObj(self, geom):
        (geomObj,) = geom._getStepsObjects()
        patches = [p._getStepsObjects()[0] for p in self.patches]
        return stepslib._py_Memb(
            self.name,
            geomObj,
            patches,
            **self._tmpParams,
        )


@nutils.FreezeAfterInit
class _DistMembrane(Membrane, nutils.ParameterizedObject):
    _FACADE_TITLE_STR = 'Only available for membranes defined in distributed tetrahedral meshes'

    def __init__(self, patches, capacitance=0, **kwargs):
        super().__init__(patches, **kwargs)
        if len(patches) != 1:
            raise ValueError(f'Membranes in distributed meshes can only have a single patch.')
        self.Capacitance = capacitance

    def _createStepsObj(self, geom):
        patches = [p._getStepsObjects()[0] for p in self.patches]
        return stepslib._py_DistMemb(self.name, geom.stepsMesh, patches)

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('F m^-2'))
    def Capacitance(self):
        """Capacitance of the compartment

        :type: Union[float, :py:class:`steps.API_2.utils.Parameter`]
        """
        return self.stepsMemb.getCapacitance()

    @Capacitance.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('F m^-2'))
    def Capacitance(self, v):
        self.stepsMemb.setCapacitance(v)


@nutils.FreezeAfterInit
class Point(numpy.ndarray):
    """Convenience class for representing 3D points

    This class inherits from :py:class:`numpy.ndarray` and can thus be used in the same way
    as a numpy array. The only difference is the possibility to access invidual coordinates
    through the *x*, *y*, and *z* properties.

    :param x: x coordinate
    :type x: float
    :param y: y coordinate
    :type y: float
    :param z: z coordinate
    :type z: float
    """

    def __new__(cls, x, y, z):
        arr = numpy.zeros(3)
        arr[0:3] = x, y, z
        return arr.view(cls)

    def __init__(self, *args, **kwargs):
        pass

    @property
    def x(self):
        """x coordinate

        :type: float, read-only
        """
        return self[0]

    @property
    def y(self):
        """y coordinate

        :type: float, read-only
        """
        return self[1]

    @property
    def z(self):
        """z coordinate

        :type: float, read-only
        """
        return self[2]

    def __hash__(self):
        return hash(tuple(self))


@nutils.FreezeAfterInit
class BoundingBox:
    """Convenience class for holding bounding box information"""

    def __init__(self, minP, maxP):
        self._min = minP
        self._max = maxP

    @property
    def min(self):
        """Minimum point of the bounding box

        :type: :py:class:`Point`, read-only
        """
        return self._min

    @property
    def max(self):
        """Maximum point of the bounding box

        :type: :py:class:`Point`, read-only
        """
        return self._max

    @property
    def center(self):
        """Center point of the bounding box

        :type: :py:class:`Point`, read-only
        """
        return (self._max - self._min) / 2 + self._min


class _BaseTetMesh(Geometry):
    """Base class for tetrahedral meshes"""

    def __init__(self, *args, _createObj=True, **kwargs):
        super().__init__(*args, _createObj=False, **kwargs)

        self.stepsMesh = self._createStepsObj() if _createObj else None

    def _createStepsObj(self):
        raise NotImplementedError()

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsMesh]

    @property
    def tets(self):
        """All tetrahedrons in the mesh

        :type: :py:class:`TetList`, read-only
        """
        return TetList(
            range(self.stepsMesh.countTets(**self._callKwargs)), self, _immutable=True, **self._lstArgs
        )

    @property
    def tris(self):
        """All triangles in the mesh

        :type: :py:class:`TriList`, read-only
        """
        return TriList(
            range(self.stepsMesh.countTris(**self._callKwargs)), self, _immutable=True, **self._lstArgs
        )

    @property
    def surface(self):
        """All triangles in the surface of the mesh

        :type: :py:class:`TriList`, read-only
        """
        return TriList(
            self.stepsMesh.getSurfTris(**self._callKwargs), self, _immutable=True, **self._lstArgs
        )

    @property
    def verts(self):
        """All vertices in the mesh

        :type: :py:class:`VertList`, read-only
        """
        return VertList(
            range(self.stepsMesh.countVertices(**self._callKwargs)), self, _immutable=True, **self._lstArgs
        )

    @property
    def bars(self):
        """All bars in the mesh

        :type: :py:class:`BarList`, read-only
        """
        return BarList(
            range(self.stepsMesh.countBars(**self._callKwargs)), self, _immutable=True, **self._lstArgs
        )

    @property
    @nutils.IgnoresGhostElements
    def bbox(self):
        """The bounding box of the mesh

        :type: :py:class:`steps.API_2.geom.BoundingBox`, read-only
        """
        return BoundingBox(
            Point(*self.stepsMesh.getBoundMin(**self._callKwargs)),
            Point(*self.stepsMesh.getBoundMax(**self._callKwargs))
        )

    @property
    @nutils.IgnoresGhostElements
    def Vol(self):
        """Total volume of the mesh

        :type: float, read-only
        """
        return self.stepsMesh.getMeshVolume(**self._callKwargs)

    def intersect(self, points, sampling=-1):
        """Computes the intersection of the current mesh and line segment(s), given their vertices

        :param points: A 2-D NumPy array (/memview) of 3D points where each element contains the 3 point
            coordinates
        :type points: numpy.ndarray
        :param sampling: if > 0, use montecarlo method with sampling points, otherwise use deterministic
            method.
        :type sampling: int

        :returns: A list of lists of tuples representing the intersected tetrahedrons, one element per line
            segment. Each tuple is made of 2 elements, a tetrahedron global identifier, and its respective
            intersection ratio.
        :rtype: List[List[Tuple[:py:class:`TetReference`, float]]]
        """
        res = self.stepsMesh.intersect(points, sampling)
        return [[(TetReference(idx, mesh=self, **self._lstArgs), rat) for idx, rat in seg] for seg in res]

    def intersectIndependentSegments(self, points, sampling=-1):
        """Similar to the intersect method but here we deal with independent segments, i.e.
        every two points we have a segment not related to previous or following ones.
        E.g. seg0 = (points[0], points[1]), seg1 = (points[2], points[3]), etc.

        :param points: A 2-D NumPy array (/memview) of 3D points where each element contains the 3 point
            coordinates
        :type points: numpy.ndarray
        :param sampling: if > 0, use montecarlo method with sampling points, otherwise use deterministic
            method.
        :type sampling: int

        :returns: A list of lists of tuples representing the intersected tetrahedrons, one element per line
            segment. Each tuple is made of 2 elements, a tetrahedron global identifier, and its respective
            intersection ratio.
        :rtype: List[List[Tuple[:py:class:`TetReference`, float]]]
        """
        res = self.stepsMesh.intersectIndependentSegments(points, sampling)
        return [[(TetReference(idx, mesh=self, **self._lstArgs), rat) for idx, rat in seg] for seg in res]


@nutils.FreezeAfterInit
class TetMesh(_BaseTetMesh):
    """Container class for static tetrahedral meshes

    This class stores the vertices points, 3D tetrahedral and 2D triangular elements that comprise
    the mesh. The indices of the elements will be stored as unsigned integers (a positive integer
    or zero) beginning at zero and incremented by 1. For example, if there are ntets number of
    tetrahedrons in the mesh, the indices of the tetrahedrons will be [0,1,2,..., (ntets-1)].

    A TetMesh object should not be created from scratch but should instead be created with one of
    the dedicated class methods:

    - :py:func:`TetMesh.Load`
    - :py:func:`TetMesh.LoadAbaqus`
    - :py:func:`TetMesh.LoadGmsh`
    - :py:func:`TetMesh.LoadVTK`
    - :py:func:`TetMesh.LoadTetGen`

    TetMesh objects can be used as a context manager in the same way as
    :py:class:`Geometry`::

        mesh = TetMesh.Load('/path/to/file')
        with mesh:
            comp1 = Comp.Create(vsys)

            ... # Declare other objects in mesh
    """

    def __init__(self, *args, _createObj=True, **kwargs):
        self._vertGroups = {}
        self._triGroups = {}
        self._tetGroups = {}
        self._vertProxy = None
        self._triProxy = None
        self._tetProxy = None
        if _createObj:
            raise Exception('Cannot create a bare TetMesh, use one of the class methods.')
        super().__init__(*args, _createObj=False, **kwargs)

    @classmethod
    def _FromStepsObject(cls, obj, comps=None, patches=None, name=None):
        """Create the interface object from a STEPS object."""
        mesh = cls(_createObj=False, name=name)
        mesh.stepsMesh = obj
        mesh.stepsGeom = obj
        comps = obj.getAllComps() if comps is None else comps
        patches = obj.getAllPatches() if patches is None else patches
        with mesh:
            for scomp in comps:
                Compartment._FromStepsObject(scomp, mesh)
            for spatch in patches:
                Patch._FromStepsObject(spatch, mesh)
            for ROIname in obj.getAllROINames():
                ROI._FromStepsObject(obj.getROI(ROIname), mesh, ROIname)
        return mesh

    @classmethod
    def FromData(cls, verts, tets, tris=[], name=None):
        """Create a mesh from geometrical data

        Construct a Tetmesh container with the folowing method: Supply a list of all
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

        :param verts: Vertex data
        :type verts: `List[int]`
        :param tets: Tetrahedron data
        :type tets: `List[int]`
        :param tris: Optional triangle data
        :type tris: `List[int]`
        :param name: Optional name for the mesh
        :type name: str

        :returns: The constructed TetMesh
        :rtype: :py:class:`TetMesh`
        """
        mesh = cls(_createObj=False, name=name)
        obj = stepslib._py_Tetmesh(verts, tets, tris)
        mesh.stepsMesh = obj
        mesh.stepsGeom = obj
        return mesh

    @classmethod
    def Load(cls, path, scale=1, strict=False, name=None):
        """Load a mesh in STEPS

        The mesh is loaded from an XML file that was previously generated by the
        :py:func:`TetMesh.Save` method.

        :param path: The root of the path where the file(s) are stored. e.g. with 'meshes/spine1'
            this function will look for the file /meshes/spine1.xml
        :type path: str
        :param scale: Optionally rescale the mesh on loading by given factor
        :type scale: float
        :param strict: Apply strict(-er) checking to the input XML
        :type strict: bool
        :param name: Optional name for the mesh
        :type name: str

        :returns: The loaded TetMesh
        :rtype: :py:class:`TetMesh`
        """
        smesh, scomps, spatches = smeshio.loadMesh(path, scale, strict)
        return TetMesh._FromStepsObject(smesh, comps=scomps, patches=spatches, name=name)

    def _loadElementProxys(self, ndprx, tetprx, triprx):
        """Load the blocks and groups from ElementProxy objects."""
        self._vertProxy, self._triProxy, self._tetProxy = ndprx, triprx, tetprx
        tmp = []
        for prx, lstCls in zip([ndprx, triprx, tetprx], [VertList, TriList, TetList]):
            # Merge blocks and groups
            dct = {n: lstCls(range(b[0], b[1] + 1), mesh=self) for n, b in prx.getBlocks().items()}
            for n, l in prx.getGroups().items():
                if n in dct:
                    warnings.warn(
                        f'Key {n} is used by a {lstCls._refCls._locStr} block and a '
                        f'group, taking the group.'
                    )
                dct[n] = lstCls(l, mesh=self)
            tmp.append(dct)
        self._vertGroups, self._triGroups, self._tetGroups = tmp

    # The shadow_mesh keyword parameter is left here for compatibility but is not documented since the
    # STEPS-CUBIT toolkit is deprecated.
    @classmethod
    def LoadAbaqus(cls, filename, scale=1, ebs=None, shadow_mesh=None, name=None):
        """Load a mesh from an ABAQUS-formated mesh file

        If blocks or groups of elements are present in the file, they will also be loaded.

        :param filename: The Abaqus filename (or path) including any suffix. Can also take a
            2-tuple containing separate paths for tetrahedron and triangle data (in this order).
        :type filename: str or `Tuple[str, str]`
        :param scale: Optionally rescale the mesh on loading by given factor
        :type scale: float
        :param ebs: Names of selected element blocks which are included in the mesh (does not apply
            if a 2-tuple is given as filename)
        :type ebs: `List[str]`
        :param name: Optional name for the mesh
        :type name: str

        :returns: The loaded TetMesh
        :rtype: :py:class:`TetMesh`
        """
        if isinstance(filename, tuple) and len(filename) == 2:
            tetfile, trifile = filename
            stepsMesh, ndprx, tetprx, triprx = smeshio.importAbaqus2(
                tetfile, trifile, scale, shadow_mesh=shadow_mesh
            )
        elif isinstance(filename, str):
            stepsMesh, ndprx, tetprx, triprx = smeshio.importAbaqus(
                filename, scale, ebs=ebs, shadow_mesh=shadow_mesh
            )
        else:
            raise TypeError(f'Expected a string or a 2-tuple of strings, got {filename} instead.')

        mesh = TetMesh._FromStepsObject(stepsMesh, name=name)
        mesh._loadElementProxys(ndprx, tetprx, triprx)
        return mesh

    @classmethod
    def LoadGmsh(cls, filename, scale=1, name=None):
        """Load a mesh from a Gmsh (2.2 ASCII)-formated mesh file

        If blocks or groups of elements are present in the file, they will also be loaded.

        :param filename: The Gmsh filename (or path) including any suffix.
        :type filename: str
        :param scale: Optionally rescale the mesh on loading by given factor
        :type scale: float
        :param name: Optional name for the mesh
        :type name: str

        :returns: The loaded TetMesh
        :rtype: :py:class:`TetMesh`
        """
        stepsMesh, ndprx, tetprx, triprx = smeshio.importGmsh(filename, scale)
        mesh = TetMesh._FromStepsObject(stepsMesh, name=name)
        mesh._loadElementProxys(ndprx, tetprx, triprx)
        return mesh

    @classmethod
    def LoadVTK(cls, filename, scale=1, name=None):
        """Load a mesh from a VTK (legacy ASCII)-formated mesh file

        If blocks or groups of elements are present in the file, they will also be loaded.

        :param filename: The VTK filename (or path) including any suffix.
        :type filename: str
        :param scale: Optionally rescale the mesh on loading by given factor
        :type scale: float
        :param name: Optional name for the mesh
        :type name: str

        :returns: The loaded TetMesh
        :rtype: :py:class:`TetMesh`
        """
        stepsMesh, ndprx, tetprx, triprx = smeshio.importVTK(filename, scale)
        mesh = TetMesh._FromStepsObject(stepsMesh, name=name)
        mesh._loadElementProxys(ndprx, tetprx, triprx)
        return mesh

    @classmethod
    def LoadTetGen(cls, pathroot, scale, name=None):
        """Load a mesh from a TetGen-formated set of files

        If blocks or groups of elements are present, they will also be loaded.

        :param pathroot: The root of the path to the mesh files. E.g. mesh/torus should be given
            to read files mesh/torus.node, mesh/torus.ele, and optionally mesh/torus.face
        :type pathroot: str
        :param scale: Optionally rescale the mesh on loading by given factor
        :type scale: float
        :param name: Optional name for the mesh
        :type name: str

        :returns: The loaded TetMesh
        :rtype: :py:class:`TetMesh`

        The TetGen file format for describing a mesh actually consists of
        multiple files with varying suffixes. Currently, this class method only
        reads meshes consisting of three files:

        * <input>.node: describing the tetrahedral mesh node points.
        * <input>.ele: describing tetrahedral elements, each of which
            consists of 4 pointers into the <input>.node file. (TetGen
            also supports 10-node elements; these 6 extra nodes are obviously
            not read by STEPS.)
        * <input>.face: describing triangular faces, each of which
            consists of 3 pointers into the <input>.node file. This file is optional.

        Other files are .vol (list of maximum volumes), .var (variant constraints)
        .neigh (list of neighbours), .smesh (simple PLC descriptions) and .edge
        (list of boundary edges) files. None of these seem relevant for our
        use cases, so we don't load them even when they are there. In particular,
        the .neigh file is computed by STEPS itself.

        Please refer to the TetGen manual (pages 31-40 in the last edition)
        for more information on these file formats

        tetgen.berlios.de/files/tetgen-manual.pdf
        """
        stepsMesh, ndprx, tetprx, triprx = smeshio.importTetGen(pathroot, scale)
        mesh = TetMesh._FromStepsObject(stepsMesh, name=name)
        mesh._loadElementProxys(ndprx, tetprx, triprx)
        return mesh

    def Save(self, path):
        """Save the mesh to an XML file

        :param path: The root of the path to store the files. e.g. 'meshes/spine1' will save
            data in /meshes/spine1.xml
        :type path: str

        This file stores the basic information about the mesh which tends to be
        common information for any software that supports tetrahedral meshes.

        * NODES are stored by cartesian coordinates.
        * TRIANGLES are stored by the indices of their 3 nodes.
        * TETRAHEDRONS are stored by the indices of their 4 nodes.

        The XML file also stores information about any compartments or
        patches created in STEPS (class :py:class:`steps.API_2.geom.Comp`
        :py:class:`steps.API_2.geom.Patch` respectively).

        * COMPARTMENT(S) are stored by:

          * their string identification.
          * a list of any volume systems added to the compartment at time of saving.
          * a list of tetrahedrons belonging to the compartment

        * PATCH(ES) are stored by:

          * their string identification.
          * a list of any surface systems added to the patch at time of saving.
          * the inner compartment id.
          * the outer compartment id (if it exists).
          * a list of triangles belonging to this patch.
        """
        smeshio.saveMesh(path, self.stepsMesh)

    def ConvertToMetis(self, path):
        """Convert the mesh data to metis connectivity data

        :param path: The path to the file in which to save the Metis data.
        :type path: str
        """
        smetis.tetmesh2metis(self.stepsMesh, path)

    @property
    def vertGroups(self):
        """Groups of vertices defined in the loaded mesh

        :type: Dict[str, :py:class:`VertList`], read-only

        .. note::
            Both blocks and groups are accessed through the same property.
        """
        return self._vertGroups

    @property
    def triGroups(self):
        """Groups of triangles defined in the loaded mesh

        :type: Dict[str, :py:class:`TriList`], read-only

        .. note::
            Both blocks and groups are accessed through the same property.
        """
        return self._triGroups

    @property
    def tetGroups(self):
        """Groups of tetrahedrons defined in the loaded mesh

        :type: Dict[str, :py:class:`TetList`], read-only

        .. note::
            Both blocks and groups are accessed through the same property.
        """
        return self._tetGroups

    @property
    def vertProxy(self):
        """Element proxy object for vertices (see :py:class:`steps.API_1.utilities.meshio.ElementProxy`)

        :type: :py:class:`steps.API_1.utilities.meshio.ElementProxy`

        .. note::
            Only available with some import methods.
        """
        return self._vertProxy

    @property
    def triProxy(self):
        """Element proxy object for triangles (see :py:class:`steps.API_1.utilities.meshio.ElementProxy`)

        :type: :py:class:`steps.API_1.utilities.meshio.ElementProxy`

        .. note::
            Only available with some import methods.
        """
        return self._triProxy

    @property
    def tetProxy(self):
        """Element proxy object for tetrahedrons (see :py:class:`steps.API_1.utilities.meshio.ElementProxy`)

        :type: :py:class:`steps.API_1.utilities.meshio.ElementProxy`

        .. note::
            Only available with some import methods.
        """
        return self._tetProxy


class _DistElemGroupsProxy:
    """Utility class for calling methods like DistMesh::getTaggedTetrahedrons"""
    def __init__(self, mesh, lstCls, method):
        self._mesh = mesh
        self._lstCls = lstCls
        self._method = method

    def keys(self):
        return self._mesh.stepsMesh.getTags(self._lstCls._dim)

    def __getitem__(self, tag):
        return self._lstCls(
            self._method(tag, **self._mesh._callKwargs),
            mesh=self._mesh, _immutable=True, **self._mesh._lstArgs
        )

    def __setitem__(self, n, v):
        raise NotImplementedError(f'Cannot add group of elements to distributed meshes.')

    def __delitem__(self, n):
        raise NotImplementedError(f'Cannot remove group of elements from distributed meshes.')


@nutils.FreezeAfterInit
class DistMesh(_BaseTetMesh):
    """Container class for distributed tetrahedral meshes

    In contrast to :py:class:`TetMesh`, only gmsh mesh files are accepted.

    :param filename: The Gmsh filename (or path) including any suffix.
    :type filename: str
    :param scale: Optionally rescale the mesh on loading by given factor
    :type scale: float
    :param mpiComm: mpi4py MPI communicator, defaults to COMM_WORLD
    :type mpiComm: :py:class:`mpi4py.MPI.Intracomm`

    ``filename`` should either be a full path to a single gmsh mesh or a path prefix for
    pre-partitioned meshes. If a single gmsh file is provided, it will be partitioned automatically.
    Pre-partitioned meshes can be loaded by providing a path prefix e.g. ``'path/to/meshName'`` and
    individual partition files should be named ``'path/to/meshName_1.msh'``, ``'path/to/meshName_2.msh'``,
    etc. The numbering of files starts at 1 and goes up to :py:attr:`steps.API_2.sim.MPI.nhosts`.

    Note that, for pre-partitioned meshes, elements are not necessarily numbered between ``0`` and ``n``
    (with ``n`` the number of elements).
    """
    def __init__(self, filename, scale=1, mpiComm=None, _createObj=True, **kwargs):
        self._filename = filename
        if _createObj:
            if mpiComm is None:
                # Only import mpi4py when necessary
                import mpi4py.MPI
                mpiComm = mpi4py.MPI.COMM_WORLD
            self._comm = mpiComm
        # In dist steps, scale == 0 means no scaling 
        self._scale = 0 if scale == 1 else scale
        
        super().__init__(_createObj=_createObj, **kwargs)

        self._callKwargs['local'] = self._local
        self._lstArgs['local'] = self._local

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsMesh]

    def _createStepsObj(self):
        library = stepslib._py_Library(self._comm)
        return stepslib._py_DistMesh(library, self._filename, self._scale)

    @contextlib.contextmanager
    def asGlobal(self):
        """Context manager method to use the mesh globally

        This should be used with the ``with`` statement in the following way::

            with mesh.asGlobal():
                lst1 = mesh.tets    # Returns a global list of tetrahedron indices
                lst2 = mesh.surface # Returns a global list of boundary triangles

        It is equivalent to using only ``with mesh:`` since this defaults to using global indices.
        ``with mesh.asGlobal():`` is however preferable since it is more explicit.

        :rtype: :py:class:`DistMesh`
        """
        try:
            self.__enter__()
            yield self
        finally:
            self.__exit__(None, None, None)

    @contextlib.contextmanager
    def asLocal(self, owned=True):
        """Context manager method to use the mesh locally

        This should be used with the ``with`` statement in the following way::

            with mesh.asLocal():
                lst1 = mesh.tets    # Returns a local list of tetrahedron indices
                lst2 = mesh.surface # Returns a local list of boundary triangles

        All properties and methods that would usually return lists of global indices will return local
        lists instead when wrapped with this context manager.

        In the above example, ``lst1`` would end up having the same value as ``mesh.tets.toLocal()``
        but the latter is less efficient since it first queries the full list of all tetrahedrons
        and then restricts it to only the local ones.

        .. note::
            This does not affect existing lists or newly created lists. For example, a tetrahedron list
            created in a ``with mesh.asLocal()`` block will still default to containing global elements
            (unless the ``local=True`` keyword argument is given at creation).
            This context manager only affects the behavior of the mesh object, and its subparts like
            compartments and patches. For example ``comp.tets`` will return a list of local indices
            in this context manager.

        :param owned: Whether only owned elements should be considered (defaults to True)
        :type owned: bool

        :rtype: :py:class:`DistMesh`
        """
        self._local = True
        self._owned = owned
        self._callKwargs['local'] = self._local
        self._lstArgs['local'] = self._local
        if not owned:
            self._callKwargs['owned'] = self._owned
        self._lstArgs['_owned'] = self._owned
        try:
            self.__enter__()
            yield self
        finally:
            self.__exit__(None, None, None)
            self._local = False
            self._owned = None
            self._callKwargs['local'] = self._local
            self._lstArgs['local'] = self._local
            if not owned:
                del self._callKwargs['owned']
            del self._lstArgs['_owned']

    def intersect(self, points, sampling=-1, raw=False, local=True):
        """Computes the intersection of the current mesh and line segment(s), given their vertices

        :param points: A 2-D NumPy array (/memview) of 3D points where each element contains the 3 point
            coordinates
        :type points: numpy.ndarray
        :param sampling: if > 0, use montecarlo method with sampling points, otherwise use deterministic
            method.
        :type sampling: int
        :param raw: If True, return raw integer tetrahedron indices.
        :type raw: bool
        :param local: if False, return global tetrahedron identifiers instead of local ones.
        :type local: bool

        This method can only be called when `mesh.asLocal()` is used, and it will return local tetrahedron
        references by default.

        :returns: A list of lists of tuples representing the intersected tetrahedrons, one element per line
            segment. Each tuple is made of 2 elements, a tetrahedron identifier (local or global depending
            on `local`), and its respective intersection ratio.
        :rtype: List[List[Tuple[Union[:py:class:`TetReference`, int], float]]]
        """
        if not self._local:
            raise Exception('Cannot use intersect method without using "with mesh.asLocal():".')
        res = self.stepsMesh.intersect(points, sampling)
        if not local:
            res = [[(self.stepsMesh.getTetGlobalIndex(idx), rat) for idx, rat in seg] for seg in res]
        if raw:
            return res
        else:
            return [[(TetReference(idx, mesh=self, local=local), rat) for idx, rat in seg] for seg in res]

    def intersectIndependentSegments(self, points, sampling=-1, raw=False, local=True):
        """Similar to the intersect method but here we deal with independent segments, i.e.
        every two points we have a segment not related to previous or following ones.
        E.g. seg0 = (points[0], points[1]), seg1 = (points[2], points[3]), etc.

        :param points: A 2-D NumPy array (/memview) of 3D points where each element contains the 3 point
            coordinates
        :type points: numpy.ndarray
        :param sampling: if > 0, use montecarlo method with sampling points, otherwise use deterministic
            method.
        :type sampling: int
        :param raw: If True, return raw integer tetrahedron indices.
        :type raw: bool
        :param local: if False, return global tetrahedron identifiers instead of local ones.
        :type local: bool

        This method can only be called when `mesh.asLocal()` is used, and it will return local tetrahedron
        references by default.

        :returns: A list of lists of tuples representing the intersected tetrahedrons, one element per line
            segment. Each tuple is made of 2 elements, a tetrahedron identifier (local or global depending
            on `local`), and its respective intersection ratio.
        :rtype: List[List[Tuple[Union[:py:class:`TetReference`, int], float]]]
        """
        if not self._local:
            raise Exception('Cannot use intersectIndependentSegments method without using "with mesh.asLocal():".')
        res = self.stepsMesh.intersectIndependentSegments(points, sampling)
        if not local:
            res = [[(self.stepsMesh.getTetGlobalIndex(idx), rat) for idx, rat in seg] for seg in res]
        if raw:
            return res
        else:
            return [[(TetReference(idx, mesh=self, local=local), rat) for idx, rat in seg] for seg in res]

    @_BaseTetMesh.tets.getter
    def tets(self):
        """All tetrahedrons in the mesh

        :type: :py:class:`TetList`, read-only
        """
        if self._local:
            return TetList(
                self.stepsMesh.getAllTetIndices(**self._callKwargs), self, _immutable=True, **self._lstArgs
            )
        else:
            return _BaseTetMesh.tets.fget(self)

    @_BaseTetMesh.tris.getter
    def tris(self):
        """All triangles in the mesh

        :type: :py:class:`TriList`, read-only
        """
        if self._local:
            return TriList(
                self.stepsMesh.getAllTriIndices(**self._callKwargs), self, _immutable=True, **self._lstArgs
            )
        else:
            return _BaseTetMesh.tris.fget(self)

    @_BaseTetMesh.bars.getter
    def bars(self):
        """All bars in the mesh

        :type: :py:class:`BarList`, read-only
        """
        if self._local:
            return BarList(
                self.stepsMesh.getAllBarIndices(**self._callKwargs), self, _immutable=True, **self._lstArgs
            )
        else:
            return _BaseTetMesh.bars.fget(self)

    @_BaseTetMesh.verts.getter
    def verts(self):
        """All vertices in the mesh

        :type: :py:class:`VertList`, read-only
        """
        if self._local:
            return VertList(
                self.stepsMesh.getAllVertIndices(**self._callKwargs), self, _immutable=True, **self._lstArgs
            )
        else:
            return _BaseTetMesh.verts.fget(self)

    @property
    def tetGroups(self):
        return _DistElemGroupsProxy(self, TetList, self.stepsMesh.getTaggedTetrahedrons)

    @property
    def triGroups(self):
        return _DistElemGroupsProxy(self, TriList, self.stepsMesh.getTaggedTriangles)

    @property
    def vertGroups(self):
        return _DistElemGroupsProxy(self, VertList, self.stepsMesh.getTaggedVertices)

    @staticmethod
    def _use_gmsh():
        return stepslib._py_DistMesh._use_gmsh()


class Reference(nutils.UsingObjects(nutils.Optional(_BaseTetMesh)), nutils.SolverPathObject, nutils.Facade):
    """Base class for all element references.

    :param idx: Index of the element
    :type idx: int
    :param mesh: Mesh object that contains the element
    :type mesh: Union[:py:class:`TetMesh`, :py:class:`DistMesh`, None]

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).

    The actual index of the element can be accessed with the idx attribute.
    """

    def __new__(cls, idx, mesh=None, _newCalled=False, **kwargs):
        # TetReference, TriReference, etc. are facades to their subclasses, the user always uses
        # the main class but the appropriate subclass is used for object instantiation.

        if cls not in Reference.__subclasses__() + _DistReference.__subclasses__() or _newCalled:
            return super(Reference, cls).__new__(cls)

        idx, mesh = cls._getIdxAndMeshFromParams(idx, mesh)

        if mesh.__class__ == DistMesh:
            return cls._distCls.__new__(cls._distCls, idx, mesh=mesh, _newCalled=True, **kwargs)
        else:
            return cls.__new__(cls, idx, mesh=mesh, _newCalled=True, **kwargs)

    def __init__(self, idx, mesh=None, anonymous=False, *args, **kwargs):
        if not anonymous:
            # Only call parent __init__ if we need the reference to be named
            super().__init__(*args, addAsElement=False, **kwargs)
        idx, mesh = self._getIdxAndMeshFromParams(idx, mesh)
        self.mesh = mesh
        self._idx = idx
        self._callKwargs = {}
        self._cloneArgs = dict(mesh=self.mesh)

    def toList(self):
        """Get a list holding only this element

        :returns: A list with only this element
        :rtype: :py:class:`RefList`
        """
        return self.__class__._lstCls([self], **self._cloneArgs)

    @property
    def idx(self):
        """Index of the element

        :Type: int, read-only
        """
        return self._idx

    @classmethod
    def _getIdxAndMeshFromParams(cls, idx, mesh):
        if mesh is None:
            # Try to infer the mesh from context-managers
            mesh = nutils.UsableObject._getUsedObjectOfClass(_BaseTetMesh)
            if mesh is None:
                # Try to infer the mesh from the idx argument
                if isinstance(idx, Reference):
                    mesh = idx.mesh
                else:
                    raise Exception('Cannot infer which mesh is associated with the reference.')
        if not isinstance(mesh, _BaseTetMesh):
            raise TypeError(f'Expected a TetMesh or DistMesh object, got {mesh} instead.')
        idx = cls._getActualIdx(idx, mesh)
        return idx, mesh

    @classmethod
    def _getActualIdx(cls, idx, mesh):
        """Return the integer idx."""
        if isinstance(idx, cls):
            if mesh != idx.mesh:
                raise Exception(
                    f'Cannot create a reference from a reference related to a different mesh.'
                )
            idx = idx._idx
        if not isinstance(idx, numbers.Integral):
            raise TypeError(f'Expected an integer index, got {idx} instead.')
        return idx

    def _getPhysicalLocation(self):
        """Return the location associated with the reference (comp for tet, patch for tri)."""
        pass

    def __hash__(self):
        return hash((id(self.mesh), self._idx))

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.mesh == other.mesh and self._idx == other._idx

    def __repr__(self):
        return f'{self.__class__._locStr}({self._idx})'

    def _solverId(self):
        """Return the id that is used to identify an object in the solver."""
        return (self._idx,)

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        mtdt = {'loc_type': self.__class__._locStr, 'loc_id': self._idx}
        # Add information about parent comp or patch
        parent = self._getPhysicalLocation()
        if parent is not None:
            for key, val in parent._simPathAutoMetaData().items():
                mtdt['parent_' + key] = val
        return mtdt


class _DistReference(Reference):
    _FACADE_TITLE_STR = 'Only available for references defined in distributed tetrahedral meshes'

    def __init__(self, idx, mesh=None, anonymous=False, *args, local=False, _owned=None, **kwargs):
        super().__init__(idx, mesh, anonymous, *args, **kwargs)

        if isinstance(idx, _DistReference):
            local = idx.isLocal()

        self._local = local
        self._callKwargs['local'] = self._local
        self._cloneArgs['local'] = self._local

        # The owned flag is just indicative of how the reference was created
        self._owned = _owned

    @classmethod
    def _getToLocalFunc(cls, mesh):
        """This is mainly used to avoid having to instantiate a VertReference when manipulating
        vert indices.
        """
        return getattr(mesh.stepsMesh, f'get{cls._locStr}LocalIndex')

    @classmethod
    def _getToGlobalFunc(cls, mesh):
        return getattr(mesh.stepsMesh, f'get{cls._locStr}GlobalIndex')

    def isLocal(self):
        """Return whether the reference uses a local index

        :rtype: bool
        """
        return self._local

    def toLocal(self, owned=True):
        """Get the local version of the reference

        If the reference is local, it returns itself. If the reference is global, and if the element
        exists in the current MPI process, a local reference to it is returned. Otherwise, 
        `None` is returned. The returned reference is thus different for each MPI process.

        :param owned: Whether only owned elements should be considered (defaults to True)
        :type owned: bool

        :rtype: :py:class:`Reference`
        """
        if self._local:
            if self._owned == owned:
                return self
            else:
                return self.toGlobal().toLocal(owned)
        method = self._getToLocalFunc(self.mesh)
        ind = method(self.idx, owned=owned)
        if ind is None:
            return None
        else:
            return self.__class__(ind, mesh=self.mesh, local=True, _owned=owned)

    def toGlobal(self, **_kwargs):
        """Get the global version of the reference

        If the reference is not local, it returns itself. Otherwise, the local index of the original
        reference is converted to global mesh index and the corresponding global reference is returned.

        :rtype: :py:class:`Reference`
        """
        if not self._local:
            return self
        method = self._getToGlobalFunc(self.mesh)
        return self.__class__(method(self.idx), mesh=self.mesh, local=False, **_kwargs)

    def _solverKeywordParams(self):
        """Return the additional keyword parameters that should be passed to the solver"""
        return dict(local=self._local)

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        mtdt = super()._simPathAutoMetaData()
        if self._local:
            # The location id in the metadata should always be global
            # Get global id without creating a new reference because it causes an issue
            # for VertReferences that querry their position upon creation
            method = getattr(self.mesh.stepsMesh, f'get{self._locStr}GlobalIndex')
            mtdt['loc_id'] = method(self.idx)
        return mtdt

    def __hash__(self):
        return hash((id(self.mesh), self._idx, self._local))

    def __eq__(self, other):
        return isinstance(other, self.__class__) and\
            (self.mesh, self._idx, self._local) == (other.mesh, other._idx, other._local)

    def __repr__(self):
        if self._local:
            return f'Local{self.__class__._locStr}({self._idx})'
        else:
            return super().__repr__()


@nutils.FreezeAfterInit
class TetReference(Reference):
    """Convenience class for accessing properties of a tetrahedron

    :param idx: Index of the tetrahedron
    :type idx: int
    :param mesh: Mesh object that contains the tetrahedron
    :type mesh: Union[:py:class:`TetMesh`, :py:class:`DistMesh`, None]

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).
    """

    _locStr = 'Tet'

    def __init__(self, idx, mesh=None, *args, **kwargs):
        super().__init__(idx, mesh=mesh, *args, **kwargs)

    @classmethod
    def _getActualIdx(cls, idx, mesh):
        """Return the integer idx."""
        try:
            idx = super()._getActualIdx(idx, mesh)
        except TypeError as ex:
            if hasattr(idx, '__iter__') and len(idx) == 3:
                idx = mesh.stepsMesh.findTetByPoint(list(idx))
                if idx == UNKNOWN_TET:
                    raise Exception(f'No tetrahedron at position {idx}.')
            else:
                raise ex
        return idx

    def __getattr__(self, name):
        """Handle the references to species, reactions, or diffusion rules."""
        if self.comp is not None:
            return getattr(self.comp, name)
        else:
            raise AttributeError

    @property
    def children(self):
        if self.comp is not None:
            return self.comp.children
        else:
            return {}

    def containsPoint(self, point):
        """Check whether a given 3D point is inside the tetrahedron

        :param point: The 3D point
        :type point: List[float]

        :rtype: bool
        """
        return self.mesh.stepsMesh.isPointInTet(point, self._idx, **self._callKwargs)

    @property
    def neighbs(self):
        """Neighboring tetrahedrons (between 0 and 4)

        :type: :py:class:`TetList`, read-only
        """
        inds = self.mesh.stepsMesh.getTetTetNeighb(self._idx, **self._callKwargs)
        return TetList([ind for ind in inds if ind != UNKNOWN_TET], **self._cloneArgs)

    @property
    def faces(self):
        """Faces of the tetrahedron

        :type: :py:class:`TriList`, read-only
        """
        return TriList(self.mesh.stepsMesh.getTetTriNeighb(self._idx, **self._callKwargs), **self._cloneArgs)

    @property
    def center(self):
        """Barycenter of the tetrahedron

        :type: :py:class:`Point`, read-only
        """
        return Point(*self.mesh.stepsMesh.getTetBarycenter(self._idx, **self._callKwargs))

    @property
    def verts(self):
        """Vertices of the tetrahedron

        :type: :py:class:`VertList`, read-only
        """
        return VertList(self.mesh.stepsMesh.getTet(self._idx, **self._callKwargs), **self._cloneArgs)

    @property
    def comp(self):
        """Compartment associated to the tetrahedron (if any)

        :type: :py:class:`Comp` or None, read-only
        """
        try:
            c = self.mesh.stepsMesh.getTetComp(self._idx, **self._callKwargs)
            if c is None:
                return None
            else:
                return getattr(self.mesh, c.getID())
        except AttributeError as ex:
            # Since comp is called in __getattr__, raising an AttributeError here can lead to
            # infinite recursion.
            raise Exception(ex)

    @property
    def Vol(self):
        """Volume of the tetrahedron

        :type: float, read-only
        """
        return self.mesh.stepsMesh.getTetVol(self._idx, **self._callKwargs)

    @property
    def qualityRER(self):
        """Radius-edge-ratio (a quality measurement) of the tetrahedron

        :type: float, read-only
        """
        return self.mesh.stepsMesh.getTetQualityRER(self._idx, **self._callKwargs)

    def _getPhysicalLocation(self):
        """Return the location associated with the reference (comp for tet, patch for tri)."""
        return self.comp

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return TetReference._locStr


class _DistTetReference(TetReference, _DistReference):
    @TetReference.neighbs.getter
    def neighbs(self):
        if self.mesh._local:
            inds = self.mesh.stepsMesh.getTetTetNeighb(self._idx, **self._callKwargs, owned=self.mesh._owned)
            return TetList([ind for ind in inds if ind != UNKNOWN_TET], **self._cloneArgs)
        else:
            return TetReference.neighbs.fget(self)


TetReference._distCls = _DistTetReference


@nutils.FreezeAfterInit
class TriReference(Reference):
    """Convenience class for accessing properties of a triangle

    :param idx: Index of the triangle
    :type idx: int
    :param mesh: Mesh object that contains the triangle
    :type mesh: Union[:py:class:`TetMesh`, :py:class:`DistMesh`, None]

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).
    """

    _locStr = 'Tri'

    def __init__(self, idx, mesh=None, *args, **kwargs):
        super().__init__(idx, mesh=mesh, *args, **kwargs)

    def __getattr__(self, name):
        """Handle the references to species or reactions."""
        if self.patch is not None:
            return getattr(self.patch, name)
        else:
            raise AttributeError

    @property
    def children(self):
        if self.patch is not None:
            return self.patch.children
        else:
            return {}

    @property
    def verts(self):
        """Vertices of the triangle

        :type: :py:class:`VertList`, read-only
        """
        return VertList(self.mesh.stepsMesh.getTri(self._idx, **self._callKwargs), **self._cloneArgs)

    @property
    def center(self):
        """Barycenter of the triangle

        :type: :py:class:`Point`, read-only
        """
        return Point(*self.mesh.stepsMesh.getTriBarycenter(self._idx, **self._callKwargs))

    @property
    def Area(self):
        """Area of the triangle

        :type: float, read-only
        """
        return self.mesh.stepsMesh.getTriArea(self._idx, **self._callKwargs)

    @property
    def tetNeighbs(self):
        """Neighboring tetrahedrons

        :type: :py:class:`TetList`, read-only
        """
        inds = self.mesh.stepsMesh.getTriTetNeighb(self._idx, **self._callKwargs)
        return TetList([ind for ind in inds if ind != UNKNOWN_TET], **self._cloneArgs)

    @property
    def triNeighbs(self):
        """Neighboring triangles

        :type: :py:class:`TriList`, read-only
        """
        inds = self.mesh.stepsMesh.getTriTriNeighbs(self._idx, **self._callKwargs)
        return TriList(inds, **self._cloneArgs)

    @property
    def patchTriNeighbs(self):
        """Neighboring triangles that are also part of the same patch
        If the triangle is not part of a :py:class:`Patch`, returns an empty :py:class:`TriList`

        :type: :py:class:`TriList`, read-only
        """
        patch = self.patch
        if patch is not None:
            inds = self.mesh.stepsMesh.getTriTriNeighb(self._idx, patch.stepsPatch, **self._callKwargs)
            inds = [i for i in inds if i != UNKNOWN_TRI]
        else:
            inds = []
        return TriList(inds, **self._cloneArgs)

    @property
    def patch(self):
        """Patch associated to the triangle (if any)

        :type: :py:class:`Patch` or None, read-only
        """
        p = self.mesh.stepsMesh.getTriPatch(self._idx, **self._callKwargs)
        if p is None:
            return None
        else:
            return getattr(self.mesh, p.getID())

    @property
    def bars(self):
        """bars of the triangle

        :type: :py:class:`VertList`, read-only
        """
        return BarList(self.mesh.stepsMesh.getTriBars(self._idx, **self._callKwargs), **self._cloneArgs)

    @property
    def norm(self):
        """Normal vector of the triangle

        :type: :py:class:`Point`, read-only
        """
        return Point(*self.mesh.stepsMesh.getTriNorm(self._idx, **self._callKwargs))

    def _getPhysicalLocation(self):
        """Return the location associated with the reference (comp for tet, patch for tri)."""
        return self.patch

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return TriReference._locStr


class _DistTriReference(TriReference, _DistReference):
    @TriReference.tetNeighbs.getter
    def tetNeighbs(self):
        if self.mesh._local:
            inds = self.mesh.stepsMesh.getTriTetNeighb(self._idx, **self._callKwargs, owned=self.mesh._owned)
            return TetList([ind for ind in inds if ind != UNKNOWN_TET], **self._cloneArgs)
        else:
            return TriReference.tetNeighbs.fget(self)


TriReference._distCls = _DistTriReference


@nutils.FreezeAfterInit
class BarReference(Reference):
    """Convenience class for accessing properties of a bar

    :param idx: Index of the bar
    :type idx: int
    :param mesh: Mesh object that contains the bar
    :type mesh: Union[:py:class:`TetMesh`, :py:class:`DistMesh`, None]

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).
    """

    _locStr = 'Bar'

    def __init__(self, idx, mesh=None, *args, **kwargs):
        super().__init__(idx, mesh=mesh, *args, **kwargs)

    @property
    def verts(self):
        """Vertices of the bar

        :type: :py:class:`VertList`, read-only
        """
        return VertList(self.mesh.stepsMesh.getBar(self._idx, **self._callKwargs), **self._cloneArgs)

    @property
    def center(self):
        v1, v2 = self.verts
        return Point(*((v1 + v2) / 2))


class _DistBarReference(BarReference, _DistReference):
    pass


BarReference._distCls = _DistBarReference


@nutils.FreezeAfterInit
class VertReference(Reference, Point):
    """Convenience class for accessing properties of a vertex

    Can be used in the same way as a :py:class:`Point`.

    :param idx: Index of the vertex
    :type idx: int
    :param mesh: Mesh object that contains the vertex
    :type mesh: Union[:py:class:`TetMesh`, :py:class:`DistMesh`, None]

    .. warning::
        When used with MPI, keep in mind that the creation of a :py:class:`VertReference`
        retrieves the coordinates of the vertex and thus potentially requires synchronization
        across MPI ranks. When using :py:class:`DistMesh`, there is no need for synchronization
        if the `local=True` keyword argument is used.

    .. note::
        In contrast to other references, VertReference can only be created by explicitely
        providing a mesh object, it cannot use the context manager syntax described in
        :py:class:`TetMesh`.
    """

    _locStr = 'Vert'

    def __new__(cls, idx, mesh, anonymous=False, *args, _newCalled=False, **kwargs):
        if not _newCalled:
            return Reference.__new__(cls, idx, mesh, anonymous=anonymous, *args, **kwargs)

        if not isinstance(mesh, _BaseTetMesh):
            raise TypeError(f'Expected a TetMesh or DistMesh object, got {mesh} instead.')
        idx = cls._getActualIdx(idx, mesh)
        return Point.__new__(cls, numpy.nan, numpy.nan, numpy.nan)

    def __init__(self, idx, mesh=None, *args, _getData=True, **kwargs):
        super().__init__(idx, mesh=mesh, *args, **kwargs)
        if _getData:
            self[0:3] = self.mesh.stepsMesh.getVertex(self._idx, **self._callKwargs)

    # Need to redefine __hash__ and __eq__ to be sure the Reference version of these methods
    # is called.
    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        return super().__eq__(other)

    def __ne__(self, other):
        return not super().__eq__(other)

    def __repr__(self):
        return super().__repr__()

    def __str__(self):
        return super().__repr__()

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return VertReference._locStr


class _DistVertReference(VertReference, _DistReference):
    def __init__(self, idx, mesh=None, *args, local=False, _fromLocal=False, **kwargs):
        if not local and _fromLocal:
            warnings.warn(
                f'Global vertex references created from local references or lists are loaded without '
                f'coordinates data (which would require global MPI synchronization).'
            )
        super().__init__(idx, mesh=mesh, *args, local=local, _getData=not _fromLocal, **kwargs)

    def toGlobal(self, **_kwargs):
        """Get the global version of the reference

        If the reference is not local, it returns itself.

        Otherwise, the local index of the original reference is converted to global mesh index and
        the corresponding global reference is returned. The returned reference will not contain coordinate
        data, which would require an MPI synchronization. In order to access coordinate data, the global
        index first needs to be synchronized across ranks.

        :rtype: :py:class:`Reference`
        """
        return super().toGlobal(_fromLocal=self._local)


VertReference._distCls = _DistVertReference


class RefList(nutils.UsingObjects(nutils.Optional(_BaseTetMesh)), nutils.SolverPathObject, nutils.Facade):
    """Base convenience class for geometrical element lists

    :param lst: The list of elements
    :type lst: Iterable[int] or range or Iterable[:py:class:`Reference`] or :py:class:`RefList`
    :param mesh: Mesh object that contains the elements
    :type mesh: Union[:py:class:`TetMesh`, :py:class:`DistMesh`, None]

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).

    Behaves like a python list but with additional functionalities.
    """

    class _OptimizationCM:
        """Context manager for handling the updating and invalidation of optimization data

        All code that modifies the list needs to be wrapped using ``with self._optimCM:``.
        Code that does not modify the list can directly use self._optimCM.
        Beware of not wrapping non-modifying code, it can lead to strange behavior.
        """

        # TODO Not urgent: Refactor to have a clear distinction (better granularity) between data 
        # setting linked to actual modification of the list and data setting for local caching
        def __init__(self, parent):
            self.parent = parent
            self._optimData = {}
            self._modifData = set()

        def __enter__(self):
            if self.parent._immutable:
                raise Exception(f'Cannot modify an immutable list.')
            self._modifData = set()
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            self._optimData = {key: val for key, val in self._optimData.items() if key in self._modifData}

        def __getitem__(self, key):
            if key in self._optimData:
                return self._optimData[key]
            return None

        def __setitem__(self, key, val):
            self._optimData[key] = val
            self._modifData.add(key)

        def __contains__(self, key):
            return key in self._optimData

    def __new__(cls, lst=None, mesh=None, _immutable=False, _newCalled=False, **kwargs):
        # TetList, TriList, etc. are facades to their subclasses, the user always uses the main class
        # but the appropriate subclass is used for object instantiation.

        if cls not in RefList.__subclasses__() or _newCalled:
            return super(RefList, cls).__new__(cls)

        if mesh is None:
            # Try to infer the mesh from context-managers
            mesh = nutils.UsableObject._getUsedObjectOfClass(_BaseTetMesh)
            # Try to infer the mesh from the lst argument
            if mesh is None:
                # If lst is a generator, do not consume it by trying to infer the mesh
                lst = None if isinstance(lst, types.GeneratorType) else lst
                lst, mesh = cls._getLstAndMeshFromParams(lst, mesh)

        if mesh.__class__ == DistMesh:
            return cls._distCls.__new__(
                cls._distCls, lst, mesh=mesh, _immutable=_immutable, _newCalled=True, **kwargs
            )
        else:
            return cls.__new__(cls, lst, mesh=mesh, _immutable=_immutable, _newCalled=True, **kwargs)

    def __init__(self, lst, mesh=None, _immutable=False, *args, **kwargs):
        super().__init__(*args, addAsElement=False, **kwargs)
        if mesh is None:
            (mesh,) = self._getUsedObjects(addAsElement=False)
        elif not isinstance(mesh, _BaseTetMesh):
            raise TypeError(f'Expected a TetMesh or DistMesh object, got a {type(mesh)} instead.')

        self.lst, self.mesh = self._getLstAndMeshFromParams(lst, mesh)

        self._immutable = _immutable
        # Optimization data wrapper, used as a context manager when modifying the list
        self._optimCM = RefList._OptimizationCM(self)
        # Keyword dict to initialize a list with the same parameters as this one
        self._cloneArgs = dict(mesh=self.mesh)

    @classmethod
    def _getIdxLst(cls, lst):
        lst = list(lst)
        if all(isinstance(e, numbers.Integral) for e in lst):
            return lst
        elif all(isinstance(e, cls._refCls) for e in lst):
            return [e._idx for e in lst]
        else:
            raise TypeError(f'Cannot use {lst} to initialize a list of {cls._refCls}.')

    @classmethod
    def _toRefList(cls, lst, mesh):
        if isinstance(lst, RefList):
            if not isinstance(lst, cls):
                raise TypeError(f'Expected a {cls}, got a {lst.__class__} instead.')
            return lst
        else:
            return cls(lst, mesh)

    @classmethod
    def _getLstAndMeshFromParams(cls, lst, mesh):
        actualLst = []
        actualMesh = mesh
        if cls in lst.__class__.__mro__:
            actualMesh = lst.mesh
            actualLst = copy.copy(lst.lst)
        elif isinstance(lst, range):
            actualLst = lst
        elif hasattr(lst, '__iter__'):
            lst = list(lst)
            actualLst = cls._getIdxLst(lst)
            if len(lst) > 0 and isinstance(lst[0], cls._refCls):
                actualMesh = lst[0].mesh
        elif lst is not None:
            raise TypeError(f'Cannot create an element list from a {type(lst)}')

        if actualMesh is None:
            raise Exception(f'Unable to identify which mesh should the list be bound to.')

        return actualLst, actualMesh

    def _splitByLocation(self):
        """Split the list in several lists that contain elements in the same comp / patch
        The order of elements is conserved. The name of the list is conserved but can be postfixed.
        """
        res = []
        loc = (None,)
        for elem in self:
            newLoc = elem._getPhysicalLocation()
            if newLoc != loc:
                loc = newLoc
                if not self._autoNamed:
                    name = f'{self.name}({loc})' if len(res) > 0 else self.name
                else:
                    name=None
                res.append(self.__class__([], **self._cloneArgs, name=name))
            res[-1].append(elem)
        return res

    def splitToROIs(self, ROIPrefix, method='grid', gridSize=None):
        """Split an element list into a square grid and define each square as a Region Of Interest

        :param ROIPrefix: A prefix for the names of the ROIs, it will be suffixed with 'x_y_z' with x, y, and
            z the integer grid coordinate of a given square ROI.
        :type ROIPrefix: str
        :param method: The method used to define the ROIs, only 'grid' is available for now.
        :type method: str
        :param gridSize: The side length of a single square from the grid, in meters.
        :type gridSize: float
        """
        rois = []
        if method == 'grid':
            gridCells = {}
            for elem in self:
                key = tuple(elem.center // gridSize)
                gridCells.setdefault(key, self.__class__([], **self._cloneArgs)).append(elem)
            for key, elems in gridCells.items():
                name = ROIPrefix + '_' + '_'.join(map(str, key))
                rois.append(ROI(elems, name=name))
        return rois


    def _modify(self):
        """Return a context manager and signal that the wrapped statements will modify the list"""
        # Simply return the _OptimizationCM for now
        return self._optimCM

    def __getitem__(self, key):
        """Access one or several element(s) from the list

        :param key: Index of the element or slice
        :type key: int or slice
        :returns: The element(s)
        :rtype: :py:class:`Reference` or :py:class:`RefList`

        Usage::

            >>> lst = RefList(range(10, 20), mesh=mesh)
            >>> lst[0] # Returns the 1st element from the list
            Tet(10)
            >>> lst[0].idx
            10
            >>> lst2 = lst[1:5] # Returns elements between the 2nd and the 5th (included).
            >>> lst2[0]
            Tet(11)

        .. note::
            The index of an element in the list is in general different from its index in the mesh.

        :meta public:
        """
        idxres = self.lst[key]
        if isinstance(idxres, (list, range)):
            return self.__class__(idxres, **self._cloneArgs)
        else:
            return self.__class__._refCls(idxres, **self._cloneArgs)

    def __iter__(self):
        """Get an iterator over the list elements

        Usage::

            >>> lst = RefList(range(3), mesh=mesh)
            >>> for ref in lst:
            ...     print(ref.idx)
            ...
            0
            1
            2
            >>> [ref.idx for ref in lst]
            [0, 1, 2]

        .. note:
            The references that it returns does not have a name, in contrast to a
            :py:class:`Reference` created by the user directly.

        :meta public:
        """
        for idx in self.lst:
            yield self.__class__._refCls(idx, **self._cloneArgs, anonymous=True)

    def append(self, e):
        """Append an element to the list

        :param e: The element
        :type e: :py:class:`Reference`

        .. note::
            The type of the element depends on the type of the list, one cannot add a
            :py:class:`TriReference` to a :py:class:`TetList`.
        """
        if not isinstance(e, self.__class__._refCls):
            raise TypeError(f'Cannot append {e}, expected a {self.__class__._refCls.__name__} object')
        self.lst = self._getLst()
        # Modifying the list
        with self._modify():
            self.lst.append(e._idx)

    def remove(self, e):
        """Remove an element from the list

        :param e: The element
        :type e: :py:class:`Reference`
        """
        if not isinstance(e, self.__class__._refCls):
            raise TypeError(f'Cannot remove {e}, expected a {self.__class__._refCls.__name__} object')
        self.lst = self._getLst()
        # Modifying the list
        with self._modify():
            self.lst.remove(e._idx)

    _ElemsAsSet_K = 'ElemsAsSet'

    def __contains__(self, elem):
        """Check whether the list contains an element

        :param elem: The element
        :type elem: :py:class:`Reference`
        :rtype: bool

        Usage::

            >>> ref1 = Reference(15, mesh=mesh)
            >>> ref2 = Reference(25, mesh=mesh)
            >>> lst = RefList(range(20), mesh=mesh)
            >>> ref1 in lst
            True
            >>> ref2 in lst
            False

        :meta public:
        """
        if isinstance(elem, self.__class__._refCls):
            # Optimization
            if RefList._ElemsAsSet_K not in self._optimCM:
                self._optimCM[RefList._ElemsAsSet_K] = set(self.lst)
            return elem._idx in self._optimCM[RefList._ElemsAsSet_K]
        else:
            return False

    def __len__(self):
        """Get the length of the list

        :rtype: int

        Usage::

            >>> lst = RefList(range(5), mesh=mesh)
            >>> len(lst)
            5

        :meta public:
        """
        return len(self.lst)

    def index(self, elem):
        """Return the index of the element elem in this list

        :param elem: The element
        :type elem: :py:class:`Reference`
        :returns: The index of the element in the class (not its index in the mesh)
        :rtype: int

        Behaves like :py:func:`list.index`.
        """
        if isinstance(elem, self.__class__._refCls):
            if self.mesh is not elem.mesh:
                raise Exception(
                    f'Cannot retrieve the index of an element that is associated to a different mesh.'
                )
            return self.lst.index(elem._idx)
        else:
            raise TypeError(f'Expected a {self.__class__._refCls}, got {elem} instead.')

    def _checkSameType(self, other):
        if other.__class__ != self.__class__:
            raise TypeError('Cannot combine lists of different types.')
        if self.mesh is not other.mesh:
            raise Exception('Cannot combine lists associated to different meshes.')

    def _getLst(self):
        return list(self.lst) if isinstance(self.lst, range) else self.lst

    def __and__(self, other):
        """Compute the intersection between two lists

        :returns: All the elements that are in both lists. Keeps the order of the left operand.
        :rtype: :py:class:`RefList`

        Usage::

            >>> lst1 = RefList(range(0, 10), mesh=mesh)
            >>> lst2 = RefList(range(5, 15), mesh=mesh)
            >>> lst3 = lst1 & lst2
            >>> [ref.idx for ref in lst3]
            [5, 6, 7, 8, 9]

        :meta public:
        """
        self._checkSameType(other)
        s2 = set(other.lst)
        res = [i for i in self.lst if i in s2]
        return self.__class__(res, **self._cloneArgs)

    def __or__(self, other):
        """Compute the union of two lists

        :returns: All the elements of the left operand followed by all the elements of the right
            operand that are not in the left one.
        :rtype: :py:class:`RefList`

        Usage::

            >>> lst1 = RefList(range(0, 10), mesh=mesh)
            >>> lst2 = RefList(range(5, 15), mesh=mesh)
            >>> lst3 = lst1 | lst2
            >>> [ref.idx for ref in lst3]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ,14]

        :meta public:
        """
        self._checkSameType(other)
        s1 = set(self.lst)
        res = self._getLst() + [j for j in other.lst if j not in s1]
        return self.__class__(res, **self._cloneArgs)

    def __sub__(self, other):
        """Compute the substraction of one list from another

        :returns: All elements that are in the left operand but are absent from the right operand.
        :rtype: :py:class:`RefList`

        Usage::

            >>> lst1 = RefList(range(0, 10), mesh=mesh)
            >>> lst2 = RefList(range(5, 15), mesh=mesh)
            >>> lst3 = lst1 - lst2
            >>> [ref.idx for ref in lst3]
            [0, 1, 2, 3, 4]
            >>> lst4 = lst2 - lst1
            >>> [ref.idx for ref in lst4]
            [10, 11, 12, 13, 14]

        :meta public:
        """
        self._checkSameType(other)
        s2 = set(other.lst)
        res = [i for i in self.lst if i not in s2]
        return self.__class__(res, **self._cloneArgs)

    def __xor__(self, other):
        """Compute the symetric difference of two lists

        :returns: All elements in the left operand that are not in the right operand followed
            by all elements of the right operand that are not in the left operand.
        :rtype: :py:class:`RefList`

        Usage::

            >>> lst1 = RefList(range(0, 10), mesh=mesh)
            >>> lst2 = RefList(range(5, 15), mesh=mesh)
            >>> lst3 = lst1 ^ lst2
            >>> [ref.idx for ref in lst3]
            [0, 1, 2, 3, 4, 10, 11, 12, 13, 14]
            >>> lst4 = lst2 ^ lst1
            >>> [ref.idx for ref in lst4]
            [10, 11, 12, 13, 14, 0, 1, 2, 3, 4]

        :meta public:
        """
        self._checkSameType(other)
        s1, s2 = set(self.lst), set(other.lst)
        res = [i for i in self.lst if i not in s2] + [j for j in other.lst if j not in s1]
        return self.__class__(res, **self._cloneArgs)

    def __add__(self, other):
        """Concatenate two lists

        :returns: All elements of the left operand followed by all elements of the right operand.
        :rtype: :py:class:`RefList`

        Usage::

            >>> lst1 = RefList(range(0, 10), mesh=mesh)
            >>> lst2 = RefList(range(5, 15), mesh=mesh)
            >>> lst3 = lst1 + lst2
            >>> [ref.idx for ref in lst3]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
            >>> lst4 = lst2 + lst1
            >>> [ref.idx for ref in lst4]
            [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        :meta public:
        """
        self._checkSameType(other)
        return self.__class__(self._getLst() + other._getLst(), **self._cloneArgs)

    def __eq__(self, other):
        """Test for the equality of two lists

        :returns: True if the lists have the same elements in the same order, False otherwise
        :rtype: bool

        Usage::

            >>> lst1 = RefList(range(0, 5), mesh=mesh)
            >>> lst2 = RefList(range(0, 10), mesh=mesh)
            >>> lst1 == lst2
            False
            >>> lst1 == lst2[0:5]
            True

        :meta public:
        """
        try:
            self._checkSameType(other)
            return len(self) == len(other) and tuple(self.lst) == tuple(other.lst)
        except Exception:
            return False

    @property
    def indices(self):
        """The indices of elements in the list

        :type: List[int], read-only
        """
        return copy.copy(self.lst)

    def __hash__(self):
        return hash((len(self), tuple(self.lst)))

    def __repr__(self):
        if self._autoNamed:
            # If the list does not have a given name, represent the list with a hash of its contents
            return f'{self.__class__.__name__}#{hash(self)}'
        else:
            return self.name

    def _solverId(self):
        """Return the id that is used to identify an object in the solver."""
        return (self._getLst(),)

    def _simPathWalkExpand(self):
        """Return an iterable of the elements that should be part of ResultSelector paths."""
        return iter(self)


class _DistRefList(RefList):
    """Base class for distributed element lists

    Distributed lists can either be constituted of local or global indices
    """

    _FACADE_TITLE_STR = 'Only available for lists defined in distributed tetrahedral meshes'

    _LocalList = 'LocalList'
    _GlobalList = 'GlobalList'

    def __init__(self, lst=None, mesh=None, _immutable=False, *args, local=False, _owned=None, **kwargs):
        super().__init__(lst, mesh, _immutable, *args, **kwargs)

        if isinstance(lst, _DistRefList):
            local = lst.isLocal()

        self._local = local
        self._cloneArgs['local'] = self._local

        # The owned flag is just indicative of how the list was created
        self._owned = _owned
        self._cloneArgs['_owned'] = self._owned

    def isLocal(self):
        """Return whether the list is composed of local element indices

        :returns: ``True`` if the list is local
        :rtype: bool
        """
        return self._local

    def toLocal(self, owned=True, _returnInds=False):
        """Get the local version of a distributed list

        If the list is local, it returns itself. Otherwise, the orginal list is filtered and only
        the elements that are hosted by the current MPI process are added to the local list, with
        their local index. The returned list is thus different for each MPI process.

        :param owned: Whether only owned elements should be considered (defaults to True)
        :type owned: bool

        :returns: A local version of the list
        :rtype: :py:class:`RefList`
        """
        # if _returnInds is True, also return the positions of the local elements in the original list
        if self._local:
            if owned:
                # Check that the returned list does not contain any non-owned elements
                return self.toGlobal().toLocal(owned=True, _returnInds=_returnInds)
            else:
                return (self, list(range(len(self)))) if _returnInds else self
        cmkey = (_DistRefList._LocalList, owned)
        if cmkey not in self._optimCM:
            method = getattr(self.mesh.stepsMesh, f'get{self._refCls._locStr}LocalIndex')
            localInds = []
            listInds = []
            for i, gInd in enumerate(self.indices):
                localInd = method(gInd, owned=owned)
                if localInd is not None:
                    localInds.append(localInd)
                    listInds.append(i)
            self._optimCM[cmkey] = (
                self.__class__(localInds, mesh=self.mesh, local=True, _owned=owned),
                listInds,
            )
        if _returnInds:
            return self._optimCM[cmkey]
        else:
            return self._optimCM[cmkey][0]

    def toGlobal(self, **_kwargs):
        """Get the global version of a distributed list

        If the list is not local, it returns itself. Otherwise, the local indices of the original
        list are converted to global mesh indices and the list of these global indices is returned.

        :returns: A global version of the list
        :rtype: :py:class:`RefList`
        """
        if not self._local:
            return self
        if _DistRefList._GlobalList not in self._optimCM:
            method = getattr(self.mesh.stepsMesh, f'get{self._refCls._locStr}GlobalIndex')
            self._optimCM[_DistRefList._GlobalList] = self.__class__(
                [ind for ind in map(method, self.indices) if ind is not None],
                mesh=self.mesh,
                local=False,
                **_kwargs
            )
        return self._optimCM[_DistRefList._GlobalList]

    def _checkSameType(self, other):
        if self._local != other._local:
            raise TypeError('Cannot combine local and global index lists.')
        return super()._checkSameType(other)

    def combineWithOperator(self, binOp):
        """Combine the local versions of this list across ranks with the given operator

        :param binOp: The binary operator used to combine the global lists
        :type binOp: Callable[[:py:class:`RefList`, :py:class:`RefList`], :py:class:`RefList`]
        :returns: The global element list resulting from the combination of local lists
        :rtype: :py:class:`RefList`
        """
        if not self._local:
            raise Exception(f'The element list is already global, cannot combine it across ranks.')

        lists = self.mesh._comm.gather(self.toGlobal().indices, root=0)

        if self.mesh._comm.Get_rank() == 0:
            allIdx = self.__class__(lists[0], mesh=self.mesh, local=False)
            for lst in lists[1:]:
                allIdx = binOp(allIdx, self.__class__(lst, mesh=self.mesh, local=False))
        else:
            allIdx = self.__class__([], mesh=self.mesh, local=False)

        lst = self.mesh._comm.bcast(allIdx.indices, root=0)

        return self.__class__(lst, mesh=self.mesh, local=False)

    def _solverKeywordParams(self):
        """Return the additional keyword parameters that should be passed to the solver"""
        return dict(local=self._local)


@nutils.FreezeAfterInit
class TetList(RefList):
    """Convenience class for tetrahedron lists

    :param lst: The list of tetrahedrons
    :type lst: Iterable[int] or range or Iterable[:py:class:`TetReference`] or :py:class:`TetList`
    :param mesh: Mesh object that contains the elements
    :type mesh: Union[:py:class:`TetMesh`, :py:class:`DistMesh`, None]

    If `lst` is omitted, the list will be initialized to an empty list.

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).

    Behaves like a :py:class:`RefList` but with additional properties specific to tetrahedrons.
    Note that :py:func:`__getitem__` and :py:func:`__contains__` are modified to accept a 3D
    point as parameter.
    """

    _refCls = TetReference
    _dim = 3

    def __init__(self, lst=None, mesh=None, _immutable=False, *args, **kwargs):
        super().__init__(lst, mesh, _immutable, *args, **kwargs)

    def __getitem__(self, key):
        """Access one or several tetrahedron(s) from the list

        :param key: Index of the tetrahedron, or slice, or 3D point
        :type key: int or slice or Tuple[float, float, float] or :py:class:`Point`
        :returns: The tetrahedron(s)
        :rtype: :py:class:`TetReference` or :py:class:`TetList`

        Usage::

            >>> lst = TetList(range(10, 20), mesh=mesh)
            >>> lst[0] # Returns the 1st element from the list
            Tet(10)
            >>> lst[0].idx
            10
            >>> lst2 = lst[1:5] # Returns elements between the 2nd and the 5th (included).
            >>> lst2[0]
            Tet(11)
            >>> mesh.tets[0, 0, 0] # Returns the tetrahedron which contains the point x=0, y=0, z=0
            Tet(123)

        .. note::
            The index of an element in the list is in general different from its index in the mesh.
            :py:class:`TetList` is the only subclass of :py:class:`RefList` that can take a 3D point
            as key of its __getitem__ special method.

        :meta public:
        """
        try:
            return super().__getitem__(key)
        except TypeError:
            if hasattr(key, '__iter__') and len(key) == 3:
                idx = self._findTetIdxByPoint(list(key))
                if idx != UNKNOWN_TET and idx in self.lst:
                    return TetReference(idx, **self._cloneArgs)
                else:
                    raise KeyError(f'No Tetrahedron exists at position {key} in this list.')
            else:
                raise

    def __contains__(self, key):
        """Check whether the list contains a tetrahedron

        :param key: The element or a 3D point
        :type key: :py:class:`TetReference` or Tuple[float, float, float] or :py:class:`Point`
        :rtype: bool

        In addition to the arguments accepted in :py:func:`RefList.__contains__`, it is possible
        to check whether a 3D point is inside one tetrahedron in the list.

        Usage::

            >>> ref1 = TetReference(15, mesh=mesh)
            >>> ref2 = TetReference(25, mesh=mesh)
            >>> lst = TetList(range(20), mesh=mesh)
            >>> ref1 in lst
            True
            >>> ref2 in lst
            False
            >>> Point(0, 0, 0) in lst
            False
            >>> Point(1e-6, 0, 5e-6) in lst
            True

        :meta public:
        """
        if isinstance(key, Reference):
            return super().__contains__(key)
        elif hasattr(key, '__iter__') and len(key) == 3:
            try:
                self[key]
                return True
            except Exception:
                return False
        return False

    _CurrShell_K = 'CurrShell'

    def dilate(self, d=1):
        """Dilates the list to tetrahedron neighbors

        One cycle of dilation corresponds to adding all tetrahedrons that are neighbors of the
        current tetrahedron list but are not part of the list. This operation is repeated ``d``
        times.

        :param d: Topological distance to grow, defaults to 1
        :type d: int
        """

        prevShell = set(self.indices)

        if TetList._CurrShell_K not in self._optimCM:
            self._optimCM[TetList._CurrShell_K] = TetList(
                [tet for tet in self if any(neighb not in self for neighb in tet.neighbs)],
                **self._cloneArgs
            )
        shell = self._optimCM[TetList._CurrShell_K]

        for i in range(d):
            newShell = set()
            for tet in shell:
                for neighb in tet.neighbs:
                    if neighb._idx not in prevShell:
                        newShell.add(neighb._idx)
            prevShell.update(newShell)
            newShell = sorted(newShell)
            self.lst = self._getLst() + newShell
            shell = TetList(newShell, **self._cloneArgs)

        with self._modify():
            self._optimCM[TetList._CurrShell_K] = shell
            self.lst = sorted(set(self.lst))

    def erode(self, d=1):
        """Erodes the list, removing surface tetrahedrons

        One cycle of erosion corresponds to removing all tetrahedrons that are on the surface of
        the list (i.e. at least one of their four neighbors is not in the list). This operation is
        repeated ``d`` times.

        :param d: Topological distance to erode, defaults to 1
        :type d: int
        """
        shell = TetList(
            [
                tet
                for tet in self
                if any(neighb not in self for neighb in tet.neighbs) or len(tet.neighbs) < 4
            ],
            **self._cloneArgs,
        )
        shell.dilate(d - 1)
        with self._modify():
            self.lst = [tet._idx for tet in self if tet not in shell]

    _CurrSurface_K = 'CurrSurface'
    _LastSurface_K = 'LastSurface'

    @property
    def surface(self):
        """Surface of the tetrahedrons in the list

        I.e. all triangles that are not shared by two tetrahedrons in the list.

        :type: :py:class:`TriList`, read-only
        """
        if TetList._CurrSurface_K not in self._optimCM:
            if TetList._LastSurface_K in self._optimCM:
                # Only compute the changes due to the addition or removal of tetrahedrons
                oldLst, oldSurf = self._optimCM[TetList._LastSurface_K]
                surf = set(oldSurf)
                for tet in TetList(oldLst, **self._cloneArgs) ^ self:
                    surf ^= set(tet.faces)
            else:
                # Compute the surface from scratch
                surf = set()
                for tet in sorted(set(self), key=lambda x:x._idx):
                    surf ^= set(tet.faces)
            # Keep the surface until next list change
            self._optimCM[TetList._CurrSurface_K] = TriList(sorted(surf, key=lambda x: x._idx), **self._cloneArgs)
            # Keep a version of the current list and surface, for diff computation later
            self._optimCM[TetList._LastSurface_K] = (copy.copy(self.lst), self._optimCM[TetList._CurrSurface_K])
        return self._optimCM[TetList._CurrSurface_K]

    @property
    def Vol(self):
        """The summed volume of all tetrahedrons in the list

        .. note::
            If the list contains duplicate elements, they will be counted only once in the
            total volume.

        :type: float, read-only
        """
        return sum(tet.Vol for tet in sorted(set(self), key=lambda t: t._idx))

    @property
    def tris(self):
        """All faces of all tetrahedrons in the list

        :type: :py:class:`TriList`, read-only
        """
        tris = TriList([], **self._cloneArgs)
        for tet in self:
            tris |= tet.faces
        return tris

    @property
    def bars(self):
        """All bars of all tetrahedrons in the list

        :type: :py:class:`BarList`, read-only
        """
        return self.tris.bars

    @property
    def verts(self):
        """All vertices of all tetrahedrons in the list

        :type: :py:class:`VertList`, read-only
        """
        return self.tris.verts

    @nutils.NamedObject.children.getter
    def children(self):
        try:
            allComps = set()
            for elem in self:
                c = elem.comp
                if c is not None:
                    allComps.add(c)
            allChildren = {}
            for c in allComps:
                allChildren.update(c.children)
            return allChildren
        except AttributeError as ex:
            # Since children is called in __getattr__, raising an AttributeError here can lead to
            # infinite recursion.
            raise Exception(ex)

    def _findTetIdxByPoint(self, pos):
        return self.mesh.stepsMesh.findTetByPoint(pos)


class _DistTetList(TetList, _DistRefList):
    """Convenience class for distributed tetrahedron list
    """

    @TetList.surface.getter
    def surface(self):
        if self._local:
            return TetList.surface.fget(self)
        else:
            return self.toLocal().surface.combineWithOperator(operator.xor)

    @TetList.tris.getter
    def tris(self):
        if self._local:
            return TetList.tris.fget(self)
        else:
            return self.toLocal().tris.combineWithOperator(operator.or_)

    @TetList.verts.getter
    def verts(self):
        if self._local:
            verts = VertList([], **self._cloneArgs)
            for tet in self:
                verts |= tet.verts
            return verts
        else:
            return self.toLocal().verts.combineWithOperator(operator.or_)

    def _findTetIdxByPoint(self, pos):
        return self.mesh.stepsMesh.findTetByPoint(pos, local=self._local)


TetList._distCls = _DistTetList


@nutils.FreezeAfterInit
class TriList(RefList):
    """Convenience class for triangle lists

    :param lst: The list of triangles
    :type lst: Iterable[int] or range or Iterable[:py:class:`TriReference`] or :py:class:`TriList`
    :param mesh: Mesh object that contains the elements
    :type mesh: Union[:py:class:`TetMesh`, :py:class:`DistMesh`, None]

    If `lst` is omitted, the list will be initialized to an empty list.

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).

    Behaves like a :py:class:`RefList` but with additional properties specific to triangles.
    """

    _refCls = TriReference
    _dim = 2

    def __init__(self, lst=None, mesh=None, _immutable=False, *args, **kwargs):
        super().__init__(lst if lst is not None else [], mesh, _immutable, *args, **kwargs)

    @property
    def edges(self):
        """Perimeter of the triangles in the list

        I.e. all bars that are not shared by two triangles in the list.

        :type: :py:class:`BarList`, read-only
        """
        edg = set()
        for tri in self:
            edg ^= set(tri.bars)
        return BarList(sorted(edg, key=lambda x: x._idx), **self._cloneArgs)

    @property
    def bars(self):
        """All bars of all triangles in the list

        :type: :py:class:`BarList`, read-only
        """
        bars = BarList([], **self._cloneArgs)
        for tri in self:
            bars |= tri.bars
        return bars

    @property
    def verts(self):
        """All vertices of all triangles in the list

        :type: :py:class:`VertList`, read-only
        """
        return self.bars.verts

    @property
    def Area(self):
        """The summed area of all triangles in the list

        .. note::
            If the list contains duplicate elements, they will be counted only once in the
            total area.

        :type: float, read-only
        """
        return sum(tri.Area for tri in sorted(set(self), key=lambda t: t._idx))

    @nutils.NamedObject.children.getter
    def children(self):
        try:
            allPatches = set()
            for elem in self:
                p = elem.patch
                if p is not None:
                    allPatches.add(p)
            allChildren = {}
            for p in allPatches:
                allChildren.update(p.children)
            return allChildren
        except AttributeError as ex:
            # Since children is called in __getattr__, raising an AttributeError here can lead to
            # infinite recursion.
            raise Exception(ex)


class _DistTriList(TriList, _DistRefList):

    @TetList.verts.getter
    def verts(self):
        if self._local:
            verts = VertList([], **self._cloneArgs)
            for tri in self:
                verts |= tri.verts
            return verts
        else:
            return self.toLocal().verts.combineWithOperator(operator.or_)


TriList._distCls = _DistTriList


@nutils.FreezeAfterInit
class BarList(RefList):
    """Convenience class for bar lists

    :param lst: The list of bars
    :type lst: Iterable[int] or range or Iterable[:py:class:`BarReference`] or :py:class:BarList`
    :param mesh: Mesh object that contains the elements
    :type mesh: Union[:py:class:`TetMesh`, :py:class:`DistMesh`, None]

    If `lst` is omitted, the list will be initialized to an empty list.

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).

    Behaves like a :py:class:`RefList` but with additional properties specific to bars.
    """

    _refCls = BarReference
    _dim = 1

    def __init__(self, lst=None, mesh=None, _immutable=False, *args, **kwargs):
        super().__init__(lst if lst is not None else [], mesh, _immutable, *args, **kwargs)

    @property
    def verts(self):
        """All vertices of all bars in the list

        :type: :py:class:`VertList`, read-only
        """
        verts = VertList([], **self._cloneArgs)
        for bar in self:
            verts |= bar.verts
        return verts


class _DistBarList(BarList, _DistRefList):
    pass


BarList._distCls = _DistBarList


@nutils.FreezeAfterInit
class VertList(RefList):
    """Convenience class for vertex lists

    :param lst: The list of vertices
    :type lst: Iterable[int] or range or Iterable[:py:class:`VertReference`] or :py:class:VertList`
    :param mesh: Mesh object that contains the elements
    :type mesh: Union[:py:class:`TetMesh`, :py:class:`DistMesh`, None]

    If `lst` is omitted, the list will be initialized to an empty list.

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).

    Behaves like a :py:class:`RefList`.
    """

    _refCls = VertReference
    _dim = 0

    def __init__(self, lst=None, mesh=None, _immutable=False, *args, **kwargs):
        super().__init__(lst if lst is not None else [], mesh, _immutable, *args, **kwargs)


class _DistVertList(VertList, _DistRefList):
    def __init__(self, *args, _fromLocal=False, **kwargs):
        super().__init__(*args, **kwargs)
        self._cloneArgs['_fromLocal'] = _fromLocal


VertList._distCls = _DistVertList


# Add a link from the reference classes to the list class
for cls in RefList.__subclasses__():
    if cls != _DistRefList:
        cls._refCls._lstCls = cls
# Do not leave cls as a global variable
del cls


@nutils.FreezeAfterInit
class ROI(nutils.UsingObjects(_BaseTetMesh), nutils.StepsWrapperObject):
    """Region of Interest class

    :param lst: Element list
    :type lst: :py:class:`TetList`, or :py:class:`TriList`, or :py:class:`VertList`

    Should be declared using the context manager syntax with :py:func:`TetMesh`::

        with mesh:
            tetROI1 = ROI.Create(mesh.tets[10:30])

    The ROI is then treated as a physical location, like patches or compartment,
    and can be accessed from the mesh with its name::

        mesh.tetROI1
    """

    _locStr = 'ROI'

    def __init__(self, lst=None, *args, _createObj=True, **kwargs):
        super().__init__(*args, **kwargs)
        (mesh,) = self._getUsedObjects()
        self.mesh = mesh
        if _createObj:
            self.stepsROI, self._elemType = self._createStepsObj(mesh, lst)
        else:
            self.stepsROI, self._elemType = None, None
        if self._elemType is not None and lst is not None:
            self.lst = self._elemType._lstCls._toRefList(lst, mesh)
        else:
            self.lst = None

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsROI]

    @classmethod
    def _FromStepsObject(cls, obj, mesh, name):
        """Create the interface object from a STEPS object."""
        revTypeMap = {
            ELEM_TET: TetReference,
            ELEM_TRI: TriReference,
            ELEM_VERTEX: VertReference,
        }
        region = cls(_createObj=False, name=name)
        region.stepsROI = obj
        region._elemType = revTypeMap[obj.type]
        region.lst = region._elemType._lstCls._toRefList(obj.indices, mesh)
        return region

    def _createStepsObj(self, mesh, lst):
        typeMap = {
            TetReference: ELEM_TET,
            TriReference: ELEM_TRI,
            VertReference: ELEM_VERTEX,
        }
        if len(lst) == 0:
            raise Exception('An empty list was given.')
        elem = next(iter(lst))
        tpe = None
        for cls, _tpe in typeMap.items():
            if isinstance(elem, cls):
                tpe = _tpe
        if tpe is None:
            raise TypeError('Unsupported element type for creating an ROI.')
        mesh.stepsMesh.addROI(self.name, tpe, lst.indices)
        return mesh.stepsMesh.getROI(self.name), elem.__class__

    def __getitem__(self, key):
        """
        Indexing can only be used with ':' or '...' to signify that we want to consider all
        elements separately instead of the full ROI.
        """
        if key == slice(None) or key is Ellipsis:
            return self.lst
        else:
            raise KeyError(
                f'ROIs can only be indexed with an empty slice ":" or the ellipsis operator "...".'
            )

    @property
    def Vol(self):
        """Total volume of the ROI, if applicable

        :type: float, read-only
        """
        return self.mesh.stepsMesh.getROIVol(self.name)

    @property
    def Area(self):
        """Total area of the ROI, if applicable

        :type: float, read-only
        """
        return self.mesh.stepsMesh.getROIArea(self.name)

    def _checkType(self, tpe):
        if self._elemType != tpe:
            raise AttributeError(f'ROI {self} is not composed of {tpe}.')

    @property
    def tets(self):
        """All tetrahedrons in the ROI

        Raises a :py:class:`TypeError` if the ROI was not created from a :py:class:`TetList`.

        :type: :py:class:`TetList`, read-only
        """
        self._checkType(TetReference)
        return self.lst

    @property
    def tris(self):
        """All triangles in the ROI

        Raises a :py:class:`TypeError` if the ROI was not created from a :py:class:`TriList`.

        :type: :py:class:`TriList`, read-only
        """
        self._checkType(TriReference)
        return self.lst

    @property
    def verts(self):
        """All vertices in the ROI

        Raises a :py:class:`TypeError` if the ROI was not created from a :py:class:`VertList`.

        :type: :py:class:`VertList`, read-only
        """
        self._checkType(VertReference)
        return self.lst

    @property
    def children(self):
        # Handle the references to species etc.
        for elem in self.lst:
            try:
                return elem._getPhysicalLocation().children
            except:
                continue
        return {}

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return ROI._locStr

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        return {'loc_type': self.__class__._locStr, 'loc_id': self.name}


@nutils.FreezeAfterInit
class DiffBoundary(nutils.UsingObjects(_BaseTetMesh), nutils.StepsWrapperObject):
    """Diffusion Boundary class

    Annotation of a group of triangles in a Tetmesh. The triangles form a boundary between two
    compartments, that may allow diffusion of some specified species.

    :param lst: List of triangles forming the boundary
    :type lst: :py:class:`TriList`
    """

    _locStr = 'DiffBoundary'
    _lstType = TriList

    _direcStr = 'direction_comp'

    def __init__(self, lst, *args, **kwargs):
        super().__init__(*args, **kwargs)
        (mesh,) = self._getUsedObjects()
        self.lst = self.__class__._lstType(lst, mesh)
        self.stepsDiffBound = self._createStepsObj(mesh)

        self._direc = None

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsDiffBound]

    def _createStepsObj(self, mesh):
        if isinstance(mesh, TetMesh):
            return stepslib._py_DiffBoundary(self.name, mesh.stepsMesh, [tri._idx for tri in self.lst])
        elif isinstance(mesh, DistMesh):
            tet1, tet2 = self.lst[0].tetNeighbs
            comp1 = tet1.comp.name
            comp2 = tet2.comp.name
            mesh._getStepsObjects()[0].addDiffusionBoundary(self.name, comp1, comp2, self.lst.indices)
            return None

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return DiffBoundary._locStr

    def _solverKeywordParams(self):
        """Return the additional keyword parameters that should be passed to the solver"""
        if self._direc is not None:
            return {DiffBoundary._direcStr: self._direc._solverId()[0]}
        else:
            return {}

    def __call__(self, direc=None):
        """Get a version of the diffusion boundary with added information

        Parentheses notation (function call) can be used on a diffusion boundary to represent a
        diffusion boundary further specified with additional information. This should never be needed
        during model declaration but can become handy for simulation control and data saving
        (see :py:class:`steps.API_2.sim.SimPath`).

        :param direc: A direction for the diffusion boundary (i.e. a target compartment for
            diffusion boundaries and a target patch for surface diffusion boundary)
        :type key: :py:class:`Comp`

        :returns: A version of the diffusion boundary that is specific to direction ``direc``.

        :meta public:
        """
        db = copy.copy(self)
        db._direc = direc
        return db

    @property
    def tris(self):
        """All triangles in the boundary, if applicable

        :type: :py:class:`TriList`, read-only
        """
        if not isinstance(self.lst, self._lstType):
            raise NotImplementedError()
        return self.lst

    @property
    def bars(self):
        """:meta private:"""
        raise NotImplementedError()

    @property
    def children(self):
        # Handle the references to species.
        for tet in self.lst[0].tetNeighbs:
            try:
                return tet.comp.children
            except:
                continue
        return {}


@nutils.FreezeAfterInit
class SDiffBoundary(DiffBoundary):
    """Surface Diffusion Boundary class

    Annotation of a group of bars in a Tetmesh. The bars form a boundary between two patches,
    that may allow diffusion of some specified species. The patches to be connected by this
    surface diffusion boundary need to be specified to avoid potential ambiguity.

    :param lst: List of bars forming the boundary
    :type lst: :py:class:`BarList`
    :param patch1: One of the patches connected to the boundary
    :type patch1: :py:class:`Patch`
    :param patch2: The other patch connected to the boundary
    :type patch2: :py:class:`Patch`
    """

    _locStr = 'SDiffBoundary'
    _lstType = BarList

    _direcStr = 'direction_patch'

    def __init__(self, lst, patch1, patch2, *args, **kwargs):
        self.patches = [patch1, patch2]
        super().__init__(lst, *args, **kwargs)

    def _createStepsObj(self, mesh):
        return stepslib._py_SDiffBoundary(
            self.name, mesh.stepsMesh, [bar._idx for bar in self.lst], [p.stepsPatch for p in self.patches]
        )

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return SDiffBoundary._locStr

    def _solverKeywordParams(self):
        """Return the additional keyword parameters that should be passed to the solver"""
        if self._direc is not None:
            return {SDiffBoundary._direcStr: self._direc._solverId()[0]}
        else:
            return {}

    @property
    def tris(self):
        """:meta private:"""
        raise NotImplementedError()

    @property
    def bars(self):
        """All bars in the boundary, if applicable

        :type: :py:class:`BarList`, read-only
        """
        if not isinstance(self.lst, self._lstType):
            raise NotImplementedError()
        return self.lst

    @property
    def children(self):
        # Handle the references to species.
        for patch in self.patches:
            try:
                return patch.children
            except:
                continue
        return {}


###################################################################################################
# Mesh partitioning


@nutils.FreezeAfterInit
class MeshPartition(nutils.Params):
    """Mesh element partition for parallel simulations

    Wrapper class that groups tetrahedrons, triangles and well-mixed partitions.

    :param mesh: Mesh being partitioned
    :type mesh: :py:class:`TetMesh`
    :param tet_hosts: List of host process indices for each tetrahedron
    :type tet_hosts: List[int]
    :param tri_hosts: Mapping between triangle index and host process index
    :type tri_hosts: Mapping[int, int]
    :param wm_hosts: List of host process indices for each well mixed compartment
    :type wm_hosts: List[int]
    """

    def __init__(self, mesh, tet_hosts=[], tri_hosts={}, wm_hosts=[]):
        super().__init__(tet_hosts=tet_hosts, tri_hosts=tri_hosts, wm_hosts=wm_hosts)
        self._mesh = mesh
        self._tet_hosts = tet_hosts
        self._tri_hosts = tri_hosts
        self._wm_hosts = wm_hosts

    @property
    def tetPart(self):
        """The partition data for tetrahedrons

        :type: List[int], read-only
        """
        return self._tet_hosts

    @property
    def triPart(self):
        """The partition data for triangles

        :type: Dict[int, int], read-only
        """
        return self._tri_hosts

    @property
    def wmPart(self):
        """The partition data for well mixed locations

        :type: List[int], read-only
        """
        return self._wm_hosts

    @property
    def tetTable(self):
        """A dictionary in the format of {host0: TetList([tet0, tet1, ...]), ...}

        :type: Dict[int, :py:class:`TetList`], read-only
        """
        return {
            h: TetList(lst, mesh=self._mesh)
            for h, lst in sgdecomp.getTetPartitionTable(self._tet_hosts).items()
        }

    @property
    def triTable(self):
        """A dictionary in the format of {host0: TriList([tri0, tri1, ...]), ...}

        :type: Dict[int, :py:class:`TriList`], read-only
        """
        return {
            h: TriList(lst, mesh=self._mesh)
            for h, lst in sgdecomp.getTriPartitionTable(self._tri_hosts).items()
        }

    def validate(self):
        """Validate the partitioning of the mesh"""
        sgdecomp.validatePartition(self._mesh.stepsMesh, self._tet_hosts, self._tri_hosts)

    def printStats(self):
        """Print out partitioning stastics

        :returns: ``(tet_stats, tri_stats, wm_stats, num_hosts, min_degree, max_degree, mean_degree)``
            [tet/tri/wm]_stats contains the number of tetrahedrons/triangles/well-mixed volumes in
            each hosting process, num_hosts provide the number of hosting processes,
            [min/max/mean]_degree provides the minimum/maximum/average connectivity degree of the
            partitioning.
        """
        return sgdecomp.printPartitionStat(
            self._tet_hosts, self._tri_hosts, self._wm_hosts, self._mesh.stepsMesh
        )


def _getTriPartitionFromTet(mesh, tet_hosts, default_tris=None):
    """Return the tri partition corresponding to the tet partition given as an argument."""
    triInds = set(default_tris.indices) if default_tris is not None else set()
    for patch in mesh._getChildrenOfType(Patch):
        triInds |= set(patch.tris.indices)
    return sgdecomp.partitionTris(mesh.stepsMesh, tet_hosts, sorted(triInds))


def LinearMeshPartition(mesh, xbin, ybin, zbin, default_tris=None):
    """Partition the mesh by linearly binning along the 3 axes

    First partition the tetrahedrons and then compute a triangle partition that
    matches the tetrahedron one.

    :param mesh: The mesh to be partitioned
    :type mesh: :py:class:`TetMesh`
    :param xbin: Number of bins on the x axis
    :type xbin: int
    :param ybin: Number of bins on the y axis
    :type ybin: int
    :param zbin: Number of bins on the z axis
    :type zbin: int
    :param default_tris: Optional list of triangles that should be partitioned even if they are
        not part of any patch
    :type default_tris: :py:class:`TriList`

    :returns: The partition object
    :rtype: :py:class:`MeshPartition`

    .. note::
        The triangle partition only takes into account triangles that are part of a
        :py:class:`Patch`.
    """
    tet_hosts = sgdecomp.linearPartition(mesh.stepsMesh, (xbin, ybin, zbin))
    tri_hosts = _getTriPartitionFromTet(mesh, tet_hosts, default_tris)
    return MeshPartition(mesh, tet_hosts=tet_hosts, tri_hosts=tri_hosts)


def MetisPartition(mesh, path, default_tris=None):
    """Partition the mesh using a Metis .epart file

    First partition the tetrahedrons according to the Metis file and then compute a triangle
    partition that matches the tetrahedron one.

    :param mesh: The mesh to be partitioned
    :type mesh: :py:class:`TetMesh`
    :param path: Path to the .epart file generated by Metis
    :type path: str
    :param default_tris: Optional list of triangles that should be partitioned even if they are
        not part of any patch
    :type default_tris: :py:class:`TriList`

    :returns: The partition object
    :rtype: :py:class:`MeshPartition`

    .. note::
        The triangle partition only takes into account triangles that are part of a
        :py:class:`Patch`.
    """
    tet_hosts = smetis.readPartition(path)
    tri_hosts = _getTriPartitionFromTet(mesh, tet_hosts, default_tris)
    return MeshPartition(mesh, tet_hosts=tet_hosts, tri_hosts=tri_hosts)

def GmshPartition(mesh, default_tris=None):
    """Retrive partition information stored in a pre-partitioned gmsh.

    The mesh needs to be pre-partitioned and stored in a single 
    Gmsh file with format version 2.2 ascii.

    :param mesh: The mesh to be partitioned
    :type mesh: :py:class:`TetMesh`
    :param default_tris: Optional list of triangles that should be partitioned even if they are
        not part of any patch
    :type default_tris: :py:class:`TriList`

    :returns: The partition object
    :rtype: :py:class:`MeshPartition`

    .. note::
        The triangle partition only takes into account triangles that are part of a
        :py:class:`Patch`.
    """
    tet_hosts = [None] * len(mesh.tets)
    for k, v in mesh.tetGroups.items():
        # partitioning tag is the forth one in gmsh v2.2 format
        # starting from 0
        if k[0] == 3:
            # gmsh partition indices start from 1
            # mpi ranks start from 0
            partition_rank = k[1] - 1
            partition_tets = v
            for tet_ref in partition_tets:
                tet_hosts[tet_ref._idx] = partition_rank
    if None in tet_hosts:
        raise Exception("Partition information not available for one or more tets.")

    tri_hosts = _getTriPartitionFromTet(mesh, tet_hosts, default_tris)
    return MeshPartition(mesh, tet_hosts=tet_hosts, tri_hosts=tri_hosts)

def MorphPartition(mesh, morph, scale=1e-6, default_tris=None):
    """Partition the mesh using morphological sectioning data

    First partition the tetrahedrons according to the morph data and then compute a triangle
    partition that matches the tetrahedron one.

    :param mesh: The mesh to be partitioned
    :type mesh: :py:class:`TetMesh`
    :param morph: Morphological sectioning data
    :type morph: :py:class:`Morph`
    :param scale: Scaling factor from morphological sectioning data to Tetmesh data,
        the default value is 1e-6, that is 1 unit (usually micron) in the morph file equals
        1e-6 unit (meter) in Tetmesh
    :type scale: float
    :param default_tris: Optional list of triangles that should be partitioned even if they are
        not part of any patch
    :type default_tris: :py:class:`TriList`

    :returns: The partition object and a mapping between partition index and morphological
        section name
    :rtype: Tuple[:py:class:`MeshPartition`, Dict[int, str]]
    """
    tet_map = smorph.mapMorphTetmesh(morph._sections, mesh.stepsMesh, scale)
    tet_hosts = []
    ind2sec = {}
    sec2ind = {}
    i = 0
    for name in tet_map:
        if name not in sec2ind:
            sec2ind[name] = i
            ind2sec[i] = name
            i += 1
        tet_hosts.append(sec2ind[name])

    tri_hosts = _getTriPartitionFromTet(mesh, tet_hosts, default_tris)
    return MeshPartition(mesh, tet_hosts=tet_hosts, tri_hosts=tri_hosts), ind2sec


###################################################################################################
# Morph sectioning


@nutils.FreezeAfterInit
class Morph:
    """Morphological sectioning data

    Contains data describing the sections of a NEURON morphology.

    This class should not be instantiated directly, it should be created using one of the dedicated
    class methods.
    """

    def __init__(self, sections):
        self._sections = sections

    @classmethod
    def Load(self, path):
        """Load morphological data from a morph file

        :param path: Path to the file
        :type path: str

        :returns: The morphological sectioning data
        :rtype: :py:class:`Morph`
        """
        with open(path, 'rb') as f:
            return Morph(pickle.load(f))

    @classmethod
    def LoadHOC(self, path):
        """Load morphological data from a NEURON hoc file

        :param path: Path to the file
        :type path: str

        :returns: The morphological sectioning data
        :rtype: :py:class:`Morph`
        """
        return Morph(smorph.hoc2morph(path))

    @classmethod
    def LoadSWC(self, path):
        """Load morphological data from an swc file

        :param path: Path to the file
        :type path: str

        :returns: The morphological sectioning data
        :rtype: :py:class:`Morph`
        """
        return Morph(smorph.swc2morph(path))

    def Save(self, path):
        """Save to a morph file

        :param path: Path to the file
        :type path: str
        """
        with open(path, 'wb') as f:
            pickle.dump(self._sections, f)
