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

import copy
import numbers
import pickle
import warnings

import numpy

import steps.API_1.geom as sgeom
from steps.API_1.geom import UNKNOWN_TET, UNKNOWN_TRI
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
    'ROI',
    'TetMesh',
    'DiffBoundary',
    'SDiffBoundary',
    'MeshPartition',
    'LinearMeshPartition',
    'MetisPartition',
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
    'Morph',
]

###################################################################################################


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
        self.stepsGeom = sgeom.Geom() if _createObj else None

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
            raise TypeError(f'Expected a VolumeSystem or a SurfaceSystem, got a {type(sys)} ' f'instead.')

    def _SetUpMdlDeps(self, mdl):
        """Set up the structures that will allow specie or reaction access from locations."""
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
        return isinstance(c, (nmodel.Reaction, nmodel.Diffusion, nmodel.Current, nmodel.Complex))

    def _createStepsObj(self, geom):
        pass

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        return {'loc_type': self.__class__._locStr, 'loc_id': self.name}


class Compartment(_PhysicalLocation):
    """Base class for compartment objects

    The same class is used to declare compartments in well-mixed models and in 3D tetrahedral
    meshes.

    Compartments in well-mixed geometry can be declared in the following ways::

        with geom:
            comp1 = Compartment.Create()
            comp2 = Compartment.Create(vsys)
            comp3 = Compartment.Create(None, vol)
            comp4 = Compartment.Create(vsys, vol)

    :param vsys: Optional, the volume system associated with this compartment.
    :type vsys: :py:class:`steps.API_2.model.VolumeSystem` or :py:class:`str`
    :param vol: Optional, volume of the compartment
    :type vol: float

    Compartments in a tetrahedral mesh can be declared in the following ways::

        with mesh:
            comp1 = Compartment.Create(tetLst)
            comp2 = Compartment.Create(tetLst, vsys)

    :param tetLst: List of tetrahedrons associated with the compartment.
    :type tetLst: :py:class:`TetList` or any argument that can be used to build one
    :param vsys: Optional, the volume system associated with this compartment.
    :type vsys: :py:class:`steps.API_2.model.VolumeSystem` or :py:class:`str`

    The volume of a compartment in a tetrahedral mesh is the total volume of the encapsulated
    tetrahedrons.
    """

    _locStr = 'Comp'

    def __init__(self, *args, _createObj=True, **kwargs):
        super().__init__(*args, **kwargs)
        (geom,) = self._getUsedObjects()
        self._isTmComp = isinstance(geom, TetMesh)

        if len(args) > 2:
            raise Exception(f'Unused parameters: {args[2:]}')
        elif len(args) < 2:
            args += (None,) * (2 - len(args))

        if self._isTmComp:
            tetLst, vsys, *args = args
            vol = None
        else:
            vsys, vol, *args = args
            tetLst = None

        self._tetLst = tetLst
        self.stepsComp = self._createStepsObj(geom) if _createObj else None
        if vsys is not None:
            self.addSystem(vsys)
        if vol is not None:
            self.Vol = vol

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsComp]

    @classmethod
    def _FromStepsObject(cls, obj, geom):
        """Create the interface object from a STEPS object."""
        comp = cls(_createObj=False, name=obj.getID())
        comp.stepsComp = obj
        # setup volume systems
        for vsysname in comp.stepsComp.getVolsys():
            comp.addSystem(vsysname)

        comp._isTmComp = isinstance(geom, TetMesh)
        comp._tetLst = TetList(obj.getAllTetIndices(), geom) if comp._isTmComp else None

        return comp

    def addSystem(self, vsys, _loc=None):
        """Add a volume system to the compartment

        :param vsys: The volume system to be added, or its name.
        :type vsys: :py:class:`steps.API_2.model.VolumeSystem` or `str`

        :returns: None
        """
        super().addSystem(vsys, _loc)
        if _loc is None:
            if isinstance(vsys, nmodel.SpaceSystem):
                vsys = vsys.name
            self.stepsComp.addVolsys(vsys)

    def _createStepsObj(self, geom):
        if self._isTmComp:
            self._tetLst = TetList._toRefList(self._tetLst, geom)
            return sgeom.TmComp(self.name, geom.stepsMesh, self._tetLst.indices)
        else:
            return sgeom.Comp(self.name, geom.stepsGeom)

    def _canBeChild(self, c):
        """Return wether c can be a child of self."""
        if not super()._canBeChild(c):
            return False
        if (isinstance(c, nmodel.Reaction) and c._isSurfaceReac()) or (
            isinstance(c, nmodel.Diffusion) and c._isSurfaceDiff()
        ):
            return False
        return True

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return Compartment._locStr

    @property
    def Vol(self):
        """Volume of the compartment

        :type: float
        """
        return self.stepsComp.getVol()

    @Vol.setter
    def Vol(self, v):
        if self._isTmComp:
            raise Exception(f'Cannot set the volume of a Tetmesh compartment.')
        self.stepsComp.setVol(v)

    @property
    def bbox(self):
        """The bounding box of the compartment

        Only available for compartments in tetrahedral meshes.

        :type: :py:class:`BoundingBox`, read-only
        """
        if not self._isTmComp:
            raise Exception(f'Well-mixed compartment do not have bounding boxes.')
        return BoundingBox(Point(*self.stepsComp.getBoundMin()), Point(*self.stepsComp.getBoundMax()))

    @property
    def tets(self):
        """All tetrahedrons in the compartment

        Only available for compartments in tetrahedral meshes.

        :type: :py:class:`TetList`, read-only
        """
        if not self._isTmComp:
            raise Exception(f'Well-mixed compartment do not have tetrahedrons.')
        return self._tetLst

    @property
    def surface(self):
        """All triangles in the surface of the compartment

        Only available for compartments in tetrahedral meshes.

        :type: :py:class:`TriList`, read-only
        """
        if not self._isTmComp:
            raise Exception(f'Well-mixed compartment do not have tetrahedrons.')
        return self._tetLst.surface


class Patch(_PhysicalLocation):
    """Base class for patch objects

    A patch is a piece of 2D surface surrounding (part of) a 3D compartment, which may be
    connected to another compartment. The same class is used to declare patches in well-mixed
    models and in 3D tetrahedral meshes.

    Patches in well-mixed geometry can be declared in the following ways::

        with geom:
            patch1 = Patch.Create(inner)
            patch2 = Patch.Create(inner, outer)
            patch3 = Patch.Create(inner, outer, vsys)
            patch4 = Patch.Create(inner, None , vsys)
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
    """

    _locStr = 'Patch'

    def __init__(self, *args, _createObj=True, **kwargs):
        super().__init__(*args, **kwargs)
        (geom,) = self._getUsedObjects()
        self._isTmPatch = isinstance(geom, TetMesh)

        if len(args) > 4:
            raise Exception(f'Unused parameters: {args[4:]}')
        elif len(args) < 4:
            args += (None, ) * (4 - len(args))

        if self._isTmPatch:
            triLst, inner, outer, ssys, *_ = args
            area = None
        else:
            inner, outer, ssys, area, *_ = args
            triLst = None

        self._triLst = triLst
        self._innerComp = inner
        self._outerComp = outer
        self.stepsPatch = self._createStepsObj(geom) if _createObj else None
        if ssys is not None:
            self.addSystem(ssys)
        if area is not None:
            self.Area = area

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsPatch]

    @classmethod
    def _FromStepsObject(cls, obj, geom):
        """Create the interface object from a STEPS object."""
        patch = cls(_createObj=False, name=obj.getID())
        patch.stepsPatch = obj
        if obj.icomp is not None:
            patch._innerComp = getattr(geom, obj.icomp.getID())
        if obj.ocomp is not None:
            patch._outerComp = getattr(geom, obj.ocomp.getID())

        # setup surface systems
        for ssysname in patch.stepsPatch.getSurfsys():
            patch.addSystem(ssysname)

        patch._isTmPatch = isinstance(geom, TetMesh)
        patch._triLst = TriList(obj.getAllTriIndices(), geom) if patch._isTmPatch else None

        return patch

    def addSystem(self, ssys):
        """Add a surface system to the patch

        :param ssys: The surface system to be added, or its name.
        :type ssys: :py:class:`steps.API_2.model.SurfaceSystem` or `str`

        :returns: None
        """
        super().addSystem(ssys, nmodel.Location.SURF)
        if isinstance(ssys, nmodel.SpaceSystem):
            ssys = ssys.name
        self.stepsPatch.addSurfsys(ssys)

    def _createStepsObj(self, geom):
        if not isinstance(self._innerComp, Compartment):
            raise TypeError(f'Expected an inner Compartment, got {self._innerComp} instead.')
        if self._outerComp is not None and not isinstance(self._outerComp, Compartment):
            raise TypeError(f'Expected a Compartment, got {self._outerComp} instead.')

        icomp = self._innerComp.stepsComp if isinstance(self._innerComp, Compartment) else None
        ocomp = self._outerComp.stepsComp if isinstance(self._outerComp, Compartment) else None
        if self._isTmPatch:
            if icomp is None and ocomp is None:
                raise Exception(
                    f'Cannot declare a Patch in a tetrahedral mesh with no neighboring compartment.'
                )
            self._triLst = TriList._toRefList(self._triLst, geom)
            return sgeom.TmPatch(self.name, geom.stepsMesh, self._triLst.indices, icomp, ocomp)
        else:
            return sgeom.Patch(self.name, geom.stepsGeom, icomp, ocomp)

    def _canBeChild(self, c):
        """Return wether c can be a child of self."""
        if not super()._canBeChild(c):
            return False
        if (isinstance(c, nmodel.Reaction) and not c._isSurfaceReac()) or (
            isinstance(c, nmodel.Diffusion) and not c._isSurfaceDiff()
        ):
            return False
        return True

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return Patch._locStr

    @property
    def Area(self):
        """Surface area of the patch

        :type: float
        """
        return self.stepsPatch.getArea()

    @Area.setter
    def Area(self, v):
        if self._isTmPatch:
            raise Exception(f'Cannot set the area of a Tetmesh patch.')
        self.stepsPatch.setArea(v)

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

    @property
    def bbox(self):
        """The bounding box of the patch

        Only available for patches in tetrahedral meshes.

        :type: :py:class:`BoundingBox`, read-only
        """
        if not self._isTmPatch:
            raise Exception(f'Well-mixed patches do not have bounding boxes.')
        return BoundingBox(Point(*self.stepsPatch.getBoundMin()), Point(*self.stepsPatch.getBoundMax()))

    @property
    def tris(self):
        """All triangles in the patch

        Only available for patches in tetrahedral meshes.

        :type: :py:class:`TriList`, read-only
        """
        if not self._isTmPatch:
            raise Exception(f'Well-mixed patches do not have triangles.')
        return self._triLst

    @property
    def edges(self):
        """Perimeter of the triangles in the patch

        I.e. all bars that are not shared by two triangles in the patch.

        Only available for patches in tetrahedral meshes.

        :type: :py:class:`BarList`, read-only
        """
        if not self._isTmPatch:
            raise Exception(f'Well-mixed patches do not have triangles.')
        return self._triLst.edges


class Membrane(_PhysicalLocation):
    """A set of patches on which membrane potential will be computed

    This class provides annotation for a group of triangles that comprise
    a surface describing a membrane in a Tetmesh. This may be the same
    description of one or several TmPatches in order for voltage-dependent
    transitions, currents and so on to be inserted in the membrane.
    A Membrane object must be available if membrane potential calculation is to be
    performed.

    :param patches: List of patches whose triangles will be added to the membrane
    :type patches: :py:class:`list`
    :param verify: Perform some checks on suitability of membrane (see below)
    :type verify: :py:class:`bool`
    :param opt_method: Optimization method (see below)
    :type opt_method: int
    :param search_percent: See below
    :type search_percent: float
    :param opt_file_name: Full path to a membrane optimization file
    :type opt_file_name: str

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

    def __init__(
        self, patches, verify=False, opt_method=1, search_percent=100.0, opt_file_name='', *args, **kwargs
    ):
        super().__init__(*args, **kwargs)
        (geom,) = self._getUsedObjects()
        self.mesh = geom
        self.patches = patches
        self.stepsMemb = self._createStepsObj(geom, verify, opt_method, search_percent, opt_file_name)

    def addSystem(self, _):
        """
        Even though a Membrane is a PhysicalLocation, surface systems should be added to the
        patches, not through the membrane.
        :meta private:
        """
        raise NotImplementedError(f'Cannot add space systems to a Membrane.')

    def _SetUpMdlDeps(self, mdl):
        """Set up the structures that will allow specie or reaction access from locations."""
        for p in self.patches:
            self.sysNames += p.sysNames
        self.sysNames = list(set(self.sysNames))

        super()._SetUpMdlDeps(mdl)

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsMemb]

    def _createStepsObj(self, geom, verify, opt_method, search_percent, opt_file_name):
        (geomObj,) = geom._getStepsObjects()
        patches = [p._getStepsObjects()[0] for p in self.patches]
        return sgeom.Memb(
            self.name,
            geomObj,
            patches,
            verify=verify,
            opt_method=opt_method,
            search_percent=search_percent,
            opt_file_name=opt_file_name,
        )

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
        return numpy.asarray([x, y, z]).view(cls)

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


class TetMesh(Geometry):
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
        self.stepsMesh = None
        self.meshInfo = {
            'vertBlocks': {},
            'tetBlocks': {},
            'triBlocks': {},
            'vertGroups': {},
            'tetGroups': {},
            'triGroups': {},
        }
        if _createObj:
            raise Exception('Cannot create a bare TetMesh, use one of the class methods.')
        super().__init__(*args, _createObj=False, **kwargs)

        self._verts = None
        self._bars = None
        self._tris = None
        self._tets = None
        self._surface = None

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsMesh]

    @classmethod
    def _FromStepsObject(cls, obj, comps=None, patches=None):
        """Create the interface object from a STEPS object."""
        mesh = cls(_createObj=False)
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
    def FromData(cls, verts, tets, tris=[]):
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

        :returns: The constructed TetMesh
        :rtype: :py:class:`TetMesh`
        """
        mesh = cls(_createObj=False)
        obj = sgeom.Tetmesh(verts, tets, tris)
        mesh.stepsMesh = obj
        mesh.stepsGeom = obj
        return mesh

    @classmethod
    def Load(cls, path, scale=1, strict=False):
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

        :returns: The loaded TetMesh
        :rtype: :py:class:`TetMesh`
        """
        smesh, scomps, spatches = smeshio.loadMesh(path, scale, strict)
        return TetMesh._FromStepsObject(smesh, comps=scomps, patches=spatches)

    def _loadElementProxys(self, ndprx, tetprx, triprx):
        """Load the blocks and groups from ElementProxy objects."""
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
    def LoadAbaqus(cls, filename, scale=1, ebs=None, shadow_mesh=None):
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

        mesh = TetMesh._FromStepsObject(stepsMesh)
        mesh._loadElementProxys(ndprx, tetprx, triprx)
        return mesh

    @classmethod
    def LoadGmsh(cls, filename, scale=1):
        """Load a mesh from a Gmsh (2.2 ASCII)-formated mesh file

        If blocks or groups of elements are present in the file, they will also be loaded.

        :param filename: The Gmsh filename (or path) including any suffix.
        :type filename: str
        :param scale: Optionally rescale the mesh on loading by given factor
        :type scale: float

        :returns: The loaded TetMesh
        :rtype: :py:class:`TetMesh`
        """
        stepsMesh, ndprx, tetprx, triprx = smeshio.importGmsh(filename, scale)
        mesh = TetMesh._FromStepsObject(stepsMesh)
        mesh._loadElementProxys(ndprx, tetprx, triprx)
        return mesh

    @classmethod
    def LoadVTK(cls, filename, scale=1):
        """Load a mesh from a VTK (legacy ASCII)-formated mesh file

        If blocks or groups of elements are present in the file, they will also be loaded.

        :param filename: The VTK filename (or path) including any suffix.
        :type filename: str
        :param scale: Optionally rescale the mesh on loading by given factor
        :type scale: float

        :returns: The loaded TetMesh
        :rtype: :py:class:`TetMesh`
        """
        stepsMesh, ndprx, tetprx, triprx = smeshio.importVTK(filename, scale)
        mesh = TetMesh._FromStepsObject(stepsMesh)
        mesh._loadElementProxys(ndprx, tetprx, triprx)
        return mesh

    @classmethod
    def LoadTetGen(cls, pathroot, scale):
        """Load a mesh from a TetGen-formated set of files

        If blocks or groups of elements are present, they will also be loaded.

        :param pathroot: The root of the path to the mesh files. E.g. mesh/torus should be given
            to read files mesh/torus.node, mesh/torus.ele, and optionally mesh/torus.face
        :type pathroot: str
        :param scale: Optionally rescale the mesh on loading by given factor
        :type scale: float

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
        mesh = TetMesh._FromStepsObject(stepsMesh)
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
    def tets(self):
        """All tetrahedrons in the mesh

        :type: :py:class:`TetList`, read-only
        """
        if self._tets is None:
            self._tets = TetList(range(self.stepsMesh.countTets()), self, _immutable=True)
        return self._tets

    @property
    def tris(self):
        """All triangles in the mesh

        :type: :py:class:`TriList`, read-only
        """
        if self._tris is None:
            self._tris = TriList(range(self.stepsMesh.countTris()), self, _immutable=True)
        return self._tris

    @property
    def surface(self):
        """All triangles in the surface of the mesh

        :type: :py:class:`TriList`, read-only
        """
        if self._surface is None:
            self._surface = TriList(self.stepsMesh.getSurfTris(), self, _immutable=True)
        return self._surface

    @property
    def verts(self):
        """All vertices in the mesh

        :type: :py:class:`VertList`, read-only
        """
        if self._verts is None:
            self._verts = VertList(range(self.stepsMesh.countVertices()), self, _immutable=True)
        return self._verts

    @property
    def bars(self):
        """All bars in the mesh

        :type: :py:class:`BarList`, read-only
        """
        if self._bars is None:
            self._bars = BarList(range(self.stepsMesh.countBars()), self, _immutable=True)
        return self._bars

    @property
    def bbox(self):
        """The bounding box of the mesh

        :type: :py:class:`BoundingBox`, read-only
        """
        return BoundingBox(Point(*self.stepsMesh.getBoundMin()), Point(*self.stepsMesh.getBoundMax()))

    @property
    def Vol(self):
        """Total volume of the mesh

        :type: float, read-only
        """
        return self.stepsMesh.getMeshVolume()

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


class Reference(nutils.UsingObjects(nutils.Optional(TetMesh)), nutils.SolverPathObject):
    """Base class for all element references.

    :param idx: Index of the element
    :type idx: int
    :param mesh: Mesh object that contains the element
    :type mesh: :py:class:`TetMesh` or None

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).

    The actual index of the element can be accessed with the idx attribute.
    """

    def __init__(self, idx, mesh=None, anonymous=False, *args, **kwargs):
        if not anonymous:
            # Only call parent __init__ if we need the reference to be named
            super().__init__(*args, addAsElement=False, **kwargs)
        if mesh is None:
            (mesh,) = self._getUsedObjects(addAsElement=False)
        if not isinstance(mesh, TetMesh):
            raise TypeError(f'Expected a TetMesh object, got {mesh} instead.')
        self.mesh = mesh
        self.idx = self._getActualIdx(idx, mesh)
        """Index of the element

        :Type: int, read-only
        """

    def toList(self):
        """Get a list holding only this element

        :returns: A list with only this element
        :rtype: :py:class:`RefList`
        """
        return self.__class__._lstCls([self], mesh=self.mesh)

    @classmethod
    def _getActualIdx(cls, idx, mesh):
        """Return the integer idx."""
        if isinstance(idx, Reference) and idx.__class__ == cls:
            if mesh != idx.mesh:
                raise Exception(
                    f'Cannot create a reference from a reference related to a different mesh.'
                )
            idx = idx.idx
        if not isinstance(idx, numbers.Integral):
            raise TypeError(f'Expected an integer index, got {idx} instead.')
        return idx

    def _getPhysicalLocation(self):
        """Return the location associated with the reference (comp for tet, patch for tri)."""
        pass

    def __hash__(self):
        return hash((id(self.mesh), self.idx))

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.mesh == other.mesh and self.idx == other.idx

    def __repr__(self):
        return f'{self.__class__._locStr}({self.idx})'

    def _solverId(self):
        """Return the id that is used to identify an object in the solver."""
        return (self.idx,)

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        mtdt = {'loc_type': self.__class__._locStr, 'loc_id': self.idx}
        # Add information about parent comp or patch
        parent = self._getPhysicalLocation()
        if parent is not None:
            for key, val in parent._simPathAutoMetaData().items():
                mtdt['parent_' + key] = val
        return mtdt


class TetReference(Reference):
    """Convenience class for accessing properties of a tetrahedron

    :param idx: Index of the tetrahedron
    :type idx: int
    :param mesh: Mesh object that contains the tetrahedron
    :type mesh: :py:class:`TetMesh` or None

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).
    """

    _locStr = 'Tet'

    def __init__(self, idx, mesh=None, *args, **kwargs):
        if mesh is None:
            (mesh,) = self._getUsedObjects(addAsElement=False)
        if hasattr(idx, '__iter__') and len(idx) == 3:
            idx = mesh.stepsMesh.findTetByPoint(list(idx))
        super().__init__(idx, mesh=mesh, *args, **kwargs)

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

    @property
    def neighbs(self):
        """Neighboring tetrahedrons (between 0 and 4)

        :type: :py:class:`TetList`, read-only
        """
        inds = self.mesh.stepsMesh.getTetTetNeighb(self.idx)
        return TetList([ind for ind in inds if ind != UNKNOWN_TET], mesh=self.mesh)

    @property
    def faces(self):
        """Faces of the tetrahedron

        :type: :py:class:`TriList`, read-only
        """
        return TriList(self.mesh.stepsMesh.getTetTriNeighb(self.idx), mesh=self.mesh)

    @property
    def center(self):
        """Barycenter of the tetrahedron

        :type: :py:class:`Point`, read-only
        """
        return Point(*self.mesh.stepsMesh.getTetBarycenter(self.idx))

    @property
    def verts(self):
        """Vertices of the tetrahedron

        :type: :py:class:`VertList`, read-only
        """
        return VertList(self.mesh.stepsMesh.getTet(self.idx), mesh=self.mesh)

    @property
    def comp(self):
        """Compartment associated to the tetrahedron (if any)

        :type: :py:class:`Comp` or None, read-only
        """
        c = self.mesh.stepsMesh.getTetComp(self.idx)
        if c is None:
            return None
        else:
            return getattr(self.mesh, c.getID())

    @property
    def Vol(self):
        """Volume of the tetrahedron

        :type: float, read-only
        """
        return self.mesh.stepsMesh.getTetVol(self.idx)

    @property
    def qualityRER(self):
        """Radius-edge-ratio (a quality measurement) of the tetrahedron

        :type: float, read-only
        """
        return self.mesh.stepsMesh.getTetQualityRER(self.idx)

    def _getPhysicalLocation(self):
        """Return the location associated with the reference (comp for tet, patch for tri)."""
        return self.comp

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return TetReference._locStr


class TriReference(Reference):
    """Convenience class for accessing properties of a triangle

    :param idx: Index of the triangle
    :type idx: int
    :param mesh: Mesh object that contains the triangle
    :type mesh: :py:class:`TetMesh` or None

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
        return VertList(self.mesh.stepsMesh.getTri(self.idx), mesh=self.mesh)

    @property
    def center(self):
        """Barycenter of the triangle

        :type: :py:class:`Point`, read-only
        """
        return Point(*self.mesh.stepsMesh.getTriBarycenter(self.idx))

    @property
    def Area(self):
        """Area of the triangle

        :type: float, read-only
        """
        return self.mesh.stepsMesh.getTriArea(self.idx)

    @property
    def tetNeighbs(self):
        """Neighboring tetrahedrons

        :type: :py:class:`TetList`, read-only
        """
        inds = self.mesh.stepsMesh.getTriTetNeighb(self.idx)
        return TetList([ind for ind in inds if ind != UNKNOWN_TRI], mesh=self.mesh)

    @property
    def patch(self):
        """Patch associated to the triangle (if any)

        :type: :py:class:`Patch` or None, read-only
        """
        p = self.mesh.stepsMesh.getTriPatch(self.idx)
        if p is None:
            return None
        else:
            return getattr(self.mesh, p.getID())

    @property
    def bars(self):
        """bars of the triangle

        :type: :py:class:`VertList`, read-only
        """
        return BarList(self.mesh.stepsMesh.getTriBars(self.idx), mesh=self.mesh)

    @property
    def norm(self):
        """Normal vector of the triangle

        :type: :py:class:`Point`, read-only
        """
        return Point(*self.mesh.stepsMesh.getTriNorm(self.idx))

    def _getPhysicalLocation(self):
        """Return the location associated with the reference (comp for tet, patch for tri)."""
        return self.patch

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return TriReference._locStr


class BarReference(Reference):
    """Convenience class for accessing properties of a bar

    :param idx: Index of the bar
    :type idx: int
    :param mesh: Mesh object that contains the bar
    :type mesh: :py:class:`TetMesh` or None

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
        return VertList(self.mesh.stepsMesh.getBar(self.idx), mesh=self.mesh)

    @property
    def center(self):
        v1, v2 = self.verts
        return Point(*((v1 + v2) / 2))


class VertReference(Reference, Point):
    """Convenience class for accessing properties of a vertex

    Can be used in the same way as a :py:class:`Point`.

    :param idx: Index of the vertex
    :type idx: int
    :param mesh: Mesh object that contains the vertex
    :type mesh: :py:class:`TetMesh`

    .. note::
        In contrast to other references, VertReference can only be created by explicitely
        providing a mesh object, it cannot use the context manager syntax described in
        :py:class:`TetMesh`.
    """

    _locStr = 'Vert'

    def __new__(cls, idx, mesh, anonymous=False):
        if not isinstance(mesh, TetMesh):
            raise TypeError(f'Expected a TetMesh object, got {mesh} instead.')
        idx = cls._getActualIdx(idx, mesh)
        return Point.__new__(cls, *mesh.stepsMesh.getVertex(idx))

    def __init__(self, idx, mesh=None, *args, **kwargs):
        super().__init__(idx, mesh=mesh, *args, **kwargs)

    # Need to redefine __hash__ and __eq__ to be sure the Reference version of these methods
    # is called.
    def __hash__(self):
        return Reference.__hash__(self)

    def __eq__(self, other):
        return Reference.__eq__(self, other)

    def __ne__(self, other):
        return not Reference.__eq__(self, other)

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return VertReference._locStr


class RefList(nutils.UsingObjects(nutils.Optional(TetMesh)), nutils.SolverPathObject):
    """Base convenience class for geometrical element lists

    :param lst: The list of elements
    :type lst: Iterable[int] or range or Iterable[:py:class:`Reference`] or :py:class:`RefList`
    :param mesh: Mesh object that contains the elements
    :type mesh: :py:class:`TetMesh` or None

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

    def __init__(self, lst, mesh=None, _immutable=False, *args, **kwargs):
        super().__init__(*args, addAsElement=False, **kwargs)
        if mesh is None:
            (mesh,) = self._getUsedObjects(addAsElement=False)
        elif not isinstance(mesh, TetMesh):
            raise TypeError(f'Expected a TetMesh object, got a {type(mesh)} instead.')
        self.mesh = mesh

        if lst.__class__ == self.__class__:
            self.lst = copy.copy(lst.lst)
            self.mesh = lst.mesh
        elif isinstance(lst, range):
            self.lst = lst
        elif hasattr(lst, '__iter__'):
            lst = list(lst)
            self.lst = self._getIdxLst(lst)
            if len(lst) > 0 and isinstance(lst[0], self.__class__._refCls):
                self.mesh = lst[0].mesh
        else:
            raise TypeError(f'Cannot create an element list from a {type(lst)}')
        self._immutable = _immutable
        # Optimization data wrapper, used as a context manager when modifying the list
        self._optimCM = RefList._OptimizationCM(self)

        if self.mesh is None:
            raise Exception(f'Unable to identify which mesh should the list be bound to.')

    @classmethod
    def _getIdxLst(cls, lst):
        lst = list(lst)
        if all(isinstance(e, numbers.Integral) for e in lst):
            return lst
        if all(isinstance(e, cls._refCls) for e in lst):
            return [e.idx for e in lst]
        else:
            raise TypeError(f'Cannot use {lst} to initialize a list of {cls._refCls}.')

    @classmethod
    def _toRefList(cls, lst, mesh):
        if isinstance(lst, RefList):
            if lst.__class__ != cls:
                raise TypeError(f'Expected a {cls}, got a {lst.__class__} instead.')
            return lst
        else:
            return cls(lst, mesh)

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
            return self.__class__(idxres, mesh=self.mesh)
        else:
            return self.__class__._refCls(idxres, mesh=self.mesh)

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
            yield self.__class__._refCls(idx, mesh=self.mesh, anonymous=True)

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
            self.lst.append(e.idx)

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
        if isinstance(elem, Reference):
            if not isinstance(elem, self.__class__._refCls):
                return False
            # Optimization
            if RefList._ElemsAsSet_K not in self._optimCM:
                self._optimCM[RefList._ElemsAsSet_K] = set(self.lst)
            return elem.idx in self._optimCM[RefList._ElemsAsSet_K]
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
            return self.lst.index(elem.idx)
        else:
            raise TypeError(f'Expected a {self.__class__._refCls}, got {elem} instead.')

    def _checkSameType(self, other):
        if not isinstance(other, self.__class__):
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
        return self.__class__(res, mesh=self.mesh)

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
        return self.__class__(res, mesh=self.mesh)

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
        return self.__class__(res, mesh=self.mesh)

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
        return self.__class__(res, mesh=self.mesh)

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
        return self.__class__(self._getLst() + other._getLst(), mesh=self.mesh)

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


class TetList(RefList):
    """Convenience class for tetrahedron lists

    :param lst: The list of tetrahedrons
    :type lst: Iterable[int] or range or Iterable[:py:class:`TetReference`] or :py:class:`TetList`
    :param mesh: Mesh object that contains the elements
    :type mesh: :py:class:`TetMesh` or None

    If `lst` is omitted, the list will be initialized to an empty list.

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).

    Behaves like a :py:class:`RefList` but with additional properties specific to tetrahedrons.
    Note that :py:func:`__getitem__` and :py:func:`__contains__` are modified to accept a 3D
    point as parameter.
    """

    _refCls = TetReference

    def __init__(self, lst=None, mesh=None, _immutable=False, *args, **kwargs):
        super().__init__(lst if lst is not None else [], mesh, _immutable, *args, **kwargs)

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
                idx = self.mesh.stepsMesh.findTetByPoint(list(key))
                if idx != UNKNOWN_TET and idx in self.lst:
                    return TetReference(idx, mesh=self.mesh)
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
                [tet for tet in self if any(neighb not in self for neighb in tet.neighbs)], mesh=self.mesh
            )
        shell = self._optimCM[TetList._CurrShell_K]

        for i in range(d):
            newShell = set()
            for tet in shell:
                for neighb in tet.neighbs:
                    if neighb.idx not in prevShell:
                        newShell.add(neighb.idx)
            prevShell.update(newShell)
            newShell = list(newShell)
            self.lst = self._getLst() + newShell
            shell = TetList(newShell, mesh=self.mesh)

        with self._modify():
            self._optimCM[TetList._CurrShell_K] = shell
            self.lst = sorted(set(self.lst))

    def erode(self, d=1):
        """Erodes the list, removing surface tetrahedrons

        One cycle of erosion corresponds to removing all tetrahedrons that are on the surface of
        the list (i.e. they have at least one neighbor tetrahedron that is not in the list). This
        operation is repeated ``d`` times.

        :param d: Topological distance to erode, defaults to 1
        :type d: int
        """
        shell = TetList(
            [
                tet
                for tet in self
                if any(neighb not in self for neighb in tet.neighbs) or len(tet.neighbs) < 4
            ],
            mesh=self.mesh,
        )
        shell.dilate(d - 1)
        with self._modify():
            self.lst = [tet.idx for tet in self if tet not in shell]

    _CurrSurface_K = 'CurrSurface'
    _LastSurface_K = 'LastSurface'

    @property
    def surface(self):
        """Surface of the tetrahedrons in the list

        I.e. all triangles that are not shared by two tetrahedrons in the list.

        :type: :py:class:`TriList`, read-only
        """
        if TetList._CurrSurface_K not in self._optimCM:
            if hasattr(self, TetList._LastSurface_K):
                # Only compute the changes due to the addition or removal of tetrahedrons
                oldLst, oldSurf = getattr(self, TetList._LastSurface_K)
                surf = set(oldSurf)
                for tet in TetList(oldLst, mesh=self.mesh) ^ self:
                    surf ^= set(tet.faces)
            else:
                # Compute the surface from scratch
                surf = set()
                for tet in set(self):
                    surf ^= set(tet.faces)
            # Keep the surface until next list change
            self._optimCM[TetList._CurrSurface_K] = TriList(sorted(surf, key=lambda x: x.idx), mesh=self.mesh)
            # Keep a version of the current list and surface, for diff computation later
            setattr(
                self, TetList._LastSurface_K, (copy.copy(self.lst), self._optimCM[TetList._CurrSurface_K])
            )
        return self._optimCM[TetList._CurrSurface_K]

    @property
    def Vol(self):
        """The summed volume of all tetrahedrons in the list

        .. note::
            If the list contains duplicate elements, they will be counted only once in the
            total volume.

        :type: float
        """
        return sum(tet.Vol for tet in set(self))

    @property
    def tris(self):
        """All faces of all tetrahedrons in the list

        :type: :py:class:`TriList`, read-only
        """
        tris = TriList([], mesh=self.mesh)
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
        allComps = set()
        for elem in self:
            c = elem.comp
            if c is not None:
                allComps.add(c)
        allChildren = {}
        for c in allComps:
            allChildren.update(c.children)
        return allChildren


class TriList(RefList):
    """Convenience class for triangle lists

    :param lst: The list of triangles
    :type lst: Iterable[int] or range or Iterable[:py:class:`TriReference`] or :py:class:`TriList`
    :param mesh: Mesh object that contains the elements
    :type mesh: :py:class:`TetMesh` or None

    If `lst` is omitted, the list will be initialized to an empty list.

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).

    Behaves like a :py:class:`RefList` but with additional properties specific to triangles.
    """

    _refCls = TriReference

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
        return BarList(sorted(edg, key=lambda x: x.idx), mesh=self.mesh)

    @property
    def bars(self):
        """All bars of all triangles in the list

        :type: :py:class:`BarList`, read-only
        """
        bars = BarList([], mesh=self.mesh)
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

        :type: float
        """
        return sum(tri.Area for tri in set(self))

    @nutils.NamedObject.children.getter
    def children(self):
        allPatches = set()
        for elem in self:
            p = elem.patch
            if p is not None:
                allPatches.add(p)
        allChildren = {}
        for p in allPatches:
            allChildren.update(p.children)
        return allChildren


class BarList(RefList):
    """Convenience class for bar lists

    :param lst: The list of bars
    :type lst: Iterable[int] or range or Iterable[:py:class:`BarReference`] or :py:class:BarList`
    :param mesh: Mesh object that contains the elements
    :type mesh: :py:class:`TetMesh` or None

    If `lst` is omitted, the list will be initialized to an empty list.

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).

    Behaves like a :py:class:`RefList` but with additional properties specific to bars.
    """

    _refCls = BarReference

    def __init__(self, lst=None, mesh=None, _immutable=False, *args, **kwargs):
        super().__init__(lst if lst is not None else [], mesh, _immutable, *args, **kwargs)

    @property
    def verts(self):
        """All vertices of all bars in the list

        :type: :py:class:`VertList`, read-only
        """
        verts = VertList([], mesh=self.mesh)
        for bar in self:
            verts |= bar.verts
        return verts


class VertList(RefList):
    """Convenience class for vertex lists

    :param lst: The list of vertices
    :type lst: Iterable[int] or range or Iterable[:py:class:`VertReference`] or :py:class:VertList`
    :param mesh: Mesh object that contains the elements
    :type mesh: :py:class:`TetMesh` or None

    If `lst` is omitted, the list will be initialized to an empty list.

    If mesh is not provided, it will be inferred, if possible, from context managers
    (see :py:class:`TetMesh`).

    Behaves like a :py:class:`RefList`.
    """

    _refCls = VertReference

    def __init__(self, lst=None, mesh=None, _immutable=False, *args, **kwargs):
        super().__init__(lst if lst is not None else [], mesh, _immutable, *args, **kwargs)


# Add a link from the reference classes to the list class
for cls in RefList.__subclasses__():
    cls._refCls._lstCls = cls
# Do not leave cls as a global variable
del cls


class ROI(nutils.UsingObjects(TetMesh), nutils.StepsWrapperObject):
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
            sgeom.ELEM_TET: TetReference,
            sgeom.ELEM_TRI: TriReference,
            sgeom.ELEM_VERTEX: VertReference,
        }
        region = cls(_createObj=False, name=name)
        region.stepsROI = obj
        region._elemType = revTypeMap[obj.type]
        region.lst = region._elemType._lstCls._toRefList(obj.indices, mesh)
        return region

    def _createStepsObj(self, mesh, lst):
        typeMap = {
            TetReference: sgeom.ELEM_TET,
            TriReference: sgeom.ELEM_TRI,
            VertReference: sgeom.ELEM_VERTEX,
        }
        if len(lst) == 0:
            raise Exception('An empty list was given.')
        elem = next(iter(lst))
        if elem.__class__ not in typeMap:
            raise TypeError('Unsupported element type for creating an ROI.')
        tpe = typeMap[elem.__class__]
        mesh.stepsMesh.addROI(self.name, tpe, [elem.idx for elem in lst])
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


class DiffBoundary(nutils.UsingObjects(TetMesh), nutils.StepsWrapperObject):
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
        return sgeom.DiffBoundary(self.name, mesh.stepsMesh, [tri.idx for tri in self.lst])

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
        return sgeom.SDiffBoundary(
            self.name, mesh.stepsMesh, [bar.idx for bar in self.lst], [p.stepsPatch for p in self.patches]
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
