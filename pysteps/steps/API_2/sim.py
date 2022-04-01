####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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

import atexit
import copy
import enum
import heapq
import importlib
import math
import numbers
import numpy
import os
import re
import sys
import warnings

from steps import stepslib

from . import model as nmodel
from . import geom as ngeom
from . import utils as nutils
from . import saving as nsaving
from . import _saving_optim as nsaving_optim
from . import simcheck as nsimcheck


__all__ = [
    'Simulation',
    'SBMLSimulation',
    'SimPath',
    'MPI',
]


###################################################################################################
# Enums


if stepslib._STEPS_USE_DIST_MESH:
    class KSPNorm(enum.IntEnum):
        DEFAULT = stepslib._py_KSPNormType.KSP_NORM_DEFAULT
        NONE = stepslib._py_KSPNormType.KSP_NORM_NONE
        PRECONDITIONED = stepslib._py_KSPNormType.KSP_NORM_PRECONDITIONED
        UNPRECONDITIONED = stepslib._py_KSPNormType.KSP_NORM_UNPRECONDITIONED
        NATURAL = stepslib._py_KSPNormType.KSP_NORM_NATURAL


    class SSAMethod(enum.IntEnum):
        SSA = stepslib._py_SSAMethod.SSA
        RSSA = stepslib._py_SSAMethod.RSSA


    class NextEventSearchMethod(enum.IntEnum):
        DIRECT = stepslib._py_SearchMethod.DIRECT
        GIBSON_BRUCK = stepslib._py_SearchMethod.GIBSON_BRUCK


    __all__ +=  ['KSPNorm', 'SSAMethod', 'NextEventSearchMethod']


###################################################################################################
# Exceptions


class SimPathExtensionError(Exception):
    """
    :meta private:
    """

    def __init__(self, elem, name=None):
        self.elem = elem
        self.name = name


class SimPathSolverMissingMethod(Exception):
    """
    :meta private:
    """

    pass


class SimPathInvalidPath(AttributeError):
    """
    :meta private:
    """

    pass


class SolverCallError(Exception):
    """
    :meta private:
    """

    pass


###################################################################################################
# Simulation


class _SimPathDescr:
    """
    Descriptor class for setting and getting values from the steps solver.
    """

    def __init__(self, name):
        self.name = name

    def __get__(self, inst, objtype):
        if inst is None:
            return self

        res = [self._getPathValue(path) for path in inst]

        return res[0] if len(res) == 1 else res

    def __set__(self, inst, val):
        paths = list(path for path in inst)
        if len(paths) == 0:
            raise SimPathInvalidPath(f'Empty path.')

        if hasattr(val, '__iter__'):
            if len(paths) == 1:
                # The value to be assigned is multidimensional
                self._setPathValue(paths[0], val)
            elif len(val) == len(paths):
                # Assign potentially different values to each path
                for path, v in zip(paths, val):
                    self._setPathValue(path, v)
            else:
                raise SimPathInvalidPath(
                    f'Path {inst} covers {len(paths)} elements while the '
                    f'values that it should be assigned contain {len(val)} '
                    f'elements.'
                )
        else:
            # Assign the same value to all paths
            for path in paths:
                self._setPathValue(path, val)

    def _getFinalPathsFunction(self, finalPaths, prefix='get'):
        """
        Return a function that, when called with no arguments, generates the values associated
        to the paths in finalPaths.
        """

        def noArgWrap(simPathComb, genFunc):
            def internal():
                return simPathComb.func([func(*args, **kwargs) for func, args, kwargs in genFunc])

            return internal

        subFuncs = []
        for path in finalPaths:
            if isinstance(path, nutils.SimPathCombiner):
                fpf = self._getFinalPathsFunction(path)
                subFuncs.append((noArgWrap(path, fpf), [], {}))
            else:
                subFuncs.append(self._getFunction(prefix, path))

        return subFuncs

    def _getPathValue(self, path):
        """Return the value corresponding to "path"."""
        if isinstance(path, nutils.SimPathCombiner):
            return path.func([self._getPathValue(p) for p in path])
        else:
            fun, params, kwparams = self._getFunction('get', path)
            try:
                return fun(*params, **kwparams)
            except Exception as e:
                raise SolverCallError(f'Exception raised during call to solver: {e}')

    def _setPathValue(self, path, val):
        """
        Set path "path" to value "val".
        "val" can be a Param object when several arguments are needed.
        """
        if isinstance(path, nutils.SimPathCombiner):
            raise Exception(f'Cannot set simulation paths that involve value combining.')

        # Parameter checks and potential conversions
        if isinstance(val, nutils.Parameter):
            param = val
        else:
            param = nutils.Parameter(val, name='')

        if self.name in SimPath._PATH_END_UNITS:
            units = SimPath._PATH_END_UNITS[self.name]
            if hasattr(units, '__call__'):
                units = units(path)

            if param._units is not None:
                # Check units and convert value to appropriate unit
                val = param.valueIn(units)
            else:
                # Set the correct unit to the Parameter and get value
                param._units = units
                val = param.value
        else:
            if param._units is not None:
                raise Exception(
                    f'Parameter {param} with units {param.units} was used to set a value with '
                    f'no units.'
                )
            else:
                # "Default" behavior
                val = param.value
        # Keep track of values that were set with a Parameter
        path[0]._parameterSimPathCallback(path + (self.name, ), param)

        fun, params, kwparams = self._getFunction('set', path)
        if isinstance(val, nutils.Params):
            args = [a.name if isinstance(a, nutils.NamedObject) else a for a in val.args]
            kwargs = {k: a.name if isinstance(a, nutils.NamedObject) else a for k, a in val.kwargs.items()}
        else:
            args = [path[-1]._solverSetValue(self.name, val)]
            kwargs = {}
        try:
            fun(*params, *args, **kwparams, **kwargs)
        except Exception as e:
            raise SolverCallError(f'Exception raised during call to solver: {e}')

    def _getFunction(self, prefix, path):
        """
        Return a pair of function and parameters necessary to get or set values in the solver
        for a specific path.
        """
        def wrap(acc, f):
            return lambda x: acc(f(x))

        sim, path = path[0], path[1:]

        funcName = prefix
        isGetter = prefix == 'get'
        params = tuple()
        kwparams = {}
        modifier = None
        for e in path:
            funcName += e._solverStr()
            params += e._solverId()
            kwparams.update(e._solverKeywordParams())

            if isGetter:
                m = e._solverModifier()
                if m is not None:
                    if modifier is None:
                        modifier = m
                    else:
                        modifier = wrap(modifier, m)

        try:
            func = getattr(sim.stepsSolver, funcName + self.name)
        except AttributeError:
            raise SimPathSolverMissingMethod(
                f'Method {funcName+self.name} is not available for solver {sim.stepsSolver}.'
            )
        if modifier is not None:
            return (lambda *x, **kw: modifier(func(*x, **kw))), params, kwparams
        else:
            return func, params, kwparams


class SimPath:
    """Describes the access to data from the solver as a path

    The root of the path is a simulation and the parts of the path are separated by dots:
    ``sim.cyt.Ca.Count`` returns the count of species ``'Ca'`` in compartment ``'cyt'`` for
    simulation ``'sim'``.

    The path usually follows the convention ``root.location.object.attribute``.
    Objects along the path can be made more specific with the square bracket syntax:
    ``sim.ERmemb.sReac1["fwd"].K`` returns the reaction constant associated to the forward part of
    reaction ``'sReac1'`` on patch ``'ERmemb'``.

    This syntax allows to either get or set the data:
    ``sim.cyt.Ca.Count = 5000`` sets the number of ``'Ca'`` in compartment ``'cyt'`` to ``5000``.

    Data can also be accessed in a grouped way:
    ``sim.ALL().Ca.Count`` returns a list of the count of species ``'Ca'`` in all locations.
    Special methods like :py:func:`ALL`, :py:func:`LIST` and :py:func:`MATCH`
    allow to select path parts based on their type or names. Special methods like :py:func:`TET`,
    :py:func:`TRI`, and :py:func:`VERT` allow the specification of single geometrical elements as
    location (e.g. ``sim.TET(tet1).Ca.Count``) while the plural version (:py:func:`TETS`, ...)
    represents grouped access to lists of geometrical elements.

    Note that if the :py:class:`SimPath` represents a list of n values, it can be set either with a
    single value or with a list of the same size n. Setting with a single value is identical to
    repeating this single value n times.

    The last part of the path needs to be a recognized attribute for the path to actually return
    a value. Recognized attributes are: , :py:attr:`Vol`, :py:attr:`Area`, :py:attr:`K`,
    :py:attr:`Active`, :py:attr:`C`, :py:attr:`H`, :py:attr:`A`, :py:attr:`D`, :py:attr:`Extent`,
    :py:attr:`Count`, :py:attr:`Conc`, :py:attr:`Amount`, :py:attr:`Clamped`, :py:attr:`Dcst`,
    :py:attr:`DiffusionActive`, :py:attr:`Potential`, :py:attr:`Capac`, :py:attr:`Res`,
    :py:attr:`VolRes`, :py:attr:`I`, :py:attr:`IClamp`, :py:attr:`V`, :py:attr:`VClamped`.

    .. note::
        Should not be directly instantiated by a user, the class is only documented for clarity.
        A :py:class:`SimPath` should only be created through a :py:class:`Simulation`.
    """

    # Internally, the paths are represented by a tree whose root is the simulation and tips are the
    # elements currently accessed.

    _endNames = {
        'Vol': """Volume, location needs to be a :py:class:`Compartment` or tetrahedron(s)""",
        'Area': """Area, location needs to be a :py:class:`Patch` or triangle(s)""",
        'K': """Reaction constant, object needs to be a :py:class:`Reaction`""",
        'Active': """Whether a reaction or diffusion rule is active, object needs to be a
            :py:class:`Reaction` or :py:class:`Diffusion`""",
        'C': """Stochastic reaction constant, object needs to be a :py:class:`Reaction`""",
        'H': """The distinct number of ways a reaction can occur h_mu, object  needs to be a
            :py:class:`Reaction`""",
        'A': """Propensity of reaction, object needs to be a :py:class:`Reaction`""",
        'D': """Diffusion constant, object needs to be a :py:class:`Diffusion`""",
        'Extent': """Number of times the reaction or diffusion rule has been executed, object
            needs to be a :py:class:`Reaction` or :py:class:`Diffusion`""",
        'Count': """Number of objects""",
        'Conc': """Concentration of objects, location needs to be a :py:class:`Compartment` or
            tetrahedron(s)""",
        'Amount': """Amount (in mol) of objects""",
        'Clamped': """Whether the number of objects is clamped""",
        'Dcst': """Diffusion constant across a diffusion boundary, object needs to be a
            :py:class:`DiffBoundary`, an additional object must be specified (e.g. a Species)""",
        'DiffusionActive': """Whether diffusion across a diffusion boundary is active object
            needs to be a :py:class:`DiffBoundary`""",
        'Potential': """Potential of a membrane, location needs to be a :py:class:`Membrane`""",
        'Capac': """Capacitance of a membrane or of (a) triangle(s), location needs to be a
            :py:class:`Membrane` or (a) triangle(s)""",
        'Res': """Electrical resistivity of a membrane, location needs to be a
            :py:class:`Membrane`""",
        'VolRes': """Bulk electrical resistivity (in ohm.m) of the volume associated with a
            membrane, location needs to be a :py:class:`Membrane`""",
        'I': """Current through the location, object needs to be a :py:class:`Current`""",
        'IClamp': """Current clamp (in A), location needs to be a vertex or a triangle, object
            needs to be :py:class:`Current`""",
        'V': """Potential of the location, location can be (a) tetrahedron(s), (a) triangle(s)
            or (a) vert(ex/ices)""",
        'VClamped': """Whether the potential is clamped at the location (same possible location
            as for ``V``)""",
    }

    # Declare units associated to some path endings
    _PATH_END_UNITS = {
        'Vol': nutils.Units('m^3'),
        'Area': nutils.Units('m^2'),
        'K': lambda path: path[2]._getRateUnits(),
        'D': nutils.Units('m^2 s^-1'),
        'Count': nutils.Units(''),
        'Conc': nutils.Units('M'),
        'Amount': nutils.Units('mol'),
        'Dcst': nutils.Units('m^2 s^-1'),
        'Potential': nutils.Units('V'),
        # TODO upon modifications: If methods for setting capacitance of a triangle are implemented,
        # check that their units are also F m^-2
        'Capac': nutils.Units('F m^-2'), 
        # 'Res': No units for Res because it corresponds to a call with resistance and reversal
        #        potential
        'VolRes': nutils.Units('ohm m'),
        'I': nutils.Units('A'),
        'IClamp': nutils.Units('A'),
        'V': nutils.Units('V'),
    }

    def __init__(self, _sim, _elems=None):
        super().__setattr__('_sim', _sim)
        super().__setattr__('_isOnlyRoot', _elems is None)
        # the tree of elements
        super().__setattr__('_elems', {_sim: {}} if _elems is None else _elems)

    def ALL(self, *cls):
        """Extend the path by selecting children

        This special method can be used in place of a location or object name when building a path.
        It will add to the path all children of the current path that are of the given classes.
        If no classes are given, it adds all possible children of the path to the path. Works in
        the same way as :py:func:`steps.API_2.utils.NamedObject.ALL`.

        :param cls: Zero or more STEPS classes
        :type cls: class

        :returns: An updated path taking into account the newly added elements.
        :rtype: :py:class:`SimPath`

        Examples::

            sim.ALL(Compartment).S1.Count        # The number of S1 Species in all compartments, it
                                                 # returns a list of length equal to the number of
                                                 # compartments in the simulation.
            sim.ALL(Patch).ALL(Species).Count    # The number of each species in each compartment,
                                                 # it returns a list with length equal to the sum
                                                 # of the number of different Species defined in
                                                 # each compartment.
            sim.comp1.ALL(Reaction).Active       # Whether each reaction in comp1 is active, it
                                                 # returns a list of booleans of length equal to
                                                 # the total number of reactions in comp1 (forward
                                                 # and backward).
            sim.ALL(Compartment, Patch).S1.Count # The number of S1 Species in all compartments and
                                                 # patches, it returns a list of length equal to the
                                                 # number of compartment and patches.
        """
        if len(cls) == 0:
            cls = [object]
        try:
            elems = self._condExtendElems(self._elems, lambda x: any(isinstance(x, c) for c in cls))
        except SimPathExtensionError as ex:
            raise SimPathInvalidPath(
                f'Element {ex.elem} has no children of type {cls}, it is not possible to extend the path.'
            )
        if len(elems) == 0:
            raise SimPathInvalidPath(
                f'None of the elements have children of type {cls}, '
                f'it is not possible to extend the path.'
            )
        return SimPath(self._sim, elems)

    def LIST(self, *lst):
        """Extend the path by selecting children in a given list

        This special method can be used in place of a location or object name when building a path.
        It will add to the path all children of the current path that are part of the user-supplied
        list. The list can contain strings corresponding to the names of the elements, or the
        elements themselves.

        :param lst: One or several elements
        :type lst: str or :py:class:`steps.API_2.utils.NamedObject`

        :returns: An updated path taking into account the newly added elements.
        :rtype: :py:class:`SimPath`

        Examples::

            sim.LIST(comp1, comp2).S1.Count            # The number of S1 Species in compartments
                                                       # comp1 and comp2 (in this order).
            sim.LIST('comp1', 'comp2').S1.Count        # Same as before, but strings are used
                                                       # instead of the elements themselves.
            sim.LIST(comp1, patch1).LIST(S1, S2).Count # Number of S1 and S2 Species in comp1 and
                                                       # patch1. The resulting order is:
                                                       # [comp1.S1, comp1.S2, patch1.S1, patch1.S2].

        .. note::
            While the order of elements with :py:func:`ALL` depends on the order of addition of the
            corresponding elements to the model, the order of elements in :py:func:`LIST` follows
            the order given by the user in `lst`.
        """
        if len(lst) == 0:
            raise SimPathInvalidPath('Cannot call LIST() without arguments.')

        lst2 = [e.name if isinstance(e, nutils.NamedObject) else e for e in lst]
        try:
            return SimPath(self._sim, self._namesExtendElems(self._elems, lst2))
        except SimPathExtensionError as ex:
            raise SimPathInvalidPath(
                f'Element {ex.elem} does not have a children named {ex.name}, '
                f'it is not possible to extend the path.'
            )

    def MATCH(self, regexp):
        """Extend the path by selecting children whose name match a regular expression

        This special method can be used in place of a location or object name when building a path.
        It will add to the path all children of the current path whose name match the user-supplied
        `regular expression <https://docs.python.org/3/howto/regex.html>`_.

        :param regexp: The regular expression to be matched to children names
        :type regexp: str`

        :returns: An updated path taking into account the newly added elements.
        :rtype: :py:class:`SimPath`

        Examples::

            sim.MATCH('comp[1-9]').S1.Count    # The number of S1 Species in compartments comp1
                                               # through comp9. Note that it only adds the
                                               # compartment that do exist in the model.
            sim.comp1.MATCH('camkII_.*').Count # The count of each element in comp1 whose name
                                               # starts with 'camkII_'.

        .. note::
            An exception will be raised if no children match the pattern.
        """
        expr = re.compile(regexp)
        try:
            elems = self._condExtendElems(self._elems, lambda n: expr.match(n.name) is not None)
        except SimPathExtensionError as ex:
            raise SimPathInvalidPath(
                f'Element {ex.elem} has no children whose name matches '
                f'"{regexp}", it is not possible to extend the path.'
            )
        if len(elems) == 0:
            raise SimPathInvalidPath(
                f'None of the elements have children whose names match'
                f'"{regexp}", it is not possible to extend the path.'
            )
        return SimPath(self._sim, elems)

    def TET(self, *tet):
        """Extend the path by selecting a given tetrahedron

        This special method can only be used in place of a location when building a path.
        It will add the given tetrahedron to the path.

        :param tet: The tetrahedron, or its index in the mesh, or a 3D point contained in a
            tetrahedron
        :type tet: Union[:py:class:`steps.API_2.geom.TetReference`, int, Tuple[float, float, float]]

        :returns: An updated path taking into account the newly added element.
        :rtype: :py:class:`SimPath`

        Example::

            sim.TET(10).S1.Count      # The number of S1 Species in the tetrahedron with index 10
                                      # in the mesh.
            sim.TET(tet1).S1.Count    # The number of S1 Species in the tet1 tetrahedron.
            sim.TET(0, 0, 0).S1.Count # The number of S1 Species in the tetrahedron that contains
                                      # the (x=0, y=0, z=0) point.
        """
        if not self._isOnlyRoot:
            raise SimPathInvalidPath('Tetrahedrons can only be selected at the root of a path.')
        if len(tet) == 1:
            tet = tet[0]
        return SimPath(self._sim, {self._sim: {ngeom.TetReference(tet, self._sim.geom): {}}})

    def TETS(self, tets=None):
        """Extend the path by selecting tetrahedrons in a list

        This special method can only be used in place of a location when building a path.
        It will add the given tetrahedrons to the path.

        :param tets: The tetrahedrons, or anything that can be used to build a
            :py:class:`steps.API_2.geom.TetList`. If no parameters are given, selects all
            tetrahedrons in the mesh.
        :type tets: :py:class:`steps.API_2.geom.TetList`

        :returns: An updated path taking into account the newly added elements.
        :rtype: :py:class:`SimPath`

        Example::

            tetLst = TetList([1, 2, 3], mesh=mesh)

            sim.TETS(tetLst).S1.Count # The number of S1 Species in the tetrahedrons from tetLst
            sim.TETS().S1.Count       # The number of S1 Species in all tetrahedrons in the mesh.
        """
        if not self._isOnlyRoot:
            raise SimPathInvalidPath('Tetrahedrons can only be selected at the root of a path.')
        if tets is None:
            tets = self._sim.geom.tets
        return SimPath(self._sim, {self._sim: {ngeom.TetList(tets, self._sim.geom): {}}})

    def TRI(self, tri):
        """Extend the path by selecting a given triangle

        This special method can only be used in place of a location when building a path.
        It will add the given triangle to the path.

        :param tri: The triangle, or its index in the mesh
        :type tri: Union[:py:class:`steps.API_2.geom.TriReference`, int]

        :returns: An updated path taking into account the newly added element.
        :rtype: :py:class:`SimPath`

        Example::

            sim.TRI(10).S1.Count      # The number of S1 Species in the triangle with index 10
                                      # in the mesh.
            sim.TRI(tri1).S1.Count    # The number of S1 Species in the tri1 triangle.
        """
        if not self._isOnlyRoot:
            raise SimPathInvalidPath('Triangles can only be selected at the root of a path.')
        return SimPath(self._sim, {self._sim: {ngeom.TriReference(tri, self._sim.geom): {}}})

    def TRIS(self, tris=None):
        """Extend the path by selecting triangles in a list

        This special method can only be used in place of a location when building a path.
        It will add the given triangles to the path.

        :param tris: The triangles, or anything that can be used to build a
            :py:class:`steps.API_2.geom.TriList`. If no parameters are given, selects all
            triangles in the mesh.
        :type tris: :py:class:`steps.API_2.geom.TriList`

        :returns: An updated path taking into account the newly added elements.
        :rtype: :py:class:`SimPath`

        Example::

            triLst = TriList([1, 2, 3], mesh=mesh)

            sim.TRIS(triLst).S1.Count # The number of S1 Species in the triangles from triLst
            sim.TRIS().S1.Count       # The number of S1 Species in all triangles in the mesh.
        """
        if not self._isOnlyRoot:
            raise SimPathInvalidPath('Triangles can only be selected at the root of a path.')
        if tris is None:
            tris = self._sim.geom.tris
        return SimPath(self._sim, {self._sim: {ngeom.TriList(tris, self._sim.geom): {}}})

    def VERT(self, vert):
        """Extend the path by selecting a given vertex

        This special method can only be used in place of a location when building a path.
        It will add the given vertex to the path.

        :param vert: The vertex, or its index in the mesh
        :type vert: Union[:py:class:`steps.API_2.geom.VertReference`, int]

        :returns: An updated path taking into account the newly added element.
        :rtype: :py:class:`SimPath`

        Example::

            sim.VERT(10).S1.Count     # The number of S1 Species in the vertex with index 10
                                      # in the mesh.
            sim.VERT(vert1).S1.Count  # The number of S1 Species in the vert1 vertex.
        """
        if not self._isOnlyRoot:
            raise SimPathInvalidPath('Vertices can only be selected at the root of a path.')
        return SimPath(self._sim, {self._sim: {ngeom.VertReference(vert, self._sim.geom): {}}})

    def VERTS(self, verts=None):
        """Extend the path by selecting vertices in a list

        This special method can only be used in place of a location when building a path.
        It will add the given vertices to the path.

        :param verts: The vertices, or anything that can be used to build a
            :py:class:`steps.API_2.geom.VertList`. If no parameters are given, selects all
            vertices in the mesh.
        :type verts: :py:class:`steps.API_2.geom.VertList`

        :returns: An updated path taking into account the newly added elements.
        :rtype: :py:class:`SimPath`

        Example::

            vertLst = VertList([1, 2, 3], mesh=mesh)

            sim.VERTS(vertLst).S1.Count # The number of S1 Species in the vertices from vertLst
            sim.VERTS().S1.Count        # The number of S1 Species in all vertices in the mesh.
        """
        if not self._isOnlyRoot:
            raise SimPathInvalidPath('Vertices can only be selected at the root of a path.')
        if verts is None:
            verts = self._sim.geom.verts
        return SimPath(self._sim, {self._sim: {ngeom.VertList(verts, self._sim.geom): {}}})

    def _condExtendElems(self, elems, pred):
        """
        Extend the element tree 'elems' by adding elements at the tips if they satisfy predicate
        'pred'. Return the new element tree.
        """
        res = {}
        for e, subs in elems.items():
            if len(subs) == 0:
                res[e] = {c: {} for _, c in e.children.items() if pred(c)}
            else:
                res[e] = self._condExtendElems(subs, pred)
            if len(res[e]) == 0:
                del res[e]
        return res

    def _namesExtendElems(self, elems, names):
        """
        Extend the element tree 'elems' by adding elements in "names" at the tips.
        Return the new element tree.
        """
        res = {}
        for e, subs in elems.items():
            if len(subs) == 0:
                res[e] = {}
                for name in names:
                    if name in e.children:
                        res[e][e.children[name]] = {}
                    elif not hasattr(e.__class__, '_SIMPATH_ONLY_CHILDREN') \
                            and name in self._sim.model.children:
                        res[e][self._sim.model.children[name]] = {}
                    else:
                        raise SimPathExtensionError(e, name)
            else:
                res[e] = self._namesExtendElems(subs, names)
        return res

    def _specifyElemsTips(self, elems, key):
        """
        Specifies further the tips of the elems tree with key.
        """
        res = {}
        for e, subs in elems.items():
            if len(subs) == 0:
                if not hasattr(e, '__getitem__'):
                    raise SimPathInvalidPath(
                        f'Element {e} cannot be further specified using square brackets.'
                    )
                newElem = e[key]
                res[newElem] = {}
            else:
                res[e] = self._specifyElemsTips(subs, key)
        return res

    def _callElemsTips(self, elems, params):
        """
        Call the tips of the elems tree with params. This is used to give additional information
        about a part of a path.
        """
        res = {}
        for e, subs in elems.items():
            if len(subs) == 0:
                if not hasattr(e, '__call__'):
                    raise SimPathInvalidPath(
                        f'Element {e} cannot be further specified using parentheses.'
                    )
                newElem = e(*params.args, **params.kwargs)
                res[newElem] = {}
            else:
                res[e] = self._callElemsTips(subs, params)
        return res

    def _getTipsAttribute(self, elems, attrName):
        """
        Return a list of attribute values from the 'elems' tree
        """
        res = []
        for e, subs in elems.items():
            if len(subs) == 0:
                try:
                    res.append(getattr(e, attrName))
                except:
                    raise SimPathInvalidPath(
                        f'Element {e} does not have any attribute named {attrName}'
                    )
            else:
                res += self._getTipsAttribute(subs, attrName)
        return res

    def _walk(self, pre=tuple(), elems=None, expand=True, combine=True):
        """
        Walk along the tree elements, yielding full paths from root to each tip.
        Expand elements that need to be expanded, like complexes or currents.
        """
        if elems is None:
            elems = self._elems

        for e, subs in elems.items():
            cls = e._simPathCombinerClass() if combine else None
            exp = e._simPathWalkExpand() if expand else [e]
            paths = []
            for e2 in exp:
                if len(subs) == 0:
                    paths.append(pre + (e2,))
                else:
                    for s in self._walk(pre + (e2,), subs, expand, combine):
                        paths.append(s)
            if cls is None:
                for path in paths:
                    yield path
            else:
                yield cls(*paths)

    def _getDescriptions(self, origDescr, pre=tuple(), elems=None):
        """
        Walk along the tree elements and yield the description of each final path.
        Descriptions are tuples of strings.
        """
        if elems is None:
            elems = self._elems

        for e, subs in elems.items():
            descriptions = []
            for e2 in e._simPathWalkExpand():
                if len(subs) == 0:
                    descriptions.append(pre + (str(e2),))
                else:
                    for s in self._getDescriptions(origDescr[1:], pre + (str(e2),), subs):
                        descriptions.append(s)
            if e._simPathCombinerClass() is None:
                for descr in descriptions:
                    yield descr
            else:
                yield pre + (str(e),) + origDescr[1:]

    def __getattr__(self, name):
        """Extend the path with a named element

        In most cases, the attribute access syntax is used to build paths, the two following lines
        are equivalent::

            sim.comp1.S1.Count
            sim.LIST('comp1').LIST('S1').Count

        Attribute access makes the paths more readable and should thus be preferred over the
        :py:func:`LIST` method when only one element is used.

        :param name: Name of the element
        :returns: An updated path taking into account the newly added element.
        :rtype: :py:class:`SimPath`

        .. note::
            This is a special python method, do not call explicitely sim.__getattr__(name).

        :meta public:
        """
        try:
            return SimPath(self._sim, self._namesExtendElems(self._elems, [name]))
        except SimPathExtensionError as ex:
            try:
                res = self._getTipsAttribute(self._elems, name)
                return res[0] if len(res) == 1 else res
            except SimPathExtensionError as ex2:
                raise SimPathInvalidPath(
                    f'{ex.elem.__class__.__name__} {ex.elem} has no children or attribute'
                    f'named {name}, it is not possible to extend the path.'
                )

    def __setattr__(self, name, value):
        """
        Prevent setting attributes except for attributes in _endNames.
        This avoids potential issues when making a typo in the endName:
        Without this, sim.comp.Ca.count = 0 (instead of Count) would not raise any
        error but would not set the count of Ca to 0 either.
        """
        if name in SimPath._endNames:
            super().__setattr__(name, value)
        else:
            raise AttributeError(
                f'Cannot set attribute {name}, only attributes in {SimPath._endNames.keys()} can be set.'
            )

    def __getitem__(self, key):
        """Specify further the last path element

        All objects that implement the `__getitem__` special method can be further specified
        in the path. The key is just passed to the corresponding object's `__getitem__` method.
        See:

        - :py:func:`steps.API_2.model.Complex.__getitem__`
        - :py:func:`steps.API_2.model.Channel.__getitem__`
        - :py:func:`steps.API_2.model.Reaction.__getitem__`
        - :py:func:`steps.API_2.model.Diffusion.__getitem__`
        - :py:func:`steps.API_2.model.Current.__getitem__`

        :param key: The specifier that will be passed to the last object(s) in the path
        :type key: Any

        :returns: An updated path taking into account the additional specification.
        :rtype: :py:class:`SimPath`

        Examples::

            sim.comp1.reac1['fwd'].K       # The reaction constant addociated with the forward
                                           # part of reac1 in comp1.
            sim.patch1.VGNaC[m3, h1].Count # The number of VGNaC in state [m3, h1] on patch1.
            sim.patch1.diff1[S1A, S1A].D   # Assuming diff1 is a diffusion rule that was declared
                                           # for a complex, this line returns the diffusion
                                           # constant for complexes in state [S1A, S1A] in patch1.

        :meta public:
        """
        return SimPath(self._sim, self._specifyElemsTips(self._elems, key))

    def __call__(self, *args, **kwargs):
        """Give additional information about the last path element

        All objects that implement the `__call__` special method can receive additional information.
        The arguments are just passed to the corresponding object's `__call__` method.
        See:

        - :py:func:`steps.API_2.model.Diffusion.__call__`

        :param \*args: The positional arguments that will be passed to the last object(s) in the
            path
        :type key: Any
        :param \*\*kwargs: The keyword arguments that will be passed to the last object(s) in the
            path
        :type key: Any

        :returns: An updated path taking into account the additional information.
        :rtype: :py:class:`SimPath`

        Examples::

            sim.TET(tet1).diff1(direc=tet2).D # Returns the diffusion constant from tet1 to tet2

        :meta public:
        """
        return SimPath(self._sim, self._callElemsTips(self._elems, nutils.Params(*args, **kwargs)))

    def __iter__(self):
        """Iterate through all elements by outputing the full path each time."""
        for path in self._walk():
            yield path

    @nutils.limitReprLength
    def __repr__(self):
        return ', '.join('.'.join(str(v) for v in path) for path in self._walk())


# Add the descriptors with their docstrings to SimPath
for _name, _docstr in SimPath._endNames.items():
    _descrobj = _SimPathDescr(_name)
    _descrobj.__doc__ = _docstr
    setattr(SimPath, _name, _descrobj)
    # Forbid the creation of objects that would have the same name because they would then be
    # inaccessible through a SimPath (descriptors have priority over __getattr__).
    nutils.NamedObject._forbiddenNames.add(_name)


class MPI:
    """Holds MPI related information and links to mpi packages

    MPI packages are only loaded if required and the user does not have to load
    them manually, changing the solver to an MPI-based one does it automatically.
    """

    EF_NONE = stepslib._py_TetAPI.EF_NONE
    """Possible value for the calcMembPot parameter of solvers that implement EField.
    Means that no EField solver is needed."""
    EF_DEFAULT = stepslib._py_TetAPI.EF_DEFAULT
    """Possible value for the calcMembPot parameter of solvers that implement EField.
    Means that serial EField simulation (Tetexact version) should be run on process 0."""
    EF_DV_BDSYS = stepslib._py_TetAPI.EF_DV_BDSYS
    """Possible value for the calcMembPot parameter of solvers that implement EField.
    Means that parallel SuperLU EField solver should be used."""
    EF_DV_PETSC = stepslib._py_TetAPI.EF_DV_PETSC
    """Possible value for the calcMembPot parameter of solvers that implement EField.
    Means that parallel PETSc EField solver should be used."""

    _usingMPI = None
    _shouldWrite = True

    _solverNameMapping = {
        'TetOpSplit': nutils._CYTHON_PREFIX + 'TetOpSplitP',
        'DistTetOpSplit': nutils._CYTHON_PREFIX + 'DistTetOpSplitP',
    }

    @classmethod
    def _loadInfos(cls):
        if cls._usingMPI is None:
            if hasattr(stepslib, 'mpiInit'):
                stepslib.mpiInit()
                atexit.register(stepslib.mpiFinish)
                cls._usingMPI = True
                cls._shouldWrite = stepslib.getRank() == 0

                # Force stderr flush when an exception is raised.
                # Without this, under some conditions, if one process raises a python exception,
                # the corresponding message is not always printed out as it should be.
                def customHook(tpe, val, bt):
                    sys.__excepthook__(tpe, val, bt)
                    sys.stderr.flush()
                    stepslib.mpiAbort()

                sys.excepthook = customHook
            else:
                raise ImportError(
                    f'Could not load cysteps_mpi.so. Please check that STEPS was built with MPI '
                    f'support.'
                )

    @nutils.classproperty
    def rank(cls):
        """Get the rank of the current process

        :type: int, read-only

        .. warning::
            Accessing this property triggers the importing of MPI related packages
        """
        cls._loadInfos()
        return stepslib.getRank()

    @nutils.classproperty
    def nhosts(cls):
        """Get the number of hosts

        :type: int, read-only

        .. warning::
            Accessing this property triggers the importing of MPI related packages
        """
        cls._loadInfos()
        return stepslib.getNHosts()

    @classmethod
    def _getSolver(cls, name):
        cls._loadInfos()
        # Solvers that require special mapping
        if name in cls._solverNameMapping:
            return getattr(stepslib, cls._solverNameMapping[name])
        else:
            # Default case
            return getattr(stepslib, nutils._CYTHON_PREFIX + name)


class _SimulationCheckpointer:
    """Class that manages the auto checkpointing

    Uses the same mechanism (i.e. the _nextSave heapq in Simulation) as result selectors.
    It only implement the ResultSelector methods called by Simulation.run() and Simulation.newRun().
    """

    def __init__(self, sim):
        self._sim = sim
        self._period = math.inf
        self._prefix = ''
        self._onlyLast = False

        self._nextTime = math.inf
        self._tind = 0
        self._startTime = 0
        self._lastName = None

    def setup(self, time, period, prefix, onlyLast):
        if period is None:
            self._period = math.inf
            self._nextTime = math.inf
            self._startTime = math.inf
        else:
            self._period = period
            self._startTime = time
            self._nextTime = time
        self._prefix = prefix
        self._onlyLast = onlyLast
        self._tind = 0
        self._lastName = None

    def _newRun(self):
        self._tind = -1
        self._updateNextSaveTime()

    def _save(self, time, _):
        newName = f'{self._prefix}_{self._sim._runId}_{time}_{self._sim._solverStr}_cp'
        self._sim.checkpoint(newName)
        if self._onlyLast and self._lastName is not None and MPI._shouldWrite:
            os.remove(self._lastName)
        self._lastName = newName
        self._updateNextSaveTime()

    def _updateNextSaveTime(self):
        self._tind += 1
        if self._period == math.inf:
            self._nextTime = math.inf
        else:
            self._nextTime = self._startTime + self._tind * self._period


@nutils.FreezeAfterInit
class Simulation(nutils.NamedObject, nutils.StepsWrapperObject, nutils.ParameterizedObject):
    """The main simulation class

    A simulation is defined from a model, a geometry, a solver and a random number generator. Once
    a :py:class:`Simulation` object is created, users can access and set simulated values through
    the simulation path syntax (see :py:class:`SimPath`). :py:func:`newRun` should be called before
    the simulation can be run and the data to be saved should be declared with
    :py:class:`steps.API_2.saving.ResultSelector` objects and added to the simulation with the
    :py:func:`toSave` method. The simulation is advanced by calling either :py:func:`step` or
    :py:func:`run`.

    :param solverName: Name of the solver to be used for the simulation (see SERIAL_SOLVERS and
        PARALLEL_SOLVERS below).
    :type solverName: str
    :param mdl: The model to be used in the simulation
    :type mdl: :py:class:`steps.API_2.model.Model`
    :param geom: The geometry to be used in the simulation
    :type geom: Union[:py:class:`steps.API_2.geom.Geometry`, :py:class:`steps.API_2.geom.TetMesh`,
        :py:class:`steps.API_2.geom.DistMesh`]
    :param rng: The random number generator to be used in the simulation
    :type rng: :py:class:`steps.API_2.rng.RNG`
    :param \*args: Positional arguments forwarded to the solver constructor, see
        :py:mod:`steps.API_1.solver` or :py:mod:`steps.API_1.mpi.solver`.
    :param \*\*kwargs: Keyword arguments forwarded to the solver constructor, see
        :py:mod:`steps.API_1.solver` or :py:mod:`steps.API_1.mpi.solver`.
    """

    SERIAL_SOLVERS = ['Wmdirect', 'Wmrssa', 'Wmrk4', 'Tetexact', 'TetODE']
    """Available serial solvers"""

    PARALLEL_SOLVERS = ['TetOpSplit', 'DistTetOpSplit']
    """Available parallel solvers"""

    def __init__(self, solverName, mdl, geom, rng=None, *args, name=None, **kwargs):
        super().__init__(name=name)

        self._model = mdl
        self._geom = geom
        self._rng = rng

        self._solverStr = solverName
        self.stepsSolver = self._createSolver(solverName, *args, **kwargs)

        self._resultSelectors = []
        self._nextSave = None
        self._runId = -1

        self._checkpointer = _SimulationCheckpointer(self)

        nsimcheck.Check(self, True)

    @property
    def model(self):
        """The model associated to the simulation

        :type: :py:class:`steps.API_2.model.Model`, read-only
        """
        return self._model

    @property
    def geom(self):
        """The geometry associated to the simulation

        :type: Union[:py:class:`steps.API_2.geom.Geometry`, :py:class:`steps.API_2.geom.TetMesh`,
            :py:class:`steps.API_2.geom.DistMesh`] read-only
        """
        return self._geom

    @property
    def rng(self):
        """The random number generator associated to the simulation

        :type: :py:class:`steps.API_2.rng.RNG`, read-only
        """
        return self._rng

    @property
    def solver(self):
        """The solver associated to the simulation

        :type: :py:class:`steps.API_1.solver._Base_Solver`, read-only

        See :py:mod:`steps.API_1.solver` for detailed API.

        .. warning::
            Direct calls to the solver should only be made if absolutely necessary. They will
            bypass data saving mechanisms implemented in :py:class:`Simulation`.
        """
        return self.stepsSolver

    def step(self):
        """Advance the simulation for one 'step'.

        In stochastic solvers this is one 'realization' of the Gillespie SSA (one reaction
        'event'). In numerical solvers (currently Wmrk4 and TetODE) this is one time-step,
        with the stepsize defined with the :py:func:`setDT` method.

        .. warning::
            If the step happens to jump over a ResultSelector saving timepoint,
            the data will be saved after the step and the associated timepoint will be
            modified to match the time after the step. If several saving timepoints are
            jumped over, the save will be performed only once. The same remark applies to
            automatic checkpointing (see :py:func:`autoCheckpoint`).
        """
        if self._nextSave is None:
            raise Exception(f'Cannot call step before calling newRun on the simulation.')

        self._discardPastSaves(currTime)

        self.stepsSolver.step()
        currTime = self.Time

        while len(self._nextSave) > 0 and self._nextSave[0][0] <= currTime:
            oldSave = self._nextSave[0]
            while self._nextSave[0][0] == oldSave[0]:
                rstime, rs = self._nextSave[0]
                rs._save(currTime, currTime)
                while rs._nextTime <= currTime:
                    rs._updateNextSaveTime()
                oldSave = heapq.heapreplace(self._nextSave, (rs._nextTime, rs))

    def run(self, t):
        """Run the simulation until a given time

        Automatically save the data that has been specified through ResultSelectors and added
        to the simulation with the :py:func:`toSave` method.
        Automatically checkpoint if the :py:func:`autoCheckpoint` method has been called.

        :param t: Run the simulation until this time (in seconds)
        :type t: float
        """
        if self._nextSave is None or self._runId < 0:
            raise Exception(f'Cannot call run before calling newRun() on the simulation.')

        currT = self.stepsSolver.getTime()

        self._discardPastSaves(currT)

        while len(self._nextSave) > 0 and self._nextSave[0][0] <= t:
            # Do not run if the time is identical to current time
            if self._nextSave[0][0] > currT:
                self.stepsSolver.run(self._nextSave[0][0])
                currT = self._nextSave[0][0]

            oldSave = self._nextSave[0]
            while self._nextSave[0][0] == oldSave[0]:
                rstime, rs = self._nextSave[0]
                rs._save(rstime, rstime)
                oldSave = heapq.heapreplace(self._nextSave, (rs._nextTime, rs))

        if t > currT:
            self.stepsSolver.run(t)

    def newRun(self):
        """Reset the solver and signal the start of a new run

        This method must be called before any call to :py:func:`run` or :py:func:`step`. Since it
        resets the solver, the initial state must be specified after the call to newRun.
        Result selectors rely on this call to separate the data from different runs.

        .. note::
            The ``'tetODE'`` solver cannot be reset.
        """
        if self.stepsSolver.__class__ != stepslib._py_TetODE:
            self.stepsSolver.reset()
        else:
            warnings.warn(
                f'Cannot reset a TetODE solver, a new run was started but the solver was not reset.'
            )

        # Update data saving structures
        self._newRun()

    def toSave(self, *selectors, dt=None, timePoints=None):
        """Add result selectors to the simulation

        :param selectors: One or several result selectors that should be saved during the
            simlulation.
        :type selectors: :py:class:`steps.API_2.saving.ResultSelector`
        :param dt: If this keyword parameter is specified, the result selectors will be saved every
            ``dt`` seconds.
        :type dt: float
        :param timePoints: A list of timepoints (in seconds) at which the result selectors should
            be saved.
        :type timePoints: List[float]
        """
        # Only optimize calls that will be called together, other selectors should be optimized independently
        if dt is not None or timePoints is not None:
            selectors = nsaving_optim.OptimizeSelectors(self, selectors)
        else:
            selectors = [nsaving_optim.OptimizeSelectors(self, [sel])[0] for sel in selectors]

        for i, rs in enumerate(selectors):
            if isinstance(rs, nsaving.ResultSelector):
                if dt is not None:
                    rs._saveWithDt(dt)
                if timePoints is not None:
                    rs._saveWithTpnts(timePoints)
                rs._addedToSimulation(i)
                self._resultSelectors.append(rs)
            else:
                raise Exception(f'Expected a ResultSelector object, got {rs} instead.')

    def toDB(self, dbh, uid, **kwargs):
        """Redirect all the added results selectors to a database

        :param dbh: The database to which the result selectors should be saved (see e.g.
            :py:class:`steps.API_2.saving.SQLiteDBHandler`).
        :type dbh: :py:class:`steps.API_2.saving.DatabaseHandler`
        :param uid: A unique identifier under which all subsequent runs should be saved
        :type uid: str
        :param kwargs: Any additional parameters that should be saved to the database along with
            the unique identifier. Values are restricted to
        :type kwargs: int, float, str or bytes

        Usage::

            sim.toSave(...) # Add the result selectors

            with SQLiteDBHandler(dbPath) as dbh:
                sim.toDB(dbh, 'MySimulation', val1=1, val2=2)
                sim.newRun()
                ... # Set initial state
                sim.run(10)
        """
        if not isinstance(dbh, nsaving.DatabaseHandler):
            raise TypeError(f'Expected a DatabaseHandler, got {dbh} instead.')
        if MPI._shouldWrite:
            dbh._newGroup(uid, self._resultSelectors, **kwargs)
            for rs in self._resultSelectors:
                rs._toDB(dbh)

    def autoCheckpoint(self, period, prefix='', onlyLast=False):
        """Activates automatic checkpointing

        After this method has been called, all subsequent calls to :py:func:`step` or
        :py:func:`run` will save a checkpoint every ``period`` seconds. This method can be called
        during simulation to change the checkpoint period or prefix. To turn off automatic
        checkpointing, use ``sim.autoCheckpoint(None)``.

        :param period: The period at which checkpointing should be performed, in seconds. None to
            disable automatic checkpointing.
        :type period: float or None
        :param prefix: An optional prefix to the checkpoint filename. Checkpoint file names follow
            the pattern ``'{prefix}_{runId}_{time}_{solver}_cp'``.
        :type prefix: str
        :param onlyLast: If True, remove previous checkpoints when a new chekpoint is being saved.
        :type onlyLast: bool
        """
        if not isinstance(period, numbers.Number) and period is not None:
            raise TypeError(f'Expected a number or None, got {period} instead.')
        if not isinstance(prefix, str):
            raise TypeError(f'Expected a string for prefix parameter, got {prefix} instead.')
        if period is not None and period <= 0:
            raise ValueError(f'The auto checkpoint period needs to be strictly positive, got {period}.')
        self._checkpointer.setup(self.Time, period, prefix, onlyLast)
        if self._nextSave is not None:
            self._initNextSave()

    def _initNextSave(self):
        """Setup the heapq that orders the save events"""
        eventLst = self._resultSelectors + [self._checkpointer]
        self._nextSave = [(rs._nextTime, rs) for rs in eventLst]
        heapq.heapify(self._nextSave)

    def _discardPastSaves(self, currT):
        """Discard saving points that are already passed."""
        while self._nextSave[0][0] < currT:
            rstime, rs = self._nextSave[0]
            rs._updateNextSaveTime()
            heapq.heapreplace(self._nextSave, (rs._nextTime, rs))

    def _newRun(self):
        """Update data saving structures and setup the next save point."""
        for rs in self._resultSelectors:
            rs._newRun()
        self._checkpointer._newRun()
        self._initNextSave()
        self._runId += 1

    def _getReferenceObject(self):
        """
        Return the object this object was derived from. Useful for getting the complex associated
        with a complex selector, etc.
        """
        # Use geom children as the simulation children
        return self.geom

    def _getSubParameterizedObjects(self):
        """Return all subobjects that can hold parameters."""
        return [self.model, self.geom]

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds."""
        return [self.stepsSolver]

    def __getattr__(self, name):
        """Redirect attribute access to a SimPath

        See :py:func:`SimPath.__getattr__`

        :meta public:
        """
        return getattr(SimPath(self), name)

    def ALL(self, *cls):
        """Redirect ALL(...) to a SimPath

        See :py:func:`SimPath.ALL`
        """
        return SimPath(self).ALL(*cls)

    def _createSolver(self, solverStr, *args, **kwargs):
        """Create and return a new solver."""
        # Set up structures depending on model
        self.model._SetUpMdlDeps()
        self.geom._SetUpMdlDeps(self.model)

        stepsGeom = self.geom._getStepsObjects()[0]

        # Get the solver class
        # Check if it is available in cysteps (not cysteps_mpi)
        try:
            if solverStr in Simulation.SERIAL_SOLVERS:
                solverCls = getattr(stepslib, nutils._CYTHON_PREFIX + solverStr)
            elif solverStr in Simulation.PARALLEL_SOLVERS:
                solverCls = MPI._getSolver(solverStr)
            else:
                raise AttributeError()
        except AttributeError:
            raise Exception(f'Unknown solver: {solverStr}')

        # Process arguments
        # Iterate through arguments and unpack nutils.Params if any is present.
        rargs = []
        rkwargs = {}
        for arg in args:
            if isinstance(arg, nutils.Params):
                rargs += arg.args
                rkwargs.update(arg.kwargs)
            else:
                rargs.append(arg)
        for key, arg in kwargs.items():
            if isinstance(arg, nutils.Params):
                rargs += arg.args
                rkwargs.update(arg.kwargs)
            else:
                rkwargs[key] = arg

        stepsRng = self.rng.stepsrng if self.rng is not None else None
        solver = solverCls(self.model.stepsModel, stepsGeom, stepsRng, *rargs, **rkwargs)

        return solver

    # Advanced parameters column names
    _ADV_PARAMS_RUN_NUMBER = 'Run #'
    _ADV_PARAMS_RUN_TIME = 'Time'
    _ADV_PARAMS_SIM_PATH = 'Simulation Path'

    def _parameterSimPathCallback(self, path, param):
        """Keep track of solver setter calls when used with Parameter objects"""
        pathStr = '.'.join(str(v) for v in path[1:])
        key = (
            (self._ADV_PARAMS_RUN_NUMBER, self._runId),
            (self._ADV_PARAMS_RUN_TIME, self.Time),
            (self._ADV_PARAMS_SIM_PATH, pathStr)
        )
        self._setAdvancedParameter(key, param)

    @classmethod
    def _getGroupedAdvParams(cls):
        """Return a list of advanced parameter property that should be grouped
        """
        return [cls._ADV_PARAMS_RUN_NUMBER, cls._ADV_PARAMS_RUN_TIME]

    #####################
    # Solver properties #
    #####################

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('s'))
    def EfieldDT(self):
        """The stepsize for membrane potential solver

        This is the time for each voltage calculation step. The SSA will run until passing this
        stepsize, so in fact each membrane potential time step will vary slightly around the dt
        so as to be aligned with the SSA.

        Default value is 1us.

        :type: float
        """
        return self.stepsSolver.getEfieldDT()

    @EfieldDT.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('s'))
    def EfieldDT(self, val):
        self.stepsSolver.setEfieldDT(val)

    @property
    @nutils.ParameterizedObject.RegisterGetter(units=nutils.Units('K'))
    def Temp(self):
        """The simulation temperature in Kelvins

        Currently, this will only influence the GHK flux rate, so will only influence simulations
        including membrane potential calculation.

        :type: float
        """

        return self.stepsSolver.getTemp()

    @Temp.setter
    @nutils.ParameterizedObject.RegisterSetter(units=nutils.Units('K'))
    def Temp(self, val):
        self.stepsSolver.setTemp(val)

    @property
    def Time(self):
        """The current simulation time in seconds

        :type: float, read-only
        """
        return self.stepsSolver.getTime()

    @property
    def A0(self):
        """The total propensity of the current simulation state

        The total propensity multiplied by an infinitesimally small time dt gives the probability
        that a reaction will occur in that dt. For Tetexact this includes the propensity from the
        extension of the SSA for diffusive flux between tetrahedral elements in the mesh.

        :type: float, read-only
        """

        return self.stepsSolver.getA0()

    @property
    def NSteps(self):
        """The number of 'realizations' of the SSA

        It corresponds to the total number of reaction (and diffusion) events in stochastic solvers.

        :type: float, read-only
        """
        return self.stepsSolver.getNSteps()

    #####################
    # Solver methods    #
    #####################

    @nutils.ParameterizedObject.SpecifyUnits(None, nutils.Units('s'))
    def setDT(self, dt):
        """Set the stepsize for numerical solvers

        Superceded by setRk4DT, but Included for backwards compatability.

        :param dt: The step size in seconds
        :type dt: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        self.stepsSolver.setDT(dt)

    @nutils.ParameterizedObject.SpecifyUnits(None, nutils.Units('s'))
    def setRk4DT(self, dt):
        """Set the stepsize for numerical solvers

        Must be called before running a simulation with these solvers (currently Wmrk4) since
        there is no default stepsize. The deterministic solver Wmrk4 implements a fixed stepsize
        (i.e. not adaptive), although the stepsize can be altered at any point during the
        simulation with this method.

        :param dt: The step size in seconds
        :type dt: Union[float, :py:class:`steps.API_2.nutils.Parameter`]
        """
        self.stepsSolver.setRk4DT(dt)

    def setEfieldTolerances(self, atol=1e-50, rtol=1e-5):
        """Set the absolute and relative tolerances for the Efield solver

        See https://petsc.org/release/docs/manual/ksp/#convergence-tests

        :param atol: Absolute tolerance
        :type atol: float
        :param rtol: Relative tolerance
        :type rtol:
        """
        self.stepsSolver.setEfieldTolerances(atol, rtol)

    def setTolerances(self, atol, rtol):
        """Set the absolute tolerance and the relative tolerance for CVODE

        Only works for ODE solvers.

        :param atol: Absolute tolerance for CVODE
        :type atol: float
        :param rtol: Relative tolerance for CVODE
        :type rtol: float
        """
        self.stepsSolver.setTolerances(atol, rtol)

    def checkpoint(self, fname):
        """Checkpoint the current simulation state to a file

        See :py:func:`autoCheckpoint` for setting up automatic checkpointing.

        :param fname: The file name / path
        :type fname: str
        """
        self.stepsSolver.checkpoint(fname)

    def restore(self, fname):
        """Restore the simulation state to a previously saved checkpoint

        :param fname: The file name / path
        :type fname: str
        """
        self.stepsSolver.restore(fname)

    def saveMembOpt(self, fname):
        """Saves the vertex optimization in the Efield structure

        :param fname: The file name / path
        :type fname: str
        """
        self.stepsSolver.saveMembOpt(fname)


# Forbid the creation of objects that would have the name of a Simulation method or property
# because they would then be inaccessible through a SimPath (methods and attributes have priority
# over __getattr__).
for _name in dir(Simulation) + dir(SimPath):
    nutils.NamedObject._forbiddenNames.add(_name)


###################################################################################################
# SBML Import


@nutils.FreezeAfterInit
class SBMLSimulation(Simulation):
    """Simulation loaded from an SBML model

    :param solverName: Name of the solver to be used for the simulation (see :py:class:`Simulation`).
    :type solverName: str
    :param filepath: Path to the SBML file
    :type filepath: str
    :param rng: The random number generator to be used in the simulation
    :type rng: :py:class:`steps.API_2.rng.RNG`
    :param maxDt: Time period between rate updates (for time-dependent rates)
    :type maxDt: float
    :param \*args: Positional arguments forwarded to the solver constructor, see
        :py:mod:`steps.API_1.solver` or :py:mod:`steps.API_1.mpi.solver`.
    :param \*\*kwargs: Keyword arguments forwarded to the solver constructor, see
        :py:mod:`steps.API_1.solver` or :py:mod:`steps.API_1.mpi.solver`.

    Optional SBML parameters:

    :param timeunits_def: Default time unit in relation to seconds (see below)
    :type timeunits_def: float
    :param volunits_def: Default volume unit in relation to cubic meters (see below)
    :type volunits_def: float
    :param subsunits_def: Default substance unit in relation to moles (see below)
    :type subsunits_def: float
    :param volume_def: Default volume (see below)
    :type volume_def: float
    :param area_def: Default area (see below)
    :type area_def: float
    :param strict_mode: False mean that all reactions will be allowed, True means that only
        fundamental reactions that can be represented in the SSA without modification will be
        allowed, the import failing if any other types of reactions are present. The vast
        majority of SBML models include some non-fundamental reactions so setting this to True
        will mean that most SBML imports will fail.
    :type strict_mode: bool
    :param solverDt: If using a numerical solver, the integration time step in seconds. If not
        specified, the default value will be set to one tenth of the ``maxDt`` value.
    :type solverDt: float

    The 'units' arguments are a means to set the default model units and unit itself shoud relate
    to s.i units. For example: to set default model units to ms, ``timeunits_def=1.0e-3``. To set
    default model substance units to micromoles, ``subsunits_def=1.0e-6``.

    .. note::
        These default units will only be used when units are not explicitly declared in the SBML
        file.

    The volume_def and area_def arguments are a way to set volume of compartments or area of
    patches with size 1.0 and no units, which are very common in SBML files. This may allow for a
    model to be run stochastically without modifying the file. Values should be given in s.i.
    units, based on metres. For example: to set default volume to 100 femtolitre:
    ``volume_def=100*1.0e-15*1.0e-3``
    """

    def __init__(
        self,
        solverName,
        filepath,
        rng,
        maxDt,
        *args,
        timeunits_def=1.0,
        volunits_def=1.0e-3,
        subsunits_def=1.0,
        volume_def=False,
        area_def=False,
        strict_mode=False,
        solverDt=None,
        **kwargs,
    ):

        self._solvStr = solverName
        self._solvArgs = args
        self._solvKwargs = kwargs
        try:
            # Not great to import inside code but the sbml module prints to stout when libsbml is
            # lacking. This was fine when it was separately imported by the user but it is not ok
            # now.
            import steps.API_1.utilities.sbml as ssbml

            self._sbml = ssbml.Interface(
                filepath,
                timeunits_def=timeunits_def,
                volunits_def=volunits_def,
                subsunits_def=subsunits_def,
                volume_def=volume_def,
                area_def=area_def,
                strict_mode=strict_mode,
            )
        except NameError:
            raise ImportError(f'Could not import libSBML, check that it is installed.')
        self._maxDt = maxDt
        self._dtInd = 0

        mdl, geom = self._getMdlGeomFromSbml()

        super().__init__(self._solvStr, mdl, geom, rng, *self._solvArgs, **self._solvKwargs)

        if hasattr(self.stepsSolver, 'setDT'):
            self.stepsSolver.setDT(self._maxDt / 10 if solverDt is None else solverDt)

    def newRun(self):
        """Reload completely the SBML model and update the ResultSelectors

        See :py:func:`Simulation.newRun` for more details.
        """
        super().newRun()
        self._sbml.setupSim(self.stepsSolver)
        self._dtInd = 0

    def run(self, t):
        """Run to time t in increments of at most maxDt

        See :py:func:`Simulation.run` for more details.
        """
        while t - self._dtInd * self._maxDt >= self._maxDt:
            self._dtInd += 1
            super().run(self._dtInd * self._maxDt)
            self._sbml.updateSim(self.stepsSolver, self._maxDt)
        super().run(t)

    def _getMdlGeomFromSbml(self):
        # Temporarily turn automatic name change for forbidden names
        nutils.NamedObject._allowForbNames = True

        mdl = nmodel.Model._FromStepsObject(self._sbml.getModel())

        sgeom = self._sbml.getGeom()
        if isinstance(sgeom, stepslib._py_Tetmesh):
            geom = ngeom.TetMesh._FromStepsObject(sgeom)
        elif isinstance(sgeom, stepslib._py_Geom):
            geom = ngeom.Geometry._FromStepsObject(sgeom)
        else:
            raise NotImplementedError()

        nutils.NamedObject._allowForbNames = False

        return mdl, geom
