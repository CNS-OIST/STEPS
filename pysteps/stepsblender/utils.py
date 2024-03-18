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

import colorsys
import enum
import importlib
import inspect
import numpy as np
import typing

####################################################################################################


class Orders(enum.IntEnum):
    GET_MESH = 0
    EXIT = 1
    GET_DATA = 2
    END = 3
    GET_MODEL = 4
    GET_ELEM_SPEC_COUNT = 5,
    GET_VES_IN_SPEC_COUNT = 6,
    GET_VES_SURF_SPEC_REL_POS = 7,
    GET_VES_POS = 8
    GET_RAFT_POS = 9
    GET_RAFT_COUNTS = 10
    GET_VES_LINKSPEC_REL_POS = 11
    GET_VES_EVENTS = 12
    GET_RAFT_EVENTS = 13
    GET_VERTS_V = 14


class Loc(enum.Enum):
    VERT = 0
    TRI = 1
    TET = 2
    VES_IN = 3
    VES_SURF = 4
    RAFT_IN = 5


class Event(enum.Enum):
    EXOCYTOSIS = 0
    ENDOCYTOSIS = 1
    RAFT_ENDOCYTOSIS = 2


####################################################################################################


def classFromStr(clsFullName):
    *modules, clsName = clsFullName.split('.')
    mod = importlib.import_module('.'.join([__package__] + modules))
    return getattr(mod, clsName)


def makeCreatableFromStr(tpe):
    origCls = typing.get_origin(tpe)
    if origCls is str:
        return origCls

    class creatableType(tpe):

        def __new__(cls, *args, **kwargs):
            if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], str):
                return origCls(eval(args[0]))
            else:
                return origCls(*args, **kwargs)

    return creatableType


colorType = typing.Annotated[
    makeCreatableFromStr(typing.Tuple[float, float, float, float]),
    dict(help='An RGBA tuple that should be supplied with quotes e.g. "(0.1, 0.5, 0.1, 1)"',
         metavar='"(r, g, b, a)"'), ]
alphaType = typing.Annotated[float, 'Alpha transparency level']
emissionType = typing.Annotated[float, 'Emission strength']

####################################################################################################


def spherical2Cartesian(spos):
    if len(spos.shape) == 1:
        return spherical2Cartesian(np.array([spos]))[0, :]
    pos = np.empty(spos.shape)
    pos[:, 0] = spos[:, 0] * np.cos(spos[:, 1]) * np.sin(spos[:, 2])
    pos[:, 1] = spos[:, 0] * np.sin(spos[:, 1]) * np.sin(spos[:, 2])
    pos[:, 2] = spos[:, 0] * np.cos(spos[:, 2])
    return pos


def cartesian2Spherical(pos):
    if len(pos.shape) == 1:
        return cartesian2Spherical(np.array([pos]))[0, :]
    spos = np.empty(pos.shape)
    spos[:, 0] = np.linalg.norm(pos, axis=1)  # Radius
    spos[:, 1] = np.arctan2(pos[:, 1], pos[:, 0])  # Theta
    spos[:, 2] = np.arccos(pos[:, 2] / spos[:, 0])  # Phi
    return spos


def sphericalInterpolation(spos1, spos2, ratio):
    if abs(spos1[1] - spos2[1]) > np.pi:
        if spos1[1] < spos2[1]:
            spos1[1] += 2 * np.pi
        else:
            spos2[1] += 2 * np.pi
    return spos1 + (spos2 - spos1) * ratio


def _noneGenerator():
    while True:
        yield None


def zipNone(*args):
    sentinel = object()
    iters = list(map(iter, [arg if arg is not None else _noneGenerator() for arg in args]))
    while iters:
        res = []
        for it in iters:
            elem = next(it, sentinel)
            if elem is sentinel:
                return
            res.append(elem)
        yield tuple(res)


####################################################################################################

_ALL_COLOR_HUES = []


def GetColor(s=1, v=1, a=1):
    global _ALL_COLOR_HUES
    if len(_ALL_COLOR_HUES) == 0:
        _ALL_COLOR_HUES += [0, 1]
        h = 0
    else:
        i, _ = max(enumerate(_ALL_COLOR_HUES[1:]), key=lambda x: x[1] - _ALL_COLOR_HUES[x[0]])
        h = (_ALL_COLOR_HUES[i + 1] + _ALL_COLOR_HUES[i]) / 2
        _ALL_COLOR_HUES = _ALL_COLOR_HUES[:i + 1] + [h] + _ALL_COLOR_HUES[i + 1:]
    return colorsys.hsv_to_rgb(h, s, v) + (a, )


####################################################################################################


class HierarchicalParameters:
    """Class for holding hierarchies of parameters

    This class holds several levels of parameter dictionaries and allows the retrieval of
    parameter values by iterating through levels, starting with the most specific.
    The hierarchy is built from nested dictionnaries, for example:

    parameters = {
        'color': (1, 0, 0, 1),         # Red
        'Species': {
            'color': (0, 1, 0, 1),     # Green
            'radius': 0.01,
            'S1': {
                'radius': 0.02,
                'color': (0, 0, 1, 1), # Blue
            },
        },
    }

    If these parameters are given to :py:class:`HDF5BlenderLoader`, all objects that declared a
    `radius` or `color` attribute in their class (see :py:class:`HierarchicalParamReader`),
    will get a value that depends on their position in the object hierarchy (abridged here):

    Loader -----> Species ---> S1 --> mesh     | radius == 0.02
           \              \       \-> material | color  == Blue
            \              \-> S2 --> mesh     | radius == 0.01
             \                    \-> material | color  == Green
              \-> Vesicles --> V1 --> mesh     |
                                  \-> material | color  == Red

    When retrieving the color value for the material of species S2, we first try to find the most
    specific value: parameters['Species']['S2']['material']['color'], if this does not exist, we then
    try parameters['Species']['S2']['color'], then parameters['Species']['color'] which exist in our
    example, so the color is set to green.
    For species S1, parameters['Species']['S1']['material']['color'] does not exist but
    parameters['Species']['S1']['color'] does, the color is thus set to blue.

    In addition to specifying hierarchies of parameter values, one can also change the class that
    will be used to instantiate any object in the hierarchy, for example, to provide a custom material
    for species of type S2, we would give:

    parameters = {
        'S2' : {
            'material': {
                '__class__': MyCustomMaterialClass,
                'myCustomParameter': 5.0,
            },
        },
    }

    With MyCustomMaterialClass inheriting from :py:class:`BlenderMaterial` and having
    `myCustomParameter` as class attribute (see :py:class:`HierarchicalParamReader`).
    Note that we did not have to specify the full hierarchy, we skipped the `Species` object, which
    means that if we have a link species called S2, the custom material will be applied to it.

    Finally, parameter objects that inherit from :py:class:`BlenderWrapper` can be loaded from the
    Blender file by giving the name of the blender object. For example, if we created a material in
    Blender called 'myCustomMaterial', we could assign it to Species S2 with:

    parameters = {
        'S2' : {
            'material': 'myCustomMaterial',
        },
    }
    """

    def __init__(self, parameters={}, _kwargs=[]):
        if isinstance(parameters, dict):
            parameters = [parameters]
        # Dictionaries ranked from most specific to least specific
        self._parameters = parameters
        self._kwargs = _kwargs

    def withKwargs(self, kwargs):
        if len(kwargs) > 0:
            return HierarchicalParameters(self._parameters, [kwargs] + self._kwargs)
        else:
            return self

    def _getClsFromDefault(self, default):
        from .objects import BlenderWrapper
        if inspect.isclass(default) and issubclass(default, BlenderWrapper):
            for cls in default.__mro__:
                if cls in BlenderWrapper.__subclasses__():
                    return cls
        return None

    def get(self, paramName, default=None, cls=None, **kwargs):
        if cls is None:
            cls = self._getClsFromDefault(default)

        objCls = None

        newParams = HierarchicalParameters([], self._kwargs)

        for dct in self._parameters + self._kwargs:
            if paramName in dct:
                val = dct[paramName]
                if not inspect.isclass(cls) or not issubclass(cls, HierarchicalParamReader):
                    return val
                else:
                    if isinstance(val, dict):
                        newParams._parameters.append(val)
                        if '__class__' in val and objCls is None:
                            objCls = val['__class__']
                            if not inspect.isclass(objCls) or not issubclass(objCls, cls):
                                raise ValueError(
                                    f'Expected a class that inherits from {cls}, got {objCls} instead.')
                    elif objCls is None:
                        newParams._parameters += self._parameters
                        if inspect.isclass(val) and issubclass(val, cls):
                            return val(parameters=newParams, **kwargs)
                        else:
                            return cls(name=val, parameters=newParams, **kwargs)

        newParams._parameters += self._parameters
        if objCls is not None:
            return objCls(parameters=newParams, **kwargs)
        else:
            # The parameter was not found in any of the dicts
            if inspect.isclass(default):
                return default(parameters=newParams, **kwargs)
            else:
                return default


class HierarchicalParamReader:
    """Base class for all classes that can be initialized from hierarchical parameters

    The parameter hierarchy mirrors the hierarchy of :py:class:`HierarchicalParamReader` objects.
    Classes that inherit from :py:class:`HierarchicalParamReader` should declare the parameters that
    they require as class attributes. The value of the class attribute will be the default
    value given to the instance attribute.

    Example::

        class classA(HierarchicalParamReader):
            val1 = 1.0
            val2 = 'str'

        a = classA(parameters={val1: 2})

    This code leads to `a.val1 == 2` and `a.val2 == 'str'`.
    
    If the default value is a class that inherits from :py:class:`HierarchicalParamReader`, an
    object from this class will be instantiated and recursively initialized with the
    HierarchicalParameters object.

    By default, all :py:class:`HierarchicalParamReader` instances have a `parent` and `nameInParent`
    attributes that will be filled in at instanciation.
    
    Example::
        
        class classB(HierarchicalParamReader):
            objA = classA
            val3 = 5.0

        objB = classB(parameters={'objA':{'val1':2}})

    This code leads to objA being automatically instantiated with:
        objB.objA.parent == objB
        objB.objA.nameInParent == 'objA'

    When creating an object, one can also pass keyword arguments to the constructors, they will be
    treated as if they were part of the hierarchical parameters but will have lower priority.

    Example::

        a = classA(parameters={val1: 2})           # with the `parameters` object only

        a = classA(val1=2)                         # Equivalent to previous call

        a = classA(parameters={'val1': 3}, val1=2) # val1 will be initialized to 3,
                                                   # because it has higher priority

    Keyword argument can thus be treated as default values given at the time of instanciation, they
    will override the default value given in the class declaration, but will be overridden by
    parameters given in the :py:class:`HierarchicalParameters` object.

    If the parameter should be of a given class that inherits from
    :py:class:`HierarchicalParamReader` but its default value should be None, the class should be
    specified through python annotations.

    Example::

        class classC(HierarchicalParamReader):
            objB: classB = None

        c1 = classC(parameters={})                     # Will lead to c1.objB == None

        c2 = classC(parameters={'objB': {'val3: 10'}}) # Will lead to c1.objB being an object of
                                                       # classB with c1.objB.val3 == 10

    Finally, additional annotations about a parameter can be supplied by using :py:class:`typing.Annotated`
    to wrap the type of the parameter.

    Example::

        class MyClass(HierarchicalParamReader):
            param1: typing.Annotated[int, 'Description of the parameter'] = 42
            param2: typing.Annotated[str, dict(help='Description', metavar='/path/to/file')] = None

        class MySubClass(MyClass):
            param1 = 123

    Note that subclasses that would want to change the default value of the parameter do not need to
    re-annotate it.

    If a string is given as a second parameter to :py:class:`typing.Annotated`, it will be interpreted
    as the description of the parameter. If a dictionary is given, it will be supplied as keyword arguments
    to :py:func:`argparse.ArgumentParser.add_argument` in the :py:module:`stepsblender.load` module to
    automatically add the parameter as a command-line argument.

    Simple types like int and float will be converted from the string provided as a command-line argument
    to the correct value, more complex types might require additional implementation, see for example
    :py:func:`stepsblender.utils.makeCreatableFromStr`.
    """
    parent = None
    nameInParent = None

    def __init__(self, parameters=None, **kwargs):
        super().__init__()
        if not isinstance(parameters, HierarchicalParameters):
            if isinstance(parameters, dict):
                parameters = HierarchicalParameters(parameters)
            else:
                raise ValueError(f'The "parameters" keyword parameter needs to be provided as a dictionary, '
                                 f'got {parameters} instead.')
        self._parameters = parameters.withKwargs(kwargs)
        self._children = []
        # Initialize attributes from class parameters
        self._loadParametersFromCls(HierarchicalParamReader)  # First load parent and nameInParent
        self._loadParametersFromCls(self.__class__)

    @staticmethod
    def getAllAnnotations(cls):
        """Return all annotations from a class and its parents"""
        annotations = {}
        for c in reversed(cls.__mro__):
            if hasattr(c, '__annotations__'):
                annotations = {**annotations, **c.__annotations__}
        return annotations

    @classmethod
    def listParameters(cls):
        annotations = HierarchicalParamReader.getAllAnnotations(cls)
        for name, defValue in inspect.getmembers(cls):
            if not name.startswith('__') and not (inspect.ismethod(defValue) or inspect.isfunction(defValue)
                                                  or inspect.isdatadescriptor(defValue)):
                yield name, defValue, annotations.get(name, None)

    @classmethod
    def listAllParameters(cls):
        res = {}
        for name, defValue, tpe in cls.listParameters():
            if inspect.isclass(defValue) and issubclass(defValue, HierarchicalParamReader):
                res[name] = defValue.listAllParameters()
            elif inspect.isclass(tpe) and issubclass(tpe, HierarchicalParamReader):
                res[name] = tpe.listAllParameters()
            else:
                res[name] = (defValue, tpe)
        return res

    def _loadParametersFromCls(self, cls):
        for name, defValue, tpe in cls.listParameters():
            if name not in self.__dict__:
                setattr(self, name, self._getParam(name, defValue, tpe))

    def _getParam(self, paramName, defValue, cls=None, **kwargs):
        value = self._parameters.get(paramName, defValue, cls, parent=self, nameInParent=paramName, **kwargs)
        if isinstance(value, HierarchicalParamReader):
            self._children.append(value)
        return value

    def _getAllChildren(self, cls=object):
        if isinstance(self, cls):
            yield self
        for child in self._children:
            if child is not self.parent:
                yield from child._getAllChildren(cls)

    @classmethod
    def using(cls, **kwargs):
        """Create a subclass of cls that uses the given keyword arguments as attributes / default value
        """

        class newClass(cls):
            pass

        for name, val in kwargs.items():
            setattr(newClass, name, val)
            if inspect.isclass(val) and issubclass(val, HierarchicalParamReader):
                newClass.__annotations__[name] = val
        return newClass


####################################################################################################


class InterpolationFunction(HierarchicalParamReader):
    n: float = 1

    def __call__(self, x):
        if self.n == 1:
            return x
        elif x <= 0.5:
            return (2 * x)**self.n / 2
        else:
            x -= 0.5
            return 1 - ((2 * (0.5 - x))**self.n) / 2


####################################################################################################


def AddBlenderDataSaving(sim, verbose=True, **kwargs):
    """Add result selectors saving all data that can be visualized in Blender

    :param sim: The simulation to which result selectors should be added
    :type sim: :py:class:`steps.API_2.sim.Simulation`
    :param verbose: Display which result selectors are being added to the simulation.
    :type verbose: bool
    :param \*\*kwargs: Keyword arguments that will be forwarded to
        :py:meth:`steps.API_2.sim.Simulation.toSave` (for example the `dt` argument)

    Usage example::

        AddBlenderDataSaving(sim, dt=1e-3)
    """
    import steps.API_2.sim as nsim
    import steps.API_2.saving as saving
    import steps.API_2.geom as geom
    import steps.API_2.model as model

    rs = saving.ResultSelector(sim)

    selectors = []

    try:
        selectors.append(rs.TETS().ALL(model.Species).Count)
    except nsim.SimPathInvalidPath:
        pass
    try:
        selectors.append(rs.TRIS().ALL(model.Species).Count)
    except nsim.SimPathInvalidPath:
        pass
    for memb in sim.geom.ALL(geom.Membrane):
        try:
            selectors.append(rs.VERTS(memb.tris.verts).V)
        except nsim.SimPathInvalidPath:
            pass
    try:
        selectors.append(rs.ALL(geom.Compartment).VESICLES().Pos)
    except nsim.SimPathInvalidPath:
        pass
    try:
        selectors.append(rs.ALL(geom.Compartment).VESICLES()('surf').POINTSPECS().PosSpherical)
    except nsim.SimPathInvalidPath:
        try:
            selectors.append(rs.ALL(geom.Compartment).VESICLES()('surf').ALL(model.Species).PosSpherical)
        except nsim.SimPathInvalidPath:
            pass
    try:
        selectors.append(rs.ALL(geom.Compartment).VESICLES()('in').ALL(model.Species).Count)
    except nsim.SimPathInvalidPath:
        pass
    try:
        selectors.append(rs.ALL(geom.Compartment).VESICLES()('surf').LINKSPECS().Pos)
    except nsim.SimPathInvalidPath:
        pass
    try:
        selectors.append(rs.ALL(geom.Compartment).VESICLES()('surf').LINKSPECS().LinkedTo)
    except nsim.SimPathInvalidPath:
        pass
    try:
        selectors.append(rs.ALL(geom.Patch).RAFTS().Pos)
    except nsim.SimPathInvalidPath:
        pass
    try:
        selectors.append(rs.ALL(geom.Patch).RAFTS().ALL(model.Species).Count)
    except nsim.SimPathInvalidPath:
        pass
    try:
        selectors.append(rs.ALL(model.Exocytosis, model.RaftEndocytosis).Events)
    except nsim.SimPathInvalidPath:
        pass
    try:
        selectors.append(rs.ALL(geom.Patch).ALL(geom.EndocyticZone).ALL(model.Endocytosis).Events)
    except nsim.SimPathInvalidPath:
        pass

    if verbose and nsim.MPI._shouldWrite:
        print('Result selectors added to the simulation:')
        for sel in selectors:
            print('\t', sel)

    sim.toSave(*selectors, **kwargs)
