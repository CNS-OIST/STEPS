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

import steps

import collections
import inspect
import linecache
import numbers
import re
import warnings

VERBOSITY = 1

__all__ = ['NamedObject', 'Params', 'SetVerbosity']

###################################################################################################
# Utility classes


class SolverPathObject:
    """Base class for all objects susceptible to be part of a SimPath."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return ''

    def _solverId(self):
        """Return the id that is used to identify an object in the solver."""
        return (self.name,)

    def _solverKeywordParams(self):
        """Return the additional keyword parameters that should be passed to the solver"""
        return {}

    def _solverModifier(self):
        """Return None or a function that will modify the output from the solver."""
        return None

    def _solverSetValue(self, v):
        """
        Return the value that should actually be set in the solver when value 'v' is given by
        the user.
        """
        return v

    def _simPathWalkExpand(self):
        """Return an iterable of the elements that should be part of ResultSelector paths."""
        return [self]

    def _simPathCombinerClass(self):
        """Return the class that needs to be used to combine expanded elements."""
        return None

    def _simPathAutoMetaData(self):
        """Return a dictionary with string keys and string or numbers values."""
        return {}

    def __hash__(self):
        return id(self)


class StepsWrapperObject:
    """Base class for all objects that are wrappers for steps objects."""

    def __init__(self, *args, **kwargs):
        super().__init__()

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds. Should be overloaded."""
        return []

    @classmethod
    def _FromStepsObject(cls, obj, *args, **kwargs):
        """Create the interface object from a STEPS object."""
        raise NotImplementedError()


class NamedObject(SolverPathObject):
    """Base class for all objects that are named and can have children

    :param name: Name of the object. If no name is provided, the object gets an automatically
        generated name based on its class.
    :type name: str

    All classes that inherit from :py:class:`NamedObject` can be built with a ``name`` keyword
    parameter. For steps objects, it corresponds to the identifiers used in
    :py:mod:`steps.API_1`. Note that some names are forbidden because they correspond to names
    of attributes or methods of classes defined in this interface. Since most objects implement
    ``__getattr__`` attribute style access to children, the names of these methods / attributes
    could clash with object names. It is thus possible that the contructor of
    :py:class:`NamedObject` raises an exception when trying to name an object with one of these
    forbidden names.

    In addition to a name, this class holds a list of children objects that can be accessed with
    :py:func:`__getattr__` and :py:func:`ALL`.

    .. note::
        This class should not be instantiated by the user, it is only documented for clarity
        since a lot of other classes inherit from it.
    """

    _nameInds = {}  # TODO Not urgent: use threading local?
    _forbiddenNames = set()
    _allowForbNames = False

    def __init__(self, *args, name=None, **kwargs):
        super().__init__(*args, **kwargs)
        if name is None:
            self.name = self.__class__._GetDefaultName()
            self._autoNamed = True
        else:
            self.name = name
            self._autoNamed = False
        if self.name in NamedObject._forbiddenNames:
            if NamedObject._allowForbNames:
                warnings.warn(
                    f"'{self.name}' is a reserved name, SimPath functionnalities might "
                    f"be impacted because of this."
                )
            else:
                raise Exception(f"Cannot call an element '{self.name}', this name is reserved")

        self._children = collections.OrderedDict()
        self._parents = {}

    def _addChildren(self, e):
        if isinstance(e, NamedObject):
            if e.name not in self.children:
                self.children[e.name] = e
            elif e is not self.children[e.name]:
                raise Exception(f'An object is already named {e.name}')
        else:
            raise Exception('Only named objects can be added to UsableObjects.')

    def __getattr__(self, name):
        """Access children of the objects as if they were attributes

        See :py:class:`NamedObject` and :py:func:`Create` for details on object naming.
        Assuming the object has a children named 'child1', one can access it as if it was an
        attribute of the object::

            obj.child1

        :meta public:
        """
        if name.startswith('__'):
            raise AttributeError
        elif name in self.children:
            return self.children[name]
        raise AttributeError(f'{self} does not have an attribute named {name}')

    def _getReferenceObject(self):
        """
        Return the object this object was derived from. Useful for getting the complex associated
        with a complex selector, etc.
        """
        return self

    @property
    def children(self):
        """
        Redirect the children to the reference object.

        :meta private:
        """
        return self._getReferenceObject()._children

    def _getChildrenOfType(self, *cls):
        """Return all children who are instances of cls."""
        for name, c in self.children.items():
            if isinstance(c, cls):
                yield c

    def ALL(self, *cls):
        """Return all children of the object, optionally filtered by class

        Takes a variable number of parameters, if no parameters are given, it returns all children
        of the object. Otherwise, if types are given, it returns the children that match at least
        one of the given types.

        :param \*cls: Variable number of classes

        :returns: A generator that iterates over all children that match the class criteria.
        :rtype: Generator[NamedObject]

        Usage::

            obj.ALL()                    # Return all children
            obj.ALL(Species)             # Return all children Species
            obj.ALL(Reaction, Diffusion) # Return all children that are either Reaction or Diffusion
        """
        return self._getChildrenOfType(*(cls if len(cls) > 0 else [object]))

    @classmethod
    def _GetDefaultName(cls):
        """Return an automatically generated name based on cls."""
        if cls not in NamedObject._nameInds:
            NamedObject._nameInds[cls] = 1
        else:
            NamedObject._nameInds[cls] += 1
        return '{}{}'.format(cls.__name__, NamedObject._nameInds[cls])

    @classmethod
    def Create(cls, *args, **kwargs):
        """Auto naming syntax for simplifying the creation of named objects

        Create one or several objects of class cls with a list of arguments and name them
        automatically based on the name of the variables they are going to be assigned to.

        Usage without arguments::

            a, b, c = Class.Create()

            # Equivalent to:

            a, b, c = Class(name = 'a'), Class(name = 'b'), Class(name = 'c')

        Usage with single arguments::

            a, b, c = Class.Create(argA, argB, argC)

            # Equivalent to:

            a, b, c = Class(argA, name = 'a'), Class(argB, name = 'b'), Class(argC, name = 'c')

        Usage with variable numbers of arguments::

            a, b, c = Class.Create(argA, Params(argB1, argB2), Params(argC1, nargC2 = val))

            # Equivalent to:

            a, b, c = Class(argA, name = 'a'),\\
                      Class(argB1, argB2, name = 'b'),\\
                      Class(argC1, nargC2 = val, name = 'c')

        Usage with global keyword argument::

            a, b, c = Class.Create(argA, argB, argC, gkwarg=val)

            # Equivalent to:

            a, b, c = Class(argA, gkwarg=val, name = 'a'),\\
                      Class(argB, gkwarg=val, name = 'b'),\\
                      Class(argC, gkwarg=val, name = 'c')

        .. warning::
            This automatic naming syntax works by reading the source file to extract the name of
            the variables. Because of this, modifying the source of a script while it is running
            is highly discouraged if the automatic naming syntax is used. Notably, this synatx
            WILL NOT work in an interactive python shell, in which there is no source file to
            parse. It WILL work as expected in a jupyter notebook.
        """
        # Automatic extraction of variable names
        frameInfo1 = inspect.stack()[0]
        frameInfo2 = inspect.stack()[1]
        fname, lineno = frameInfo2.filename, frameInfo2.lineno
        callLine = linecache.getline(fname, lineno).strip(' \t\n\\')
        # Discard lines that correspond only to arguments to Create()
        while lineno > 0 and '.Create(' not in callLine:
            lineno -= 1
            callLine = linecache.getline(fname, lineno)
        lineno -= 1
        currLine = linecache.getline(fname, lineno)
        # Get the full line of variable names, in case line continuation is used
        # We do not need to worry about comments since line continuation does not allow the use of
        # comments before or after it.
        while lineno > 1 and currLine.endswith('\\\n'):
            callLine = currLine.strip(' \t\n\\') + callLine
            lineno -= 1
            currLine = linecache.getline(fname, lineno)
        varNameExpr = '[_a-zA-Z][_a-zA-Z0-9]*'
        p = re.compile(
            rf'^\s*({varNameExpr}(,\s*{varNameExpr})*)\s*=\s*(?:{varNameExpr}\.)*'
            rf'{cls.__name__}\.{frameInfo1.function}\(.*$'
        )
        m = p.match(callLine)
        if m is not None:
            names = [s.strip() for s in m.group(1).split(',')]
            if len(args) == 0:
                if len(names) > 1:
                    res = [cls(name=name, **kwargs) for name in names]
                else:
                    res = cls(name=names[0], **kwargs)
            elif len(names) == 1:
                res = cls(*args, **kwargs, name=names[0])
            elif len(args) == len(names):
                res = []
                for name, arg in zip(names, args):
                    if isinstance(arg, Params):
                        res.append(cls(*arg.args, **arg.kwargs, **kwargs, name=name))
                    else:
                        res.append(cls(arg, **kwargs, name=name))
            else:
                raise Exception(
                    f'The number of arguments ({len(args)}) does not match the number '
                    f'of objects to be created ({len(names)}).'
                )
            return res
        else:
            raise Exception(
                f'The line {callLine} does not match the expected format for automatic assignement.'
            )

    def __repr__(self):
        return self.name

    def _exitCallback(self, parent):
        """Method to be called when we get out of the context manager that used the parent object."""
        pass


class Params:
    """Container class for grouping arguments"""

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs


class UsableObject(NamedObject):
    """Base class for steps objects that can be used as context managers ('with' keyword)."""

    _currUsed = {}  # TODO Not urgent: use threading local?

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.__class__ not in UsableObject._currUsed:
            UsableObject._currUsed[self.__class__] = self._getDefaultCurrUsedVal()

    def _getDefaultCurrUsedVal(self):
        return None

    def __enter__(self):
        if UsableObject._currUsed[self.__class__] is None:
            UsableObject._currUsed[self.__class__] = self
        else:
            raise Exception(f'Cannot use two {self.__class__.__name__} objects simultaneously.')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        UsableObject._currUsed[self.__class__] = None
        if (exc_type, exc_val, exc_tb) == (None, None, None):
            for name, c in self.children.items():
                c._exitCallback(self)


class MultiUsable(UsableObject):
    """Base class for steps objects that can be used several times simultaneously."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _getDefaultCurrUsedVal(self):
        return []

    def __enter__(self):
        UsableObject._currUsed[self.__class__].append(self)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        UsableObject._currUsed[self.__class__].remove(self)
        if (exc_type, exc_val, exc_tb) == (None, None, None):
            for name, c in self.children.items():
                c._exitCallback(self)


class Optional:
    """Wrapper to tag objects as optional."""

    def __init__(self, elem):
        self.elem = elem


def UsingObjects(*usedCls):
    """
    Return a base class for steps objects that depends on other steps objects for their creation.
    Take a list of UsableObject subclasses as argument.
    For example, Species requires a reference to a Model object to be created, so Species will
    inherit from UsingObjects(Model) and can then only be created inside a 'with mdl:' block.
    Usage:
        - StepsObj requires an A:

            Interface code:

                class StepsObj(UsingObjects(A)):
                    def __init__(...):
                        a, = self._getUsedObjects()

            User code:

                a = A()
                with a:
                    s = StepsObj()

        - StepsObj requires an A and either a B or a C:

            Interface code:

                class StepsObj(UsingObjects(A, (B, C))):
                    def __init__(...):
                        a, b, c = self._getUsedObjects()

            User code:

                a = A()
                b = B()
                with a:
                    with b:
                        s = StepsObj() # c will be set to None
    """
    # Check that the used classes inherit from UsableObject
    allCls = []
    for uc in usedCls:
        if isinstance(uc, Optional):
            uc = uc.elem
        if isinstance(uc, tuple):
            allCls += list(uc)
        else:
            allCls.append(uc)
    for uc in allCls:
        if UsableObject not in uc.__mro__:
            raise TypeError(f'{uc} is not a {UsableObject}.')

    class ObjectUser(NamedObject):
        def __init__(self, *args, addAsElement=True, **kwargs):
            super().__init__(*args, **kwargs)
            # Register the object as child of the used objects
            list(self._getUsedObjects(addAsElement=addAsElement, **kwargs))

        def _getObjOfClass(self, cls, addAsElement=True):
            allCls = collections.deque([cls])
            while len(allCls) > 0:
                uc = allCls.popleft()
                allCls.extend(uc.__subclasses__())
                if uc in UsableObject._currUsed and UsableObject._currUsed[uc] is not None:
                    if not isinstance(UsableObject._currUsed[uc], list) and addAsElement:
                        UsableObject._currUsed[uc]._addChildren(self)
                        self._parents[cls] = UsableObject._currUsed[uc]
                        for chld in self._getAdditionalChildren():
                            UsableObject._currUsed[uc]._addChildren(chld)
                            chld._parents[cls] = UsableObject._currUsed[uc]

                    return UsableObject._currUsed[uc]
            return [] if MultiUsable in cls.__mro__ else None

        def _getUsedObjects(self, addAsElement=True, **kwargs):
            """
            Return a generator that gives the used object instances, in the same order as given
            to UsingObjects.
            """
            for uc in usedCls:
                if isinstance(uc, Optional):
                    uc = uc.elem
                    opt = True
                else:
                    opt = False
                if isinstance(uc, tuple):
                    ok = False
                    for uc2 in uc:
                        obj = self._getObjOfClass(uc2, addAsElement)
                        ok |= obj is not None and obj != []
                        yield obj
                    if not ok and not opt:
                        ucNames = [uc2.__name__ for uc2 in uc]
                        clsName = self.__class__.__name__
                        ucNamesStr = ' or a '.join(ucNames)
                        raise Exception(f'Cannot declare a {clsName} out of a {ucNamesStr}.')
                else:
                    obj = self._getObjOfClass(uc, addAsElement)
                    if obj is None and not opt:
                        raise Exception(f'Cannot declare a {self.__class__.__name__} out of a {uc.__name__}.')
                    yield obj

        def _getAdditionalChildren(self):
            """
            Return a list of NamedObjects that should be added as children of the used object
            along with the current object.
            """
            return []

        def _getParentOfType(self, cls):
            return self._parents[cls]

    return ObjectUser


class SimPathCombiner:
    """
    Grouping class for combining the values of subpaths.
    """

    def __init__(self, *paths):
        self.paths = paths

    def __iter__(self):
        return iter(self.paths)

    def func(self, valgen):
        """Combining function, takes an iterable as argument."""
        pass


class SumSimPathComb(SimPathCombiner):
    """A Simpath value combiner that simply sums the values."""

    def __init__(self, *paths):
        super().__init__(*paths)

    def func(self, valgen):
        """Combining function, takes an iterable as argument."""
        return sum(valgen)


class classproperty:
    """Like property but for classes instead of objects"""

    def __init__(self, func):
        self.func = func
        self.__doc__ = func.__doc__

    def __get__(self, obj, objtype):
        return self.func(objtype)


def limitReprLength(func):
    """Decorator for limiting the length of __repr__ methods."""
    REPR_MAX_LENGTH = 200
    ellips = '...'

    def repr(self):
        res = func(self)
        if len(res) > REPR_MAX_LENGTH:
            return res[0 : min(len(res), REPR_MAX_LENGTH - len(ellips))] + ellips
        else:
            return res

    return repr


def formatKey(key, sz, forceSz=False):
    """
    Format a __getitem__ key into a tuple of minimum size sz, transforming '...' into the
    appropriate number of ':'.
    """
    if not isinstance(key, tuple):
        key = (key,)
    # If ellipsis is used, expand it first
    if any(s is Ellipsis for s in key):
        if key.count(Ellipsis) > 1:
            raise KeyError('Cannot use the ellipsis operator ("...") more than once.')
        # If the ellipsis operator was used correctly, fill the missing slots with ':'
        i = key.index(Ellipsis)
        m = max(0, len(key) - sz + 1)
        key = key[:i] + (slice(None),) * m + key[i + 1 :]
    if forceSz:
        key += (slice(None),) * (sz - len(key))
        if len(key) > sz:
            raise Exception(f'Too many dimensions in key indexing.')
    return key


def getSliceIds(s, sz):
    """Return the indices coresponding to slice s of a structure that has length sz."""
    if isinstance(s, numbers.Integral):
        if s >= sz:
            raise IndexError()
        if s >= 0:
            s = slice(s, s + 1, None)
        else:
            s = slice(s % sz, (s % sz) + 1, None)
    if s.start is not None and s.start >= sz:
        raise IndexError()
    return range(*s.indices(sz))


def key2str(key):
    """Return a string describing the key."""
    if not isinstance(key, tuple):
        key = (key,)
    res = []
    for k in key:
        if k is Ellipsis:
            res.append('...')
        elif isinstance(k, slice):
            st = k.start if k.start is not None else ''
            st += f':{k.stop}' if k.stop is not None else ':'
            st += f':{k.step}' if k.step is not None else ''
            res.append(st)
        else:
            res.append(str(k))
    return ', '.join(res)


def args2str(*args, **kwargs):
    """Return a string representation of the arguments."""
    lst = []
    for arg in args:
        if inspect.isclass(arg):
            lst.append(arg.__name__)
        else:
            lst.append(str(arg))
    lst += [f'{name}={val}' for name, val in kwargs.items()]
    return ', '.join(lst)


def SetVerbosity(v):
    """
    Set verbosity to the specified level.

    :param v: Verbosity level
    :type v: int

    The higher the vebosity, the more messages are displayed on standard output.
    """
    global VERBOSITY
    VERBOSITY = v
    # Compatibility with API_1
    steps._suppress_greet = (v == 0)
    steps._quiet = (v == 0)


def _print(msg, prio):
    """
    Print a message if its priority permits it and if we are in rank one in case of an MPI
    simulation.
    """
    from . import sim as nsim
    if prio <= VERBOSITY and nsim.MPI._shouldWrite:
        print(msg)
