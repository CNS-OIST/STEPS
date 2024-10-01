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

import steps

import collections
import copy
from enum import Enum
import functools
import inspect
import itertools
import linecache
import numbers
import numpy
import os
import re
import shutil
import subprocess
import sys
import tempfile
import types
import warnings

VERBOSITY = 1

__all__ = ['Parameter', 'ExportParameters', 'NamedObject', 'Params', 'SetVerbosity']

###################################################################################################
# Parameter and related utility classes


class Units:
    """Represent a physical unit along with a scale

    An object of this class basically consists of a vector of integers of length 7. Each value corresponds
    to a physical dimension in the SI unit system.
    A unit object is also associated with a scale, representing a multiplying factor.
    """

    _SI_PREFIXES = {
        'Y' : (24 , 'yotta-'),
        'Z' : (21 , 'zetta-'),
        'E' : (18 , 'exa-'  ),
        'P' : (15 , 'peta-' ),
        'T' : (12 , 'tera-' ),
        'G' : (9  , 'giga-' ),
        'M' : (6  , 'mega-' ),
        'k' : (3  , 'kilo-' ),
        'h' : (2  , 'hecto-'),
        'da': (1  , 'deca-' ),
        ''  : (0  , ''      ),
        'd' : (-1 , 'deci-' ),
        'c' : (-2 , 'centi-'),
        'm' : (-3 , 'milli-'),
        'u' : (-6 , 'micro-'),
        'n' : (-9 , 'nano-' ),
        'p' : (-12, 'pico-' ),
        'f' : (-15, 'femto-'),
        'a' : (-18, 'atto-' ),
        'z' : (-21, 'zepto-'),
        'y' : (-24, 'yocto-'),
    }
    _SI_UNITS = {
        # Base units
        's': (
            numpy.array([1, 0, 0, 0, 0, 0, 0], numpy.intc), # Exponent
            0, # Scale in the form 10^x so 0 means 1, 1 means 10, etc
            'second',
            'time',
        ),
        'm': (
            numpy.array([0, 1, 0, 0, 0, 0, 0], numpy.intc),
            0,
            'meter',
            'length',
        ),
        'g': (
            numpy.array([0, 0,  1, 0, 0, 0, 0], numpy.intc),
            -3, # Should be kg but 'k' is already part of the prefixes
            'gram',
            'mass',
        ),
        'A': (
            numpy.array([0, 0, 0, 1, 0, 0, 0], numpy.intc),
            0,
            'ampere',
            'electric current',
        ),
        'K': (
            numpy.array([0, 0, 0, 0, 1, 0, 0], numpy.intc),
            0,
            'kelvin',
            'thermodynamic temperature',
        ),
        'mol': (
            numpy.array([0, 0, 0, 0, 0, 1, 0], numpy.intc),
            0,
            'mole',
            'amount of substance',
        ),
        'cd': (
            numpy.array([0, 0, 0, 0, 0, 0, 1], numpy.intc),
            0,
            'candela',
            'luminous intensity',
        ),
        # Derived units
        'V': (
            numpy.array([-3, 2, 1, -1, 0, 0, 0], numpy.intc),
            0,
            'volt',
            'electrical potential difference',
        ),
        'F': (
            numpy.array([4, -2, -1, 2, 0, 0, 0], numpy.intc),
            0,
            'farad',
            'capacitance',
        ),
        'ohm': (
            numpy.array([-3, 2, 1, -2, 0, 0, 0], numpy.intc),
            0,
            'ohm',
            'electrical resistance',
        ),
        'S': (
            numpy.array([3, -2, -1, 2, 0, 0, 0], numpy.intc),
            0,
            'siemens',
            'electrical conductance'
        ),
        'L': (
            numpy.array([0, 3, 0, 0, 0, 0, 0], numpy.intc),
            -3,
            'litre',
            'volume',
        ),
        'M': (
            numpy.array([0, -3, 0, 0, 0, 1, 0], numpy.intc),
            3,
            'molar',
            r'concentration (mol L\ :sup:`-1`)',
        ),
        'C': (
            numpy.array([1, 0, 0, 1, 0, 0, 0], numpy.intc),
            0,
            'coulomb',
            'electric charge'
        ),
        'J': (
            numpy.array([-2, 2, 1, 0, 0, 0, 0], numpy.intc),
            0,
            'joule',
            'energy, work, heat'
        )
    }

    _SPACER_RE = re.compile(r'\s+')
    _SCALED_UNIT_RE = re.compile(
        r'({pref})({dim})([\s\)\^]|$)'.format(
            pref='|'.join(sorted(_SI_PREFIXES.keys(), key=lambda x:-len(x))),
            dim='|'.join(sorted(_SI_UNITS.keys(), key=lambda x:-len(x))),
        )
    )
    _SCALED_GROUP_RE = re.compile(
        r'({pref})\('.format(
            pref='|'.join(sorted(_SI_PREFIXES.keys(), key=lambda x:-len(x))),
        )
    )
    _EXPONENT_RE = re.compile(r'\^(-?\d+)')

    _SPACE_SIMPLIF_1 = re.compile(r'\s+')
    _SPACE_SIMPLIF_2 = re.compile(r'\(\s+')
    _SPACE_SIMPLIF_3 = re.compile(r'\s+\)')

    __slots__ = ['_str', '_scale', '_exponents']

    @staticmethod
    def _parseUnitString(s):
        exponents = numpy.array([0] * 7, numpy.intc)
        scale = 0

        i = 0
        while i < len(s):
            # Match one of the three tokens
            m_spacer = Units._SPACER_RE.match(s[i:])
            m_unit = Units._SCALED_UNIT_RE.match(s[i:])
            m_group = Units._SCALED_GROUP_RE.match(s[i:])
            if m_spacer is not None:
                # Spaces
                i += m_spacer.end()
                continue
            elif m_unit is not None:
                # SI Unit with prefix
                pref, unit, _ = m_unit.groups()

                tmp_exp, tmp_scale = Units._SI_UNITS[unit][:2]
                tmp_scale += Units._SI_PREFIXES[pref][0]

                i += m_unit.end(2)
            elif m_group is not None:
                # Parentheses group with prefix
                cnt = 1
                for j in range(i + m_group.end(), len(s)):
                    cnt += 1 if s[j] == '(' else (-1 if s[j] == ')' else 0)
                    if cnt == 0:
                        break
                if cnt > 0:
                    raise Exception(f'Unmatched parenthesis in {s}')

                tmp_exp, tmp_scale = Units._parseUnitString(s[i + m_group.end():j])
                tmp_scale += Units._SI_PREFIXES[m_group.group(1)][0]

                i = j+1
            else:
                raise Exception(f'Could not parse "{s[i:]}" in "{s}".')

            # Check for exponents
            expo = 1
            if i < len(s):
                m_expo = Units._EXPONENT_RE.match(s[i:])
                if m_expo is not None:
                    expo = int(m_expo.group(1))

                    i+= m_expo.end()
            
            # Add to current values
            exponents += tmp_exp * expo
            scale += tmp_scale * expo

        return exponents, scale

    def __init__(self, units):
        """Create the Units from a formatted string."""

        if not isinstance(units, str):
            raise TypeError(f'Expected a string, got {units} instead.')

        self._str = units.strip()
        if len(self._str) > 0:
            # Remove extra spaces
            self._str = Units._SPACE_SIMPLIF_1.sub(' ', self._str)
            self._str = Units._SPACE_SIMPLIF_2.sub('(', self._str)
            self._str = Units._SPACE_SIMPLIF_3.sub(')', self._str)

            self._exponents, self._scale = Units._parseUnitString(self._str)
        else:
            self._exponents = numpy.array([0] * 7, numpy.intc)
            self._scale = 0

    def _multiplyWith(self, other):
        """Combine self with an other Units object by multiplication"""
        res = Units('')
        res._exponents += self._exponents
        res._exponents += other._exponents
        res._scale = self._scale + other._scale
        res._str = f'{self._str} {other._str}'.strip()
        return res

    def _raiseToPower(self, other):
        """Raise self to an integer power"""
        res = Units('')
        res._exponents = self._exponents * other
        res._scale = self._scale * other
        res._str = Units._strRaiseToPower(self._str, other)
        return res

    def _divideWith(self, other):
        """Combine self with an other Units object by division"""
        res = Units('')
        res._exponents += self._exponents
        res._exponents -= other._exponents
        res._scale = self._scale - other._scale
        res._str = f'{self._str} {Units._strRaiseToPower(other._str, -1)}'.strip()
        return res

    @staticmethod
    def _strRaiseToPower(sstr, power):
        if len(sstr) > 0:
            if ' ' in sstr or '(' in sstr:
                return f'({sstr})^{power}'
            else:
                if '^' in sstr:
                    powind = sstr.index('^')
                    newpow = int(sstr[powind+1:]) * power
                    return f'{sstr[:powind]}^{newpow}'
                else:
                    return f'{sstr}^{power}'
        else:
            return ''

    def __eq__(self, other):
        return (
            isinstance(other, Units) and
            (self._exponents == other._exponents).all() and
            self._scale == other._scale
        )

    def __hash__(self):
        return hash((self._str, self._scale, tuple(self._exponents)))

    def __repr__(self):
        return self._str

    def _compatibleWith(self, other):
        if not isinstance(other, Units):
            raise TypeError(f'Expected a Units object, got {other} instead.')
        return (self._exponents == other._exponents).all()

    def _isDimensionless(self):
        return (self._exponents == 0).all() and self._scale == 0

    def _toUnicode(self):
        """Return a unicode string representing the unit."""
        _SUPER_CHAR_MAP = {str(i): c for i, c in enumerate('⁰¹²³⁴⁵⁶⁷⁸⁹')}
        _SUPER_CHAR_MAP['-'] = '⁻'

        def prefix(val):
            if val == 'u':
                return 'μ'
            else:
                return val

        res = Units._SCALED_UNIT_RE.sub(
            lambda m: prefix(m.group(1)) + m.group(2) + m.group(3),
            self._str,
        )
        res = Units._EXPONENT_RE.sub(
            lambda m: ''.join(_SUPER_CHAR_MAP[c] for c in m.group(0)[1:]),
            res,
        )
        return res

    def _toLatex(self):
        """Return a LaTeX formated string representing the unit.
        Uses the siunitx latex package"""

        def prefix(val):
            if val == 'u':
                return r'\micro '
            else:
                return val

        res = Units._EXPONENT_RE.sub(
            lambda m: '^{' + m.group(1) + '}',
            self._str,
        )
        res = Units._SCALED_UNIT_RE.sub(
            lambda m: prefix(m.group(1)) + m.group(2) + m.group(3),
            res,
        )
        res = r'\si{' + res + '}'
        res = res.replace(' ', '.')
        res = res.replace(r'\micro.', r'\micro ')
        return res


class ParameterizedObject:
    """Base class for all objects holding Parameter objects

    Classes that inherit from ParameterizedObject can declare parameters that will then be used to
    construct parameter tables. The simplest way to declare parameters is to decorate properties getter
    and setter with RegisterGetter and RegisterSetter respectively (see examples in geom.Compartment
    e.g.). This allows custom code to be executed during the calls to getter and setters.

    It is also possible to declare parameters that are not associated to custom code by using
    cls._registerParameter(...). It will simply add the corresponding property along with its getter
    and setter.

    Note that if the setter is never called, the parameter is never added to self._parameters and will
    thus not be part of parameter tables. This is the desired behavior as we do not want to consider that
    e.g. 'Vol' is a parameter for tetrahedral compartments. Since it is never set, it will not be treated as
    a parameter but simply as a property.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Simple static parameters
        self._parameters = {}

    def _getAllParams(self):
        """Return all parameters"""
        return list(self._parameters.values())

    def _getSubParameterizedObjects(self):
        """Return all subobjects that can hold parameters."""
        return []

    def _getParameter(self, key):
        """Get a parameter"""
        return self._parameters[key] if key in self._parameters else None

    def _setParameter(self, key, param, units=None):
        """Set a parameter"""
        if not isinstance(param, Parameter):
            param = Parameter(param, units, name='')
        self._parameters[key] = param

    def _includeinParamTables(self):
        """Whether the object should be included in parameter tables"""
        return True

    @classmethod
    def _getDisplayName(cls):
        """Name that will be used as column title during parameter export
        Defaults to class name.
        """
        return cls.__name__

    @staticmethod
    def RegisterGetter(units=Units('')):
        """
        Decorator for wrapping the getter of a property so that the correct parameter value is
        returned. Since parameters are kept on the python side, this means that the cython bindings
        do not need to implement the corresponding getters.

        :meta private:
        """
        def wrapper(func, name=None):
            def getter(self):
                propName = func.__name__ if name is None else name
                param = self._getParameter(propName)
                if param is not None:
                    # Get expected units
                    # It needs to be done here because the unit can depend on the object.
                    if hasattr(units, '__call__'):
                        _units = units(self)
                    else:
                        _units = units

                    if _units is None:
                        return param.value
                    else:
                        return param.valueIn(_units)
                else:
                    # If the property was not set, call the getter. This means that properties
                    # that were not explicitely set will not be listed as parameters in parameter
                    # tables.
                    return func(self)
            getter.__doc__ = func.__doc__
            getter._units = units
            return getter
        return wrapper

    @staticmethod
    def RegisterSetter(units=Units('')):
        """
        Decorator for wrapping the setter of a property so that parameter objects are
        automatically added to self._parameters dict.

        :meta private:
        """
        def wrapper(func, name=None):
            def setter(self, val):
                # Get expected units
                # It needs to be done here because the unit can depend on the object.
                if hasattr(units, '__call__'):
                    _units = units(self)
                else:
                    _units = units

                # Handle parameter
                if not isinstance(val, Parameter):
                    val = Parameter(val, name='', units=_units)
                elif val._units is None:
                    val._units = _units

                if val._units is not None and not val._units._compatibleWith(_units):
                    raise Exception(
                        f'Expected a value in "{_units}", got a value in "{val._units}" instead'
                    )

                # Call the actual setter (if RegisterSetter was used as decorator)
                if func is not None:
                    func(self, val.valueIn(_units))

                # Only add the property if the setter did not raise an exception
                propName = func.__name__ if name is None else name
                self._setParameter(propName, val)
            setter.__doc__ = func.__doc__
            setter._units = units
            return setter
        return wrapper

    @classmethod
    def RegisterParameter(cls, name, units, defVal=None, addSetter=True):
        """
        Register a simple parameter that does not require custom code execution for setting or getting.
        If custom code execution is needed, RegisterGetter and RegisterSetter should instead be used as
        decorators of properties.
        Optionally supply a default value that will be returned if the parameter is not set.

        :meta private:
        """
        fget = ParameterizedObject.RegisterGetter(units=units)(lambda x: defVal, name=name)
        fset = ParameterizedObject.RegisterSetter(units=units)(None, name=name) if addSetter else None

        # Add the parameter as a property
        setattr(cls, name, property(fget, fset))

    @staticmethod
    def SpecifyUnits(*units, **kwunits):
        """
        Wrapper decorator for specifying units of parameters.

        :meta private:
        """
        def wrapper(func):
            def newFunc(*args, **kwargs):
                newArgs = []
                for arg, unit in itertools.zip_longest(args, units):
                    if isinstance(arg, Parameter):
                        if unit is not None:
                            newArgs.append(arg.valueIn(unit))
                        elif arg._units is None or arg._units._isDimensionless():
                            newArgs.append(arg._value)
                        else:
                            raise Exception(
                                f'Expected a non-dimensional value, got a value in {arg.units} '
                                f'instead.'
                            )
                    else:
                        newArgs.append(arg)

                newKwargs = {}
                for key, arg in kwargs.items():
                    unit = kwunits.get(key, None)
                    if isinstance(arg, Parameter):
                        if unit is not None:
                            newKwargs[key] = arg.valueIn(units)
                        elif arg._units is None or arg._units._isDimensionless():
                            newKwargs[key] = arg._value
                        else:
                            raise Exception(
                                f'Expected a non-dimensional value, got a value in {arg.units} '
                                f'instead.'
                            )
                    else:
                        newKwargs[key] = arg

                return func(*newArgs, **newKwargs)

            return newFunc
        return wrapper


class AdvancedParameterizedObject(ParameterizedObject):
    """Base class for parameterized objects with more complex key for parameters

    This is usefull for objects that have parameter values that can change several times during
    a single run, e.g. simulation path values.
    
    Parameters held in this class have a description that is more complex than just a string name
    The key in the self._parameters dict should be a tuple with the following format:
    (('prop1', val1), ('prop2', 'val2'), ...)
    """
    _ADV_PARAMS_DEFAULT_NAME = 'Name'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _getParameter(self, name):
        """Get a parameter"""
        return super()._getParameter(self._getCurrentParamKey() + ((self._ADV_PARAMS_DEFAULT_NAME, name),))

    def _setParameter(self, name, param, units=None):
        """Set a parameter"""
        key = self._getCurrentParamKey()
        super()._setParameter(key + ((self._ADV_PARAMS_DEFAULT_NAME, name),), param, units=units)

    def _getCurrentParamKey(self):
        """Return a key representing the current state of the object"""
        raise NotImplementedError()

    @classmethod
    def _getGroupedAdvParams(cls):
        """Return a list of advanced parameter property that should be grouped
        """
        return []

class Parameter:
    r"""Class to describe a parameter in a STEPS script

    This class allows the user to declare values as parameters of the script. Parameter objects
    can have units (see example below) and are used to create STEPS objects or to set some of
    their properties. Declaring Parameter objects makes it possible for the user to do automatic
    unit conversions and to give a name to a parameter that is used with several STEPS objects.
    When parameter tables are automatically generated using :py:func:`ExportParameters`, the
    informations that the user supplied to the ``Parameter`` objects will be displayed in the tables.

    :param value: The value of the parameter
    :type value: Any but will mostly be used with float
    :param units: A string describing the unit of the parameter (see below for more explanations).
        If it is not supplied, the parameter will infer its units when used with a STEPS object.
    :type units: str
    :param name: A user-supplied name for the parameter. If it is not given, the Parameter object
        will try to take its name from the name of the variable it is assigned to.
    :type name: str
    :param \*\*kwargs: Keyword arguments that will be displayed in parameter tables.

    Usage::

        k_on  = Parameter(2, 'uM^-1 s^-1')
        k_off = Parameter(1, 's^-1', comment='Comment about the k_off parameter')

        ...

        with mdl, vsys:
            SA + SB <r[1]> SC
            r[1].K = k_on, k_off

    In the above example, ``k_on`` will internally be converted to the appropriate STEPS unit (M^-1 s^-1)
    and both ``k_on`` and ``k_off`` will appear in the parameter tables generated from
    :py:func:`ExportParameters` with their name and associated units. In addition, any other keyword
    arguments (``comment=...`` here) will also be shown in the parameter tables.

    When given to a STEPS object, if a Parameter is declared with a unit, a test will be performed to
    check whether the given unit is of the expected physical dimension. If not, an exception will be
    raised. It is thus good practise to specify the units of parameters as a sanity check.

    Note that Parameters can be combined using standard arithmetic operators ``+``, ``-``, ``*``, ``/``,
    ``**`` (see :py:func:`Parameter.__add__`, etc.) to create new Parameters. These operations will then
    be visible in the exported parameter tables (see :py:func:`ExportParameters`).

    .. warning::
        Converting between units should be done by calling the :py:func:`convertTo` method. For example,
        to convert a Parameter in V to mV, one should use::

            pot = Parameter(-0.07, 'V')
            potInmV = pot.convertTo('mV')

        One should never use arithmetic operations like ``pot * 1e3`` which would yield a parameter of
        -70 V instead.

    Finally, note that when a Parameter is explicitely converted to a number, it will yield its value in SI
    units::

        >>> pot = Parameter(-70, 'mV')
        >>> float(pot)
        -0.07

    Units specification:

        The string that specifies the units should be composed of physical units (see available units below)
        optionally raised to some integer power and optionally prefixed with a scale (``'u'`` for micro,
        ``'k'`` for kilo, etc. see list below).

        Examples::
            
            'mm^2'               # Square millimeters
            'uM^-1 s^-1'         # Per micromolar per second
            '(mol m^-2)^-1 s^-1' # Per (moles per square meters) per second
            'mV'                 # Millivolts

        Note that groups of units can be wrapped in parentheses and raised to a power. This enhances
        the readability of some units. Reaction rates for surface-surface reactions of order 2 thus have
        units of ``'(mol m^-2)^-1 s^-1'`` which is equivalent to ``'mol^-1 m^2 s^-1'`` but the former
        makes it clear that a surface 'concentration' is involved.

        Available prefixes:

        {prefixTable}

        Available units:

        {unitTable}


    """

    _NAME_EXTRACT_RE = re.compile(r'(?:(?:^\s*(\w+)\s*=\s*)|.*)Parameter(\(.*$)')

    __slots__ = ['_defLoc', '_name', '_fullname', '_value', '_units', '_kwargs', '_composedPrio', '_dependencies']

    _PARAMETER_USAGE_RECORD = None

    def __init__(self, value, units=None, name=None, _composedPrio=0, _dependencies=[], **kwargs):
        self._defLoc = None
        self._fullname = None

        if name is None:
            # Call from user, try to infer name from variable name
            frameInfo = inspect.stack()[1]
            fname, lineno = frameInfo.filename, frameInfo.lineno
            fullLine = ''
            line = linecache.getline(fname, lineno)
            match = None
            # Discard lines that correspond only to arguments to Parameter()
            while lineno >= 0 and (line.endswith('\\\n') or match is None):
                fullLine = line.rstrip('\\\n') + fullLine
                match = Parameter._NAME_EXTRACT_RE.match(fullLine)
                lineno -= 1
                line = linecache.getline(fname, lineno)
            if match is not None and match.group(1) is not None:
                # Check that there is only one parameter object on the rhs
                pcount = 0
                ind = 0
                parStr = match.group(2)
                for i, c in enumerate(parStr):
                    if c == '(':
                        pcount += 1
                    elif c == ')':
                        pcount -= 1
                    if pcount == 0:
                        ind = i
                        break

                if re.match(r'^(\s*[;#].*|\s*)$', parStr[ind+1:]) is not None:
                    # If there is indeed only one parameter on the rhs
                    self._name = match.group(1)
                    self._defLoc = (fname, lineno + 1)
                else:
                    self._name = None
            else:
                # If it fails, keep None
                self._name = None
        elif isinstance(name, str):
            self._name = name if len(name) > 0 else None
        else:
            raise TypeError(f'Expected a string, got {name} instead')

        if isinstance(units, str):
            units = Units(units)
        elif units is not None and not isinstance(units, Units):
            raise TypeError(f'Expected a string or None, got {units} instead.')

        if isinstance(value, Parameter):
            # If we are creating a new parameter from some combination of parameters
            self._fullname = value._name
            _dependencies = _dependencies + [value]
            if units is not None:
                value = value.valueIn(units)
            else:
                units = value._units
                value = value._value

        self._value = value
        self._units = units

        self._kwargs = {}
        for key in kwargs.keys():
            if hasattr(self, key):
                raise Exception(f'Cannot use a keyword parameter named {key}.')
        self._kwargs = kwargs

        self._composedPrio = _composedPrio
        self._dependencies = _dependencies

    @property
    def name(self):
        """Name of the parameter

        :type: str, read-only
        """
        return self._name

    @property
    def value(self):
        """Value of the parameter (in its units)

        :type: Any (usually float), read-only
        """
        return self._value

    @property
    def units(self):
        """Units of the parameter

        :type: Union[str, None], read-only
        """
        return str(self._units) if self._units is not None else None

    def convertTo(self, unit):
        """Convert a parameter to the specified unit

        :param unit: The unit to be converted to, it needs to be compatible with the units
            of the parameter.
        :type unit: str

        :returns: The converted parameter.
        :rtype: :py:class:`Parameter`
        """
        return self._convertTo(Units(unit))

    def valueIn(self, unit):
        """Return the value of the parameter in the specified unit.

        If the parameter is itself a parameterized object, do not attempt any conversion.
        :meta private:
        """
        if self.value is None:
            return None
        elif isinstance(self.value, ParameterizedObject):
            return self.value
        elif self._units is not None:
            if not self._units._compatibleWith(unit):
                raise Exception(f'Expected "{unit}", got "{self._units}" instead.')
            multiplier = 10**int(self._units._scale - unit._scale)
            # Avoid extending lists and tuples, apply the multiplication to each element
            if isinstance(self.value, (list, tuple)):
                return self.value.__class__(v * multiplier for v in self.value)
            else:
                return self.value * multiplier
        elif unit is None:
            return self.value
        else:
            raise Exception(f'Expected a nondimensional unit, got {unit} instead.')

    def _valueInSI(self):
        """Return the value of the parameter in SI unit system

        Beware that STEPS uses M (molar) for concentration while SI uses mol m^-3
        """
        if self._units is not None:
            return self._value * 10**self._units._scale
        return self._value

    def _convertTo(self, unit):
        """Return a new parameter with the specified unit
        Take a Unit object instead of a string.
        """
        return Parameter(
            self.valueIn(unit),
            unit,
            name=self._name if self._isNamed() else '',
            _composedPrio=self._composedPrio,
            _dependencies=self._dependencies,
        )

    def _isNamed(self):
        return self._name is not None

    def _isUserDefined(self):
        return self._composedPrio == 0

    def _checkUsableInOp(self):
        """Return whether self can be used in an arithmetic operation"""
        if self._units is None:
            raise Exception(f'Cannot use a parameter with undefined units in arithmetic operations.')

        # Record parameter usage, useful for knowing which parameters are used in e.g. VDepRate.
        if Parameter._PARAMETER_USAGE_RECORD is not None:
            lst, depth = Parameter._PARAMETER_USAGE_RECORD
            lst[depth].append(self)

    def _getAllDependencies(self, alreadyUsed=None):
        """Recursively return all Parameters dependencies"""
        if alreadyUsed is None:
            alreadyUsed = set()

        for param in self._dependencies:
            if param not in alreadyUsed:
                alreadyUsed.add(param)
                yield param
                yield from param._getAllDependencies(alreadyUsed=alreadyUsed)

    def __getattr__(self, name):
        """Access properties of the parameter as if they were attributes

        Keywords arguments passed to the constructor of Parameter can then be accessed as if they were
        attributes::

            >>> k_on = Parameter(2, 'uM^-1 s^-1', comment='Comment about k_on')
            >>> k_on.comment
            'Comment about k_on'

        Note that the properties can only be defined at object creation and are all read-only.

        :meta public:
        """
        if not name.startswith('__') and name in self._kwargs:
            return self._kwargs[name]
        else:
            raise AttributeError(f'Attribute "{name}" was not declared for the Parameter.')

    def __eq__(self, other):
        return (
            isinstance(other, Parameter) and
            (self._name, self._value, self._units, tuple(self._kwargs.items())) == 
            (other._name, other._value, other._units, tuple(other._kwargs.items()))
        )

    def __hash__(self):
        return hash((self._name, self._value, self._units, tuple(self._kwargs.items())))

    def __add__(self, other, _inv=False):
        """Add Parameters with the ``+`` operator

        :param other: The other Parameter or a number
        :type other: Union[:py:class:`Parameter`, float]
        :returns: The Parameter resulting from the addition of both operands. If both
            operands are Parameters, the units of the result matches the units of the leftmost
            named parameter. If one operand is a number, it implicitely takes the units of the other
            operand. The name of the returned Parameter describes the operation.
        :rtype: :py:class:`Parameter`

        Usage::

            >>> param1 = Parameter(1, 'uM s^-1')
            >>> param2 = Parameter(2, 'nM ms^-1')
            >>> param3 = param1 + param2
            >>> param3.units
            'uM s^-1'
            >>> param3.value
            3
            >>> param3.name
            'param1 + param2'

        :meta public:
        """
        if isinstance(other, numbers.Number):
            other = Parameter(other, self._units, name='')
        elif not isinstance(other, Parameter):
            raise TypeError(f'Expected a Parameter or a number, got {other} instead.')

        self._checkUsableInOp()
        other._checkUsableInOp()

        deps = [obj for obj in [self, other] if obj._isNamed()]

        if not self._isNamed() and other._isNamed():
            units = other._units
            self = self._convertTo(units)
        else:
            units = self._units
            other = other._convertTo(units)

        val = self._value + other._value

        if not self._isNamed() and not other._isNamed():
            name = ''
        else:
            name = f'{other} + {self}' if _inv else f'{self} + {other}'

        return Parameter(val, units, name=name, _composedPrio=3, _dependencies=deps)

    def __sub__(self, other, _inv=False):
        """Subtract Parameters with the ``-`` operator

        :param other: The other Parameter or a number
        :type other: Union[:py:class:`Parameter`, float]
        :returns: The Parameter resulting from the subtraction of both operands. If both
            operands are Parameters, the units of the result matches the units of the leftmost
            named parameter. If one operand is a number, it implicitely takes the units of the other
            operand. The name of the returned Parameter describes the operation.
        :rtype: :py:class:`Parameter`

        Usage::

            >>> param1 = Parameter(1, 'uM')
            >>> param2 = Parameter(100, 'nM')
            >>> param3 = param1 - param2
            >>> param3.units
            'uM'
            >>> param3.value
            0.9
            >>> param3.name
            'param1 - param2'

        :meta public:
        """
        if isinstance(other, numbers.Number):
            other = Parameter(other, self._units, name='')
        elif not isinstance(other, Parameter):
            raise TypeError(f'Expected a Parameter or a number, got {other} instead.')

        self._checkUsableInOp()
        other._checkUsableInOp()

        deps = [obj for obj in [self, other] if obj._isNamed()]

        if not self._isNamed() and other._isNamed():
            units = other._units
            self = self._convertTo(units)
        else:
            units = self._units
            other = other._convertTo(units)

        val = (other._value - self._value) if _inv else (self._value - other._value)

        if not self._isNamed() and not other._isNamed():
            name = ''
        else:
            selfStr = f'({self})' if self._composedPrio >= 3 else str(self)
            otherStr = f'({other})' if other._composedPrio >= 3 else str(other)
            name = f'{other} - {selfStr}' if _inv else f'{self} - {otherStr}'

        return Parameter(val, units, name=name, _composedPrio=3, _dependencies=deps)

    def __mul__(self, other, _inv=False):
        """Multiply Parameters with the ``*`` operator

        :param other: The other Parameter or a number
        :type other: Union[:py:class:`Parameter`, float]
        :returns: The Parameter resulting from the multiplication of both operands. If both
            operands are Parameters, the units of the result consists in the product of the units
            of the operands. If one operand is a number, it is implicitely treated as dimensionless.
            The name of the returned Parameter describes the operation.
        :rtype: :py:class:`Parameter`

        Usage::

            >>> param1 = Parameter(1, 'uM')
            >>> param2 = Parameter(5, 's^-1')
            >>> param3 = param1 * param2
            >>> param3.units
            'uM s^-1'
            >>> param3.value
            5
            >>> param3.name
            'param1 * param2'

        :meta public:
        """
        if isinstance(other, numbers.Number):
            other = Parameter(other, '', name='')
        elif not isinstance(other, Parameter):
            raise TypeError(f'Expected a Parameter or a number, got {other} instead.')

        self._checkUsableInOp()
        other._checkUsableInOp()

        deps = [obj for obj in [self, other] if obj._isNamed()]

        if _inv:
            units = other._units._multiplyWith(self._units)
        else:
            units = self._units._multiplyWith(other._units)

        val = self._value * other.value

        if not self._isNamed() and not other._isNamed():
            name = ''
        else:
            selfStr = f'({self})' if self._composedPrio > 2 else str(self)
            if not self._isNamed() and self._units is not None and len(self.units) > 0:
                selfStr = f'({selfStr} {self.units})'
            otherStr = f'({other})' if other._composedPrio > 2 else str(other)
            if not other._isNamed() and other._units is not None and len(other.units) > 0:
                otherStr = f'({otherStr} {other.units})'
            name = f'{otherStr} * {selfStr}' if _inv else f'{selfStr} * {otherStr}'

        return Parameter(val, units, name=name, _composedPrio=2, _dependencies=deps)

    def __truediv__(self, other, _inv=False):
        """Divide Parameters with the ``/`` operator

        :param other: The other Parameter or a number
        :type other: Union[:py:class:`Parameter`, float]
        :returns: The Parameter resulting from the division of both operands. If both
            operands are Parameters, the units of the result consists in the division of the units
            of the operands. If one operand is a number, it is implicitely treated as dimensionless.
            The name of the returned Parameter describes the operation.
        :rtype: :py:class:`Parameter`

        Usage::

            >>> param1 = Parameter(1, 'uM')
            >>> param2 = Parameter(5, 's')
            >>> param3 = param1 / param2
            >>> param3.units
            'uM s^-1'
            >>> param3.value
            0.2
            >>> param3.name
            'param1 / param2'

        :meta public:
        """
        if isinstance(other, numbers.Number):
            other = Parameter(other, '', name='')
        elif not isinstance(other, Parameter):
            raise TypeError(f'Expected a Parameter or a number, got {other} instead.')

        self._checkUsableInOp()
        other._checkUsableInOp()

        deps = [obj for obj in [self, other] if obj._isNamed()]

        if _inv:
            units = other._units._divideWith(self._units)
            val = other._value / self.value
        else:
            units = self._units._divideWith(other._units)
            val = self._value / other.value

        if not self._isNamed() and not other._isNamed():
            name = ''
        else:
            selfStr = f'({self})' if self._composedPrio > 2 else str(self)
            if not self._isNamed() and self._units is not None and len(self.units) > 0:
                selfStr = f'({selfStr} {self.units})'
            otherStr = f'({other})' if other._composedPrio > 2 else str(other)
            if not other._isNamed() and other._units is not None and len(other.units) > 0:
                otherStr = f'({otherStr} {other.units})'
            name = f'{otherStr} / {selfStr}' if _inv else f'{selfStr} / {otherStr}'

        return Parameter(val, units, name=name, _composedPrio=2, _dependencies=deps)

    def __pow__(self, power):
        """Raise to an integer power with the ``**`` operator

        :param power: The integer power
        :type power: int
        :returns: The Parameter resulting from raising self to the given power. The units of the
            result are the units of self raised to power. The name of the returned Parameter
            describes the operation.
        :rtype: :py:class:`Parameter`

        Usage::

            >>> param1 = Parameter(2, 'um')
            >>> param2 = param1 ** 3
            >>> param2.units
            'um^3'
            >>> param2.value
            8
            >>> param2.name
            'param1 ** 3'

        :meta public:
        """
        if not isinstance(power, int):
            raise TypeError(f'Expected an integer, got {power} instead.')

        self._checkUsableInOp()

        deps = [self] if self._isNamed() else []

        units = self._units._raiseToPower(power)
        val = self._value ** power

        if self._isNamed():
            name = f'({self}) ** {power}' if self._composedPrio > 1 else f'{self} ** {power}'
        else:
            name = ''

        return Parameter(val, units, name=name, _composedPrio=1, _dependencies=deps)

    def __radd__(self, other):
        return self.__add__(other, _inv=True)

    def __rsub__(self, other):
        return self.__sub__(other, _inv=True)

    def __rmul__(self, other):
        return self.__mul__(other, _inv=True)

    def __rtruediv__(self, other):
        return self.__truediv__(other, _inv=True)

    def __round__(self, ndigits=None):
        """Round parameter value in SI

        :param ndigits: Optional, number of digits after the decimal point, defaults to None (i.e.
            rounds to nearest integer).
        :type ndigits: int
        :returns: The Parameter resulting from the rounding to ndigits digits. When rounding a parameter,
            its value is first converted to SI and the SI value is rounded. The returned parameter keeps
            the same units. The name of the returned Parameter describes the operation.
        :rtype: :py:class:`Parameter`

        Usage::

            >>> param1 = Parameter(12.7, 'dm')
            >>> param2 = round(param1)
            >>> param2.units # The units are unchanged
            'dm'
            >>> param2.value # 12.7 dm is first converted to 1.27 m and rounded to 1 m = 10 dm
            10
            >>> param2.name
            'round(param1)'
            >>> round(param1, 1).value # Keep more digits: 1.27 m rounded to 1.3 m = 13 dm
            13.0

        :meta public:
        """
        self._checkUsableInOp()

        deps = [self] if self._isNamed() else []

        units = self._units

        if self._units is not None:
            scl = self._units._scale
            val = round(self._value * 10 ** scl, ndigits) * 10 ** -scl
        else:
            val = round(self._value, ndigits)

        if self._isNamed():
            name = f'round({self}' + (f', {ndigits}' if ndigits is not None else '') + ')'
        else:
            name = ''

        return Parameter(val, units, name=name, _composedPrio=0, _dependencies=deps)

    def __neg__(self):
        """Unary negative operator

        :returns: The Parameter resulting from multiplication of its value by -1.
        :rtype: :py:class:`Parameter`

        Usage::

            >>> param1 = Parameter(2, 'mV')
            >>> param2 = -param1
            >>> param2.units
            'mV'
            >>> param2.value
            -2
            >>> param2.name
            '-param1'

        :meta public:
        """
        self._checkUsableInOp()

        deps = [self] if self._isNamed() else []

        units = self._units
        val = -1 * self.value

        if self._isNamed():
            name = f'-({self})' if self._composedPrio > 2 else f'-{self}'
        else:
            name = ''

        return Parameter(val, units, name=name, _composedPrio=2, _dependencies=deps)

    def __float__(self):
        return float(self._valueInSI())

    def __int__(self):
        return int(self._valueInSI())

    def __bool__(self):
        return bool(self._value)

    def __repr__(self):
        if self._isNamed():
            return self._name
        else:
            return f'{self.value}'

    @classmethod
    def _startRecording(cls):
        if cls._PARAMETER_USAGE_RECORD is None:
            cls._PARAMETER_USAGE_RECORD = [[[]], 0]
        else:
            cls._PARAMETER_USAGE_RECORD[0].append([])
            cls._PARAMETER_USAGE_RECORD[1] += 1

    @classmethod
    def _endRecordingAndGetUsage(cls):
        lst, depth = cls._PARAMETER_USAGE_RECORD

        res = []
        for i in range(depth, len(lst)):
            res += lst[i]

        if depth == 0:
            cls._PARAMETER_USAGE_RECORD = None
        else:
            cls._PARAMETER_USAGE_RECORD[1] -= 1

        return res


def _extractParameterized(obj, pre=tuple(), exploredObjs=None):
    """Yield tuples that describe a path to parameterized objects."""
    if exploredObjs is None:
        exploredObjs = set()

    if isinstance(obj, ParameterizedObject):
        # Only yield obj if it should be included in param tables
        if obj._includeinParamTables():
            yield pre + (obj,)
        
        if obj not in exploredObjs:
            # Check subobjects
            for subObj in obj._getSubParameterizedObjects():
                yield from _extractParameterized(subObj, pre, exploredObjs)

            # Check Parameters that are themselves Parameterized
            for param in obj._getAllParams():
                if isinstance(param, Parameter) and isinstance(param._value, ParameterizedObject):
                    yield from _extractParameterized(param._value, pre, exploredObjs)

    if isinstance(obj, NamedObject) and obj not in exploredObjs:
        # Handle children
        for subObj in obj._getChildrenOfType(ParameterizedObject, NamedObject):
            yield from _extractParameterized(subObj, pre + (obj,), exploredObjs)

    exploredObjs.add(obj)


def _dictRowsToTable(rowList):
    """
    Take a list of rows in which each row is a list of key-value pair and output a
    list of column name and a list of rows in which rows just contain the value associated
    with the column or None if it is not defined for this row.
    Discard columns that only contain None.
    """
    colNames = {}
    ind = 0
    for row in rowList:
        for cname, val in row.items():
            if cname not in colNames and val is not None:
                colNames[cname] = ind
                ind += 1

    table = []
    for row in rowList:
        values = [None] * len(colNames)
        for cname, val in row.items():
            if val is not None:
                values[colNames[cname]] = val
        table.append(tuple(values))

    return colNames, table


def _tableToCSV(rowList, filename, sep='\t', **kwargs):
    """Export a single table to a CSV file
    Return the path to the exported file.
    """

    colNames, table = _dictRowsToTable(rowList)
    filePath = f'{filename}.csv'

    with open(filePath, 'w') as f:
        f.write(sep.join(_formatValue(cn, bold=True, **kwargs) for cn in colNames) + '\n')
        for row in table:
            f.write(sep.join(map(str, [val if val is not None else '' for val in row])) + '\n')

    return filePath


_TEXT_FORMAT_UNICODE = 'unicode'
_TEXT_FORMAT_LATEX = 'latex'

def _formatValue(value, **kwargs):
    """Return a formated string to be used in parameter tables"""

    textFormat = kwargs.get('textFormat', _TEXT_FORMAT_UNICODE)
    numPrecision = kwargs.get('numPrecision', 10)
    isComputation = kwargs.get('isComputation', False)
    funcRelFilePath = kwargs.get('funcRelFilePath', True)
    latexPrefix = kwargs.get('latexPrefix', None)
    bold = kwargs.get('bold', False)

    if isinstance(value, Parameter):
        # Do not modify Parameter objects, they will be expanded later
        return value
    elif hasattr(value, '__code__'):
        filename = value.__code__.co_filename
        m = re.match(r'<ipython-input-(\d)+-', filename)
        if m is not None:
            filename = f'ipython-{m.group(1)}'
        if textFormat == _TEXT_FORMAT_LATEX:
            qualName = _formatValue(value.__qualname__, **kwargs)
            if funcRelFilePath:
                filename = os.path.relpath(filename)
            filename = _formatValue(filename, **kwargs)
            return f'\\texttt{{{qualName}(...)}} at \\texttt{{{filename}}}:{value.__code__.co_firstlineno}'
        else:
            return f'{value.__qualname__}(...) at {filename}:{value.__code__.co_firstlineno}'
    elif isinstance(value, _ParamList):
        # Try to match simple patterns
        if len(value) > 2:
            if isinstance(value[0], numbers.Number):
                srtList = sorted(value)
                step = srtList[1] - srtList[0]
                if srtList == list(numpy.arange(srtList[0], srtList[-1] + step, step)):
                    start = _formatValue(srtList[0], **kwargs)
                    finish = _formatValue(srtList[-1], **kwargs)
                    stp = _formatValue(step, **kwargs)
                    return f'{start} to {finish}' + (f' step {stp}' if step != 1 else '')
        fvalues = [_formatValue(v, **kwargs) for v in value]
        return ', '.join(fv for fv in fvalues if fv is not None)
    elif isinstance(value, steps.API_2.model._SubReactionList):
        if textFormat == _TEXT_FORMAT_LATEX:
            if value._isFwd:
                return f'\\detokenize{{{value._parent.lhs}}} $\\rightarrow$~\\detokenize{{{value._parent.rhs}}}'
            else:
                return f'\\detokenize{{{value._parent.rhs}}} $\\rightarrow$~\\detokenize{{{value._parent.lhs}}}'
        else:
            return str(value)
    elif isinstance(value, Units):
        if textFormat == _TEXT_FORMAT_UNICODE:
            return value._toUnicode()
        elif textFormat == _TEXT_FORMAT_LATEX:
            return value._toLatex()
        else:
            return value._str
    elif isinstance(value, numbers.Number):
        if textFormat == _TEXT_FORMAT_LATEX:
            return f'\\num{{{value:.{numPrecision}g}}}'
        else:
            return f'{value:.{numPrecision}g}'
    elif value is None:
        return None
    else:
        if textFormat == _TEXT_FORMAT_LATEX:
            ret = str(value)
            if isComputation:
                ret = re.sub(r' \*\* ([0-9]+)', lambda m: f'^{{{m.group(1)}}}', ret)
                ret = re.sub(r'([a-zA-Z]\w+)', lambda m: f'\\mathrm{{\\detokenize{{{m.group(1)}}}}}', ret)
                ret = re.sub(r'([0-9]+\.[0-9]+)', lambda m: f'{float(m.group(1)):{numPrecision}g}', ret)
                ret = re.sub(r'([^\s]+) / ([^\s]+)', lambda m: f'\\frac{{{m.group(1)}}}{{{m.group(2)}}}', ret)
                ret = ret.replace('*', '\\times')
                ret = f'${ret}$'
            elif latexPrefix is None or not ret.startswith(latexPrefix):
                toEscape = '#_^&%${}'
                for c in toEscape:
                    ret = ret.replace(c, f'\\{c}')
            else:
                ret = ret[len(latexPrefix):]
            if bold:
                ret = f'\\textbf{{{ret}}}'
            return ret
        else:
            return str(value)


class _ParamList:
    """Utility class for parameter grouping"""
    def __init__(self, lst):
        self.lst = lst

    def __iter__(self):
        return iter(self.lst)

    def __getitem__(self, i):
        return self.lst[i]

    def __hash__(self):
        return hash(frozenset(self.lst))

    def __eq__(self, other):
        return isinstance(other, _ParamList) and frozenset(self.lst) == frozenset(other.lst)

    def __repr__(self):
        return 'PL' + str(self.lst)

    def __len__(self):
        return len(self.lst)


def _groupParams(params, groupingInds=None):
    """
    Take a list of parameter tuples as input and outputs a list of tuples of _ParamLists
    whose cartesian product yields the orginal parameter tuples.
    """
    n = len(params[0])
    if groupingInds is None:
        groupingInds = list(range(n))

    groupings = [params]
    for i in range(n):
        if i in groupingInds:
            # group by n-1 column
            val2Lst = {}
            grouped = False
            for row in params:
                val = (row[:i], row[i + 1 :])
                val2Lst.setdefault(val, []).append(row[i])
                if len(val2Lst[val]) > 1:
                    grouped = True
            if grouped:
                newParams = [val[0] + (_ParamList(lst),) + val[1] for val, lst in val2Lst.items()]
                groupings.append(_groupParams(newParams))
    return min(groupings, key=lambda x: len(x))


def _exportToCSV(cls2Tables, cls2NamedParams, filename, **kwargs):
    """Export parameter data to a series of CSV files
    One table per file, separate file for named parameters.
    Return a list of paths to exported files
    """

    allPaths = []

    for cls, rowList in cls2Tables.items():
        allPaths.append(_tableToCSV(rowList, f'{filename}_{cls._getDisplayName()}', **kwargs))

    rowList = []
    names = set()
    for cls, rows in cls2NamedParams.items():
        for row in rows:
            if row['Name'] not in names:
                rowList.append(row)
                names.add(row['Name'])
    if len(rowList) > 0:
        allPaths.append(_tableToCSV(rowList, f'{filename}_Parameters', **kwargs))

    return allPaths


def _exportToTEX(
        cls2Tables, cls2NamedParams, filename,
        csv2latexArgs = '-l 99999 -s t -x -r 2 -z -c 0.9 -e', **kwargs
    ):
    """Export parameter data to a series of TEX files
    One table per file, separate file for named parameters.
    Requires csv2latex to be installed.
    Return a list of paths to exported files.
    """
    allPaths = []

    tmpPath = os.path.join(tempfile.gettempdir(), os.path.basename(filename))
    csvPaths = _exportToCSV(cls2Tables, cls2NamedParams, tmpPath, **kwargs)
    for path in csvPaths:
        out = subprocess.check_output(
            f'csv2latex {csv2latexArgs} {path}',
            shell=True, universal_newlines=True
        )

        texPath = filename + path[len(tmpPath):-4] + '.tex'
        with open(texPath, 'w') as f:
            for line in out.split('\n'):
                if 'documentclass' in line:
                    f.write('\\documentclass{standalone}\n')
                    continue
                if line == '\\begin{document}':
                    f.write('\\usepackage{siunitx}\n')
                f.write(line+'\n')

        allPaths.append(texPath)
        os.remove(path)

    return allPaths


def _exportToPDF(
        cls2Tables, cls2NamedParams, filename,
        **kwargs
    ):
    """Export parameter data to a series of PDF files
    One table per file, separate file for named parameters.
    Requires csv2latex and pdflatex to be installed.
    Return a list of paths to exported files.
    """
    allPaths = []

    tmpPath = os.path.join(tempfile.gettempdir(), os.path.basename(filename))
    texPaths = _exportToTEX(cls2Tables, cls2NamedParams, tmpPath, **kwargs)
    for path in texPaths:
        tmpDirPath = os.path.dirname(os.path.abspath(path))
        outDir = os.path.dirname(os.path.abspath(filename))
        subprocess.call(f'pdflatex {path}', shell=True, cwd=tmpDirPath)
        shutil.copy(path[:-4] + '.pdf', outDir)

        allPaths.append(os.path.join(outDir, os.path.basename(path[:-4]) + '.pdf'))
        os.remove(path)

    return allPaths


def ExportParameters(container, filename, method='csv', hideComputations=False, hideColumns=[], unitsToSimplify=[], **kwargs):
    r"""Export container parameters to tables

    Extracts all parameters that are used in the container and its contained objects. Parameter
    tables are then created for each object type.

    :param container: The container whose parameters are going to be exported. Most of the time,
        this will be the :py:class:`steps.API_2.sim.Simulation` object.
    :type container: Any STEPS object
    :param filename: Path and prefix for the exported table files (e.g. './parameter/Sim1' would
        yield files like './parameter/Sim1_Parameters.csv', etc.).
    :type filename: str
    :param method: Specify which file format should be used. Available options are 'csv', 'tex',
        and 'pdf'.
    :type method: str
    :param hideComputations: If true, arithmetic combinations of parameters will not be displayed
        in the tables.
    :type hideComputations: bool
    :param hideColumns: List of column names that should not be displayed in the tables.
    :type hideColumns: List[str]
    :param unitsToSimplify: List of units that should always be displayed in this specific way if
        equivalent units appear in the table. For example, some surface reaction rate constants are
        better expressed as '(mol m^-2)^-1 s^-1' so we would give this string in the list to avoid
        seeing e.g. 'mol^-1 m^2 s^-1' (which is equivalent but harder to read) in the parameter tables.
    :type unitsToSimplify: List[str]
    :param \*\*kwargs: Additional keyword arguments specified below.

    Keyword arguements depend on the `method` used but some of them are common to all methods.

    Common additional keyword arguments:

    :param numPrecision: Number of significant digits to be displayed for floating point numbers.
    :type numPrecision: int

    :param funcRelFilePath: Use relative file path for specifying the source file of a function that
        is used in a parameter. Defaults to `True`.
    :type funcRelFilePath: bool

    The 'csv' export method will export parameters as tab-separated values in a text file. Since some
    cells can contain commas, the separator has been chosen to be tabs instead of the more common comma.

    Additional keyword arguments for the 'csv' format:

    :param textFormat: Either 'unicode' or 'latex', defaults to 'unicode'. 'latex' will introduce
        latex formating in the cells, in case other csv to latex conversion programs need to be used.
    :type textFormat: str

    Both the 'tex' and 'pdf' export methods are experimental. The 'tex' method will export the parameters
    to independent .tex files each containing one table. It requires the installation of the `csv2latex`
    conversion program. The 'pdf' method will simply call `pdflatex` on the tex files resulting from the
    'tex' method.

    Additional keyword arguments for the 'tex' and 'pdf' formats:

    :param csv2latexArgs: A string of command line arguments to be passed to `csv2latex`, see the 
        `manual page <http://manpages.ubuntu.com/manpages/bionic/man1/csv2latex.1.html>`_ for details.
        Defaults to ``-l 99999 -s t -x -r 2 -z -c 0.9 -e``
    :type csv2latexArgs: str

    :param latexPrefix: If provided, special latex character will not be escaped for strings that are
        prefixed with `latexPrefix`. For example, if we specified a `Source = '\cite{ArticleRef}'`
        keyword argument to a `Parameter` object, the `{` and `}` characters will be escaped by default.
        This can be prevented by using e.g. `latexPrefix='__'` and `Source = '__\cite{ArticleRef}'`.
    :type latexPrefix: str
    """
    # TODO improvement: Add param to restrict value setting to t=0 or maybe a function that takes time
    # and returns whether it should be considered ?

    # Force latex text format if needed
    if method in ['tex', 'pdf']:
        kwargs['textFormat'] = _TEXT_FORMAT_LATEX

    # Insert the simplification of non-dimensional units
    unitsToSimplify = [Units('')] + [Units(unitstr) for unitstr in unitsToSimplify]

    # Extract all parameterized objects
    allParameterizedPaths = list(_extractParameterized(container))

    # If several paths end on the same object, only keep the longest
    # If several paths have the same maximum length, keep track of all of them
    uniquePaths = {}
    for path in allParameterizedPaths:
        key = id(path[-1])
        if key in uniquePaths:
            currLen = len(uniquePaths[key][0])
            if len(path) < currLen:
                continue
            if len(path) > currLen:
                del uniquePaths[key]
        uniquePaths.setdefault(key, []).append(path)

    # Group paths by class to establish tables
    cls2Path = {}
    for _, paths in uniquePaths.items():
        cls2Path.setdefault(paths[0][-1].__class__, []).append(paths)

    # Generate tables per class
    cls2Tables = {}
    cls2NamedParams = {}
    for cls, allPathLists in cls2Path.items():
        rowList = cls2Tables.setdefault(cls, [])
        for pathsList in allPathLists:
            # Get the corresponding object, they all end on the same
            *_, obj = pathsList[0]

            row = {}
            row[cls._getDisplayName()] = _formatValue(obj, **kwargs)

            if any(len(path) > 0 for *path, _ in pathsList):
                defList = _ParamList(
                    [path[-2] for path in pathsList if len(path) > 1 and path[-2] is not container]
                )
                row['Defined in'] = _formatValue(defList, **kwargs)

            if isinstance(obj, AdvancedParameterizedObject):
                param2Tuples = {}
                for tuples, param in obj._parameters.items():
                    param2Tuples.setdefault(param, []).append(tuples)

                newRows = []
                for param, tuples in param2Tuples.items():
                    propNames, propTable = _dictRowsToTable([{key:val for key, val in tples} for tples in  tuples])

                    groupingInds = [i for i, pn in enumerate(propNames) if pn in cls._getGroupedAdvParams()]
                    groups = _groupParams(propTable, groupingInds=groupingInds)
                    for prod in groups:
                        newProd = {**row}
                        for lst, prop in zip(prod, propNames.keys()):
                            if not isinstance(lst, _ParamList):
                                lst = _ParamList([lst])
                            newProd[prop] = _formatValue(lst, **kwargs)
                        newRows.append({**newProd, 'Value': _formatValue(param, **kwargs)})

                rowList += newRows
            else:
                rowList.append({
                    **row,
                    **{key: _formatValue(val, **kwargs) for key, val in obj._parameters.items()},
                })

        # Expand Parameters into more columns and keep track of named parameters
        namedParams = []
        for i, row in enumerate(rowList):
            newRow = {}
            for name, param in row.items():
                if isinstance(param, Parameter):
                    if param._isNamed():
                        namedParams.append(param)
                    # Also add named parameters in dependency tree
                    namedParams += [dep for dep in param._getAllDependencies() if dep._isNamed()]

                    if not hideComputations or (param._isNamed() and param._isUserDefined()):
                        newRow[f'{name} Parameter'] = _formatValue(param._name, isComputation=True, **kwargs)
                    # Convert to standard units if available
                    value = param._value
                    units = param._units
                    if param._units is not None:
                        for u in unitsToSimplify:
                            if param._units._compatibleWith(u):
                                value = param.valueIn(u)
                                units = u
                                break

                    if value is not None:
                        newRow[f'{name}'] = _formatValue(value, **kwargs)
                        if units is not None:
                            newRow[f'{name} Units'] = _formatValue(units, **kwargs)
                        if not param._isNamed():
                            for pname, pval in param._kwargs.items():
                                newRow[f'{name} {pname}'] = _formatValue(pval, **kwargs)
                else:
                    newRow[name] = param
            rowList[i] = newRow

        # Only keep named params that are not dependent on other params
        # and add the missing params
        expand = True
        dontExpand = []
        while expand:
            newNamedParams = []
            expand = False
            for param in namedParams:
                if len(param._dependencies) > 0 and param not in dontExpand:
                    for dep in param._dependencies:
                        if dep not in newNamedParams:
                            newNamedParams.append(dep)

                    if param._fullname is not None:
                        # Add the parameter itself if it was declared as a parameter explicitely
                        newNamedParams.append(param)
                        # But prevent it from being expanded in subsequent iterations
                        dontExpand.append(param)

                    expand = True
                elif param not in newNamedParams:
                    newNamedParams.append(param)
            namedParams = newNamedParams

        cls2NamedParams[cls] = [
            {
                'Name':_formatValue(np._name, **kwargs),
                'Value':_formatValue(np._value, **kwargs),
                'Units':_formatValue(np._units, **kwargs) if np._units is not None else '',
                'Computed From':_formatValue(np._fullname, isComputation=True, **kwargs),
                **{key: _formatValue(val, **kwargs) for key, val in np._kwargs.items()}
            } for np in namedParams
        ]

    # Hide columns
    for hc in hideColumns:
        for clsmap in [cls2Tables, cls2NamedParams]:
            for cls, table in cls2Tables.items():
                for row in table:
                    try:
                        del row[hc]
                    except:
                        continue

    if method == 'csv':
        return _exportToCSV(cls2Tables, cls2NamedParams, filename, **kwargs)
    elif method == 'tex':
        return _exportToTEX(cls2Tables, cls2NamedParams, filename, **kwargs)
    elif method == 'pdf':
        return _exportToPDF(cls2Tables, cls2NamedParams, filename, **kwargs)
    else:
        raise NotImplementedError(f'Unknown export method "{method}".')



###################################################################################################

# The prefix used for naming classes in the cython bindings
_CYTHON_PREFIX = '_py_'

###################################################################################################
# STEPS objects utility classes


class SolverPathObject:
    """Base class for all objects susceptible to be part of a SimPath."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _solverStr(self):
        """Return the string that is used as part of method names for this specific object."""
        return self.__class__.__name__

    def _solverId(self):
        """Return the id that is used to identify an object in the solver."""
        try:
            return (self.name,)
        except AttributeError:
            return tuple()

    def _solverKeywordParams(self):
        """Return the additional keyword parameters that should be passed to the solver"""
        return {}

    def _solverModifier(self):
        """Return None or a function that will modify the output from the solver."""
        return None

    def _solverSetValue(self, valName, v):
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

    def _simPathCheckParent(self):
        """
        Determines whether the object needs to be in the parent's children in order to be added to
        a simulation path
        """
        return True

    def __hash__(self):
        return id(self)


class SolverRunTimeObject:
    """Base class for object whose value will only be evaluated at runtime"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _runTimeFunc(self, opType, sim, prefix, suffix, descriptor):
        """Return a function that will be called during runtime evaluation of the path."""
        return None


class StepsWrapperObject:
    """Base class for all objects that are wrappers for steps objects."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _getStepsObjects(self):
        """Return a list of the steps objects that this named object holds. Should be overloaded."""
        return []

    @classmethod
    def _FromStepsObject(cls, obj, *args, **kwargs):
        """Create the interface object from a STEPS object."""
        raise NotImplementedError()


class KeyOrderedDict:
    """Utility class: dict in which iteration is ordered based on key values"""
    def __init__(self):
        self._dict = {}
        self._keys = []
        self._sorted = True

    def keys(self):
        return iter(self)

    def items(self):
        for key in self:
            yield (key, self._dict[key])

    def __contains__(self, key):
        return key in self._dict

    def __getitem__(self, key):
        return self._dict[key]

    def __setitem__(self, key, v):
        if key not in self._dict:
            self._keys.append(key)
            self._sorted = False
        self._dict[key] = v

    def __delitem__(self, key):
        del self._dict[key]
        self._keys.remove(key)

    def __iter__(self):
        if not self._sorted:
            self._keys.sort()
            self._sorted = True
        return iter(self._keys)


class NamedObject(SolverPathObject):
    """Base class for all objects that are named and can have children

    :param name: Name of the object. If no name is provided, the object gets an automatically
        generated name based on its class.
    :type name: str

    All classes that inherit from :py:class:`NamedObject` can be built with a ``name`` keyword
    parameter. For steps objects, it corresponds to the identifiers used in STEPS. Note that some
    names are forbidden because they correspond to names of attributes or methods of classes
    defined in this interface. Since most objects implement ``__getattr__`` attribute style access
    to children, the names of these methods / attributes could clash with object names. It is thus
    possible that the contructor of :py:class:`NamedObject` raises an exception when trying to
    name an object with one of these forbidden names.

    In addition to a name, this class holds a list of children objects that can be accessed with
    :py:func:`__getattr__` and :py:func:`ALL`.

    .. note::
        This class should not be instantiated by the user, it is only documented for clarity
        since a lot of other classes inherit from it.
    """

    _nameInds = {}  # TODO Not urgent: use threading local?
    _forbiddenNames = set()
    _allowForbNames = False
    _children = {} # Use default empty value to avoid infinite recursion in __getattr__

    def __init__(self, *args, name=None, **kwargs):
        super().__init__(*args, **kwargs)
        if name is None:
            self.name = self.__class__._GetDefaultName()
            self._autoNamed = True
        else:
            self.name = name
            self._autoNamed = False
        if len(NamedObject._forbiddenNames) == 0:
            NamedObject._loadForbiddenNames()
        if self.name in NamedObject._forbiddenNames:
            if NamedObject._allowForbNames:
                warnings.warn(
                    f"'{self.name}' is a reserved name, SimPath functionnalities might "
                    f"be impacted because of this."
                )
            else:
                raise Exception(f"Cannot call an element '{self.name}', this name is reserved")

        self._children = KeyOrderedDict()
        self._parents = {}

    @classmethod
    def _loadForbiddenNames(cls):
        """Load names that will not be allowed as NamedObject names.

        This prevents cases in which a STEPS object would be inaccessible through SimPath, Simulation, or
        ResultSelector, because its name is the same as the name of an attribute, a method or a property of
        these classes (they have priority over __getattr__).
        """
        from steps.API_2 import sim as nsim
        from steps.API_2 import saving as nsaving

        classes = [nsim.SimPath, nsim.Simulation, nsaving.ResultSelector, nsaving._ResultPath,
                nsaving._ResultList, nsaving._ResultCombiner]
        for cl in classes:
            cls._forbiddenNames |= set(dir(cl))

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
        # TODO later release: revert this temporary change for split meshes
        # if name.startswith('__'):
        if name.startswith('__') and not name.startswith('__MESH'):
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
        r"""Return all children of the object, optionally filtered by class

        Takes a variable number of parameters, if no parameters are given, it returns all children
        of the object. Otherwise, if types are given, it returns the children that match at least
        one of the given types.

        :param \*cls: Variable number of classes

        :returns: A generator that iterates over all children that match the class criteria.
        :rtype: Generator[NamedObject]

        Usage::

            obj.ALL()                    # Return all children
            obj.ALL(Species)             # Return all children Species
            obj.ALL(Reaction, Diffusion) # Return all children that are either Reaction or
                                         # Diffusion
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

        Note that if a `name` keyword argument is provided to `Params(...)`, it will override the
        name inferred from the destination variable.

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
                        if 'name' in arg.kwargs:
                            res.append(cls(*arg.args, **arg.kwargs, **kwargs))
                        else:
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
        return self.name if self._children is not NamedObject._children else super().__repr__()


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

        self._currentUsers = []

    def _getDefaultCurrUsedVal(self):
        return None

    def _registerUser(self, cls, user, addAsElement):
        self._currentUsers.append(user)
        if addAsElement:
            self._addChildren(user)
            user._parents[cls] = self
            for chld in user._getAdditionalChildren():
                self._addChildren(chld)
                chld._parents[cls] = self

    def __enter__(self):
        if UsableObject._currUsed[self.__class__] is None:
            UsableObject._currUsed[self.__class__] = self
        else:
            raise Exception(f'Cannot use two {self.__class__.__name__} objects simultaneously.')
        self._currentUsers = []
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        UsableObject._currUsed[self.__class__] = None
        if (exc_type, exc_val, exc_tb) == (None, None, None):
            for user in self._currentUsers:
                user._exiting(self)

    @staticmethod
    def _getUsedObjectOfClass(cls):
        isMulti = MultiUsable in cls.__mro__
        allCls = collections.deque([cls])
        while len(allCls) > 0:
            uc = allCls.popleft()
            allCls.extend(uc.__subclasses__())
            if uc in UsableObject._currUsed and UsableObject._currUsed[uc] is not None:
                if isMulti:
                    # Return a copy of the list, it would otherwise create problems
                    # when objects keep track of the values returned by _getUsedObjects()
                    return copy.copy(UsableObject._currUsed[uc])
                else:
                    return UsableObject._currUsed[uc]
        return [] if isMulti else None


class MultiUsable(UsableObject):
    """Base class for steps objects that can be used several times simultaneously."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _getDefaultCurrUsedVal(self):
        return []

    def __enter__(self):
        UsableObject._currUsed[self.__class__].append(self)
        self._currentUsers = []
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        UsableObject._currUsed[self.__class__].remove(self)
        if (exc_type, exc_val, exc_tb) == (None, None, None):
            for user in self._currentUsers:
                user._exiting(self)


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
            self._exited = set()
            # Register the object as child of the used objects
            list(self._getUsedObjects(addAsElement=addAsElement, **kwargs))

        def _getObjOfClass(self, cls, addAsElement=True):
            obj = UsableObject._getUsedObjectOfClass(cls)
            if obj is not None and not isinstance(obj, list):
                obj._registerUser(cls, self, addAsElement)
            return obj

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

        def _exiting(self, parent):
            """The callback should only be called once"""
            if parent not in self._exited:
                self._exited.add(parent)
                self._exitCallback(parent)

        def _exitCallback(self, parent):
            """
            Method to be called the first time we get out of a context manager in which the object as been
            declared.
            """
            pass

        def _decl(self):
            """
            Return a string that describes the declaration of the object. It returns its name by default
            but it can contain e.g. the file and line at which it was declared.
            """
            return self.name

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


def FreezeAfterInit(cls):
    """Class decorator for preventing dynamical attribute creation outside __init__"""

    def _setattr(self, name, val):
        if self.__freezeCounter > 0 or name in self.__dict__ or hasattr(self.__class__, name):
            object.__setattr__(self, name, val)
        else:
            raise AttributeError(f'There is no attribute named {name} in {self}.')

    oldInit = cls.__init__
    def newInit(self, *args, **kwargs):
        self.__freezeCounter += 1
        oldInit(self, *args, **kwargs)
        self.__freezeCounter -= 1

    cls.__freezeCounter = 0
    # Wrap existing __setattr__, if any
    if '__setattr__' in cls.__dict__:
        oldSetAttr = cls.__setattr__
        def _newSetAttr(self, name, val):
            try:
                _setattr(self, name, val)
            except AttributeError:
                oldSetAttr(self, name, val)
        cls.__setattr__ = _newSetAttr
    else:
        cls.__setattr__ = _setattr

    cls.__init__ = functools.wraps(cls.__init__)(newInit)

    return cls


def IgnoresGhostElements(func):
    """Method decorator to raise a warning if the method is used with non-owned elements"""

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if self._local and not self._owned:
            warnings.warn(f'Local {func.__name__} accesses do not take non-owned elements into account.')
        if 'owned' in self._callKwargs:
            owned = self._callKwargs['owned']
            del self._callKwargs['owned']
            try:
                val = func(self, *args, **kwargs)
            finally:
                self._callKwargs['owned'] = owned
            return val
        else:
            return func(self, *args, **kwargs)
    return wrapper


###################################################################################################
# File version management utilities

class Versioned:
    """Base class for objects that can come in more than one version

    The main use case for this is when behavior changes with different STEPS versions but we still need
    to have code that is compatible with previous STEPS versions.

    Several methods sharing the same name can be defined and the `_versionRange` decorator should be
    used to associate a version range to each of them. Once the `_setVersion` method is called on an
    object of this class, the methods corresponding to the supplied version will be used by the object.

    The `_setVersion` method is called at construction with STEPS current version.
    """

    _methods = {}

    def __init__(self, *args, version=steps.__version__, **kwargs):
        super().__init__(*args, **kwargs)

        self._setVersion(version)

    @staticmethod
    def _versionRange(above=None, belowOrEq=None):
        """
        Decorator that helps selecting the correct version of a method

        When a method is decorated with this, it will only be used if the version passed to
        `_setVersion` method is inside the range specified here (above < version <= belowOrEq).
        """
        minVersion = Versioned._parseVersion(above)
        maxVersion = Versioned._parseVersion(belowOrEq)
        class decorator:
            def __init__(self, func):
                Versioned._methods.setdefault(func.__name__, []).append((minVersion, maxVersion, func))
                self.func = func
            def __set_name__(self, owner, name):
                methods = {}
                for cls in owner.__mro__:
                    for _name, st in cls.__dict__.get('_versionedMethods', {}).items():
                        methods.setdefault(_name, set())
                        methods[_name] |= st
                setattr(owner, '_versionedMethods', methods)

                if name in Versioned._methods:
                    for funcTuple in Versioned._methods[name]:
                        owner._versionedMethods.setdefault(name, set()).add(funcTuple)
                    del Versioned._methods[name]
                setattr(owner, name, self.func)

        return decorator

    def _setVersion(self, version):
        """Set the version of a versioned object

        It will check all methods that were registered with the `_versionRange` decorator.
        If the version parameter is inside the version range of a specific method, this method will be
        monkey patched to the object.
        """
        version = Versioned._parseVersion(version)
        self._version = version
        cls = self.__class__
        if hasattr(cls, '_versionedMethods'):
            for name, versions in cls._versionedMethods.items():
                res = None
                for minv, maxv, func in versions:
                    if (minv is None or minv < version) and (maxv is None or version <= maxv):
                        if res is not None:
                            raise Exception(
                                f'Several methods could be used for method {name} with version {version}.'
                            )
                        res = func
                if res is None:
                    raise Exception(f'No methods were found for {name} with version {version}.')
                setattr(self, name, types.MethodType(res, self))

    @staticmethod
    def _parseVersion(versionStr):
        """Parse a version string to a tuple of integers

        `packaging.version.parse` would be preferable but in order to not add another dependency, we will use
        this simple parsing method.
        """
        if versionStr is None:
            return None
        if isinstance(versionStr, tuple):
            return versionStr
        try:
            return tuple(map(int, versionStr.split('.')))
        except:
            raise ValueError(f'Invalid version string "{versionStr}".')


###################################################################################################
# Various utilities

class ReadOnlyDictInterface:
    """Base class for objects that should behave like read-only dicts
    __getitem__() and keys() should be implemented.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        raise NotImplementedError('Cannot write data, the dictionary is read-only.')

    def __getitem__(self, key):
        raise NotImplementedError()

    def keys(self):
        raise NotImplementedError()

    def __iter__(self):
        for key in self.keys():
            yield key

    def items(self):
        for key in self.keys():
            yield key, self[key]

    def values(self):
        for key in self.keys():
            yield self[key]

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        return key in self.keys()

    def __eq__(self, val):
        return sorted(self.items()) == sorted(val.items())


class MutableDictInterface(ReadOnlyDictInterface):
    """Base class for objects that should behave like dicts
    __getitem__(), __setitem__() and keys() should be implemented.
    """
    def __setitem__(self, key, value):
        raise NotImplementedError()


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
    if isinstance(s, list):
        return s
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


def nparray(data):
    """Return a numpy array with the appropriate dtype"""
    try:
        return numpy.array(data)
    except ValueError:
        return numpy.array(data, dtype=object)


def getValueIfAllIdentical(lst):
    """If all values in lst are identical, return a list with just this value;
    otherwise return None"""
    try:
        val = next(iter(lst))
    except StopIteration:
        return None
    return None if any(v != val for v in lst) else [val]

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


def extractArgs(args, kwargs, lst, raiseIfMoreArgs=True):
    """
    Extract arguments from args and kwargs based on name and type.

    This function also removes arguments from args and kwargs if they got extracted.

    :param args: Tuple of arguments
    :type args: Tuple[Any, ...]
    :param kwargs: Keyword arguments
    :type kwargs: Dict[str, Any]
    :param lst: The list of name, type and default value
    :type lst: List[Tuple[str, Union[Type, Tuple[Type, ...]], Any]]
    :param raiseIfMoreArgs: If True, raise an exception if some positional arguments were not extracted
    :type raiseIfMoreArgs: bool

    :return: Tuple of values of size n + 2 with n the length of the `lst` parameter. The two last
        elements are the updated `args` and `kwargs` after extraction
    :rtype: Tuple[Any, ..., Tuple[Any, ...], Dict[str, Any]]
    """
    args = list(args)
    kwargs = dict(kwargs)
    res = []
    for name, tpe, default in lst:
        val = default
        if name in kwargs:
            val = kwargs[name]
            del kwargs[name]
        else:
            for i, arg in enumerate(args):
                if isinstance(arg, tpe) or arg == default:
                    val = arg
                    del args[i]
                    break
        res.append(val)
    if raiseIfMoreArgs and len(args) > 0:
        raise Exception(f'Unused arguments: {args}')
    return tuple(res) + (args, kwargs)


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


class MessagesTypes(str, Enum):
    WARNING = '\033[93m'
    ERROR = '\033[91m'
    SUCCESS = '\033[92m'
    UNDERLINE = '\033[4m'
    BOLD = '\033[1m'
    HEADER = BOLD + UNDERLINE
    _END_COLOR = '\033[0m'


def _print(msg, prio, tpe=None, indent=0):
    """
    Print a message if its priority permits it and if we are in rank one in case of an MPI
    simulation.
    """
    from steps.API_2 import sim as nsim
    if prio <= VERBOSITY and nsim.MPI._shouldWrite:
        msg, *sublines = msg.split('\n')
        if tpe is not None and hasattr(sys.stdout, 'isatty') and sys.stdout.isatty():
            msg = tpe + msg + MessagesTypes._END_COLOR
        if indent > 0:
            msg = '    ' * indent + msg
        print(msg)
        for sl in sublines:
            _print(sl, prio, tpe=tpe, indent=indent + 1)

###################################################################################################
# Docstrings utilities


def _paddedStringsfromRow(row, colLengths):
    return ' '.join(val + ' '*(cl - len(val)) for val, cl in zip(row, colLengths))


def GetDocstringTable(colNames, rows, indentStr=''):
    lines = []
    colLengths = []
    for c, cname in enumerate(colNames):
        colLengths.append(max(len(cname), max(len(row[c]) for row in rows)))
    sepLine = ' '.join('='*cl for cl in colLengths)
    lines.append(sepLine)
    lines.append(indentStr + _paddedStringsfromRow(colNames, colLengths))
    lines.append(indentStr + sepLine)
    for row in rows:
        lines.append(indentStr + _paddedStringsfromRow(row, colLengths))
    lines.append(indentStr + sepLine)

    return '\n'.join(lines)


class Facade:
    """Base class for classes that are exposed to the user but can instantiate objects from their
    subclasses. For example, ``Compartment`` is a facade to _DistCompartment, etc.
    Subclasses can define the ``_FACADE_TITLE_STR`` class attribute to give a section title in the
    documentation.
    """
    _FACADE_ATTR_NAME = '_FACADE_TITLE_STR'

###################################################################################################
# Docstring editing


# Insert Prefix and Units information into the Parameter docstring
prefixTable = GetDocstringTable(
    ['Prefix', 'Name', 'Base 10'], 
    [(f"``'{pref}'``", name, fr'10\ :sup:`{exp}`') for pref, (exp, name) in Units._SI_PREFIXES.items()],
    '        '
)
unitTable = GetDocstringTable(
    ['Unit', 'Name', 'Quantity'],
    [(f"``'{unit}'``", name, quant) for unit, (_, _, name, quant) in Units._SI_UNITS.items()],
    '        '
)
Parameter.__doc__ = Parameter.__doc__.format(
    prefixTable=prefixTable,
    unitTable=unitTable,
)

