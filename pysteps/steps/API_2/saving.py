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
import datetime
import functools
import math
import numbers
import operator
import os
import pickle
import re
import sqlite3
import struct
import sys
import warnings

from xml.etree import ElementTree

import numpy
import steps

from . import geom as ngeom
from . import _saving_optim as nsaving_optim
from . import sim as nsim
from . import utils as nutils

__all__ = [
    'ResultSelector',
    'DatabaseHandler',
    'DatabaseGroup',
    'SQLiteDBHandler',
    'SQLiteGroup',
    'HDF5Handler',
    'HDF5Group',
    'XDMFHandler',
]

###################################################################################################
# Result selectors


class _MetaData:
    """
    Small utility class for handling metadata setting and getting through __getitem__ and
    __setitem__. Behaves like a dict but performs additional checks when setting values.
    """

    def __init__(self, parent):
        self._parent = parent
        self._dict = {}

    def _checkKey(self, key):
        """Check that the given key is valid."""
        if not isinstance(key, str):
            raise TypeError(
                f'MetaData can only be accessed by first specifying a string key, '
                f'got a {type(key)} instead.'
            )

    def _clear(self):
        self._dict = {}

    def __getitem__(self, key):
        """Return the metadata corresponding to key."""
        self._checkKey(key)
        if key not in self._dict:
            raise KeyError(f'No metadata corresponding to key {key}.')
        return nutils.nparray(self._dict[key])

    def __setitem__(self, key, val, _internal=False):
        """Set the metadata corresponding to key."""
        self._checkKey(key)
        if self._parent._savingStarted() and not _internal:
            raise Exception(f'Cannot save metadata once sim.newRun() has been called.')

        # Convert to list in case val is a generator
        lst = list(val)

        if len(lst) != self._parent._getEvalLen():
            raise Exception(
                f'Expected a list of length {self._parent._getEvalLen()}, got a list '
                f'of length {len(lst)}.'
            )
        if not all(isinstance(v, (numbers.Number, str)) or v is None for v in lst):
            raise TypeError(f'Metadata can only be composed of numbers and / or strings.')

        if key in self._dict and self._dict[key] != lst and not _internal:
            warnings.warn(
                f'The metadata associated with key {key} was already set, replacing with new values.'
            )

        self._dict[key] = lst

    def __iter__(self):
        return iter(self._dict)

    def keys(self):
        return self._dict.keys()

    def items(self):
        return self._dict.items()

    def get(self, key, default):
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        return key in self._dict


class ResultSelector:
    """Class to describe which data should be saved during simulation

    :param sim: The simulation for which we want to select data
    :type sim: :py:class:`steps.API_2.sim.Simulation`

    This class works in a way that is very similar to :py:class:`steps.API_2.sim.SimPath`, paths to
    the data that should be saved are built in the same way, using dot syntax. For
    :py:class:`steps.API_2.sim.SimPath`, the root of the path is the simulation itself and when a
    path is completed with a property (e.g. `Count`), it returns the actual value in the
    simulation. Since :py:class:`ResultSelector` aims at describing the data to be saved, we have
    to use a different root for our paths::

        >>> sim.comp1.S1.Count
        13
        >>> rs = ResultSelector(sim)
        >>> rs.comp1.S1.Count
        comp1.S1.Count

    While the path whose root is the actual simulation returns a number, the path whose root is the
    result selector object does not.

    Any methods defined in :py:class:`steps.API_2.sim.SimPath` can be used to build result selector
    paths. In addition, result selectors can be combined using standard arithmetic operators (see
    :py:func:`ResultSelector.__add__`, etc.) and can be concatenated with ``<<`` (see
    :py:func:`ResultSelector.__lshift__`)::

        rs1 = rs.comp1.S1.Count + rs.comp1.S2.Count  # This result selector will save a single
                                                     # value that corresponds to the sum of S1
                                                     # and S2 counts in comp1.

        rs2 = rs.comp1.S1.Count << rs.comp1.S2.Count # This one will save 2 values, the count of
                                                     # S1 in comp1 and the count of S2 in comp1.

    Result selectors can also transform data and only save the result of the transformation::

        rs3 = rs.SUM(rs.TETS(tetlst).S1.Count) # This will save only one value: the total number
                                               # of S1 in all the tetrahedrons in tetlst.

    Once we defined all our result selectors, we need to add them to the
    :py:class:`steps.API_2.sim.Simulation` so that  the corresponding data gets saved during
    simulation runs. This is done with e.g.::

        sim.toSave(rs1, rs2, rs3, dt=0.01) # Save the three result selectors every 0.01 seconds.

    After simulations have been run, results can be accessed with the same result selector
    objects::

        rs1.data[0] # Accessing the data saved during run 0
        rs1.time[0] # The time points associated to each saving for run 0

    Usage of result selectors is presented in more details in the user guide.
    """

    def __init__(self, sim, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if not isinstance(sim, nsim.Simulation):
            raise TypeError(f'Expected a Simulation object, got {sim} instead.')
        self.sim = sim
        self._dataHandler = _MemoryDataHandler(self)

        self._saveDt = None
        self._saveTpnts = None
        self._saveTind = None
        self._nextTime = math.inf

        self._addedToSim = False
        self._selectorInd = None
        self._optimGroupInd = None

        self._labels = None
        self._metaData = _MetaData(self)

        # When the result selector gets distributed across MPI rank, it keeps a list of the indexes
        # of its values in the original result selector
        self._distrInds = None
        self._fullLen = None

    @classmethod
    def FromFile(cls, path):
        """Load data that has been saved to a file

        :param path: Path to the file
        :type path: str

        Result selectors that have been saved to file (with the :py:func:`toFile` method), can
        then be loaded in a different python process and be used in the same way as in the
        simulaiton process.

        Usage::

            rs1 = ResultSelector.FromFile('path/to/file')

            plt.plot(rs1.time[0], rs1.data[0])
            plt.legend(rs1.labels)
            ...
        """
        return _ReadOnlyResultSelector(_FileDataHandler(None, path))

    @property
    def time(self):
        """Get the time points at which data saving was done for this result selector

        An accessor to the timepoints data that should then be indexed with square
        brackets notation. The underlying data it two dimensional; the first dimension
        corresponds to runs and the second to time.

        :type: Data accessor, read-only

        Usage assuming 5 runs of 1s with data saving every 10ms::

            >>> rs1.time[0]      # Time points of first run
            array([0., 0.01, 0.02, ..., 0.98, 0.99, 1.])
            >>> rs1.time[0, -1]  # Last time point of first run
            array(1.)
            >>> rs1.time[0][-1]  # Same as above
            array(1.)
            >>> rs1.time[:, -1]  # Last time point of all 5 runs
            array([1., 1., 1., 1., 1.])
            >>> rs1.time[1:3, 0] # First time point of 2nd and 3rd runs
            array([0, 0])
            >>> rs1.time[...]    # All time points of all runs
            array([[0., 0.01, 0.02, ..., 0.98, 0.99, 1.],
                   [0., 0.01, 0.02, ..., 0.98, 0.99, 1.],
                   [0., 0.01, 0.02, ..., 0.98, 0.99, 1.],
                   [0., 0.01, 0.02, ..., 0.98, 0.99, 1.],
                   [0., 0.01, 0.02, ..., 0.98, 0.99, 1.]])

        .. warning::
            Although the type of this property implements square bracket element access, it is
            not a list or an array itself and does not directly contain the data. The data is only
            really accessed upon using the square bracket notation. To force the retrieval of all
            the data, it is possible to use the ellipsis notation in square brackets:
            ``rs.time[...]``.
        """
        self._checkAddedToSim()
        return self._dataHandler.time()

    @property
    def data(self):
        """Get the data that was saved by this result selector

        An accessor to the data that should then be indexed with square brackets notation
        The underlying data it three dimensional; the first dimension corresponds to runs, the
        second to time, and the third to saved paths.

        :type: Data accessor, read-only

        Usage assuming 5 runs of 3s, saving 3 values every 1 ms::

            >>> rs1.data[0]      # Data from the first run
            array([[312., 221.,   0.],
                   [310., 219.,   2.],
                   [308., 217.,   4.],
                   ...
                   [206., 115., 106.],
                   [205., 114., 107.],
                   [205., 114., 107.]])
            >>> rs1.data[0, -1]  # Data corresponding to the last time point of first run
            array([205., 114., 107.])
            >>> rs1.data[0][-1]  # Same as above
            array([205., 114., 107.])
            >>> rs1.data[:, -1]  # Data corresponding to the last time point of all 5 runs
            array([[205., 114., 107.],
                   [189.,  98., 123.],
                   [188.,  97., 124.],
                   [185.,  95., 127.],
                   [198., 107., 114.]])
            >>> rs1.data[0, :, 0] # First saved value for all time points of first run
            array([312., 310, 308, ..., 206, 205, 205])
            >>> rs1.data[...]    # All data from all runs
            array([[[312., 221.,   0.],
                    [310., 219.,   2.],
                    [308., 217.,   4.],
                    ...,
                    [206., 115., 106.],
                    [205., 114., 107.],
                    [205., 114., 107.]],
            ...
                   [[312., 221.,   0.],
                    [309., 218.,   3.],
                    [305., 214.,   7.],
                    ...,
                    [199., 108., 113.],
                    [199., 108., 113.],
                    [198., 107., 114.]]])

        .. warning::
            Although the type of this property implements square bracket element access, it is
            not a list or an array itself and does not directly contain the data. The data is only
            really accessed upon using the square bracket notation. To force the retrieval of all
            the data, it is possible to use the ellipsis notation in square brackets:
            ``rs.data[...]``.
        """
        self._checkAddedToSim()
        return self._dataHandler.data()

    @property
    def labels(self):
        """A list of descriptions of the values saved by the result selector

        :type: List[str]

        By default labels are automatically generated from the result selector. Assuming 3 saved
        values, one can access their values with::

            >>> rs1.labels # Default values, built from the simulation paths used for saving
            ['comp.molA.Count', 'comp.molB.Count', 'comp.molC.Count']

        The labels can also be set by the user but it needs to be done before
        :py:func:`steps.API_2.sim.Simulation.newRun` has been called. Assuming 3 saved value, one
        would write::

            >>> rs1.labels = ['custom1', 'custom2', 'custom3']

        Labels be saved to whichever support the result selector is being saved to (memory, file,
        database, etc.).
        """
        return self._labels

    @labels.setter
    def labels(self, lbls):
        """Set custom labels."""
        lbls = list(lbls)
        if len(lbls) != self._getEvalLen():
            raise Exception(
                f'Expected a list of length {self._getEvalLen()}, got a list of length {len(lbls)}.'
            )
        if self._dataHandler._savingStarted():
            raise Exception(f'Cannot modify the labels once sim.newRun() has been called.')
        self._labels = lbls

    @property
    def metaData(self):
        """Meta data relative to the values saved by the result selector

        :type: Mapping[str, List[Union[str, int, float, None]]]

        This property allows the user to save additional static (i.e. not time-dependent) data
        about the values being saved by the result selector. It works as a mapping between
        arbitrary string keys and lists of values that have the same length as the number of values
        saved by the result selector.

        The meta data needs to be set before :py:func:`steps.API_2.sim.Simulation.newRun` has been
        called. Assuming 3 values saved, one could write::

            >>> rs1.metaData['key1'] = ['str1', 'str2', 'str3']
            >>> rs1.metaData['key2'] = [1, 2, 3]
            >>> rs1.metaData['key1']
            array(['str1', 'str2', 'str3'], dtype='<U4')
            >>> 'key2' in rs1.metaData
            True
            >>> 'key3' in rs1.metaData
            False

        Like labels, meta data will be saved to whichever support the result selector is being
        saved to (memory, file, database, etc.).

        .. note::
            Some path elements automatically define their own meta data, one can always check which
            meta data is already declared with e.g. ``print(rs1.metaData.keys())``

        .. warning::
            Although the type of this property implements square bracket key access, it is
            not a dict itself and does not directly contain the data. The data is only
            really accessed upon using the square bracket notation. However,  it does implement
            ``keys()``, ``items()``, ``__iter__()`` and ``__contains__()`` so it can be used like
            a dict to some extent.
        """
        return self._metaData

    def toFile(self, path, buffering=-1):
        """Specify that the data should be saved to a file

        :param path: The path to the file
        :type path: str
        :param buffering: The buffering parameter passed to the ``open()`` function, see
            https://docs.python.org/3/library/functions.html#open for details
        :type buffering: int

        This method should be called before :py:func:`steps.API_2.sim.Simulation.newRun`
        has been called. The file is written in a custom binary format and can be read in a
        different python process by creating a result selector from file with
        :py:func:`ResultSelector.FromFile`.

        .. warning::
            After all simulations are finished, depending on the buffering policy, it is possible
            that the file does not contain all the data. The data will be flushed to the file upon
            destruction of the result selector (when the python process ends for example). This
            should not create any issues for using the result selector in the process in which it
            was created (because the data that might not be written to file is kept in memory) but
            it could create issues when trying to load the file from another python process while
            the first one is still running.
        """
        self._checkComplete()
        self._dataHandler = _FileDataHandler(self, path, self._getEvalLen(), buffering)

    def _newRun(self):
        """Signal that a new run of the simulation started."""
        self._dataHandler._newRun()

        # Initialize time save points
        if self._saveDt is not None:
            self._saveTind = 0
            self._nextTime = 0
        elif self._saveTpnts is not None and len(self._saveTpnts) > 0:
            self._saveTind = 0
            self._nextTime = self._saveTpnts[0]
        else:
            self._saveTind = None
            self._nextTime = math.inf

    def save(self):
        """Trigger saving of the result selector at the current simulation time

        Most saving should be done automatically by providing a ``dt`` or a ``timePoints`` argument
        to the :py:func:`steps.API_2.sim.Simulation.toSave` method but it is possible to manually
        decide when to save data by calling ``save()`` on a result selector during simulation.

        Usage::

            for r in range(NBRUNS):
                for t in timePoints:
                    sim.run(t)
                    rs1.save() # Saving values manually
        """
        self._checkAddedToSim()
        self._save(self.sim.Time, (self.sim.Time, self.sim._runId))

    def _toDB(self, dbhanlder):
        """
        Specify that the data should be saved to a database. 
        """
        self._checkComplete()
        self._dataHandler = dbhanlder._getDataHandler(self)

    def _saveWithDt(self, dt):
        """Specify that the data needs to be saved every dt seconds."""
        self._checkComplete()
        self._saveDt = dt
        self._saveTpnts = []
        self._saveTind = 0
        self._nextTime = 0

    def _saveWithTpnts(self, tpnts):
        """Specify at which time points the data should be saved."""
        self._saveTpnts = tpnts
        self._saveDt = None
        self._saveDtStart = None
        self._saveTind = 0
        self._nextTime = self._saveTpnts[0]

    def _addedToSimulation(self, ind, rsGroupInd):
        """Specify that the result selector was added to a simulation with index ind."""
        self._checkComplete()
        self._addedToSim = True
        self._selectorInd = ind
        self._optimGroupInd = rsGroupInd

    def _concat(self, other):
        """Concatenate two result selectors into a _ResultList."""
        return _ResultList([self, other], self.sim)

    def _save(self, t, solvStateId=None):
        """Save the data using self._dataHandler."""
        self._dataHandler.save(t, self._evaluate(solvStateId))
        self._updateNextSaveTime()

    def _updateNextSaveTime(self):
        """Update the time of the next save."""
        if self._saveTind is not None:
            self._saveTind += 1
            if self._saveDt is not None:
                self._nextTime = self._saveTind * self._saveDt
            elif self._saveTind < len(self._saveTpnts):
                self._nextTime = self._saveTpnts[self._saveTind]
            else:
                self._nextTime = math.inf

    def _distribute(self):
        """Distribute the path across MPI ranks if it involves mesh elements of a distributed meshes."""
        return self, False

    def _evaluate(self, solvStateId=None):
        """
        Return a list of the values to save. An optional integer can be given to uniquely
        identify a solver state, this is useful for optimizing solver calls (i.e. not calling
        several times the same thing if the solver state did not change).
        """
        pass

    def _getEvalLen(self):
        """Return the number of values that _evaluate() will return."""
        pass

    def __getattr__(self, name):
        """Redirect attribute access to a SimPath

        See :py:func:`steps.API_2.sim.SimPath.__getattr__`.

        .. note::
            This method should not be called explicitely, it is only documented for clarity.

        :meta public:
        """
        try:
            return super().__getattr__(name)
        except AttributeError:
            return getattr(_ResultPath(self.sim), name)

    def _checkAddedToSim(self):
        """Check that the ResultSelector was added to the Simulation."""
        if not self._addedToSim:
            raise Exception(
                f'Cannot access data from a ResultSelector that was not added to a '
                f'simulation with the "toSave" method.'
            )

    def _checkCompatible(self, other):
        """
        Check that 'other' is a ResultSelector that is associated to the same simulation as self
        """
        if not isinstance(other, ResultSelector):
            raise TypeError(f'Cannot combine a ResultSelector with {other}.')
        if self.sim != other.sim:
            raise Exception(f'Cannot combine ResultSelectors associated to different simulations.')
        self._checkComplete()
        other._checkComplete()

    def _checkComplete(self):
        """Raise an exception if the result selector is not complete."""
        raise Exception(f'{self} is not a complete ResultSelector.')

    def _savingStarted(self):
        """Return whether data started being saved."""
        return self._dataHandler._savingStarted()

    def _binaryOp(self, other, op, symetric=False, opStr='{0} {1}'):
        """Return a _ResultCombiner that represents the binary operation op."""
        if isinstance(other, numbers.Number):

            def opFunc(x):
                return [op(v, other) for v in x]

            def lblFunc(i, children):
                allLbl = []
                for c in children:
                    allLbl += c.labels
                return opStr.format(allLbl[i], other)

            return _ResultCombiner(
                opFunc,
                lambda x: x,
                [self],
                self.sim,
                labelFunc=lblFunc,
                strDescr=opStr.format(self._strDescr(), other),
            )
        elif isinstance(other, ResultSelector):
            self._checkCompatible(other)
            if other._getEvalLen() == 1:

                def opFunc(x):
                    return [op(v, x[-1]) for v in x[:-1]]

                otherLbl = other.labels[0]

                def lblFunc(i, children):
                    allLbl = []
                    for c in children:
                        allLbl += c.labels
                    return opStr.format(allLbl[i], otherLbl)

                return _ResultCombiner(
                    opFunc,
                    lambda x: x - 1,
                    [self, other],
                    self.sim,
                    labelFunc=lblFunc,
                    strDescr=opStr.format(self._strDescr(), other._strDescr()),
                )
            elif symetric and self._getEvalLen() == 1:
                return other._binaryOp(self, op, True, opStr)
            elif other._getEvalLen() == self._getEvalLen():
                n = self._getEvalLen()

                def opFunc(x):
                    return [op(a, b) for a, b in zip(x[:n], x[n:])]

                otherLbls = other.labels

                def lblFunc(i, children):
                    allLbl = []
                    for c in children:
                        allLbl += c.labels
                    return opStr.format(allLbl[i], otherLbls[i])

                return _ResultCombiner(
                    opFunc,
                    lambda x: x // 2,
                    [self, other],
                    self.sim,
                    labelFunc=lblFunc,
                    strDescr=opStr.format(self._strDescr(), other._strDescr()),
                )
            else:
                raise Exception(
                    f'Cannot apply binary operation {opStr.format("","")}, '
                    f'incompatible output lengths: "{self}" has an output '
                    f'length of {self._getEvalLen()} while "{other}" has an '
                    f'output length of {other._getEvalLen()}.'
                )

        else:
            raise TypeError(f'Cannot combine a resultSelector with {other} using {op}.')

    def __lshift__(self, other):
        """Concatenate two result selectors with the ``<<`` operator

        :param other: The other result selector
        :type other: :py:class:`ResultSelector`
        :returns: The result selector resulting from the concatenation of both operands. Its
            length is thus the sum of both of the operands' lengths.
        :rtype: :py:class:`ResultSelector`

        Usage::

            rs2 = rs.comp1.S1.Count << rs.comp1.S2.Count # rs2 will save 2 values, the count of S1
                                                         # in comp1 and the count of S2 in comp1.

        :meta public:
        """
        self._checkCompatible(other)
        return self._concat(other)

    def __mul__(self, other):
        """Multiply result selectors with the * operator

        :param other: The other result selector or a number
        :type other: Union[:py:class:`ResultSelector`, float]
        :returns: The result selector resulting from the multiplication of both operands. If both
            operands are result selectors and have the same size, this corresponds to the
            element-wise product of values. If one of the operand is a number or a result selector
            of length 1, all values of the result selector are multiplied with this single value.
        :rtype: :py:class:`ResultSelector`

        Usage::

            rs3 = 10 * rs.TETS(tetlst).S1.Count                       # rs3 will save the number
                                                                      # of S1 in each tetrahedron
                                                                      # in tetlst, multiplied by
                                                                      # 10.

            rs4 = rs.TETS(tetlst).S1.Count * rs.TETS(tetlst).S2.Count # rs4 will save the product
                                                                      # of the number of S1 and
                                                                      # the number of S2 in each
                                                                      # tetrahedron in tetLst.

        :meta public:
        """
        return self._binaryOp(other, operator.mul, symetric=True, opStr='({0} * {1})')

    def __rmul__(self, other):
        return self._binaryOp(other, operator.mul, symetric=True, opStr='({1} * {0})')

    def __truediv__(self, other):
        """Divide result selectors with the ``/`` operator

        :param other: The other result selector or a number
        :type other: Union[:py:class:`ResultSelector`, float]
        :returns: The result selector resulting from the division of both operands. If both
            operands are result selectors and have the same size, this corresponds to the
            element-wise division of values. If one of the operand is a number or a result selector
            of length 1, all values of the result selectors are divided by this single value (or
            this single value is divided by all values from the result selector, depending on
            order).
        :rtype: :py:class:`ResultSelector`

        Usage::

            rs3 = rs.TETS(tetlst).S1.Count / 10                       # rs3 will save the number
                                                                      # of S1 in each tetrahedron
                                                                      # in tetlst, divided by 10.

            rs4 = 1 / rs.TETS(tetlst).S1.Count                        # rs4 will save the inverse
                                                                      # of the number of S1 in
                                                                      # each tetrahedron in
                                                                      # tetlst, divided by 10.

            rs5 = rs.TETS(tetlst).S1.Count / rs.TETS(tetlst).S2.Count # rs5 will save the ratio of
                                                                      # S1 to S2 in each
                                                                      # tetrahedron in tetLst.

        :meta public:
        """
        return self._binaryOp(other, operator.truediv, symetric=False, opStr='({0} / {1})')

    def __rtruediv__(self, other):
        return self._binaryOp(other, lambda a, b: b / a, symetric=False, opStr='({1} / {0})')

    def __add__(self, other):
        """Add result selectors with the ``+`` operator

        :param other: The other result selector or a number
        :type other: Union[:py:class:`ResultSelector`, float]
        :returns: The result selector resulting from the addition of both operands. If both
            operands are result selectors and have the same size, this corresponds to the
            element-wise addition of values. If one of the operand is a number or a result selector
            of length 1, this single value is added to all values of the result selector.
        :rtype: :py:class:`ResultSelector`

        Usage::

            rs3 = 10 + rs.TETS(tetlst).S1.Count                       # rs3 will save the number
                                                                      # of S1 in each tetrahedron
                                                                      # in tetlst, increased by
                                                                      # 10.

            rs4 = rs.TETS(tetlst).S1.Count + rs.TETS(tetlst).S2.Count # rs4 will save the sum
                                                                      # of the number of S1 and
                                                                      # the number of S2 in each
                                                                      # tetrahedron in tetLst.

        :meta public:
        """
        return self._binaryOp(other, operator.add, symetric=True, opStr='({0} + {1})')

    def __radd__(self, other):
        return self._binaryOp(other, operator.add, symetric=True, opStr='({1} + {0})')

    def __sub__(self, other):
        """Subtract result selectors with the ``-`` operator

        :param other: The other result selector or a number
        :type other: Union[:py:class:`ResultSelector`, float]
        :returns: The result selector resulting from the subtraction of both operands. If both
            operands are result selectors and have the same size, this corresponds to the
            element-wise subtraction of values. If one of the operand is a number or a result selector
            of length 1, this single value is subtracted from all values of the result selectors
            (or each value from the result selector is subtracted from the single value, depending
            on order).
        :rtype: :py:class:`ResultSelector`

        Usage::

            rs3 = rs.TETS(tetlst).S1.Count - 10                       # rs3 will save the number
                                                                      # of S1 in each tetrahedron
                                                                      # in tetlst, minus 10.

            rs4 = 10 - rs.TETS(tetlst).S1.Count                       # rs4 will save 10 minus the
                                                                      # number of S1 for each
                                                                      # tetrahedron in tetlst.

            rs5 = rs.TETS(tetlst).S1.Count - rs.TETS(tetlst).S2.Count # rs5 will save the number
                                                                      # of S1 minus the number of
                                                                      # S2 in each tetrahedron in
                                                                      # tetLst.

        :meta public:
        """
        return self._binaryOp(other, operator.sub, symetric=False, opStr='({0} - {1})')

    def __rsub__(self, other):
        return self._binaryOp(other, lambda a, b: b - a, symetric=False, opStr='({1} - {0})')

    def __pow__(self, other):
        """Exponentiate result selectors with the ** operator

        :param other: The other result selector or a number
        :type other: Union[:py:class:`ResultSelector`, float]
        :returns: The result selector resulting from the exponentiation of both operands. If both
            operands are result selectors and have the same size, this corresponds to the
            element-wise exponentiation of values. If one of the operand is a number or a result
            selector of length 1, this single value is exponentiated by each value of the result
            selector (or each value in the result selector is exponentiated by the single value,
            depending on order).
        :rtype: :py:class:`ResultSelector`

        Usage::

            rs3 = rs.TETS(tetlst).S1.Count ** 2                       # rs3 will save the square
                                                                      # of the number of S1 in
                                                                      # each tetrahedron in
                                                                      # tetlst.

        :meta public:
        """
        return self._binaryOp(other, operator.pow, symetric=False, opStr='({0} ** {1})')

    # Needed for the heapq ordering in Simulation
    def __lt__(self, other):
        return True

    @classmethod
    def SUM(cls, sel):
        """Sum of all values from a result selector

        :param sel: Result selector whose values should be summed
        :type sel: :py:class:`ResultSelector`

        :returns: A result selector with a single value that corresponds to the sum of the values
            in ``sel``.
        :rtype: :py:class:`ResultSelector`

        Usage::

            rs3 = rs.SUM(rs.TETS(tetLst).S1.Count) # The total number of S1 in tetLst
        """
        return _ResultCombiner(
            lambda x: [sum(x)],
            lambda x: 1,
            [sel],
            sel.sim,
            labelFunc=lambda _, d: f"SUM({' + '.join(c._strDescr() for c in d)})",
            strDescr=f'SUM({sel._strDescr()})',
        )

    @classmethod
    def MIN(cls, sel):
        """Minimum of all values from a result selector

        :param sel: Result selector whose values should be used
        :type sel: :py:class:`ResultSelector`

        :returns: A result selector with a single value that corresponds to the minimum of the values
            in ``sel``.
        :rtype: :py:class:`ResultSelector`

        Usage::

            rs3 = rs.MIN(rs.TETS(tetLst).S1.Count) # The minimum number of S1 in tetLst
        """
        return _ResultCombiner(
            lambda x: [min(x)],
            lambda x: 1,
            [sel],
            sel.sim,
            labelFunc=lambda _, d: f"MIN({' + '.join(c._strDescr() for c in d)})",
            strDescr=f'MIN({sel._strDescr()})',
        )

    @classmethod
    def MAX(cls, sel):
        """Maximum of all values from a result selector

        :param sel: Result selector whose values should be used
        :type sel: :py:class:`ResultSelector`

        :returns: A result selector with a single value that corresponds to the maximum of the values
            in ``sel``.
        :rtype: :py:class:`ResultSelector`

        Usage::

            rs3 = rs.MAX(rs.TETS(tetLst).S1.Count) # The maximum number of S1 in tetLst
        """
        return _ResultCombiner(
            lambda x: [max(x)],
            lambda x: 1,
            [sel],
            sel.sim,
            labelFunc=lambda _, d: f"MAX({' + '.join(c._strDescr() for c in d)})",
            strDescr=f'MAX({sel._strDescr()})',
        )


class _ResultPath(ResultSelector):
    """
    Represents a SimPath to be saved during simulation.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.simpath = nsim.SimPath(self.sim)
        self._descriptor = None
        self._compiledFuncs = None
        self._len = None
        self._descr = [str(self.sim)]

        # When distributing result paths, the simpath can cover more data than needs to be saved.
        # We need to call non-distributed methods on all ranks but we only need to save the results
        # on rank 0.
        # _simPathMask defines a boolean mask determining which parts of the simpath that should be
        # saved locally, the other parts are not saved and the corresponding labels and metadata are
        # not saved either.
        self._simPathMask = None

    def _strDescr(self):
        """Return a generic description of the ResultSelector."""
        return '.'.join(self._descr[1:])

    def _evaluate(self, solvStateId=None):
        """Return a list of the values to save."""
        if self._compiledFuncs is None:
            self._compiledFuncs = self._descriptor._getFinalPathsFunction(self.simpath)
        if self._simPathMask is None:
            return [f(*args, **kwargs) for f, args, kwargs in self._compiledFuncs]
        else:
            res = []
            for i, (f, args, kwargs) in enumerate(self._compiledFuncs):
                val = f(*args, **kwargs)
                if self._simPathMask[i]:
                    res.append(val)
            return res

    def _getEvalLen(self):
        """Return the number of values that _evaluate() will return."""
        return self._len

    def _checkComplete(self):
        """Raise an exception if the path is not complete."""
        if self._descriptor is None:
            raise Exception(f'{self} is incomplete.')

    def _concat(self, other):
        """Concatenate two result selectors into a _ResultList."""
        self._checkComplete()
        return super()._concat(other)

    def _binaryOp(self, other, op, symetric=False, opStr='{0} {1}'):
        """Return a _ResultCombiner that represents the binary operation op."""
        self._checkComplete()
        return super()._binaryOp(other, op, symetric, opStr)

    def _distribute(self):
        """Distribute the path across MPI ranks if it involves mesh elements of a distributed meshes."""
        if self._distrInds is not None:
            return self, False

        self._fullLen = self._getEvalLen()
        self.simpath, self._distrInds, spMask, changed = self.simpath._distribute()
        if changed:
            if self.simpath is None:
                return None, changed
            self._simPathMask = numpy.array(spMask)
            # Recompute length, labels and metadata
            self._finalize()
        return self, changed

    def _finalize(self):
        """Finalize the result path by computing length, labels, and metadata"""
        self._len = len([p for p in self.simpath])

        # Compute labels
        self._labels = []
        # for descr in self.simpath._getDescriptions(tuple(self._descr)):
        *_descr, endName = self._descr
        for descr in self.simpath._getDescriptions(tuple(_descr)):
            self._labels.append('.'.join(descr[1:] + (endName,)))

        # Compute automatic metadata
        mtdt = {}
        for i, path in enumerate(self.simpath._walk()):
            for key, lst in mtdt.items():
                lst.append(None)
            if isinstance(path, nutils.SimPathCombiner):
                # If the path is a combination of paths, only consider metadata that is
                # common to all paths
                dct = {}
                for elem in path.paths[0]:
                    dct.update(elem._simPathAutoMetaData())
                for p in path.paths[1:]:
                    dct2 = {}
                    for elem in p:
                        dct2.update(elem._simPathAutoMetaData().items())
                    dct = {key: dct[key] for key in dct.keys() & dct2.keys() if dct[key] == dct2[key]}
                for key, val in dct.items():
                    if key not in mtdt:
                        mtdt[key] = [None] * (i + 1)
                    mtdt[key][i] = val
            else:
                # Otherwise save all metadata for the path
                for elem in path:
                    for key, val in elem._simPathAutoMetaData().items():
                        if key not in mtdt:
                            mtdt[key] = [None] * (i + 1)
                        mtdt[key][i] = val

        # Restrict data to simPathMask, if applicable
        if self._simPathMask is not None:
            self._len = sum(self._simPathMask)
            self._labels = [lbl for lbl, msk in zip(self._labels, self._simPathMask) if msk]
            mtdt = {key:[v for v, msk in zip(lst, self._simPathMask) if msk] for key, lst in mtdt.items()}

        # Set the metadata
        for key, lst in mtdt.items():
            self.metaData.__setitem__(key, lst, _internal=True)

    def __getattr__(self, name):
        self._descr.append(name)
        if name not in nsim.SimPath._endNames:
            self.simpath = getattr(self.simpath, name)
        else:
            self._descriptor = getattr(nsim.SimPath, name)
            self._endName = name
            self._finalize()
        return self

    def __call__(self, *args, **kwargs):
        self.simpath = self.simpath(*args, **kwargs)
        self._descr[-1] += f'({nutils.args2str(*args, **kwargs)})'
        return self

    def __getitem__(self, key):
        self.simpath = self.simpath[key]
        self._descr[-1] += f'[{nutils.key2str(key)}]'
        return self

    def __repr__(self):
        return self._strDescr()


class _ResultList(ResultSelector):
    """Represents the concatenation of several ResultSelectors."""

    def __init__(self, lst, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.children = lst

        self._finalize()

        # Do not initialize metadata directly
        self._metaData = None

    def _strDescr(self):
        """Return a generic description of the ResultSelector."""
        return ', '.join(c._strDescr() for c in self.children)

    def _evaluate(self, solvStateId=None):
        """Return a list of the values to save."""
        res = []
        for c in self.children:
            res += list(c._evaluate(solvStateId))
        return res

    def _getEvalLen(self):
        """Return the number of values that _evaluate() will return."""
        return self._evalLen

    def _concat(self, other):
        """Concatenate two result selectors into a _ResultList."""
        if other.__class__ is _ResultList:
            return _ResultList(self.children + other.children, self.sim)
        else:
            return _ResultList(self.children + [other], self.sim)

    def _checkComplete(self):
        """Raise an exception if the result selector is not complete."""
        for c in self.children:
            c._checkComplete()

    def _finalize(self):
        """Compute labels and evel length"""
        self._evalLen = 0
        for c in self.children:
            self._evalLen += c._getEvalLen()

        # concatenate labels
        self._labels = []
        for c in self.children:
            self._labels += c.labels

    def _distribute(self):
        """Distribute the path across MPI ranks if it involves mesh elements of a distributed meshes."""
        if self._distrInds is not None:
            return self, False

        self._fullLen = self._getEvalLen()
        newLst = []
        globalIdx = 0
        changed = False
        self._distrInds = []
        for c in self.children:
            totLen = c._getEvalLen()
            nc, cChanged = c._distribute()

            if nc is not None:
                newLst.append(nc)
                distrInds = c._distrInds if c._distrInds is not None else range(c._getEvalLen())
                self._distrInds += [globalIdx + idx for idx in distrInds]

            changed |= cChanged
            globalIdx += totLen

        self.children = newLst
        self._finalize()
        # Update metaData
        for key, vals in self.metaData.items():
            if len(vals) != self._getEvalLen():
                self._metaData.__setitem__(key, [vals[inds] for ind in self._distrInds], _internal=True)

        return self, changed

    def _computeMetaData(self):
        """Compute the concatenation of metadata."""
        metaDataKeys = set()
        for c in self.children:
            metaDataKeys.update(c.metaData.keys())
        mtdt = {}
        for key in metaDataKeys:
            lst = []
            for c in self.children:
                try:
                    lst += list(c.metaData[key])
                except KeyError:
                    lst += [None] * c._getEvalLen()
            mtdt[key] = lst
        return mtdt

    def _getAllTerminalChildren(self):
        """Return all children that are not ResultLists recursively"""
        for c in self.children:
            if isinstance(c, _ResultList):
                yield from c._getAllTerminalChildren()
            else:
                yield c

    @ResultSelector.metaData.getter
    def metaData(self):
        if self._metaData is None:
            self._metaData = _MetaData(self)
            mtdt = self._computeMetaData()
            for key, lst in mtdt.items():
                self._metaData._dict[key] = lst
        return self._metaData

    def __repr__(self):
        return self._strDescr()


class _ResultCombiner(_ResultList):
    """
    Transforms results using function func that takes an iterable and outputs a list.
    function lenFunc should take the length of children output as an argument and return the
    length of the combiner output.
    """

    def __init__(self, func, lenFunc, *args, labelFunc=None, strDescr=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.func = func
        self._len = lenFunc(super()._getEvalLen())
        self._labelFunc = labelFunc if labelFunc is not None else lambda a, b: self.func

        if strDescr is not None:
            self._descr = strDescr
        else:
            self._descr = f"{self.func}({super()._strDescr()})"

        self._labels = [self._labelFunc(i, self.children) for i in range(self._len)]

    def _concat(self, other):
        """Concatenate two result selectors into a _ResultList."""
        return _ResultList([self, other], self.sim)

    def _strDescr(self):
        """Return a generic description of the ResultSelector."""
        return self._descr

    def _evaluate(self, solvStateId=None):
        """Return a list of the values to save."""
        return self.func(super()._evaluate(solvStateId))

    def _getEvalLen(self):
        """Return the number of values that _evaluate() will return."""
        return self._len

    def _distribute(self):
        """Distribute the path across MPI ranks if it involves mesh elements of a distributed meshes."""
        if self._distrInds is not None:
            return self, False

        self._fullLen = self._getEvalLen()
        if not nsim.MPI._shouldWrite:
            self._len = 0
            self._labels = []
            self._distrInds = []
            if self._metaData is not None:
                self._metaData._clear()
            return self, True
        else:
            self._distrInds = list(range(self._getEvalLen()))
            return self, False

    @ResultSelector.metaData.getter
    def metaData(self):
        if self._metaData is None:
            self._metaData = _MetaData(self)
            mtdt = self._computeMetaData()
            # Only keep common metadata
            for key, lst in mtdt.items():
                if len(set(lst)) == 1:
                    self._metaData._dict[key] = [lst[0]] * self._getEvalLen()
        return self._metaData

    def __repr__(self):
        return self._strDescr()


###################################################################################################
# Read only ResultSelectors


class _ReadOnlyResultSelector:
    """
    Only implement data access methods of ResultSelector
    """

    def __init__(self, handler):
        self._dataHandler = handler

    @property
    def time(self):
        """Return an accessor to the timepoints data."""
        return self._dataHandler.time()

    @property
    def data(self):
        """Return an accessor to the saved data."""
        return self._dataHandler.data()

    @property
    def labels(self):
        """Return a list of strings describing the things being saved."""
        return self._dataHandler.labels()

    @property
    def metaData(self):
        """Return the metadata associated to the ResultSelector."""
        return self._dataHandler.metaData()


###################################################################################################
# Data handlers


class _DataHandler:
    """
    Interface for data saving classes.
    """

    def __init__(self, parent, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._parent = parent
        self._runId = -1

    def time(self):
        """Return an accessor to the timepoints data."""
        pass

    def data(self):
        """Return an accessor to the saved data."""
        pass

    def labels(self):
        """Return the labels of the data being saved."""
        raise NotImplementedError()

    def _newRun(self):
        """Signal that a new run of the simulation started."""

        self._runId += 1

    def save(self, t, row):
        """Save the data."""
        pass

    def _savingStarted(self):
        """Return whether data started being saved."""
        return self._runId >= 0

    @classmethod
    def _checkCanAccess(cls):
        if not nsim.MPI._shouldWrite:
            raise Exception(f'Cannot access ResultSelector data out of the rank 0 process while using MPI.')


class _MemoryDataHandler(_DataHandler):
    """
    Data handler for saving data in memory.
    """

    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self.saveData = []
        self.saveTime = []

    def time(self):
        """Return an accessor to the timepoints data."""
        return _MemoryDataAccessor(self.saveTime, 2)

    def data(self):
        """Return an accessor to the saved data."""
        return _MemoryDataAccessor(self.saveData, 3)

    def _newRun(self):
        """Signal that a new run of the simulation started."""
        super()._newRun()
        self.saveData.append([])
        self.saveTime.append([])

    def save(self, t, row):
        """Save the data."""
        self.saveTime[-1].append(t)
        self.saveData[-1].append(copy.copy(row))


class _FileDataHandler(_DataHandler, nutils.Versioned):
    """
    Data handler for saving data to files.
    """

    HEADER_FORMAT = '>QQQQ'
    DATA_FORMAT = '>d'
    DEFAULT_BUFFER_SIZE = 4096
    INT_SIZE = 4

    FILE_FORMAT_STR = '__steps_version__'
    FILE_FORMAT_OLDEST_VERSION = '3.6.0'

    def __init__(self, parent, path, evalLen=None, buffering=-1, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self._savePath = path
        self._readOnly = evalLen is None
        self._evalLen = evalLen if evalLen is not None else 1
        self._saveFile = None
        self._saveBuffering = buffering
        self._fileHeaderInfo = None
        self._filePrevPos = None

        self.saveData = collections.deque([], self._getDequeMaxSize())
        self.saveTime = []

        self._labels = None
        self._metaData = None
        # TODO Not urgent: make labels and metadata readonly
        self._labelEndPos = None

        # If we are reading from a file, we need to set the version
        if self._readOnly:
            version = self.metaData()._dict.get(
                _FileDataHandler.FILE_FORMAT_STR,
                _FileDataHandler.FILE_FORMAT_OLDEST_VERSION
            )
            self._setVersion(version)

    def __del__(self):
        if hasattr(self, '_saveFile') and self._saveFile is not None:
            self._finalizeFile()

    def time(self):
        """Return an accessor to the timepoints data."""
        self._checkCanAccess()
        self._finalizeFile()
        return _FileDataAccessor(self._savePath, parent=self, saveTime=True)

    def data(self):
        """Return an accessor to the saved data."""
        self._checkCanAccess()
        self._finalizeFile()
        return _FileDataAccessor(self._savePath, parent=self, saveTime=False)

    def labels(self):
        """Return the labels of the data being saved."""
        self._checkCanAccess()
        if self._labels is None:
            self._readLabelsAndMetaData()
        return self._labels

    def metaData(self):
        """Return the metaData of the data being saved."""
        self._checkCanAccess()
        if self._metaData is None:
            self._readLabelsAndMetaData()
        md = _MetaData(None)
        md._dict = self._metaData
        return md

    @property
    def _dataStartPos(self):
        if self._labelEndPos is None:
            self.labels()
        return self._labelEndPos

    def _newRun(self):
        """Signal that a new run of the simulation started."""
        super()._newRun()
        self.saveData.clear()
        self.saveTime.append([])
        if nsim.MPI._shouldWrite:
            self._writeRunHeader(self._runId, 0, 1 + self._evalLen)

    def save(self, t, row):
        """Save the data."""
        self.saveTime[-1].append(t)
        self.saveData.append(list(row))
        if nsim.MPI._shouldWrite:
            self._writeToFile(t, self.saveData[-1])

    def _openFile(self):
        """Open the file in the correct mode."""
        if self._saveFile is None:
            if self._fileHeaderInfo is None:
                self._saveFile = open(self._savePath, 'wb', buffering=self._saveBuffering)
            else:
                self._saveFile = open(self._savePath, 'r+b', buffering=self._saveBuffering)
                self._saveFile.seek(0, 2)
        return self._saveFile

    def _writeRunHeader(self, runId, nbRows, nbCols, writeNext=True):
        """Write the header line of a run."""
        self._openFile()
        if self._fileHeaderInfo is not None:
            nxtPos = self._saveFile.seek(0, 1)
            self._saveFile.seek(self._filePrevPos, 0)
            if writeNext:
                self._fileHeaderInfo[3] = nxtPos
            self._saveFile.write(struct.pack(_FileDataHandler.HEADER_FORMAT, *self._fileHeaderInfo))
            self._saveFile.seek(nxtPos, 0)
            self._saveFile.flush()
        else:
            self._writeLabelsAndMetaData()

        if writeNext:
            self._fileHeaderInfo = [runId, nbRows, nbCols, 0]
            self._filePrevPos = self._saveFile.seek(0, 1)
            self._saveFile.write(struct.pack(_FileDataHandler.HEADER_FORMAT, *self._fileHeaderInfo))

    def _writeLabelsAndMetaData(self):
        """Write the labels and the metadata to the file header."""
        # Labels
        lbls = self._parent.labels
        self._writeInt(len(lbls))
        for l in lbls:
            self._writeStr(l)

        # MetaData
        mtdt = self._parent.metaData._dict
        if _FileDataHandler.FILE_FORMAT_STR in mtdt:
            raise Exception(
                f'The metadata contains the reserved key name "{_FileDataHandler.FILE_FORMAT_STR}"'
            )
        mtdt[_FileDataHandler.FILE_FORMAT_STR] = steps.__version__

        data = pickle.dumps(mtdt)
        self._writeInt(len(data))
        self._saveFile.write(data)

    def _writeInt(self, i):
        """Write int i to the binary file."""
        self._saveFile.write(i.to_bytes(_FileDataHandler.INT_SIZE, byteorder='big'))

    def _writeStr(self, s):
        """Write string s to the binary file."""
        bs = s.encode('ascii')
        self._writeInt(len(bs))
        self._saveFile.write(bs)

    @staticmethod
    def _readInt(f):
        """Read an int from binary file f."""
        return int.from_bytes(f.read(_FileDataHandler.INT_SIZE), byteorder='big')

    @staticmethod
    def _readStr(f):
        """Read a string from binary file f."""
        strLen = _FileDataHandler._readInt(f)
        return f.read(strLen).decode('ascii')

    def _readLabelsAndMetaData(self):
        """Open the file and read labels."""
        with open(self._savePath, 'rb') as f:
            # Labels
            nbLbls = _FileDataHandler._readInt(f)
            self._labels = []
            for i in range(nbLbls):
                self._labels.append(_FileDataHandler._readStr(f))

            mtdtSz = _FileDataHandler._readInt(f)
            # TODO Not urgent: make the dict readonly
            self._metaData = pickle.loads(f.read(mtdtSz))

            self._labelEndPos = f.seek(0, 1)

    @nutils.Versioned._versionRange(belowOrEq=FILE_FORMAT_OLDEST_VERSION)
    def _writeToFile(self, t, vals):
        """Write the data to file."""
        self._openFile()
        self._saveFile.write(struct.pack('>d' + 'd' * len(vals), t, *vals))
        self._fileHeaderInfo[1] += 1

    @nutils.Versioned._versionRange(above=FILE_FORMAT_OLDEST_VERSION)
    def _writeToFile(self, t, vals):
        """Write the data to file."""
        self._openFile()
        pickle.dump((t, vals), self._saveFile)
        self._fileHeaderInfo[1] += 1

    def _finalizeFile(self):
        """Flush the file buffer and close the file."""
        # Only write things if the result selector was created from a simulation, not a file path
        if nsim.MPI._shouldWrite and not self._readOnly:
            self._writeRunHeader(None, None, None, writeNext=False)
            self._saveFile.close()
            self._saveFile = None

    def _getDequeMaxSize(self):
        """Return the length of the buffer deque."""
        if self._saveBuffering != -1:
            buf = self._saveBuffering
        else:
            buf = _FileDataHandler.DEFAULT_BUFFER_SIZE
        return max(1, buf // self._evalLen)


class _DBDataHandler(_DataHandler):
    pass


class _SQLiteDataHandler(_DBDataHandler):
    """
    Data handler for saving to sqlite db file.
    """

    TABLE_NAME_TEMPLATE = 'Group_{}_Selector_{}'
    COLUMN_NAME_TEMPLATE = 'Col_{} real'

    def __init__(
        self, parent, dbh, commitFreq, *args, groupId=None, rsid=None, tableName=None, nbCols=None, **kwargs
    ):
        super().__init__(parent, *args, **kwargs)
        self._dbh = dbh
        self._conn = dbh._conn
        self._commitFreq = commitFreq
        self._commitInd = 0

        self._initialized = False

        self._groupId = groupId
        self._rsid = rsid
        self._tableName = tableName
        self._nbCols = nbCols

        self._labels = None
        self._metaData = None

    def time(self):
        """Return an accessor to the timepoints data."""
        self._checkCanAccess()
        return _SQLiteDataAccessor(
            self._dbh, self._groupId, self._rsid, self._tableName, self._nbCols, saveTime=True
        )

    def data(self):
        """Return an accessor to the saved data."""
        self._checkCanAccess()
        return _SQLiteDataAccessor(
            self._dbh, self._groupId, self._rsid, self._tableName, self._nbCols, saveTime=False
        )

    def labels(self):
        """Return the labels of the saved data."""
        self._checkCanAccess()
        if self._labels is None:
            self._labels = self._dbh._labelsQuerry(self._groupId, self._rsid)
        return self._labels

    def metaData(self):
        """Return the metadata of the saved data."""
        self._checkCanAccess()
        if self._metaData is None:
            self._metaData = self._dbh._metaDataQuerry(self._groupId, self._rsid)
        return self._metaData

    def _initialize(self):
        """
        Create the table and initialize everything. Should only be called after the first newRun.
        """
        lbls = self._parent.labels
        self._groupId = self._dbh._groupId
        self._rsid = self._parent._selectorInd
        colStr = ', '.join(_SQLiteDataHandler.COLUMN_NAME_TEMPLATE.format(i) for i in range(len(lbls)))

        self._nbCols = len(lbls)
        self._tableName = _SQLiteDataHandler.TABLE_NAME_TEMPLATE.format(self._groupId, self._rsid)
        self._insertStr = f"INSERT INTO {self._tableName} VALUES ({','.join('?'*(2+self._nbCols))});"

        # Check if the table already exists
        rows = self._conn.execute(
            f"SELECT name FROM sqlite_master WHERE type='table' AND name='{self._tableName}'"
        ).fetchall()
        if len(rows) == 0:
            # Create table
            self._conn.execute(f'CREATE TABLE {self._tableName} (runid int, time real, {colStr});')
            # Add table info to main table
            self._conn.execute(
                f'INSERT INTO {SQLiteDBHandler._RS_MAIN_TABLE_NAME} VALUES (?,?,?,?,?);',
                (self._groupId, self._rsid, self._parent._strDescr(), self._tableName, self._nbCols),
            )
            # Add labels
            self._conn.executemany(
                f'INSERT INTO {SQLiteDBHandler._RS_LABEL_TABLE_NAME} VALUES (?,?,?,?);',
                [(self._groupId, self._rsid, i, lbl) for i, lbl in enumerate(lbls)],
            )
            # Add MetaData
            self._conn.execute(
                f'INSERT INTO {SQLiteDBHandler._RS_META_DATA_TABLE_NAME} VALUES (?,?,?);',
                (self._groupId, self._rsid, pickle.dumps(self._parent.metaData._dict)),
            )
        else:
            # Initialize the runId to the last recorded one
            rid = self._conn.execute(f'SELECT MAX(runid) FROM {self._tableName}').fetchone()[0]
            if rid is not None:
                self._runId = rid

        self._conn.commit()
        self._cursor = self._conn.cursor()
        self._initialized = True

    def _newRun(self):
        """Signal that a new run of the simulation started."""
        if not self._initialized and nsim.MPI._shouldWrite:
            self._initialize()
        super()._newRun()

    def save(self, t, row):
        """Save the data."""
        if nsim.MPI._shouldWrite:
            self._cursor.execute(self._insertStr, (self._runId, t) + tuple(row))
            self._commitInd += 1
            if self._commitInd % self._commitFreq == 0:
                self._conn.commit()


class _HDF5DataHandler(_DBDataHandler, nutils.Versioned):
    """
    Data handler for saving to HDF5 file.
    """

    _RS_COLREMAPPING_NAME = 'ColumnRemapping'
    _LABELS_DSET_NAME = 'labels'
    _METADATA_GROUP_NAME = 'metaData'
    _RUNS_GROUP_NAME = 'runs'
    _DATA_DSET_NAME = 'data'
    _TIME_DSET_NAME = 'time'
    _RUN_GROUP_TEMPLATE = 'Run_{}'

    def __init__(self, dbh, parent, group, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self._dbh = dbh
        self._group = group
        self._initialized = False
        self._timeInd = None

        self._time = None
        self._data = None

        # Vector representing the permutation that should be applied before saving the data to file
        # It is used by XDMF data handler to ensure contiguous blocks of data.
        self._colRemapping = None
        self._revColRemapping = None

        if parent is None:
            # Load column remapping if we are reading data
            self._loadColumnRemap()

    def time(self):
        """Return an accessor to the timepoints data."""
        self._checkCanAccess()
        return _HDF5DataAccessor(self, True)

    def data(self):
        """Return an accessor to the saved data."""
        self._checkCanAccess()
        return _HDF5DataAccessor(self, False)

    def labels(self):
        """Return the labels of the data being saved."""
        self._checkCanAccess()
        return [lbl.decode('utf-8') for lbl in self._group[_HDF5DataHandler._LABELS_DSET_NAME]]

    def metaData(self):
        """Return the metadata of the saved data."""
        self._checkCanAccess()
        res = {}
        for name, dset in self._group[_HDF5DataHandler._METADATA_GROUP_NAME].items():
            if isinstance(dset[0], bytes):
                res[name] = [v.decode('utf-8') for v in dset]
                # Try to convert strings to numbers, if possible
                for i, v in enumerate(res[name]):
                    for tpe in [int, float]:
                        try:
                            res[name][i] = tpe(v)
                            break
                        except ValueError:
                            pass
                # Restore Nones
                res[name] = [None if v == 'None' else v for v in res[name]]
            else:
                res[name] = [None if numpy.isnan(v) else v for v in dset]
        return res

    def _checkCanAccess(self):
        if not self._dbh._shouldWrite:
            raise Exception(
                f'Cannot access HDF5 data out of the rank 0 process while using non-distributed '
                f'simulation with MPI.'
            )

    def _initialize(self):
        """
        Create the subgroups and initialize everything. Should only be called after the first newRun.
        """
        dskwargs = self._dbh._dataSetKWargs
        if _HDF5DataHandler._RUNS_GROUP_NAME not in self._group:
            self._group.create_group(_HDF5DataHandler._RUNS_GROUP_NAME, track_order=True)
        else:
            # Initialize the runId to the last recorded one
            self._runId = len(self._group[_HDF5DataHandler._RUNS_GROUP_NAME]) - 1
        if _HDF5DataHandler._LABELS_DSET_NAME not in self._group:
            self._group.create_dataset(
                _HDF5DataHandler._LABELS_DSET_NAME,
                data=[lbl.encode('utf-8') for lbl in self._parent.labels],
                **dskwargs
            )
        if _HDF5DataHandler._METADATA_GROUP_NAME not in self._group:
            mtdtGroup = self._group.create_group(_HDF5DataHandler._METADATA_GROUP_NAME)
            for key, vals in self._parent.metaData._dict.items():
                # Try to handle different types
                if any(isinstance(v, str) for v in vals):
                    vals = [v.encode('utf-8') for v in map(str, vals)]
                elif any(isinstance(v, numbers.Number) or v is None for v in vals):
                    vals = [numpy.nan if v is None else v for v in vals]
                elif len(vals) > 0:
                    raise TypeError(
                        f'Metadata contains values that are not strings or numbers, they cannot '
                        f'be saved to HDF5 format.'
                    )
                if len(vals) > 0:
                    mtdtGroup.create_dataset(key, data=vals, **dskwargs)
        self._loadColumnRemap()

        self._initialized = True

    def _loadColumnRemap(self):
        """Load Column remapping, if available"""
        if _HDF5DataHandler._RS_COLREMAPPING_NAME in self._group:
            # Load column remapping if it was already saved to file
            colRemapping = numpy.array(self._group[_HDF5DataHandler._RS_COLREMAPPING_NAME])
            if self._colRemapping is not None and any(list(colRemapping != self._colRemapping)):
                raise Exception(
                    f'The column remapping saved in the HDF5 file for result selector '
                    f'{self._parent._selectorInd} ({self._parent}) is different from the one that was computed '
                    f'for this simulation. Try to save your data to a different HDF5 file.'
                )
            self._colRemapping = colRemapping
            # Compute reverse mapping for data reading
            self._revColRemapping = numpy.array([0] * len(self._colRemapping))
            for i, v in enumerate(self._colRemapping):
                self._revColRemapping[v] = i

    def _newRun(self):
        """Signal that a new run of the simulation started."""
        dskwargs = self._dbh._dataSetKWargs
        if not self._initialized and self._group is not None:
            self._initialize()
        self._timeInd = -1
        super()._newRun()

        if self._group is not None:
            n = self._parent._getEvalLen()
            runGroup = self._group[_HDF5DataHandler._RUNS_GROUP_NAME].create_group(
                _HDF5DataHandler._RUN_GROUP_TEMPLATE.format(self._runId)
            )
            self._data = runGroup.create_dataset(
                _HDF5DataHandler._DATA_DSET_NAME, (1, n), maxshape=(None, n), dtype='f', **dskwargs
            )
            self._time = runGroup.create_dataset(
                _HDF5DataHandler._TIME_DSET_NAME, (1,), maxshape=(None,), dtype='f', **dskwargs
            )

    def save(self, t, row):
        """Save the data."""
        if self._group is not None:
            self._timeInd += 1
            if self._timeInd >= self._time.shape[0]:
                self._time.resize(self._timeInd + 1, axis=0)
                self._data.resize(self._timeInd + 1, axis=0)
            self._time[self._timeInd] = t
            if self._colRemapping is None:
                self._data[self._timeInd, :] = row 
            else:
                self._data[self._timeInd, :] = numpy.array(row)[self._colRemapping]


class _HDF5DistribDataHandler(_HDF5DataHandler):
    """
    Data handler for loading several HDF5 files that have been saved in a distributed way
    """
    def __init__(self, hdfGroup, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._hdfGroup = hdfGroup
        self._colMap = numpy.array(self._group[HDF5Handler._RS_DIST_IND_MAP_NAME])
        self._nbCols = self._colMap.shape[1]
        self._rsInd = self._group.attrs[HDF5Handler._RS_INDEX_ATTR]
        self._dbUid = self._hdfGroup.name

        self._lbls = None
        self._mtdt = None
        self._usedRanks = None

        self._setUpRankFiles()

    def _setUpRankFiles(self):
        self._usedRanks = set(self._colMap[0,:])
        for rnk in self._usedRanks:
            # Open HDF5 files
            if rnk not in self._dbh._distribRankDBHs:
                if rnk == nsim.MPI.rank:
                    self._dbh._distribRankDBHs[rnk] = self._dbh
                else:
                    self._dbh._distribRankDBHs[rnk] = HDF5Handler(
                        HDF5Handler._DISTRIBUTED_HDF_SUFFIX.format(self._dbh._pathPrefix, rnk),
                        hdf5FileKwArgs=self._dbh._fileKwArgs,
                        version=self._version
                    )
            # Load read-only result selector
            rnkDBH = self._dbh._distribRankDBHs[rnk]
            rsDict = rnkDBH._distribRS.setdefault(self._dbUid, {})
            if self._rsInd not in rsDict:
                rnkDBH._checkOpenFile()
                rsGroup = rnkDBH._file[self._dbUid][HDF5Handler._RS_GROUP_NAME.format(self._rsInd)]
                rsDict[self._rsInd] = _ReadOnlyResultSelector(
                    _HDF5DataHandler(rnkDBH, None, rsGroup, version=self._version)
                ) 

    def time(self):
        """Return an accessor to the timepoints data."""
        # All times should be the same, can return the first one
        rnk = self._colMap[0, 0]
        return self._dbh._distribRankDBHs[rnk]._distribRS[self._dbUid][self._rsInd].time

    def data(self):
        """Return an accessor to the saved data."""
        return _HDF5DistribDataAccessor(self)

    def labels(self):
        """Return the labels of the data being saved."""
        if self._lbls is None:
            self._lbls = [''] * self._nbCols
            rnk2Lbls = {}
            for rnk in self._usedRanks:
                rnk2Lbls[rnk] = self._dbh._distribRankDBHs[rnk]._distribRS[self._dbUid][self._rsInd].labels
            for i, (rnk, locIdx) in enumerate(self._colMap.T):
                self._lbls[i] = rnk2Lbls[rnk][locIdx]
        return self._lbls

    def metaData(self):
        """Return the metadata of the saved data."""
        if self._mtdt is None:
            self._mtdt = {}
            rnk2Mtdt = {}
            for rnk in self._usedRanks:
                rnk2Mtdt[rnk] = self._dbh._distribRankDBHs[rnk]._distribRS[self._dbUid][self._rsInd].metaData
            for i, (rnk, locIdx) in enumerate(self._colMap.T):
                for key, vals in rnk2Mtdt[rnk].items():
                    self._mtdt.setdefault(key, [None] * self._nbCols)[i] = vals[locIdx]
        return self._mtdt


class _XDMFDataHandler(_HDF5DataHandler):
    """
    Data handler for writing xdmf files for each run while saving data to HDF5 file.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _newRun(self):
        """Signal that a new run of the simulation started."""
        super()._newRun()
        self._dbh._newRun(self._runId)

    def save(self, t, row):
        """Save the data."""
        super().save(t, row)
        self._dbh._newTimeStep(t, self._parent, self._timeInd)


###################################################################################################
# Data accessors


def _sliceData(data, key):
    """
    Slice multidimentional data in nested lists according to key. Use the numpy slicing
    conventions. 'key' should be a tuple that can only contain integers or slices.
    Return the data in nested lists.
    """
    if len(key) == 1:
        return data[key[0]]
    k = key[0]
    if isinstance(k, slice):
        res = []
        for sub in data[k]:
            res.append(_sliceData(sub, key[1:]))
        return res
    else:
        return _sliceData(data[k], key[1:])


class _MemoryDataAccessor:
    """
    Data accessor for _MemoryDataHandler
    """

    def __init__(self, data, nbDims):
        self._data = data
        self._nbDims = nbDims

    def __getitem__(self, key):
        key = nutils.formatKey(key, self._nbDims, forceSz=True)
        return nutils.nparray(_sliceData(self._data, key))

    def __array__(self):
        return nutils.nparray(self._data)

    def __len__(self):
        return len(self._data)


class _FileDataAccessor(nutils.Versioned):
    """
    Data accessor for _FileDataHandler
    """

    HEADER_SIZE = struct.calcsize(_FileDataHandler.HEADER_FORMAT)
    DATA_SIZE = struct.calcsize(_FileDataHandler.DATA_FORMAT)
    DEFAULT_MAXRUNID = sys.maxsize

    class UnexpectedEnd(Exception):
        pass

    # TODO Optimization: save the number of runs and related data and only update it if the file
    # was changed

    def __init__(self, fp, parent, saveTime=False):
        self._fp = fp
        self._saveTime = saveTime
        self._parentHandler = parent
        self._dataStartPos = parent._dataStartPos

        self._file = open(self._fp, 'rb')
        self._fileInfo = {}
        self._nbDims = 2 if saveTime else 3

        self._setVersion(self._parentHandler._version)

    def __del__(self):
        if hasattr(self, '_file') and self._file is not None:
            self._file.close()

    def __len__(self):
        if 'len' not in self._fileInfo:
            pos = self._file.seek(0, 1)
            self._file.seek(self._dataStartPos)
            nb = 0
            try:
                runId, nbRows, nbCols, nxt = struct.unpack(
                    _FileDataHandler.HEADER_FORMAT, self._file.read(_FileDataAccessor.HEADER_SIZE)
                )
                nb += 1
                while nxt != 0:
                    self._file.seek(nxt)
                    runId, nbRows, nbCols, nxt = struct.unpack(
                        _FileDataHandler.HEADER_FORMAT, self._file.read(_FileDataAccessor.HEADER_SIZE)
                    )
                    nb += 1
            except struct.error:
                pass

            self._file.seek(pos)
            self._fileInfo['len'] = nb
        return self._fileInfo['len']

    def __getitem__(self, key, forceArray=False):
        key = nutils.formatKey(key, self._nbDims, forceSz=True)

        # If possible, try to access from memory
        if not self._parentHandler._readOnly:
            if self._saveTime:
                return nutils.nparray(_sliceData(self._parentHandler.saveTime, key))
            elif self._parentHandler._fileHeaderInfo is not None:
                idxs = nutils.getSliceIds(key[0], sz=self._parentHandler._fileHeaderInfo[0] + 1)
                if len(idxs) == 1 and idxs[0] == self._parentHandler._fileHeaderInfo[0]:
                    nbRows = self._parentHandler._fileHeaderInfo[1]
                    lenDeque = len(self._parentHandler.saveData)

                    inds = nutils.getSliceIds(key[1], sz=nbRows)
                    if all(nbRows - lenDeque <= ti < nbRows for ti in inds):
                        res = [self._parentHandler.saveData[ti - (nbRows - lenDeque)] for ti in inds]
                        if forceArray:
                            if isinstance(key[2], slice):
                                return nutils.nparray([[row[key[2]] for row in res]])
                            else:
                                return nutils.nparray([[[row[key[2]]] for row in res]])
                        else:
                            mk = (slice(None) if isinstance(key[1], slice) else 0, key[2])
                            return nutils.nparray(_sliceData(res, mk))

        # Otherwise, read from file
        res = []
        # Find the number of runs first
        nbRuns = len(self)

        if nbRuns == 0:
            raise IndexError(f'Cannot access data, nothing has been written to the file.')

        # Read header
        self._file.seek(self._dataStartPos)
        runId, nbRows, nbCols, nxt = struct.unpack(
            _FileDataHandler.HEADER_FORMAT, self._file.read(_FileDataAccessor.HEADER_SIZE)
        )

        # Iterate through runs
        for ind in nutils.getSliceIds(key[0], sz=nbRuns):
            while runId != ind:
                if nxt == 0:
                    break
                pos = self._file.seek(nxt)
                try:
                    runId, nbRows, nbCols, nxt = struct.unpack(
                        _FileDataHandler.HEADER_FORMAT, self._file.read(_FileDataAccessor.HEADER_SIZE)
                    )
                except struct.error:
                    break
            if runId != ind:
                if isinstance(key[0], numbers.Integral) or key[0].stop is not None:
                    raise IndexError(f'Run {ind} is not in the file.')
                else:
                    break
            # handle the cases in which the file was only partially written and nbRows == 0
            if nxt == 0 and nbRows == 0:
                warnings.warn(
                    f'Run {ind} from file {self._fp} was not correctly written to file, the '
                    f'corresponding data will be partial.'
                )
                nbRows = None
                if (isinstance(key[1], numbers.Integral) and key[1] < 0) or (
                    isinstance(key[1], slice)
                    and (
                        (key[1].start is not None and key[1].start < 0)
                        or (key[1].stop is not None and key[1].stop < 0)
                    )
                ):
                    raise IndexError('Cannot access partially written data using negative indices.')

            if nbRows is None:
                nbRows = _FileDataAccessor.DEFAULT_MAXRUNID
            rowInds = nutils.getSliceIds(key[1], sz=nbRows)
            res.append([])
            # Read actual data
            try:
                for t, line in self._readRows(rowInds, nbCols):
                    if self._saveTime:
                        res[-1].append(t)
                    else:
                        line = line[key[2]] if isinstance(key[2], slice) else [line[key[2]]]
                        res[-1].append(line)
            except _FileDataAccessor.UnexpectedEnd:
                if nxt == 0:
                    break
                else:
                    raise IndexError(
                        f'Could not load time slice {ti} of run {ind} from {self._fp}.'
                        f' The file might be corrupted.'
                    )

        if forceArray:
            return nutils.nparray(res)
        mk = tuple(slice(None) if isinstance(k, slice) else 0 for k in key)
        return nutils.nparray(_sliceData(res, mk))

    @nutils.Versioned._versionRange(belowOrEq=_FileDataHandler.FILE_FORMAT_OLDEST_VERSION)
    def _readRows(self, rowInds, nbCols):
        datFormat = _FileDataHandler.DATA_FORMAT[0] + _FileDataHandler.DATA_FORMAT[1] * nbCols
        pos = self._file.seek(0, 1)
        try:
            for ti in rowInds:
                self._file.seek(pos + ti * nbCols * _FileDataAccessor.DATA_SIZE)
                if self._saveTime:
                    t, *line = struct.unpack(
                        _FileDataHandler.DATA_FORMAT, self._file.read(_FileDataAccessor.DATA_SIZE)
                    )
                else:
                    t, *line = struct.unpack(
                        datFormat, self._file.read(_FileDataAccessor.DATA_SIZE * nbCols)
                    )

                yield t, line
        except (EOFError, struct.error):
            raise _FileDataAccessor.UnexpectedEnd()

    @nutils.Versioned._versionRange(above=_FileDataHandler.FILE_FORMAT_OLDEST_VERSION)
    def _readRows(self, rowInds, nbCols):
        currti = -1
        try:
            for ti in rowInds:
                # Find desired row
                while currti < ti:
                    t, line = pickle.load(self._file)
                    currti += 1
                yield t, line

        except (EOFError, pickle.UnpicklingError):
            raise _FileDataAccessor.UnexpectedEnd()

    def __array__(self):
        return self.__getitem__(slice(None, None, None), forceArray=True)


class _SQLiteDataAccessor:
    """
    Data accessor for SQLite database
    """

    def __init__(self, dbh, groupid, rsid, tabName, nbCols, saveTime=False):
        self._dbh = dbh
        self._groupid = groupid
        self._rsid = rsid
        self._tabName = tabName
        self._nbCols = nbCols
        self._saveTime = saveTime
        self._nbDims = 2 if saveTime else 3

        self._colLst = [_SQLiteDataHandler.COLUMN_NAME_TEMPLATE.format(ci) for ci in range(self._nbCols)]

    def __getitem__(self, key, forceArray=False):
        key = nutils.formatKey(key, self._nbDims, forceSz=True)
        res = []
        for ri in nutils.getSliceIds(key[0], sz=len(self)):
            if self._saveTime:
                timeDat = self._dbh._conn.execute(
                    f'SELECT time FROM {self._tabName} WHERE runid={ri} ORDER BY time'
                ).fetchall()
                res.append([timeDat[i][0] for i in nutils.getSliceIds(key[1], len(timeDat))])
            else:
                colStr = ','.join(self._colLst[i] for i in nutils.getSliceIds(key[2], self._nbCols))
                allDat = self._dbh._conn.execute(
                    f'SELECT {colStr} FROM {self._tabName} WHERE runid={ri} ORDER BY time'
                ).fetchall()
                res.append([list(allDat[i]) for i in nutils.getSliceIds(key[1], len(allDat))])

        if forceArray:
            return nutils.nparray(res)
        mk = tuple(slice(None) if isinstance(k, slice) else 0 for k in key)
        return nutils.nparray(_sliceData(res, mk))

    def __len__(self):
        return self._dbh._conn.execute(f'SELECT MAX(runid) FROM {self._tabName}').fetchone()[0] + 1

    def __array__(self):
        return self.__getitem__(slice(None, None, None), forceArray=True)


class _HDF5DataAccessor(nutils.Versioned):
    """
    Data accessor for HDF5 files
    """

    def __init__(self, handler, saveTime=False, **kwargs):
        super().__init__(**kwargs)
        self._handler = handler
        self._saveTime = saveTime
        self._nbDims = 2 if saveTime else 3

    def __getitem__(self, key, forceArray=False):
        key = nutils.formatKey(key, self._nbDims, forceSz=True)
        res = []
        runs = self._handler._group[_HDF5DataHandler._RUNS_GROUP_NAME]
        for ri in nutils.getSliceIds(key[0], sz=len(self)):
            res.append([])
            runGrp = runs[_HDF5DataHandler._RUN_GROUP_TEMPLATE.format(ri)]
            runTime = runGrp[_HDF5DataHandler._TIME_DSET_NAME]
            runData = runGrp[_HDF5DataHandler._DATA_DSET_NAME]
            for i in nutils.getSliceIds(key[1], runTime.shape[0]):
                if self._saveTime:
                    res[-1].append(runTime[i])
                else:
                    remapKey = key[2]
                    if self._handler._revColRemapping is not None:
                        remapKey = self._handler._revColRemapping[key[2]]
                    if isinstance(remapKey, slice):
                        res[-1].append(runData[i, remapKey])
                    elif hasattr(remapKey, '__iter__'):
                        res[-1].append([runData[i, j] for j in remapKey])
                    else:
                        res[-1].append([runData[i, remapKey]])
        if forceArray:
            return nutils.nparray(res)
        mk = tuple(slice(None) if isinstance(k, slice) else 0 for k in key)
        return nutils.nparray(_sliceData(res, mk))

    def __len__(self):
        return len(self._handler._group[_HDF5DataHandler._RUNS_GROUP_NAME])

    def __array__(self):
        return self.__getitem__(slice(None, None, None), forceArray=True)


class _HDF5DistribDataAccessor(_HDF5DataAccessor):
    """
    Data accessor for HDF5 files
    """

    def __init__(self, handler, **kwargs):
        super().__init__(handler, saveTime=False, **kwargs)

    def __getitem__(self, key, forceArray=False):
        key = nutils.formatKey(key, self._nbDims, forceSz=True)
        runKey, timeKey, colKey = key
        rsInd = self._handler._rsInd
        rnk2RsAndInds = {}
        allColInds = nutils.getSliceIds(colKey, sz=self._handler._nbCols)
        for i, ci in enumerate(allColInds):
            rnk, li = self._handler._colMap[:, ci]
            if rnk not in rnk2RsAndInds:
                rs = self._handler._dbh._distribRankDBHs[rnk]._distribRS[self._handler._dbUid][rsInd]
                rnk2RsAndInds[rnk] = (rs, [], [])
            rnk2RsAndInds[rnk][1].append(i)
            rnk2RsAndInds[rnk][2].append(li)
        res = None
        for rnk, (rs, resInds, locInds) in rnk2RsAndInds.items():
            locData = rs.data.__getitem__((runKey, timeKey, locInds), forceArray=True)
            if res is None:
                nbRuns, nbTpts, _ = locData.shape
                res = numpy.zeros((nbRuns, nbTpts, len(allColInds)))
            res[:, :, resInds] = locData

        if forceArray:
            return nutils.nparray(res)
        mk = tuple(slice(None) if isinstance(k, slice) else 0 for k in key)
        return nutils.nparray(_sliceData(res, mk))

    def __len__(self):
        # Number of runs should be the same in all files, return the first
        rnk = self._handler._colMap[0, 0]
        rs = self._handler._dbh._distribRankDBHs[rnk]._distribRS[self._handler._dbUid][self._handler._rsInd]
        return len(rs.data)


###################################################################################################
# Database handlers


class DatabaseHandler:
    """Base class for all database handlers."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._parameters = None
        self._param2Groups = None

    def _getDataHandler(self, rs):
        """Return a _DBDataHandler for ResultSelector rs."""
        pass

    def _getFilePaths(self):
        """Return a list of file paths managed by this rank"""
        return []

    def __getitem__(self, key):
        """Access a run group from its unique identifier"""
        raise NotImplementedError()

    def __iter__(self):
        """Iterate over run groups in the database"""
        raise NotImplementedError()

    def _initializeParameters(self):
        """Setup parameters and param2Groups data structures"""
        if self._parameters is None or self._param2Groups is None:
            self._parameters = {}
            self._param2Groups = {}
            for group in self:
                for name, val in group.parameters.items():
                    self._parameters.setdefault(name, set()).add(val)
                    self._param2Groups.setdefault((name, val), set()).add(group)

    @property
    def parameters(self):
        """All parameter values from all run groups

        A dictionary whose keys are parameter names and values are sets of possible values

        :rtype: Mapping[str, Set[Any]], read-only
        """
        self._initializeParameters()
        return copy.deepcopy(self._parameters)

    def filter(self, **kwargs):
        r"""Return all run groups that match the given parameter values

        :param \*\*kwargs: Keyword arguments specifying the values of parameters that the filtered
            run groups must match.

        :return: A set of run groups whose parameter values match the parameters supplied by keyword
            arguments
        :rtype: Set[DatabaseGroup]
        """
        self._initializeParameters()
        res = None
        for keyVal in kwargs.items():
            try:
                if res is None:
                    res = copy.copy(self._param2Groups[keyVal])
                else:
                    res = res & self._param2Groups[keyVal]
            except KeyError:
                raise KeyError(f'Could not find any run group with {keyVal[0]} == {keyVal[1]}')
        return res if res is not None else set(self)

    def get(self, **kwargs):
        r"""Get the run group that matches the given parameter values

        :param \*\*kwargs: Keyword arguments specifying the values of parameters that the run group
            must match.

        :return: The single run group that matches the given parameter values
        :rtype: DatabaseGroup

        If several or none of the groups match these values, an exception will be raised.
        """
        groups = self.filter(**kwargs)
        if len(groups) != 1:
            paramVals = ' and '.join(f'{name} == {val}' for name, val in kwargs.items())
            raise Exception(
                f'Expected a single run group to match ({paramVals}), got {len(groups)} run groups instead.'
            )
        else:
            return next(iter(groups))


class DatabaseGroup:
    """Base class for all database run groups"""

    def __init__(self, dbh, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._dbh = dbh

    @property
    def name(self):
        """The unique identifier of the group

        :type: str, read-only
        """
        raise NotImplementedError()

    def __hash__(self):
        return hash((self._dbh, self.name))

    def __eq__(self, other):
        return isinstance(other, self.__class__) and (self._dbh, self.name) == (other._dbh, other.name)


class SQLiteDBHandler(DatabaseHandler):
    r"""SQLite database handler

    :param path: The path to the SQLite database file
    :type path: str
    :param \*args: Transmitted to :py:func:`sqlite3.connect`, see
        `documentation <https://docs.python.org/3/library/sqlite3.html#sqlite3.connect>`__ for
        details
    :param commitFreq: How frequently the data should be committed to the database. For example,
        this value is set to 10 by default which means that every 10 saving events, the data will
        be committed. If a result selector is saved every 10ms, it means the data will be committed
        to database every 100ms.
    :type commitFreq: int
    :param \*\*kwargs: Transmitted to :py:func:`sqlite3.connect`, see
        `documentation <https://docs.python.org/3/library/sqlite3.html#sqlite3.connect>`__ for
        details

    Handles the connection to a SQLite database and enables the saving of result selectors to that
    database. In contrast to the regular saving of result selectors (to memory or to file), it is
    possible to define groups of runs identified by a unique string so that the same database file
    can be used for several (sequential) runs of scripts.

    The database handler should be used as a context manager that wraps all simulation code. Inside
    this wrapped block, the user should call the :py:func:`steps.API_2.sim.Simulation.toDB` method
    to indicate that all results selectors associated to the simulation should be saved in the
    database. In this call, the user should provide the unique simulation group identifier as well
    as optional parameters that will also be saved to the database.

    Usage when saving::

        sim.toSave(rs1, rs2, rs3, dt=0.01)                # Add the result selectors to the
                                                          # simulation.

        with SQLiteDBHandler(dbPath) as dbh:              # Create database handler.

            sim.toDB(dbh, 'MySimulation', val1=1, val2=2) # Create a new group of runs in the
                                                          # database with identifier
                                                          # 'MySimulation' and save additional
                                                          # parameters val1 and val2.

            for i in range(NBRUNS):                       # Run a series of runs, all of them
                sim.newRun()                              # being associated to the
                ...                                       # 'MySimulation' group.
                sim.run(...)

    Note that after calling `sim.toDB(...)` it is still possible to force the saving of some result
    selectors to files by calling ``toFile(...)`` on them. Result selectors that contain a high
    number of values to save are probably better saved to a file. The name of the file can be added
    as a keyword parameter to the ``simtoDB(...)`` call to simplify loading.

    Usage when accessing data from the database::

        with SQLiteDBHandler(dbPath) as dbh:              # Create database handler.

            val1 = dbh['MySimulation'].val1               # Querying a parameter value from the
                                                          # 'MySimulation' group.

            rs1, rs2, rs3 = dbh['MySimulation'].results   # Querying the result selectors that
                                                          # were saved for the 'MySimulation'
                                                          # group. They are returned in the same
                                                          # order ad they were added to the
                                                          # simulation.

            plt.plot(rs1.time[0], rs1.data[0])            # The results selectors can be used as
                                                          # if they had been declared in the same
                                                          # process.

    """

    _RS_MAIN_TABLE_NAME = 'ResultSelectors'
    _RS_LABEL_TABLE_NAME = 'Labels'
    _RS_META_DATA_TABLE_NAME = 'MetaData'
    _GROUP_TABLE_NAME = 'SimGroups'
    _DEFAULT_COMMIT_FREQ = 10
    _GROUP_TABLE_KEYS = ['groupid', 'timestamp', 'uniqueid', 'nbselectors']

    def __init__(self, path, *args, commitFreq=-1, **kwargs):
        super().__init__(*args, **kwargs)
        self._path = path
        # Only rank 0 should actually connect to the database
        if nsim.MPI._shouldWrite:
            self._conn = sqlite3.connect(path, *args, **kwargs)
            self._conn.row_factory = sqlite3.Row
            self._connected = True
            self._createTables()
        else:
            self._conn = None
            self._connected = False

        self._commitFreq = commitFreq if commitFreq > 0 else SQLiteDBHandler._DEFAULT_COMMIT_FREQ

        self._groupId = None
        # self._simInd = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._close()

    def __del__(self):
        self._close()

    def _close(self):
        """Commit and close the connection."""
        if self._connected:
            self._conn.commit()
            self._conn.close()
            self._connected = False

    def _checkConnection(self):
        if not self._connected:
            if nsim.MPI._shouldWrite:
                raise Exception(f'The connection to the database has been closed.')
            else:
                raise Exception(f'Cannot access the database out of the rank 0 process while using MPI.')

    def _getFilePaths(self):
        """Return a list of file paths managed by this rank"""
        return [self._path] if self._conn is not None else []

    def _getDataHandler(self, rs):
        """Return a _DBDataHandler for ResultSelector rs."""
        return _SQLiteDataHandler(rs, self, self._commitFreq)

    def _newGroup(self, sim, uid, selectors, **kwargs):
        """Initialize the database and add a new run group."""
        if nsim.MPI._shouldWrite:
            self._checkConnection()

            # Check if the group already exists
            rows = self._conn.execute(
                f"SELECT * FROM {SQLiteDBHandler._GROUP_TABLE_NAME} WHERE uniqueid = '{uid}'"
            ).fetchall()
            if len(rows) == 0:
                # If it doesn't, create it
                typeMap = {int: 'int', float: 'real', str: 'text', bytes: 'BLOB'}
                colNames = ['timestamp', 'uniqueid', 'nbselectors']
                values = [datetime.datetime.now(), uid, len(selectors)]
                for colName, val in kwargs.items():
                    if type(val) not in typeMap:
                        raise TypeError(
                            f'Cannot process {colName}={val} because val is not from one of '
                            f'these types: {typeMap.keys()}'
                        )
                    try:
                        self._conn.execute(
                            f'ALTER TABLE {SQLiteDBHandler._GROUP_TABLE_NAME} ADD COLUMN '
                            f'{colName} {typeMap[type(val)]}'
                        )
                    except sqlite3.OperationalError:
                        pass
                    colNames.append(colName)
                    values.append(val)

                c = self._conn.cursor()
                c.execute(
                    f"INSERT INTO {SQLiteDBHandler._GROUP_TABLE_NAME}({','.join(colNames)}) "
                    f"VALUES ({','.join('?'*len(values))})",
                    values,
                )
                self._groupId = c.lastrowid
            else:
                row = rows[0]
                # get existing group id
                self._groupId = row[SQLiteDBHandler._GROUP_TABLE_KEYS.index('groupid')]
                # Checks parameters
                params = {
                    k: v
                    for k, v in zip(row.keys(), row)
                    if k not in SQLiteDBHandler._GROUP_TABLE_KEYS and v is not None
                }
                if kwargs != params:
                    raise Exception(
                        f'The keyword arguments provided to the toDB method ({kwargs}) '
                        f'do not match with the keyword arguments in the database for '
                        f'the same unique identifier ({params}).'
                    )
                # Check selectors
                if row[SQLiteDBHandler._GROUP_TABLE_KEYS.index('nbselectors')] != len(selectors):
                    idnbs = SQLiteDBHandler._GROUP_TABLE_KEYS.index('nbselectors')
                    raise Exception(
                        f'The {uid} run group saved in the database is associated with '
                        f'{row[idnbs]} resultSelectors while the current simulation is '
                        f'associated with {len(selectors)}'
                    )
                allRs = self._conn.execute(
                    f'SELECT * FROM {SQLiteDBHandler._RS_MAIN_TABLE_NAME} '
                    f'WHERE groupid={self._groupId} ORDER BY rsid'
                ).fetchall()
                for dbrs, simrs in zip(allRs, selectors):
                    groupId, rsid, descr, tabName, nbCols = dbrs
                    if simrs._strDescr() != descr:
                        raise Exception(
                            f'The result selector that was previously used for this '
                            f'unique identifier ({descr}) differs from the one being '
                            f'currently used ({simrs._strDescr()}).'
                        )
                    if simrs._getEvalLen() != nbCols:
                        raise Exception(
                            f'The result selector that was previously used for this '
                            f'unique identifier had {nbCols} columns while the current '
                            f' one has {simrs._getEvalLen()} columns.'
                        )
                    # check labels
                    dblbls = self._labelsQuerry(groupId, rsid)
                    if simrs.labels != dblbls:
                        raise Exception(
                            f'The result selector that was previously used for this '
                            f'unique identifier had different column labels. Expected '
                            f'{dblbls} but got {simrs.labels} instead.'
                        )
                    # check metadata
                    dbmd = self._metaDataQuerry(groupId, rsid)
                    simmd = simrs.metaData._dict
                    if simmd != dbmd:
                        raise Exception(
                            f'The result selector that was previously used for this '
                            f'unique identifier had different metadata. Expected '
                            f'{dbmd} but got {simmd} instead.'
                        )

            self._conn.commit()

        return selectors

    def _labelsQuerry(self, groupId, rsid):
        """Return labels for ResultSelector rsid in group groupid."""
        self._checkConnection()
        rows = self._conn.execute(
            f'SELECT label FROM '
            f'{SQLiteDBHandler._RS_LABEL_TABLE_NAME} '
            f'WHERE groupid={groupId} AND rsid={rsid} '
            f'ORDER BY colind'
        ).fetchall()
        return [row[0] for row in rows]

    def _metaDataQuerry(self, groupId, rsid):
        """Return metadata for ResultSelector rsid in group groupid."""
        self._checkConnection()
        dat = self._conn.execute(
            f'SELECT data FROM '
            f'{SQLiteDBHandler._RS_META_DATA_TABLE_NAME} '
            f'WHERE groupid={groupId} AND rsid={rsid} '
        ).fetchone()[0]
        return pickle.loads(dat)

    def _createTables(self):
        """Create the tables if they do not exist."""
        self._checkConnection()
        self._conn.execute(
            f'CREATE TABLE IF NOT EXISTS {SQLiteDBHandler._GROUP_TABLE_NAME} '
            f'(groupid INTEGER PRIMARY KEY AUTOINCREMENT, timestamp date, '
            f'uniqueid text UNIQUE, nbselectors int);'
        )
        self._conn.execute(
            f'CREATE TABLE IF NOT EXISTS {SQLiteDBHandler._RS_MAIN_TABLE_NAME} '
            f'(groupid int, rsid int, descr text, tabName text, nbcols int);'
        )
        self._conn.execute(
            f'CREATE TABLE IF NOT EXISTS {SQLiteDBHandler._RS_LABEL_TABLE_NAME} '
            f'(groupid int, rsid int, colind int, label text);'
        )
        self._conn.execute(
            f'CREATE TABLE IF NOT EXISTS {SQLiteDBHandler._RS_META_DATA_TABLE_NAME} '
            f'(groupid int, rsid int, data blob);'
        )
        self._conn.commit()

    def __getitem__(self, key):
        """Access a SQLite group from its unique identifier

        :param key: Unique identifier to the group
        :type key: str
        :returns: The associated SQLite group
        :rtype: :py:class:`SQLiteGroup`

        See :py:class:`SQLiteDBHandler` for usage examples.

        Raises a ``KeyError`` if the key is not in the database.
        :meta public:
        """
        self._checkConnection()
        if not isinstance(key, str):
            raise TypeError(f'Expected a unique identifier string, got {key} instead.')
        rows = self._conn.execute(
            f"SELECT * FROM {SQLiteDBHandler._GROUP_TABLE_NAME} WHERE uniqueid == '{key}'"
        ).fetchall()
        if len(rows) == 0:
            raise KeyError(f'{key} does not exist in {self._path}.')
        return SQLiteGroup(self, rows[0])

    def __iter__(self):
        """Iterate over SQLite groups in the database

        Usage::

            with SQLiteDBHandler(dbPath) as dbh:              # Create database handler.

                for group in dbh:                             # Iterate over all groups

                    val1 = group.val1                         # Access group data

        :meta public:
        """
        self._checkConnection()
        rows = self._conn.execute(
            f'SELECT * FROM {SQLiteDBHandler._GROUP_TABLE_NAME} ORDER BY groupid'
        ).fetchall()
        for row in rows:
            yield SQLiteGroup(self, row)


class SQLiteGroup(DatabaseGroup):
    """A class representing a group of runs in a SQLite database

    .. note::
        This class should never be instantiated by the user, it is obtained through
        :py:class:`SQLiteDBHandler` instead.
    """

    def __init__(self, dbh, row, *args, **kwargs):
        super().__init__(*args, dbh=dbh, **kwargs)
        self._row = row
        self._dict = {}
        for k, v in zip(row.keys(), row):
            if v is not None:
                self._dict[k] = v

    def __getattr__(self, name):
        """Attribute access for parameters of the group

        :param name: Name of the parameter, as defined in the original call to ``sim.toDB(...)``
        :type name: str

        :returns: The corresponding parameter value

        See :py:class:`SQLiteDBHandler` for usage examples.

        :meta public:
        """
        if name not in self._dict:
            raise AttributeError(f'{name} is not an attribute of {self}.')
        return self._dict[name]

    @DatabaseGroup.name.getter
    def name(self):
        """The unique identifier of the group

        :type: str, read-only
        """
        return self._dict['uniqueid']

    @property
    def results(self):
        """A list of all result selectors that were saved

        :type: List[:py:class:`ResultSelector`], read-only

        The result selectors are returned in the same order as they were added to the simulation
        with the :py:func:`steps.API_2.sim.Simulation.toSave` method.

        See :py:class:`SQLiteDBHandler` for usage examples.
        """
        res = [None] * self.nbselectors

        rows = self._dbh._conn.execute(
            f'SELECT * FROM {SQLiteDBHandler._RS_MAIN_TABLE_NAME} WHERE groupid={self.groupid}'
        ).fetchall()
        for groupId, rsid, descr, tableName, nbCols in rows:
            res[rsid] = _ReadOnlyResultSelector(
                _SQLiteDataHandler(
                    None, self._dbh, None, groupId=groupId, rsid=rsid, tableName=tableName, nbCols=nbCols
                )
            )

        return res

    @property
    def parameters(self):
        """A dictionary of all parameters defined for this group

        :type: Mapping[str, Any], read-only

        Usage::

            >>> with SQLiteDBHandler(dbPath) as dbh:
            ...     dbh['MySimulation'].parameters
            {'val1': 1, 'val2': 2}
        """
        return {k: v for k, v in self._dict.items() if k not in SQLiteDBHandler._GROUP_TABLE_KEYS}


class HDF5Handler(DatabaseHandler, nutils.Versioned):
    """HDF5 File handler

    :param pathPrefix: Path and prefix for the HDF5 file(s) (e.g. './data/HDF5Data' would
        yield one file named './data/HDF5Data.h5' when the simulation is not distributed and several
        files named './data/HDF5Data_rank0.h5', './data/HDF5Data_rank1.h5', etc. when the simulation
        is distributed.
    :type pathPrefix: str
    :param hdf5FileKwArgs: Keyword arguments transmitted to :py:func:`h5py.File`, see
        `documentation <https://docs.h5py.org/en/stable/high/file.html#h5py.File>`__ for
        details
    :type hdf5FileKwArgs: dict
    :param hdf5DatasetKwArgs: Keyword arguments transmitted to :py:func:`h5py.Group.create_dataset`, see
        `documentation <https://docs.h5py.org/en/stable/high/group.html#h5py.Group.create_dataset>`__ for
        details. Most notably, compression-related argument can be set there.
    :type hdf5FileKwArgs: dict

    Handles reading and writing to an HDF5 file and enables the saving of result selectors to that
    file. In contrast to the regular saving of result selectors (to memory or to file), it is
    possible to define groups of runs identified by a unique string so that the same HDF5 file
    can be used for several (sequential) runs of scripts.

    The HDF5Handler should be used as a context manager that wraps all simulation code. Inside
    this wrapped block, the user should call the :py:func:`steps.API_2.sim.Simulation.toDB` method
    to indicate that all results selectors associated to the simulation should be saved in the
    HDF5 file. In this call, the user should provide the unique simulation group identifier as well
    as optional parameters that will also be saved to the file.

    Usage when saving::

        sim.toSave(rs1, rs2, rs3, dt=0.01)                # Add the result selectors to the
                                                          # simulation.

        with HDF5Handler('./path/to/Prefix') as hdf:      # Create database handler.

            sim.toDB(hdf, 'MySimulation', val1=1, val2=2) # Create a new group of runs in the
                                                          # HDF5 file with identifier
                                                          # 'MySimulation' and save additional
                                                          # parameters val1 and val2.

            for i in range(NBRUNS):                       # Run a series of runs, all of them
                sim.newRun()                              # being associated to the
                ...                                       # 'MySimulation' group.
                sim.run(...)

    Note that, in contrast with :py:class:`SQLiteDBHandler`, there is no use for forcing the saving of
    some result selectors to files by calling ``toFile(...)`` on them. HDF5 files can contain high
    amounts of data.

    Usage when accessing data from the database::

        with HDF5Handler('./path/to/Prefix') as hdf:      # Create database handler.

            val1 = hdf['MySimulation'].val1               # Querying a parameter value from the
                                                          # 'MySimulation' group.

            rs1, rs2, rs3 = hdf['MySimulation'].results   # Querying the result selectors that
                                                          # were saved for the 'MySimulation'
                                                          # group. They are returned in the same
                                                          # order as they were added to the
                                                          # simulation.

            plt.plot(rs1.time[0], rs1.data[0])            # The results selectors can be used as
                                                          # if they had been declared in the same
                                                          # process.

    Note that :py:class:`XDMFHandler` inherits from :py:class:`HDF5Handler` and generates `.xmf` files
    that point to the HDF5 files and can be read by data visualization software such as
    `Paraview <https://www.paraview.org/>`_. 
    """

    _TIMESTAMP_ATTR_NAME = 'timestamp'
    _STEPS_VERSION_ATTR_NAME = 'steps_version'
    _NB_DISTR_RANKS_ATTR_NAME = 'nb_distributed_ranks'
    _GROUP_DEFAULT_ATTRS = [_TIMESTAMP_ATTR_NAME, _STEPS_VERSION_ATTR_NAME, _NB_DISTR_RANKS_ATTR_NAME]

    _RS_DESCRIPTION_ATTR = 'Description'
    _RS_INDEX_ATTR = 'RSIndex'
    _RS_GROUP_NAME = 'ResultSelector{}'
    _RS_DISTGROUP_NAME = 'DistributedResultSelector{}'
    _RS_DIST_IND_MAP_NAME = 'DistributedColumnMap'
    _DISTRIBUTED_HDF_SUFFIX = '{}_rank{}'
    _HDF_EXTENSION = '.h5'

    def __init__(self, pathPrefix, hdf5FileKwArgs={}, hdf5DatasetKwArgs={}, **kwargs):
        super().__init__(**kwargs)

        self._pathPrefix = pathPrefix
        self._path = None
        self._file = None

        self._currGroup = None
        self._shouldWrite = nsim.MPI._shouldWrite
        self._fileKwArgs = hdf5FileKwArgs
        self._dataSetKWargs = hdf5DatasetKwArgs

        self._nbSavingRanks = None

        # HDF5 database handlers needed when loading distributed data
        self._distribRankDBHs = {} # rank -> dbh
        # Local read-only result selectors needed when loading distributed data
        self._distribRS = {} # dbUID -> { rsInd -> rs }

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._close()

    def __del__(self):
        self._close()

    def _close(self):
        """Close the file"""
        if hasattr(self, '_file') and self._file is not None:
            self._file.close()
        for rnk, dbh in self._distribRankDBHs.items():
            if rnk != nsim.MPI.rank:
                dbh._close()

    def _checkOpenFile(self, sim=None):
        if self._file is None:
            import h5py
            if sim is None:
                self._path = self._pathPrefix + HDF5Handler._HDF_EXTENSION
                if not os.path.isfile(self._path):
                    self._path = (
                        HDF5Handler._DISTRIBUTED_HDF_SUFFIX.format(self._pathPrefix, 0) +
                        HDF5Handler._HDF_EXTENSION
                    )
                if not os.path.isfile(self._path):
                    raise Exception(
                        f'Cannot load any HDF files with prefix {self._pathPrefix}.'
                    )
                self._file = h5py.File(self._path, 'r', **self._fileKwArgs)
            else:
                if sim._isDistributed():
                    self._path = HDF5Handler._DISTRIBUTED_HDF_SUFFIX.format(self._pathPrefix, nsim.MPI.rank)
                elif nsim.MPI._shouldWrite:
                    self._path = self._pathPrefix
                else:
                    raise Exception(
                        f'Cannot access the HDF5 file out of the rank 0 process while using MPI.'
                    )
                self._path += HDF5Handler._HDF_EXTENSION
                self._file = h5py.File(self._path, 'a', **self._fileKwArgs)
        elif not self._file:
            raise Exception(f'The HDF5 was closed.')

    def _getFilePaths(self):
        """Return a list of file paths managed by this rank"""
        return [self._path] if self._path is not None else []

    def _getRsHDFGroup(self, rs, groupNamePattern=None):
        if groupNamePattern is None:
            groupNamePattern = HDF5Handler._RS_GROUP_NAME
        if self._shouldWrite and rs._getEvalLen() > 0:
            rsName = groupNamePattern.format(rs._selectorInd)
            if rsName not in self._currGroup:
                rsgroup = self._currGroup.create_group(rsName)
                rsgroup.attrs[HDF5Handler._RS_DESCRIPTION_ATTR] = rs._strDescr()
                rsgroup.attrs[HDF5Handler._RS_INDEX_ATTR] = rs._selectorInd
            return self._currGroup[rsName]
        else:
            return None

    def _getDataHandler(self, rs, groupNamePattern=None):
        """Return a _DBDataHandler for ResultSelector rs."""
        if groupNamePattern is None:
            groupNamePattern = HDF5Handler._RS_GROUP_NAME
        return _HDF5DataHandler(self, rs, self._getRsHDFGroup(rs, groupNamePattern), version=self._version)

    def _checkSelectors(self, uid, selectors):
        """Check that selectors in the HDF5 file match the simulation selectors"""
        selectorNames = [
            n for n in self._currGroup if n.startswith(HDF5Handler._RS_GROUP_NAME.format(''))
        ]
        nonEmptySelectors = [rs for rs in selectors if rs._getEvalLen() > 0]
        if len(selectorNames) != len(nonEmptySelectors):
            raise Exception(
                f'The {uid} run group saved in the database is associated with '
                f'{len(selectorNames)} resultSelectors while the current simulation is '
                f'associated with {len(nonEmptySelectors)}.'
            )
        for rsName, simrs in zip(selectorNames, nonEmptySelectors):
            handler = self._getDataHandler(simrs)
            grouprs = self._currGroup[rsName]
            descr = grouprs.attrs[HDF5Handler._RS_DESCRIPTION_ATTR]
            rsInd = grouprs.attrs[HDF5Handler._RS_INDEX_ATTR]
            if simrs._strDescr() != descr:
                raise Exception(
                    f'The result selector that was previously used for this '
                    f'unique identifier ({descr}) differs from the one being '
                    f'currently used ({simrs._strDescr()}).'
                )
            if simrs._selectorInd != rsInd:
                raise Exception(
                    f'The index result selector that was previously used for this '
                    f'unique identifier ({rsInd}) differs from the one being '
                    f'currently used ({simrs._selectorInd}).'
                )
            if simrs._getEvalLen() != len(handler.data()[0, 0, :]):
                raise Exception(
                    f'The result selector that was previously used for this unique identifier '
                    f'had {grouprs["data"].shape[2]} columns while the current one has '
                    f'{simrs._getEvalLen()} columns.'
                )
            # check labels
            if handler.labels() != simrs.labels:
                raise Exception(
                    f'The result selector that was previously used for this '
                    f'unique identifier had different column labels. Expected '
                    f'{grouprs["labels"]} but got {simrs.labels} instead.'
                )
            # check metadata
            filemd = handler.metaData()
            simmd = simrs.metaData._dict
            if simmd != filemd:
                raise Exception(
                    f'The result selector that was previously used for this '
                    f'unique identifier had different metadata. Expected '
                    f'{filemd} but got {simmd} instead.'
                )

    def _checkDistributedSelectors(self, uid, rsInd2DistColMap):
        """Check that the distributed column map of each distributed result selector in the HDF5 file matches
        the simulation values"""
        selectorNames = [
            n for n in self._currGroup if n.startswith(HDF5Handler._RS_DISTGROUP_NAME.format(''))
        ]
        if len(selectorNames) != len(rsInd2DistColMap):
            if nsim.MPI._shouldWrite:
                raise Exception(
                    f'The {uid} run group saved in the database is associated with '
                    f'{len(selectorNames)} resultSelectors while the current simulation is '
                    f'associated with {len(rsInd2DistColMap)}.'
                )
            else:
                raise Exception(
                    f'The HDF5 file associated with MPI rank {nsim.MPI.rank} contains '
                    f'{len(selectorNames)} distributed column remapping dataset while it should '
                    f'contain {len(rsInd2DistColMap)}.'
                )
        for rsInd, rsName in enumerate(selectorNames):
            grouprs = self._currGroup[rsName]
            simColMap = rsInd2DistColMap[rsInd]
            rsColMap = numpy.array(grouprs[HDF5Handler._RS_DIST_IND_MAP_NAME])
            if simColMap.shape[1] != rsColMap.shape[1]:
                raise Exception(
                    f'The result selector that was previously used for this unique identifier '
                    f'had {rsColMap.shape[1]} columns while the current one has '
                    f'{simColMap.shape[1]} columns.'
                )
            for i, ((rsRnk, rsIdx), (simRnk, simIdx)) in enumerate(zip(rsColMap.T, simColMap.T)):
                if (rsRnk, rsIdx) != (simRnk, simIdx):
                    raise Exception(
                        f'In Result selector {rsInd}, according to the data in the HDF5 file, column '
                        f'{i} should be mapped to rank {rsRnk} and local column {rsIdx}, but in the'
                        f'current simulation, it is mapped to rank {simRnk} and local column {simIdx}.'
                    )

    def _distributeSelectors(self, sim, selectors):
        """Distribute result selectors
        Return a list of distributed selectors and, on rank 0, a map between result selector index
        and distributed column mapping: for each column in the non-distributed selector, it contains
        the rank in which the value will be saved and the local index in the distributed selector.
        """
        optimGroups = [[]]
        rs2FullLen = {}
        for rs in selectors:
            # Distribute the result selector
            rs, changed = rs._distribute()
            if rs is not None:
                # Re-create optimization groups
                if len(optimGroups[-1]) == 0 or optimGroups[-1][-1][0]._optimGroupInd == rs._optimGroupInd:
                    optimGroups[-1].append((rs, changed))
                else:
                    optimGroups.append([])

        distribSelectors = []
        for optGrp in optimGroups:
            if len(optGrp) > 0:
                distRs, distChanged = zip(*optGrp)
                # Reoptimize the calls if something changed after distribution
                if any(distChanged):
                    distRs = nsaving_optim.OptimizeSelectors(sim, distRs)
                distribSelectors += distRs

        import mpi4py
        rsInd2DistColMap = {}
        # Build the mapping between full selector and local distributed selectors on rank 0
        localDistrInds = {rs._selectorInd: (rs._fullLen, rs._distrInds) for rs in distribSelectors}
        allDistribInds = mpi4py.MPI.COMM_WORLD.gather(localDistrInds, root=0)
        if nsim.MPI._shouldWrite:
            allRsInds = set()
            rsFullLens = {}
            for dct in allDistribInds:
                for rsInd, (fullLen, _) in dct.items():
                    allRsInds.add(rsInd)
                    rsFullLens[rsInd] = max(rsFullLens.get(rsInd, 0), fullLen)
            for rsInd in allRsInds:
                dcm = -1 * numpy.ones((2, rsFullLens[rsInd]), dtype=numpy.int64)
                for rnk, distInds in enumerate(allDistribInds):
                    if rsInd in distInds:
                        for localInd, ind in enumerate(distInds[rsInd][1]):
                            dcm[0, ind] = rnk
                            dcm[1, ind] = localInd
                rsInd2DistColMap[rsInd] = dcm

        return distribSelectors, rsInd2DistColMap

    def _newGroup(self, sim, uid, selectors, **kwargs):
        """Initialize the file and add a new run group."""
        self._shouldWrite = nsim.MPI._shouldWrite or sim._isDistributed()
        self._nbSavingRanks = nsim.MPI.nhosts if sim._isDistributed() else 1
        if self._shouldWrite:
            self._checkOpenFile(sim)

            # Check if the group already exists
            if uid not in self._file:
                # If it doesn't, create it
                self._currGroup = self._file.create_group(uid, track_order=True)
                self._currGroup.attrs[HDF5Handler._TIMESTAMP_ATTR_NAME] = str(datetime.datetime.now())
                self._currGroup.attrs[HDF5Handler._STEPS_VERSION_ATTR_NAME] = steps.__version__
                self._currGroup.attrs[HDF5Handler._NB_DISTR_RANKS_ATTR_NAME] = self._nbSavingRanks
                try:
                    for argName, val in kwargs.items():
                        self._currGroup.attrs[argName] = val
                except Exception as ex:
                    # If attributes could not be set, delete the group before raising the exception
                    self._currGroup = None
                    del self._file[uid]
                    raise ex
                # Pre-create distributed result selector groups if needed
                if sim._isDistributed() and self._nbSavingRanks > 1:
                    # All ranks distribute their selectors
                    distribSelectors, rsInd2DistColMap = self._distributeSelectors(sim, selectors)

                    if nsim.MPI._shouldWrite:
                        # Only rank 0 writes the full result selectors
                        for rs in selectors:
                            handler = self._getDataHandler(rs, groupNamePattern=HDF5Handler._RS_DISTGROUP_NAME)
                            handler._group.create_dataset(
                                HDF5Handler._RS_DIST_IND_MAP_NAME,
                                data=rsInd2DistColMap[rs._selectorInd],
                                **self._dataSetKWargs
                            )

                    selectors = distribSelectors
                # Pre-create local result selector groups to be sure to have the correct order
                for rs in selectors:
                    self._getRsHDFGroup(rs)
            else:
                self._currGroup = self._file[uid]
                # Load version and check it matches the currently used version
                version = self._parseVersion(self._currGroup.attrs[HDF5Handler._STEPS_VERSION_ATTR_NAME])
                if version != self._version:
                    raise Exception(
                        f'Cannot add results to a group that was created with a different STEPS version. '
                        f'Group version: {version}, Current version: {self._version}'
                    )
                # Check that the number of distributed ranks is the same
                fileNbRanks = self._currGroup.attrs[HDF5Handler._NB_DISTR_RANKS_ATTR_NAME]
                if fileNbRanks != self._nbSavingRanks:
                    raise Exception(
                        f'Cannot add results to a group that was created with a different number of MPI '
                        f'processes. The {uid} group was created with {fileNbRanks} while the current '
                        f'simulation is being run with {self._nbSavingRanks} processes.'
                    )
                # Checks parameters
                params = {
                    k: v for k, v in self._currGroup.attrs.items() if k not in HDF5Handler._GROUP_DEFAULT_ATTRS
                }
                if kwargs != params:
                    raise Exception(
                        f'The keyword arguments provided to the toDB method ({kwargs}) '
                        f'do not match with the keyword arguments in the HDF5 file for '
                        f'the same unique identifier ({params}).'
                    )
                # Check distributed result selectors
                if sim._isDistributed() and self._nbSavingRanks > 1:
                    # All ranks distribute their selectors
                    selectors, rsInd2DistColMap = self._distributeSelectors(sim, selectors)
                    self._checkDistributedSelectors(uid, rsInd2DistColMap)
                # Check local result selectors
                self._checkSelectors(uid, selectors)
        else:
            self._close()

        return selectors

    def __getitem__(self, key):
        """Access an HDF5 group from its unique identifier

        :param key: Unique identifier to the group
        :type key: str
        :returns: The associated HDF5 group
        :rtype: :py:class:`HDF5Group`

        See :py:class:`HDF5Handler` for usage examples.

        Raises a ``KeyError`` if the key is not in the file.

        :meta public:
        """
        self._checkOpenFile()
        if not isinstance(key, str):
            raise TypeError(f'Expected a unique identifier string, got {key} instead.')
        if key not in self._file:
            raise KeyError(f'{key} does not exist in {self._path}.')
        group = self._file[key]
        return HDF5Group(self, group, version=group.attrs[HDF5Handler._STEPS_VERSION_ATTR_NAME])

    def __iter__(self):
        """Iterate over STEPS groups in the file

        Usage::

            with HDF5Handler(filePath) as hdf:    # Create database handler.

                for group in hdf:                 # Iterate over all groups

                    val1 = group.val1             # Access group data

        Note that not all HDF5 groups in the file will be iterated on, only the ones that were added
        by STEPS. Other groups will be ignored by STEPS.

        :meta public:
        """
        self._checkOpenFile()
        res = []

        def visit(name, group):
            if all(attr in group.attrs for attr in HDF5Handler._GROUP_DEFAULT_ATTRS):
                res.append(HDF5Group(self, group, version=group.attrs[HDF5Handler._STEPS_VERSION_ATTR_NAME]))
        self._file.visititems(visit)
        for gr in res:
            yield gr


class HDF5Group(DatabaseGroup, nutils.Versioned):
    """A class representing a group of runs in an HDF5 file

    .. note::
        This class should never be instantiated by the user, it is obtained through
        :py:class:`HDF5Handler` instead.
    """

    def __init__(self, dbh, group, *args, **kwargs):
        super().__init__(*args, dbh=dbh, **kwargs)
        self._group = group

    def __getattr__(self, name):
        """Attribute access for parameters of the group

        :param name: Name of the parameter, as defined in the original call to ``sim.toDB(...)``
        :type name: str

        :returns: The corresponding parameter value

        See :py:class:`SQLiteDBHandler` for usage examples.

        :meta public:
        """
        if name not in self._group.attrs:
            raise AttributeError(f'{name} is not an attribute of {self}.')
        return self._group.attrs[name]

    @DatabaseGroup.name.getter
    def name(self):
        """The unique identifier of the group

        :type: str, read-only
        """
        path = self._group.name.split('/')
        return path[-1]

    @property
    def results(self):
        """A list of all result selectors that were saved

        :type: List[:py:class:`ResultSelector`], read-only

        The result selectors are returned in the same order as they were added to the simulation
        with the :py:func:`steps.API_2.sim.Simulation.toSave` method.

        See :py:class:`HDF5Handler` for usage examples.
        """
        # First check if distributed data available
        distrGroupNames = [
            gn for gn in self._group if re.match(HDF5Handler._RS_DISTGROUP_NAME.format(r'\d+'), gn) is not None
        ]
        if len(distrGroupNames) > 0:
            return [
                _ReadOnlyResultSelector(
                    _HDF5DistribDataHandler(self, self._dbh, None, self._group[gn], version=self._version)
                ) for gn in distrGroupNames
            ]
        else:
            return [
                _ReadOnlyResultSelector(
                    _HDF5DataHandler(self._dbh, None, self._group[gn], version=self._version)
                ) for gn in self._group if re.match(HDF5Handler._RS_GROUP_NAME.format(r'\d+'), gn) is not None
            ]

    @property
    def parameters(self):
        """A dictionary of all parameters defined for this group

        :type: Mapping[str, Any], read-only

        Usage::

            >>> with HDF5Handler(filePath) as hdf:
            ...     hdf['MySimulation'].parameters
            {'val1': 1, 'val2': 2}
        """
        return {k: v for k, v in self._group.attrs.items() if k not in HDF5Handler._GROUP_DEFAULT_ATTRS}


class XDMFHandler(HDF5Handler):
    """XDMF / HDF5 File handler

    :param pathPrefix: Path and prefix for the HDF5 file(s), see :py:class:`HDF5Handler`.
    :type pathPrefix: str
    :param hdf5FileKwArgs: see :py:class:`HDF5Handler`.
    :type hdf5FileKwArgs: dict
    :param hdf5DatasetKwArgs: see :py:class:`HDF5Handler`.
    :type hdf5FileKwArgs: dict
    :param xdmfFolder: Path to the folder to which XDMF files should be written. If `None`, it uses
        the folder in which HDF5 files will be saved.
    :type xdmfFolder: Union[str, None]

    The `XDMF file format <https://www.xdmf.org/>`_ uses XML files with the `.xmf` extension to describe
    data saved in an HDF5 file. Scientific visualization tools like `Paraview <https://www.paraview.org/>`_
    can read `.xmf` files and access the corresponding data in HDF5 files to display the mesh and the data
    associated with it.

    The :py:class:`XDMFHandler` database handler inherits from :py:class:`HDF5Handler` and thus behaves
    like an HDF5 database handler. The main difference is that :py:class:`XDMFHandler` also saves mesh
    information to the HDF5 file and generates `.xmf` XDMF files that describe all the data that can
    be visualized on meshes. Since it works in the same way as :py:class:`HDF5Handler`, usage examples
    can be seen there.

    Regardless on whether the data saved by a :py:class:`ResulSelector` is specific to mesh elements, it
    will be saved in the same way as if :py:class:`HDF5Handler` was used. Data that is specific to mesh
    elements (e.g. count of species in a tetrahedron, concentration of species in a region of interest, etc.)
    will be described in the `.xmf` file.

    The following mesh locations (see :py:class:`steps.API_2.sim.SimPath`) are supported:

    +---------------------+--------------------------+----------------+
    | Location            | Result selector          | XDMF data type |
    +=====================+==========================+================+
    | Tetrahedrons        | ``rs.TETS(tetLst)...``   | Cell data      |
    +---------------------+--------------------------+----------------+
    | Triangles           | ``rs.TRIS(triLst)...``   | Cell data      |
    +---------------------+--------------------------+----------------+
    | Vertices            | ``rs.VERTS(vertLst)...`` | Node data      |
    +---------------------+--------------------------+----------------+
    | Regions of Interest | ``rs.ROIname...``        | Grid data      |
    +---------------------+--------------------------+----------------+
    | Compartments        | ``rs.compName...``       | Grid data      |
    +---------------------+--------------------------+----------------+
    | Patches             | ``rs.patchName...``      | Grid data      |
    +---------------------+--------------------------+----------------+

    Technical note: the result selectors that involve mesh data (data that will be described in the `.xmf`
    file) might be stored in the HDF5 file in a way that is different from how it is normally stored with
    :py:class:`HDF5Handler`. Notably, the order of :py:class:`ResultSelector` columns might be changed to
    allow contiguous data access when loading the data into scientific visualization softwares. These
    differences do not affect users as long as the data is read using :py:func:`HDF5Handler.__getitem__`,
    columns will be correctly reordered to match the original :py:class:`ResultSelector` order.
    """

    _MESH_GROUP_NAME = 'mesh'
    _FILE_NAME_PATTERN = '{uid}_Run{run}_rank{rank}.xmf'
    _FULL_FILE_NAME_PATTERN = '{uid}_Run{run}_Full.xmf'

    def __init__(self, pathPrefix, hdf5FileKwArgs={}, hdf5DatasetKwArgs={}, xdmfFolder=None, **kwargs):
        super().__init__(pathPrefix, hdf5FileKwArgs, hdf5DatasetKwArgs, **kwargs)
        if xdmfFolder is None:
            xdmfFolder = os.path.dirname(pathPrefix)
        self._xdmfFolder = xdmfFolder

        self._xdmfTree = None
        self._fullXdmf = None
        self._currUID = None
        self._sim = None
        self._temporalGrids = None
        self._spatialGrids = None
        self._currTime = None
        self._currRun = None
        self._savedFilePaths = []

    def _close(self):
        """Close the file"""
        super()._close()
        if self._xdmfTree is not None:
            self._writeXMLTree()
            self._xdmfTree = None

    def _getFilePaths(self):
        """Return a list of file paths managed by this rank"""
        # To be saved soon:
        tobeSaved = []
        for fp in [self._getCurrFullXDMFFilePath(), self._getCurrXDMFFilePath()]:
            if fp not in self._savedFilePaths:
                tobeSaved.append(fp)
        return super()._getFilePaths() + self._savedFilePaths + tobeSaved

    def _getDataHandler(self, rs, groupNamePattern=None):
        """Return a _DBDataHandler for ResultSelector rs."""
        return _XDMFDataHandler(self, rs, self._getRsHDFGroup(rs, groupNamePattern))

    def _newGroup(self, sim, uid, selectors, **kwargs):
        """Initialize the file and add a new run group."""
        selectors = super()._newGroup(sim, uid, selectors, **kwargs)
        if self._xdmfTree is not None:
            self._writeXMLTree()
        self._currUID = uid
        self._sim = sim
        self._xdmfTree = None
        self._fullXdmf = None
        self._currRun = None

        if not isinstance(self._sim.geom, ngeom._BaseTetMesh):
            raise TypeError(f'XDMF data saving only works with tetrahedral meshes.')

        self._setUpSelectors(selectors)

        return selectors

    def _getAttributeName(self, obj, val, loctpe):
        name = f'{obj}.{val}' if obj is not None else val
        if loctpe is not None:
            name += f' ({loctpe})'
        return name

    def _extractAndAddSpatialGridInfo(self, elem2infos, refCls, rsColRemaps, ROIgrids):
        elem2infos = elem2infos[refCls._locStr]
        lstCls = refCls._lstCls
        # Group by set of saved values
        elemMap = {}
        for elemIdx, infos in elem2infos.items():
            key = frozenset(info[0:2] for info in infos)
            for info in infos:
                elemMap.setdefault(key, []).append((elemIdx, ) + info)

        lstKwArgs = dict(mesh=self._sim.geom)
        if self._sim._isDistributed():
            lstKwArgs['local'] = True
        coveredElems = refCls._lstCls([], **lstKwArgs)
        # Treat set of saved values independently
        for rsSet, elemInfo in elemMap.items():
            # Remove duplicate values from different result selectors
            elemInfo = {(idx, obj, val): (rsId, rsPos) for idx, obj, val, rsId, rsPos in elemInfo}

            # Get element list for each result selector
            rs2Elems = {}
            for (idx, obj, val), (rsId, rsPos) in elemInfo.items():
                rs2Elems.setdefault(rsId, set()).add(idx)

            # Compute the intersections of element lists
            grids = []
            for rsId, elemInds in rs2Elems.items():
                elemLst = lstCls(sorted(elemInds), **lstKwArgs)
                # Add intersection with ROI grids
                for roitets, roiRsVals, loc in ROIgrids[refCls]:
                    inter = (roitets & elemLst) - coveredElems
                    if len(inter) > 0:
                        grids.append((inter, copy.copy(roiRsVals), loc))
                        coveredElems |= inter
                newGrids = []
                for g, rsVals, loc in grids:
                    if len(elemLst) > 0:
                        inter = g & elemLst
                        if len(inter) > 0:
                            newGrids.append((inter, copy.copy(rsVals), loc))
                            g -= inter
                            elemLst -= inter
                    if len(g) > 0:
                        newGrids.append((g, copy.copy(rsVals), loc))
                if len(elemLst) > 0:
                    coveredElems |= elemLst
                    newGrids.append((elemLst, [], (refCls,)))
                grids = newGrids

            # Add The resulting element lists as our distinct spatial grids
            elem2Grid = {}
            for grid, rsVals, loc in grids:
                gridPos = len(self._spatialGrids[refCls])
                self._spatialGrids[refCls].append((grid, rsVals, loc))
                for elem in grid:
                    elem2Grid[elem.idx] = gridPos

            # Group data by result selectors and grid
            rs2Infos = {}
            for (idx, obj, val), (rsId, rsPos) in elemInfo.items():
                rs2Infos.setdefault(
                    (rsId, elem2Grid[idx]), {}
                ).setdefault(
                    (obj, val), []
                ).append((idx, rsPos))

            # Add mapping between result selectors and spatial grids
            for (rsId, gridPos), val2elems in rs2Infos.items():
                for (obj, val), idxLst in val2elems.items():
                    idxLst.sort(key=lambda x: x[0])
                    elemIdxs, rsPoss = zip(*idxLst)
                    start = len(rsColRemaps.get(rsId, []))
                    center = 'Node' if refCls == ngeom.VertReference else 'Cell'
                    self._spatialGrids[refCls][gridPos][1].append(
                        (obj, val, None, rsId, start, 1, len(elemIdxs), center)
                    )
                    rsColRemaps.setdefault(rsId, [])
                    rsColRemaps[rsId] += list(rsPoss)

        # Add portions of ROI grid that were not added previously
        for roiElems, roiRsVals, loc in ROIgrids[refCls]:
            remaining = roiElems - coveredElems
            if len(remaining) > 0:
                self._spatialGrids[refCls].append((remaining, roiRsVals, loc))

    def _setUpSelectors(self, selectors):
        self._spatialGrids = {
            ngeom.TetReference: [],
            ngeom.TriReference: [],
            ngeom.VertReference: [],
        }
        self._rs2Grids = {(nsim.MPI.rank, rs._selectorInd): [] for rs in selectors}

        locStr2RefCls = {
            cls._locStr: cls for cls in [ngeom.TetReference, ngeom.TriReference, ngeom.VertReference]
        }
        # Extract info from result selectors
        # elem2infos will contain data using global ids in the case of distributed meshes
        distrElem2infos = {
            ngeom.TetReference._locStr: {},
            ngeom.TriReference._locStr: {},
            ngeom.VertReference._locStr: {},
        }
        nonDistrElem2infos = {
            ngeom.ROI._locStr: {},
            ngeom.Compartment._locStr: {},
            ngeom.Patch._locStr: {},
        }
        for rs in selectors:
            rsId = (nsim.MPI.rank, rs._selectorInd)
            if 'loc_type' in rs.metaData and 'loc_id' in rs.metaData:
                types = rs.metaData['loc_type']
                inds = rs.metaData['loc_id']
                obj_ids = rs.metaData.get('obj_id', None)
                parent_obj_ids = rs.metaData.get('parent_obj_id', None)
                for i, (tpe, ind) in enumerate(zip(types, inds)):
                    val = rs._labels[i].split('.')[-1]
                    if obj_ids is not None:
                        obj = obj_ids[i]
                    elif parent_obj_ids is not None:
                        obj = parent_obj_ids[i]
                    else:
                        obj = None

                    # If the mesh is distributed, we need to use local inds
                    if self._sim._isDistributed() and tpe in locStr2RefCls:
                        ind = locStr2RefCls[tpe]._distCls._getToLocalFunc(self._sim.geom)(ind)

                    if tpe in distrElem2infos:
                        distrElem2infos[tpe].setdefault(ind, []).append((obj, val, rsId, i))
                    elif tpe in nonDistrElem2infos:
                        nonDistrElem2infos[tpe].setdefault(ind, []).append((obj, val, rsId, i))

        if self._sim._isDistributed():
            # Update local nonDistrElem2infos with the ones from other ranks
            import mpi4py
            allNDElem2Infos = mpi4py.MPI.COMM_WORLD.allgather(nonDistrElem2infos)
            for rnk, e2i in enumerate(allNDElem2Infos):
                for locStr, dct in e2i.items():
                    for ind, lst in dct.items():
                        nonDistrElem2infos[locStr].setdefault(ind, [])
                        nonDistrElem2infos[locStr][ind] += lst

        # Merge distributable and non-distributable elem2infos
        elem2infos = {**distrElem2infos, **nonDistrElem2infos}

        rsColRemaps = {}

        # Regions of interest
        # ROIgrids will contain data using local ids in the case of distributed meshes
        ROIgrids = {
            ngeom.TetReference: [],
            ngeom.TriReference: [],
            ngeom.VertReference: [],
        }
        # Pre-split the grids according to compartments and patches
        # self._hierarchicalGrids = {}
        for comp in self._sim.geom.ALL(ngeom.Compartment):
            if self._sim._isDistributed():
                compTets = comp.tets.toLocal()
            else:
                compTets = comp.tets
            if len(compTets) > 0:
                loc = (ngeom.TetReference, comp.name)
                ROIgrids[ngeom.TetReference].append((compTets, [], loc))
        for patch in self._sim.geom.ALL(ngeom.Patch):
            if self._sim._isDistributed():
                patchTris = patch.tris.toLocal()
            else:
                patchTris = patch.tris
            if len(patchTris) > 0:
                loc = (ngeom.TriReference, patch.name)
                ROIgrids[ngeom.TriReference].append((patchTris, [], loc))
        # Treat compartments and patches data as ROI data
        for roiCls in [ngeom.ROI, ngeom.Compartment, ngeom.Patch]:
            for name, valsLst in elem2infos[roiCls._locStr].items():
                zone = getattr(self._sim.geom, name)
                if roiCls == ngeom.ROI:
                    elems = zone[:]
                elif roiCls == ngeom.Compartment:
                    elems = zone.tets
                elif roiCls == ngeom.Patch:
                    elems = zone.tris
                if self._sim._isDistributed():
                    elems = elems.toLocal()
                rsVals = [
                    (obj, val, roiCls._locStr, rsId, rsPos, 1, 1, 'Grid')
                    for obj, val, rsId, rsPos in valsLst
                ]
                newGrids = []
                for elems2, rsVals2, loc2 in ROIgrids[elems._refCls]:
                    inter = elems & elems2
                    elems -= inter
                    elems2 -= inter
                    if len(inter) > 0:
                        newGrids.append((inter, rsVals + rsVals2, loc2))
                    if len(elems2) > 0:
                        newGrids.append((elems2, rsVals2, loc2))
                if len(elems) > 0:
                    newGrids.append((elems, rsVals, (elems._refCls,)))
                ROIgrids[elems._refCls] = newGrids

        # Tetrahedron grids
        self._extractAndAddSpatialGridInfo(elem2infos, ngeom.TetReference, rsColRemaps, ROIgrids)
        # Add empty tet grid at the end
        if self._sim._isDistributed():
            tets = self._sim.geom.tets.toLocal()
        else:
            tets = self._sim.geom.tets
        for tetInds, *_ in self._spatialGrids[ngeom.TetReference]:
            tets -= tetInds
        if len(tets) > 0:
            self._spatialGrids[ngeom.TetReference].append((tets, [], (ngeom.TetReference,)))

        # Triangle grids
        self._extractAndAddSpatialGridInfo(elem2infos, ngeom.TriReference, rsColRemaps, ROIgrids)

        # Vertices
        self._extractAndAddSpatialGridInfo(elem2infos, ngeom.VertReference, rsColRemaps, ROIgrids)

        # Fill rs2Grids map
        for refCls, grids in self._spatialGrids.items():
            for i, (elems, rsVals, _) in enumerate(grids):
                for obj, val, loctpe, rsId, start, step, nVals, center in rsVals:
                    self._rs2Grids.setdefault(rsId, []).append((obj, val, loctpe, start, step, nVals, refCls, i, center))

        if self._shouldWrite:
            # Write column remapping of result selectors, if needed
            for rs in selectors:
                n = rs._getEvalLen()
                colRemap = rsColRemaps.get((nsim.MPI.rank, rs._selectorInd), [])
                if len(colRemap) > 0:
                    if len(colRemap) < n:
                        colRemap += sorted(set(range(n)) - set(colRemap))
                    colRemap = numpy.array(colRemap, dtype=numpy.int64)
                    # Write the column remap if different from neutral remap
                    if any(a != b for a, b in zip(colRemap, range(len(colRemap)))):
                        rsgroup = self._getRsHDFGroup(rs)
                        if _HDF5DataHandler._RS_COLREMAPPING_NAME in rsgroup:
                            if any(a != b for a, b in zip(rsgroup[_HDF5DataHandler._RS_COLREMAPPING_NAME], colRemap)):
                                raise Exception(
                                    f'Column remapping was different for previous runs. Try saving to an '
                                    f'empty HDF5 file.'
                                )
                        else:
                            rsgroup.create_dataset(_HDF5DataHandler._RS_COLREMAPPING_NAME,
                                                   data=colRemap, **self._dataSetKWargs)
        
        if self._sim._isDistributed():
            # Synchronize non-distributed result selector data
            import mpi4py
            allRS2Len = mpi4py.MPI.COMM_WORLD.allgather({rs._selectorInd: rs._getEvalLen() for rs in
                selectors})
            allColRemaps = mpi4py.MPI.COMM_WORLD.allgather(rsColRemaps)

            self._grid2NonDistrRS = {}
            for (rnk, rsIdx), gridVals in self._rs2Grids.items():
                if rnk != nsim.MPI.rank:
                    for obj, val, loctpe, start, step, nVals, gridCls, gridInd, center in gridVals:
                        assert nVals == 1
                        localRemap = allColRemaps[nsim.MPI.rank].get((rnk, rsIdx), None)
                        remoteRemap = allColRemaps[rnk].get((rnk, rsIdx), None)
                        if localRemap is not None:
                            start = localRemap[start]
                        if remoteRemap is not None:
                            start = remoteRemap.index(start)
                        remoteRsLen = allRS2Len[rnk][rsIdx]
                        self._grid2NonDistrRS.setdefault((gridCls, gridInd), []).append(
                            (obj, val, loctpe, start, center, rnk, rsIdx, remoteRsLen)
                        )

    def _getCurrXDMFFilePath(self):
        fileName = XDMFHandler._FILE_NAME_PATTERN.format(
            uid=self._currUID, run=self._currRun, rank=nsim.MPI.rank
        )
        return os.path.join(self._xdmfFolder, fileName)

    def _getCurrFullXDMFFilePath(self):
        fileName = XDMFHandler._FULL_FILE_NAME_PATTERN.format(
            uid=self._currUID, run=self._currRun
        )
        return os.path.join(self._xdmfFolder, fileName)

    def _writeXMLTree(self):
        if self._shouldWrite:
            tree = ElementTree.ElementTree(self._xdmfTree)
            filePath = self._getCurrXDMFFilePath()
            self._savedFilePaths.append(filePath)
            tree.write(filePath)

            if self._fullXdmf is not None:
                tree = ElementTree.ElementTree(self._fullXdmf)
                filePath = self._getCurrFullXDMFFilePath()
                self._savedFilePaths.append(filePath)
                tree.write(filePath)

    @staticmethod
    def _getHierarchicalParent(hierarchicalGrids, path):
        xmlPath = '/Xdmf/Domain/Grid'
        for i in range(len(path)):
            gridName = f'{path[i]._locStr}Grids' if i == 0 else path[i]
            xmlPath += f"/Grid[@Name='{gridName}']"
            if path[:i+1] not in hierarchicalGrids:
                # First element is a reference class, not a string, requires special treatment.
                hierarchicalGrids[path[:i+1]] = ElementTree.SubElement(
                    hierarchicalGrids[path[:i]], 'Grid', Name=gridName, GridType='Collection',
                    CollectionType='Spatial'
                )
        return hierarchicalGrids[path], xmlPath

    @staticmethod
    def _createXMLDocument():
        xdmfTree = ElementTree.Element('Xdmf', Version='2.0')
        xdmfTree.set('xmlns:xi', 'http://www.w3.org/2001/XInclude')
        dom = ElementTree.SubElement(xdmfTree, 'Domain')
        # ElementTree.register_namespace('xi', 'https://www.w3.org/2001/XInclude/')

        hierarchicalGrids = {}
        # Add root
        hierarchicalGrids[tuple()] = ElementTree.SubElement(
            dom, 'Grid', Name=f'SpatialGrids', GridType='Collection', CollectionType='Spatial'
        )
        return xdmfTree, hierarchicalGrids

    def _newRun(self, rid):
        if self._currRun != rid:
            self._temporalGrids = {
                ngeom.TetReference: [],
                ngeom.TriReference: [],
                ngeom.VertReference: [],
            }

            # Write previous XML tree
            if self._xdmfTree is not None:
                self._writeXMLTree()

            # Create XML document
            self._xdmfTree, self._hierarchicalGrids = self._createXMLDocument()

            # Add temporal grids
            allLocations = []
            for refCls, grids in self._spatialGrids.items():
                for i, (elems, _, loc) in enumerate(grids):
                    if self._shouldWrite:
                        if self._sim._isDistributed():
                            gridName = f'{refCls._locStr}Grid{i}_rank{nsim.MPI.rank}'
                        else:
                            gridName = f'{refCls._locStr}Grid{i}'
                        parentGrid, xmlPath = self._getHierarchicalParent(self._hierarchicalGrids, loc)
                        xmlPath += f"/Grid[@Name='{gridName}']"
                        # Add the grid
                        tempGrid = ElementTree.SubElement(
                            parentGrid, 'Grid', Name=gridName, GridType='Collection', CollectionType='Temporal'
                        )
                        allLocations.append((loc, xmlPath))
                    else:
                        tempGrid = None
                        xmlPath = None
                    self._temporalGrids[refCls].append((tempGrid, xmlPath, None, None))
                    # Write the grids for t=0, this way we are sure that all grids are written at least once
                    # even if no data is associated to it.
                    self._writeGrid(refCls, i, 0, 0)

            self._currTime = None
            self._currRun = rid

            if self._sim._isDistributed():
                # Gather allLocations to rank 0
                import mpi4py
                elemStr2ElemCls = {elemCls._locStr: elemCls for elemCls in self._temporalGrids.keys()}
                allLocations = [((loc[0]._locStr,) + loc[1:], path) for loc, path in allLocations]
                fileName = XDMFHandler._FILE_NAME_PATTERN.format(
                    uid=self._currUID, run=self._currRun, rank=nsim.MPI.rank
                )
                allInfos = mpi4py.MPI.COMM_WORLD.gather((fileName, allLocations), root=0)

                if nsim.MPI.rank == 0:
                    # Write a common xmf file on rank 0
                    self._fullXdmf, fullhierarchicalGrids = self._createXMLDocument()
                    for fileName, locations in allInfos:
                        for loc, xmlPath in locations:
                            loc = (elemStr2ElemCls[loc[0]],) + loc[1:]
                            parent, _ = self._getHierarchicalParent(fullhierarchicalGrids, loc)
                            ElementTree.SubElement(parent, 'xi:include', href=fileName,
                                                   xpointer=f'xpointer({xmlPath})')

    def _getHyperSlab(self, parent, fileName, rsIdx, rsLen, tind, start, step, nVals):
        hyperslab = ElementTree.SubElement(
            parent, 'DataItem', ItemType='HyperSlab', Dimensions=f'{nVals}'
        )
        dims = ElementTree.SubElement(
            hyperslab, 'DataItem', Dimensions='3 2', NumberType='Int', Format='XML'
        )
        dims.text = f'{tind} {start} 1 {step} 1 {nVals}'
        data_item = ElementTree.SubElement(
            hyperslab, 'DataItem', DataType='Float', Dimensions=f'{tind+1} {rsLen}',
            Format='HDF', Precision='8'
        )
        rsName = HDF5Handler._RS_GROUP_NAME.format(rsIdx)
        hdffn = os.path.basename(fileName)
        data_item.text = f'{hdffn}:{self._currUID}/{rsName}/runs/Run_{self._currRun}/data'
        return hyperslab

    def _newTimeStep(self, t, rs, tind):
        rsId = (nsim.MPI.rank, rs._selectorInd)
        for obj, val, loctpe, start, step, nVals, gridCls, gridInd, center in self._rs2Grids[rsId]:
            tempGrid, xmlPath, currTime, currGrid = self._temporalGrids[gridCls][gridInd]
            if currTime != t:
                currGrid = self._writeGrid(gridCls, gridInd, t, tind)

                # Only triggered once per grid and per timestep
                if self._sim._isDistributed():
                    # Add attributes that are not in this rank
                    if (gridCls, gridInd) in self._grid2NonDistrRS:
                        for _obj, _val, _loctpe, _start, _center, rnk, rsIdx, rsLen in self._grid2NonDistrRS[(gridCls, gridInd)]:
                            att = ElementTree.SubElement(
                                currGrid, 'Attribute', Name=self._getAttributeName(_obj, _val, _loctpe), AttributeType='Scalar',
                                Center=_center
                            )
                            fileName = HDF5Handler._DISTRIBUTED_HDF_SUFFIX.format(
                                self._pathPrefix, rnk) + HDF5Handler._HDF_EXTENSION
                            hyperslab = self._getHyperSlab(att, fileName, rsIdx, rsLen, tind, _start, 1, 1)

            if self._shouldWrite:
                att = ElementTree.SubElement(
                    currGrid, 'Attribute', Name=self._getAttributeName(obj, val, loctpe), AttributeType='Scalar',
                    Center=center
                )
                if center == 'Grid':
                    # Save the value in the xdmf file directly instead of using a hyperslab
                    data_item = ElementTree.SubElement(
                        att, 'DataItem', DataType='Float', Dimensions=f'1', Format='XML'
                    )
                    data_item.text = str(rs.data[self._currRun, tind, start])
                else:
                    hyperslab = self._getHyperSlab(
                        att, self._path, rs._selectorInd, rs._getEvalLen(), tind, start, step, nVals
                    )

    def _writeGrid(self, gridCls, gridInd, t, tind):
        if self._shouldWrite:
            if XDMFHandler._MESH_GROUP_NAME not in self._currGroup:
                # Write mesh to HDF5 file
                meshGroup = self._currGroup.create_group(XDMFHandler._MESH_GROUP_NAME)
            else:
                meshGroup = self._currGroup[XDMFHandler._MESH_GROUP_NAME]
        else:
            meshGroup = None

        # Write xdmf description of mesh data
        tpeMap = {
            ngeom.TetReference: (6, 5),
            ngeom.TriReference: (4, 4),
            ngeom.VertReference: (1, 1),
        }
        if self._sim._isDistributed():
            gridName = f'{gridCls._locStr}Grid{gridInd}_rank{nsim.MPI.rank}'
        else:
            gridName = f'{gridCls._locStr}Grid{gridInd}'
        elemCode, elemColNb = tpeMap[gridCls]

        cond = f'{gridName}/XYZ' not in meshGroup if meshGroup is not None else None
        if not self._sim._isDistributed():
            # If the mesh is not distributed, all ranks need to be involved in the calls to get mesh data
            import mpi4py
            cond = mpi4py.MPI.COMM_WORLD.bcast(cond, root=0)

        elems, _, loc = self._spatialGrids[gridCls][gridInd]
        if cond:
            allVerts = []
            topo = []
            if gridCls == ngeom.VertReference:
                allVerts = [numpy.array(vert) for vert in elems]
            else:
                vertInds = {}
                for elem in elems:
                    localInds = []
                    for vert in elem.verts:
                        if vert not in vertInds:
                            vertInds[vert] = len(allVerts)
                            vertPos = numpy.array(vert)
                            if self._shouldWrite:
                                allVerts.append(vertPos)
                        localInds.append(vertInds[vert])
                    if self._shouldWrite:
                        topo.append(numpy.array([elemCode] + localInds))

            if self._shouldWrite:
                meshGroup.create_dataset(f'{gridName}/XYZ', data=numpy.array(allVerts), **self._dataSetKWargs)
                if len(topo) > 0:
                    meshGroup.create_dataset(f'{gridName}/topology',
                                             data=numpy.array(topo), **self._dataSetKWargs)

        if self._shouldWrite:
            tempGrid, xmlPath, _, currGrid = self._temporalGrids[gridCls][gridInd]

            grid = ElementTree.SubElement(
                tempGrid, 'Grid', Name=f'{gridName}_{tind}', GridType='Uniform'
            )
            ElementTree.SubElement(grid, 'Time', Value=f'{t}')

            if currGrid is None:
                nverts = meshGroup[f'{gridName}/XYZ'].shape[0]
                if f'{gridName}/topology' in meshGroup:
                    nelems = meshGroup[f'{gridName}/topology'].shape[0]
                else:
                    nelems = nverts
                if nelems != len(elems):
                    raise Exception(
                        f'Previous simulations were run with a different XDMF mesh splitting. '
                        f'Try saving to a different HDF5 file.'
                    )

                hdfFileName = os.path.basename(self._path)
                if gridCls == ngeom.VertReference:
                    topo = ElementTree.SubElement(
                        grid, 'Topology', TopologyType='PolyVertex', NumberOfElements=f'{nelems}'
                    )
                else:
                    topo = ElementTree.SubElement(
                        grid, 'Topology', TopologyType='Mixed', Dimensions=f'{nelems}'
                    )
                    data_item = ElementTree.SubElement(
                        topo, 'DataItem', DataType='Int', Dimensions=f'{nelems * elemColNb}', Format='HDF',
                        Precision='8'
                    )
                    data_item.text = f'{hdfFileName}:{self._currUID}/{XDMFHandler._MESH_GROUP_NAME}/{gridName}/topology'

                geo = ElementTree.SubElement(grid, 'Geometry', GeometryType='XYZ')
                data_item = ElementTree.SubElement(
                    geo, 'DataItem', DataType='Float', Dimensions=f'{nverts * 3}', Format='HDF', Precision='8'
                )
                data_item.text = f'{hdfFileName}:{self._currUID}/{XDMFHandler._MESH_GROUP_NAME}/{gridName}/XYZ'
            else:
                topo = ElementTree.SubElement(grid, 'Topology', Reference='XML')
                topo.text = f'{xmlPath}/Grid/Topology'
                geo = ElementTree.SubElement(grid, 'Geometry', Reference='XML')
                geo.text = f'{xmlPath}/Grid/Geometry'
        else:
            tempGrid, xmlPath, grid = None, None, None

        self._temporalGrids[gridCls][gridInd] = (tempGrid, xmlPath, t, grid)
        return grid
