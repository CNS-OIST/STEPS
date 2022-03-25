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

import collections
import copy
import datetime
import math
import numbers
import operator
import pickle
import sqlite3
import struct
import sys
import warnings

import numpy
import steps

from . import utils as nutils
from . import sim as nsim


__all__ = [
    'ResultSelector',
    'SQLiteDBHandler',
    'SQLiteGroup',
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

    def __getitem__(self, key):
        """Return the metadata corresponding to key."""
        self._checkKey(key)
        if key not in self._dict:
            raise KeyError(f'No metadata corresponding to key {key}.')
        return numpy.array(self._dict[key])

    def __setitem__(self, key, val):
        """Set the metadata corresponding to key."""
        self._checkKey(key)
        if self._parent._savingStarted():
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

        if key in self._dict and self._dict[key] != lst:
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
                                                     # value that corresponds to the sum of S1 and
                                                     # S2 counts in comp1.

        rs2 = rs.comp1.S1.Count << rs.comp1.S2.Count # This one will save 2 values, the count of S1
                                                     # in comp1 and the count of S2 in comp1.

    Result selectors can also transform data and only save the result of the transformation::

        rs3 = rs.SUM(rs.TETS(tetlst).S1.Count) # This will save only one value: the total number of
                                               # S1 in all the tetrahedrons in tetlst.

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

        self._labels = None
        self._metaData = _MetaData(self)

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

            >>> rs1.labels # Default values, established from the simulation paths used for saving
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

    def _addedToSimulation(self, ind):
        """Specify that the result selector was added to a simulation with index ind."""
        self._checkComplete()
        self._addedToSim = True
        self._selectorInd = ind

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

            rs3 = 10 * rs.TETS(tetlst).S1.Count                       # rs3 will save the number of
                                                                      # S1 in each tetrahedron in
                                                                      # tetlst, multiplied by 10.

            rs4 = rs.TETS(tetlst).S1.Count * rs.TETS(tetlst).S2.Count # rs4 will save the product
                                                                      # of the number of S1 and the
                                                                      # number of S2 in each
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

            rs3 = rs.TETS(tetlst).S1.Count / 10                       # rs3 will save the number of
                                                                      # S1 in each tetrahedron in
                                                                      # tetlst, divided by 10.

            rs4 = 1 / rs.TETS(tetlst).S1.Count                        # rs4 will save the inverse
                                                                      # of the number of S1 in each
                                                                      # tetrahedron in tetlst,
                                                                      # divided by 10.

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

            rs3 = 10 + rs.TETS(tetlst).S1.Count                       # rs3 will save the number of
                                                                      # S1 in each tetrahedron in
                                                                      # tetlst, increased by 10.

            rs4 = rs.TETS(tetlst).S1.Count + rs.TETS(tetlst).S2.Count # rs4 will save the sum
                                                                      # of the number of S1 and the
                                                                      # number of S2 in each
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

            rs3 = rs.TETS(tetlst).S1.Count - 10                       # rs3 will save the number of
                                                                      # S1 in each tetrahedron in
                                                                      # tetlst, minus 10.

            rs4 = 10 - rs.TETS(tetlst).S1.Count                       # rs4 will save 10 minus the
                                                                      # number of S1 for each
                                                                      # tetrahedron in tetlst.

            rs5 = rs.TETS(tetlst).S1.Count - rs.TETS(tetlst).S2.Count # rs5 will save the number of
                                                                      # S1 minus the number of S2
                                                                      # in each tetrahedron in
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

            rs3 = rs.TETS(tetlst).S1.Count ** 2                       # rs3 will save the square of
                                                                      # the number of S1 in each
                                                                      # tetrahedron in tetlst.

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
        self._compiledFuncs = None
        self._len = None
        self._descr = [str(self.sim)]

        # labels are only computed once the path is fully defined

    def _strDescr(self):
        """Return a generic description of the ResultSelector."""
        return '.'.join(self._descr[1:])

    def _evaluate(self, solvStateId=None):
        """Return a list of the values to save."""
        return [f(*args, **kwargs) for f, args, kwargs in self._compiledFuncs]

    def _getEvalLen(self):
        """Return the number of values that _evaluate() will return."""
        return self._len

    def _checkComplete(self):
        """Raise an exception if the path is not complete."""
        if self._compiledFuncs is None:
            raise Exception(f'{self} is incomplete.')

    def _concat(self, other):
        """Concatenate two result selectors into a _ResultList."""
        self._checkComplete()
        return super()._concat(other)

    def _binaryOp(self, other, op, symetric=False, opStr='{0} {1}'):
        """Return a _ResultCombiner that represents the binary operation op."""
        self._checkComplete()
        return super()._binaryOp(other, op, symetric, opStr)

    def __getattr__(self, name):
        if name not in nsim.SimPath._endNames:
            self.simpath = getattr(self.simpath, name)
        else:
            descr = getattr(nsim.SimPath, name)
            self._compiledFuncs = descr._getFinalPathsFunction(self.simpath)
            self._len = len([p for p in self.simpath])
            self._endName = name

            # Compute labels
            self._labels = []
            for descr in self.simpath._getDescriptions(tuple(self._descr)):
                self._labels.append('.'.join(descr[1:] + (name,)))

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

            for key, lst in mtdt.items():
                self.metaData[key] = lst

        self._descr.append(name)
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

        self._evalLen = 0
        for c in self.children:
            self._evalLen += c._getEvalLen()

        # concatenate labels
        self._labels = []
        for c in self.children:
            self._labels += c.labels

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
                    # children metadata should already be computed, so we can access _dict directly
                    lst += c.metaData._dict[key]
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


# Forbid the creation of objects that would have the name of a ResultSelector method or property
# because they would then be inaccessible through a SimPath (methods and attributes have priority
# over __getattr__).
for _name in dir(ResultSelector) + dir(_ResultPath) + dir(_ResultList) + dir(_ResultCombiner):
    nutils.NamedObject._forbiddenNames.add(_name)


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


class _DBDatahandler(_DataHandler):
    pass


class _SQLiteDataHandler(_DBDatahandler):
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
        return numpy.array(_sliceData(self._data, key))

    def __array__(self):
        return numpy.array(self._data)

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
                return numpy.array(_sliceData(self._parentHandler.saveTime, key))
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
                                return numpy.array([[row[key[2]] for row in res]])
                            else:
                                return numpy.array([[[row[key[2]]] for row in res]])
                        else:
                            mk = (slice(None) if isinstance(key[1], slice) else 0, key[2])
                            return numpy.array(_sliceData(res, mk))

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
            except UnexpectedEnd:
                if nxt == 0:
                    break
                else:
                    raise IndexError(
                        f'Could not load time slice {ti} of run {ind} from {self._fp}.'
                        f' The file might be corrupted.'
                    )

        if forceArray:
            return numpy.array(res)
        mk = tuple(slice(None) if isinstance(k, slice) else 0 for k in key)
        return numpy.array(_sliceData(res, mk))

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
            raise UnexpectedEnd()

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
            raise UnexpectedEnd()

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
            return numpy.array(res)
        mk = tuple(slice(None) if isinstance(k, slice) else 0 for k in key)
        return numpy.array(_sliceData(res, mk))

    def __len__(self):
        return self._dbh._conn.execute(f'SELECT MAX(runid) FROM {self._tabName}').fetchone()[0] + 1

    def __array__(self):
        return self.__getitem__(slice(None, None, None), forceArray=True)


###################################################################################################
# Database handlers


class DatabaseHandler:
    """Base class for all database handlers."""

    def __init__(self, *args, **kwargs):
        pass

    def _getDataHandler(self, rs):
        """Return a _DBDatahandler for ResultSelector rs."""
        pass


class SQLiteDBHandler(DatabaseHandler):
    """SQLite database handler

    :param path: The path to the SQLite database file
    :type path: str
    :param \*args: Transmitted to :py:func:`sqlite3.connect`, see
        `documentation <https://docs.python.org/3/library/sqlite3.html#sqlite3.connect>`_ for
        details
    :param commitFreq: How frequently the data should be committed to the database. For example,
        this value is set to 10 by default which means that every 10 saving events, the data will
        be committed. If a result selector is saved every 10ms, it means the data will be committed
        to database every 100ms.
    :type commitFreq: int
    :param \*\*kwargs: Transmitted to :py:func:`sqlite3.connect`, see
        `documentation <https://docs.python.org/3/library/sqlite3.html#sqlite3.connect>`_ for
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
                                                          # database with identifier 'MySimulation'
                                                          # and save additional parameters val1 and
                                                          # val2.

            for i in range(NBRUNS):                       # Run a series of runs, all of them being
                sim.newRun()                              # associated to the 'MySimulation' group.
                ...
                sim.run(...)

    Note that after calling `sim.toDB(...)` it is still possible to force the saving of some result
    selectors to files by calling ``toFile(...)`` on them. Result selectors that contain a high
    number of values to save are probably better saved to a file. The name of the file can be added
    as a keyword parameter to the ``simtoDB(...)`` call to simplify loading.

    Usage when accessing data from the database::

        with SQLiteDBHandler(dbPath) as dbh:              # Create database handler.

            val1 = dbh['MySimulation'].val1               # Querying a parameter value from the
                                                          # 'MySimulation' group.

            rs1, rs2, rs3 = dbh['MySimulation'].results   # Querying the result selectors that were
                                                          # saved for the 'MySimulation' group.
                                                          # They are returned in the same order as
                                                          # they were added to the simulation.

            plt.plot(rs1.time[0], rs1.data[0])            # The results selectors can be used as if
                                                          # they had been declared in the same
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

    def _getDataHandler(self, rs):
        """Return a _DBDatahandler for ResultSelector rs."""
        return _SQLiteDataHandler(rs, self, self._commitFreq)

    def _newGroup(self, uid, selectors, **kwargs):
        """Initialize the database and add a new run group."""
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


class SQLiteGroup:
    """A class representing a group of runs in a SQLite database

    .. note::
        This class should never be instantiated by the user, it is obtained through
        :py:class:`SQLiteDBHandler` instead.
    """

    def __init__(self, dbh, row):
        self._dbh = dbh
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
        """
        if name not in self._dict:
            raise AttributeError(f'{name} is not an attribute of {self}.')
        return self._dict[name]

    @property
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
