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

"""Unit tests for the steps.API_2.saving._HDF5CompoundObjHandler class."""

import importlib.util
import numpy as np
import os
import tempfile
import unittest

from steps import interface

from steps.saving import _HDF5CompoundObjHandler


# Stub classes
class DbhStub(object):
    def __init__(self, group):
        self._dataSetKWargs = {}


class HDF5CompoundObjects(unittest.TestCase):
    """Test saving compound objects to HDF5 files"""
    def setUp(self):
        self.createdFiles = set()

    def tearDown(self):
        super().tearDown()
        for path in self.createdFiles:
            if os.path.isfile(path):
                os.remove(path)

    def checkWriteRead(self, obj, callback=None, multi=False, **kwargs):
        import h5py

        _, path = tempfile.mkstemp(prefix=f'{self.__class__.__name__}testFile', suffix='.h5')
        self.createdFiles.add(path)

        groupname = 'testgroup'
        with h5py.File(path, 'w') as f:
            grp = f.create_group(groupname)

            compObjHandler = _HDF5CompoundObjHandler(grp, DbhStub(grp), **kwargs)

            if multi:
                inds = []
                for ob in obj:
                    inds.append(compObjHandler.write(ob))
            else:
                ind = compObjHandler.write(obj)

            if callback is not None:
                callback(compObjHandler)

        with h5py.File(path, 'r') as f:
            grp = f[groupname]

            compObjHandler = _HDF5CompoundObjHandler(grp, DbhStub(grp), **kwargs)

            if multi:
                obj2 = []
                for ind in inds:
                    obj2.append(compObjHandler.read(ind))
            else:
                obj2 = compObjHandler.read(ind)

        if isinstance(obj2, np.ndarray):
            obj2 = list(obj2)
        self.assertEqual(obj, obj2)

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testSingleNumbers(self):
        self.checkWriteRead(2)
        self.checkWriteRead(6.25)

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testNumbersList(self):
        self.checkWriteRead(list(range(10)))
        self.checkWriteRead([1.5, 2.5, 8.5, 2.3])
        self.checkWriteRead([1.5, 2.5, 8.5, 2.3] + list(range(5)))

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testStrings(self):
        self.checkWriteRead('test')
        self.checkWriteRead('string with unicode characters 2.0±0.1')
        self.checkWriteRead('test\nline break')

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testListOfLists(self):
        self.checkWriteRead([f'Elem{i}' for i in range(5)])

        self.checkWriteRead([list(range(10)) for i in range(5)])
        self.checkWriteRead([[list(range(10)) for j in range(5)] for i in range(5)])
        self.checkWriteRead([[i * j / 7 for j in range(5, 10+i)] for i in range(5)])
        self.checkWriteRead([[list(range(i * j  + 2)) for j in range(5, 10+i)] for i in range(5)])

        # Mixed types
        self.checkWriteRead(list(range(5)) + [list(range(2 * i + 1)) for i in range(5)])
        self.checkWriteRead(list(range(5)) + [list(range(2 * i + 1)) for i in range(5)] + [1.2 + 3.2])
        self.checkWriteRead(list(range(5)) + [list(range(2 * i + 1)) for i in range(5)] + [1.2 + 3.2] + ['test'])

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testDict(self):
        self.checkWriteRead({i: 2 * i for i in range(5)})
        self.checkWriteRead({i: i / 7 for i in range(5)})
        self.checkWriteRead({i: list(range(i)) for i in range(5)})
        self.checkWriteRead({i: [1.2, 3.2] + list(range(i)) for i in range(5)})

        # String keys
        self.checkWriteRead({f'{i}': 2 * i for i in range(5)})

        # Nested dicts
        self.checkWriteRead({i: {j / 7: list(range(j)) for j in range(2*i + 5)} for i in range(5)})
        self.checkWriteRead({i: {j / 7: [{k: 3*k if k % 2 == 0 else list(range(k))} for k in range(j)] for j in range(2*i + 5)} for i in range(5)})

        # Empty dict
        self.checkWriteRead({})

    def _getCallBack(self, tpes, expectedSizes, expectedNumbers):
        def checkCaching(handler):
            for tpe, exps, expn in zip(tpes, expectedSizes, expectedNumbers):
                if tpe in handler._dsets:
                    dset = handler._dsets[tpe]
                    self.assertEqual(len(dset), exps)
                nb = sum(handler._compDset[:,0] == tpe)
                self.assertEqual(nb, expn)
        return checkCaching

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testCaching(self):
        L, N, M = 3, 7, 10
        data = [{j: [1 / (k + 1) for k in range(L)] for j in range(N)} for i in range(M)]

        values = [
            (
                [_HDF5CompoundObjHandler._DATA_TYPE.INT],
                (
                    [N * M],
                    [M]
                ),
                (
                    [N],
                    [1]
                ),
            ),
            (
                [_HDF5CompoundObjHandler._DATA_TYPE.FLOAT],
                (
                    [L * N * M],
                    [N * M]
                ),
                (
                    [L],
                    [1]
                ),
            ),
            (
                [
                    _HDF5CompoundObjHandler._DATA_TYPE.INT,
                    _HDF5CompoundObjHandler._DATA_TYPE.FLOAT,
                    _HDF5CompoundObjHandler._DATA_TYPE.LIST,
                ],
                (
                    [N * M, L * N * M, 3 * M + N * M],
                    [M, N * M, 1 + M]
                ),
                (
                    [N, L, 3 * M + N],
                    [1, 1, 2]
                ),
            ),
            (
                [
                    _HDF5CompoundObjHandler._DATA_TYPE.INT,
                    _HDF5CompoundObjHandler._DATA_TYPE.FLOAT,
                    _HDF5CompoundObjHandler._DATA_TYPE.LIST,
                    _HDF5CompoundObjHandler._DATA_TYPE.DICT,
                ],
                (
                    [N * M, L * N * M, 3 * M + N * M, M],
                    [M, N * M, 1 + M, M],
                ),
                (
                    [N, L, M + 2 + N, 1],
                    [1, 1, 2, 1],
                ),
            ),
        ]
        for tpes, nocache, cache in values:
            self.checkWriteRead(
                data,
                callback=self._getCallBack(tpes, *nocache), cachedTypes=[]
            )
            self.checkWriteRead(
                data,
                callback=self._getCallBack(tpes, *cache), cachedTypes=tpes
            )

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testStringCaching(self):
        L, N, M = 3, 7, 10
        string = 'test'
        slen = len(string)
        data = [{j: [string for k in range(L)] for j in range(N)} for i in range(M)]

        self.checkWriteRead(
            data,
            callback=self._getCallBack(
                [_HDF5CompoundObjHandler._DATA_TYPE.STRING],
                [L * N * M * slen],
                [L * N * M]
            ), cachedTypes=[]
        )
        self.checkWriteRead(
            data,
            callback=self._getCallBack(
                [_HDF5CompoundObjHandler._DATA_TYPE.STRING],
                [slen],
                [1]
            ), cachedTypes=[_HDF5CompoundObjHandler._DATA_TYPE.STRING]
        )

    @unittest.skipIf(importlib.util.find_spec('h5py') is None, 'h5py not available')
    def testSeveralObjects(self):
        data = [
            'string',
            list(range(10)),
            [{j: [1 / (k + 1) for k in range(1 + 2*j)] for j in range(i + 5)} for i in range(10)]
        ]
        self.checkWriteRead(data, multi=True)


def suite():
    all_tests = []
    all_tests.append(unittest.makeSuite(HDF5CompoundObjects, "test"))
    return unittest.TestSuite(all_tests)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
