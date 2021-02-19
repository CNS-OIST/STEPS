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
 
import itertools

import numpy

from steps.API_1.geom import INDEX_DTYPE

from . import saving as nsaving
from . import utils as nutils


# Nothing to import, this module is not supposed to be imported by users
__all__ = []

###################################################################################################
# Result Path optimizations


def getGroupedBatchCallNodes(allNodes):
    solvCalls = []

    npCall = None
    for i, node in enumerate(allNodes):
        if isinstance(node, _OptimCallNode) and node.hasBatchCall():
            call, ind = node.call
            if call == npCall:
                groupCall, currSlice = solvCalls[-1]
                if isinstance(groupCall._slc, slice):
                    # If the indices cannot be expressed as a slice, convert to list
                    if (groupCall._slc.stop - 1 >= ind) or (
                        groupCall._slc.step is not None
                        and (ind - (groupCall._slc.stop - 1)) != groupCall._slc.step
                    ):
                        rangeArgs = (groupCall._slc.start, groupCall._slc.stop)
                        if groupCall._slc.step is not None:
                            rangeArgs += (groupCall._slc.step,)
                        groupCall._slc = list(range(*rangeArgs)) + [ind]
                    else:
                        if groupCall._slc.step is None:
                            groupCall._slc = slice(
                                groupCall._slc.start, ind + 1, ind - (groupCall._slc.stop - 1)
                            )
                        else:
                            groupCall._slc = slice(groupCall._slc.start, ind + 1, groupCall._slc.step)
                else:
                    groupCall._slc.append(ind)
                solvCalls[-1][1] = slice(currSlice.start, currSlice.stop + 1)
            else:
                npCall = call
                solvCalls.append([_OptimGroupedCallNodes(npCall, slice(ind, ind + 1)), slice(i, i + 1)])
        else:
            npCall = None
            solvCalls.append((node, slice(i, i + 1)))

    return solvCalls


class _OptimizedResultPath:
    """
    Wrapper class around a ResultPath that redirects part or all of the value evaluations
    to an optimized solver call.
    """

    def __init__(self, allNodes):
        self._len = len(allNodes)
        self._solvCalls = getGroupedBatchCallNodes(allNodes)

        self._array = numpy.zeros(self._len, dtype=numpy.float64)

    def _evaluate(self, solvStateId=None):
        for call, rslice in self._solvCalls:
            self._array[rslice] = call._evaluate(solvStateId)
        return self._array

    # TODO Not urgent: If depends on a single other, just return call._evaluate(solvStateId)


##########


class _OptimizedNPSolverCall:
    """
    Class handling the numpy array associated to a specific optimized solver call.
    Makes sure that the call is only executed once even if several Optimized result paths depend
    on it.
    """

    def __init__(self, func, args, kwargs, sz):
        self._func = func
        self._args = self._processArgs(args)
        self._kwargs = kwargs

        self._stateId = None
        self._array = numpy.zeros(sz, dtype=numpy.float64)

    def _processArgs(self, args):
        new_args = []
        for arg in args:
            if isinstance(arg, _TrueParamList) and isinstance(arg[0], int):
                new_args.append(numpy.array(arg, dtype=INDEX_DTYPE))
            else:
                new_args.append(arg)
        return new_args

    def _evaluate(self, solvStateId):
        if solvStateId != self._stateId:
            self._func(*self._args, self._array, **self._kwargs)
            self._stateId = solvStateId
        return self._array


class _StandardSolverCall:
    def __init__(self, func, args, kwargs):
        self._func = func
        self._args = args
        self._kwargs = kwargs

        self._stateId = None
        # Assume that standard solver calls always output a single float
        self._array = numpy.array([0], dtype=numpy.float64)

    def _evaluate(self, solvStateId):
        if solvStateId != self._stateId:
            self._array[0] = self._func(*self._args, **self._kwargs)
            self._stateId = solvStateId
        return self._array


class _FuncSolverCall:
    def __init__(self, func, optimFuncNodes):
        self._func = func
        self._optimFuncNodes = optimFuncNodes

        self._stateId = None
        self._paramVals = numpy.zeros(self._optimFuncNodes[-1][1].stop, dtype=numpy.float64)
        self._array = numpy.array([0], dtype=numpy.float64)

    def _evaluate(self, solvStateId):
        if solvStateId != self._stateId:
            for call, rslice in self._optimFuncNodes:
                self._paramVals[rslice] = call._evaluate(solvStateId)
            self._array[0] = self._func(self._paramVals)
            self._stateId = solvStateId
        return self._array


##########
# Nodes that can be part of allValues are _OptimCallNode and _OptimFuncNode,
# grouped node is for batch NP calls


class _OptimCallNode:
    def __init__(self, name, args, kwargs, endName):
        self.name = name
        self.args = args
        self.kwargs = kwargs
        self.endName = endName
        self.call = (
            None,
            None,
        )  # Tuple[_OptimizedNPSolverCall, relative ind] the index is None if the call is a standard call

    def __hash__(self):
        return hash((self.name, self.args, frozenset(self.kwargs.items()), self.endName))

    def __eq__(self, other):
        return isinstance(other, _OptimCallNode) and (
            self.name,
            self.args,
            frozenset(self.kwargs.items()),
            self.endName,
        ) == (other.name, other.args, frozenset(other.kwargs.items()), other.endName)

    def hasBatchCall(self):
        return self.call[1] is not None

    def _evaluate(self, solvStateId):
        return self.call[0]._evaluate(solvStateId)


class _OptimGroupedCallNodes:
    def __init__(self, optimNPcall, slc):
        self._call = optimNPcall
        self._slc = slc

    def _evaluate(self, solvStateId):
        return self._call._evaluate(solvStateId)[self._slc]

    # def __init__(self, optimNPcall, start, end, step):
    # self._call = optimNPcall
    # self._start = start
    # self._end = end
    # self._step = step

    # def _evaluate(self, solvStateId):
    # return self._call._evaluate(solvStateId)[self._start:self._end:self._step]


class _OptimFuncNode:
    def __init__(self, func, indices, allValues):
        self.func = func
        self.indices = indices
        self.allValues = allValues

        self.call = None

    def __hash__(self):
        return hash((self.func, self.indices))

    def __eq__(self, other):
        return isinstance(other, _OptimFuncNode) and (self.func, self.indices) == (other.func, other.indices)

    def hasBatchCall(self):
        return any(self.allValues[i].hasBatchCall() for i in self.indices)

    def _setUp(self):
        self.call = _FuncSolverCall(
            self.func, getGroupedBatchCallNodes([self.allValues[i] for i in self.indices])
        )

    def _evaluate(self, solvStateId):
        return self.call._evaluate(solvStateId)


##########


class _ParamList:
    def __init__(self, lst):
        self.lst = lst

    def __iter__(self):
        return iter(self.lst)

    def __hash__(self):
        return hash(frozenset(self.lst))

    def __eq__(self, other):
        return isinstance(other, _ParamList) and frozenset(self.lst) == frozenset(other.lst)

    def __repr__(self):
        return 'PL' + str(self.lst)


class _TrueParamList(list):
    def __init__(self, lst):
        super().__init__(lst)


class _OptimValuesList:
    def __init__(self, sim):
        self.lst = []
        self.node2Ind = {}
        self.solver = sim.solver

    def addNode(self, node):
        if node not in self.node2Ind:
            self.node2Ind[node] = len(self.lst)
            self.lst.append(node)
        return self.node2Ind[node]

    def __getitem__(self, ind):
        return self.lst[ind]

    def setUpCallNodes(self):
        # group by func name and kwargs
        key2Nodes = {}
        for node in self.lst:
            if isinstance(node, _OptimCallNode):
                key = (node.name, node.endName, frozenset(node.kwargs.items()))
                key2Nodes.setdefault(key, []).append(node)

        call2Batch = {}
        for (name, endName, kwargs), nodes in key2Nodes.items():
            kwargs = dict(kwargs)
            params = [node.args for node in nodes]
            groupings = self.groupParams(params)
            for group in groupings:
                for funcName, args, corresCalls in self.getFunctions(group, name, endName):
                    # Retrieve the indexes of the corresponding single solver calls in allValues
                    func = getattr(self.solver, funcName)
                    if len(corresCalls) > 1:
                        # If we are calling a batch function, map all original call to the optimized one
                        optCall = _OptimizedNPSolverCall(func, args, kwargs, len(corresCalls))
                        for i, corrArgs in enumerate(corresCalls):
                            node = _OptimCallNode(name, corrArgs, kwargs, endName)
                            self.lst[self.node2Ind[node]].call = (optCall, i)
                    else:
                        node = _OptimCallNode(name, args, kwargs, endName)
                        call = _StandardSolverCall(func, args, kwargs)
                        self.lst[self.node2Ind[node]].call = (call, None)

        # Set up optimizations for func nodes
        for node in self.lst:
            if isinstance(node, _OptimFuncNode):
                node._setUp()
                # TODO Not urgent: further optimizations?

    def getFunctions(self, group, name, endName):
        """Take a group as parameter and output a list of function and their arguments"""
        maskGenerator = itertools.product(
            *[[True, False] if isinstance(val, _ParamList) else [False] for val in group]
        )

        allCalls = []
        for mask in maskGenerator:
            funcName = 'get'
            funcArgs = tuple()
            nbCalls = 1
            batchCall = False
            for subName, useBatch, val in zip(name, mask, group):
                if isinstance(val, _ParamList):
                    if useBatch:
                        funcName += 'Batch' + subName
                        funcArgs += (_TrueParamList(val.lst),)
                        batchCall = True
                    else:
                        funcName += subName
                        funcArgs += (val,)
                        nbCalls *= len(val.lst)
                else:
                    funcName += subName
                    funcArgs += (val,)

            funcName += endName
            if batchCall:
                funcName += 'sNP'
            if hasattr(self.solver, funcName):
                allCalls.append((nbCalls, funcName, funcArgs))

        if len(allCalls) == 0:
            raise Exception(f'The solver does not implement get{"".join(name)}{endName}(...).')

        _, funcName, funcArgs = min(allCalls, key=lambda x: x[0])
        combargs = [arg.lst if isinstance(arg, _ParamList) else [arg] for arg in funcArgs]
        for args in itertools.product(*combargs):
            combSingleCall = [a if isinstance(a, _TrueParamList) else [a] for a in args]
            correspondingCallArgs = list(itertools.product(*combSingleCall))
            yield funcName, args, correspondingCallArgs

    @staticmethod
    def groupParams(params):
        """
        Take a list of parameter tuples as input and outputs a list of tuples of _ParamLists
        whose cartesian product yields the orginal parameter tuples.
        """
        n = len(params[0])
        groupings = [params]
        for i in range(n):
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
                groupings.append(_OptimValuesList.groupParams(newParams))
        return min(groupings, key=lambda x: len(x))


def _ExtractOptimValuesFromPath(path, rs, allValues, rs2IndLst):
    def wrap(acc, f):
        return lambda x: acc(f(x))

    if isinstance(path, nutils.SimPathCombiner):
        # If the path is a SimPathCOmbiner, add all dependent values and add the combiner as an _OptimFuncNode
        indices = []
        for path2 in path:
            _ExtractOptimValuesFromPath(path2, rs, allValues, indices)
        ind = allValues.addNode(_OptimFuncNode(path.func, tuple(indices), allValues))
    else:
        funcName = tuple()
        funcArgs = tuple()
        funcKwArgs = {}
        modifier = None
        for e in path[1:]:
            funcName += (e._solverStr(),)
            arg = e._solverId()
            if len(arg) > 1:
                funcArgs += (arg,)
            else:
                funcArgs += arg
            funcKwArgs.update(e._solverKeywordParams())
            m = e._solverModifier()
            if m is not None:
                if modifier is None:
                    modifier = m
                else:
                    modifier = wrap(modifier, m)

        ind = allValues.addNode(_OptimCallNode(funcName, funcArgs, funcKwArgs, rs._endName))
        if modifier is not None:
            ind = allValues.addNode(_OptimFuncNode(modifier, (ind,), allValues))

    # Add the path index
    rs2IndLst.append(ind)


def OptimizeSelectors(sim, selectors):
    """
    Take a list of ResultSelectors and returns a list containing result selectors and
    optimized result selectors.
    """

    # Get all result paths
    rs2Ind = {}
    for sel in selectors:
        if isinstance(sel, nsaving._ResultList):
            for sel2 in sel._getAllTerminalChildren():
                if isinstance(sel2, nsaving._ResultPath):
                    rs2Ind[sel2] = []
        elif isinstance(sel, nsaving._ResultPath):
            rs2Ind[sel] = []

    # Extract calls that could be optimized
    allValues = _OptimValuesList(sim)

    for rs, lst in rs2Ind.items():
        for path in rs.simpath:
            _ExtractOptimValuesFromPath(path, rs, allValues, lst)

    allValues.setUpCallNodes()

    # Monkey patch result selectors
    for rs, indLst in rs2Ind.items():
        if any(allValues[i].hasBatchCall() for i in indLst):
            orp = _OptimizedResultPath([allValues[i] for i in indLst])
            rs._evaluate = orp._evaluate

    return selectors
