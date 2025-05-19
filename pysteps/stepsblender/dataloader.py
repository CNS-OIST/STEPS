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

import argparse
from multiprocessing.managers import BaseManager
import numpy as np
from queue import Queue
import re
import warnings

from .utils import Orders, Loc, Event, zipNone, cartesian2Spherical


class _HDF5BlenderDataLoader:

    def __init__(self, HDFPath, queue_snd, queue_rcv):
        self._hdfPath = HDFPath
        self._queue_snd = queue_snd
        self._queue_rcv = queue_rcv

    def isIncluded(self, name):
        return any(reg.match(name)
                   for reg in self._include) or not any(reg.match(name) for reg in self._exclude)

    def serve(self):
        import steps.API_2.saving as saving
        import steps.API_2.utils as utils

        self._allElems = {}
        self._tri2Tets = {}
        self._spec2Rs = {}
        self._ves2SpecsRs = {}
        self._ves2Rs = {}
        self._maxTind = None

        self._currVesPos = {}
        self._currRaftPos = {}
        self._ves2Rad = {}
        self._vesEvents = {}
        self._raftEvents = {}

        self._dbInd, self._rInd, include, exclude, envInfo = self._queue_rcv.get()
        # Check that Blender's numpy version is compatible with the locally installed one
        servNpVer = utils.Versioned._parseVersion(np.__version__)
        blendNpVer = utils.Versioned._parseVersion(envInfo.get('numpy_version', '0.0'))
        if servNpVer >= (2, 0) and blendNpVer < (1, 26) or blendNpVer >= (2,0) and servNpVer < (1, 26):
            self._queue_snd.put(Orders.EXIT)
            raise Exception(
                f'The numpy version installed for the data loading server {servNpVer} is not '
                f'compatible with the numpy version installed in Blender {blendNpVer}. Try '
                f'{"down" if servNpVer >= (2,0) else "up"}grading the numpy version installed '
                f'for the data server. See '
                f'https://numpy.org/doc/stable//numpy_2_0_migration_guide.html#note-about-pickled-files '
                f'for details on which versions are compatible.'
            )
        self._queue_snd.put(Orders.OK)

        self._include = [re.compile(reg) for reg in include.split(',')] if include != '' else []
        self._exclude = [re.compile(reg) for reg in exclude.split(',')] if exclude != '' else []

        utils.SetVerbosity(3)
        with saving.HDF5Handler(self._hdfPath, internalKwArgs=dict(maxFullLoadSize=1024**3)) as hdf:
            if self._dbInd is None:
                group, *rem = hdf
                if len(rem) > 0:
                    raise ValueError('Several run groups exist in the HDF5 file, but no run group identifier '
                                     'was provided')
                self._dbInd = group.name
            else:
                group = hdf[self._dbInd]
            results = group.results
            meshGroup = group._group[saving.XDMFHandler._MESH_GROUP_NAME]

            while True:
                order, *args = self._queue_rcv.get()
                ret = None
                print('Received order', Orders(order))
                if order is Orders.GET_MESH:
                    ret = self._getMesh(hdf, group, meshGroup, *args)
                elif order is Orders.GET_MODEL:
                    ret = self._getModelObjects(hdf, group, results, *args)
                elif order == Orders.GET_ELEM_SPEC_COUNT:
                    ret = self._getElemSpecCounts(*args)
                elif order == Orders.GET_VES_IN_SPEC_COUNT:
                    ret = self._getVesInSpecCounts(*args)
                elif order == Orders.GET_VES_SURF_SPEC_REL_POS:
                    ret = self._getVesSurfSpecRelPos(*args)
                elif order == Orders.GET_VES_LINKSPEC_REL_POS:
                    ret = self._getVesLinkSpecRelPos(*args)
                elif order == Orders.GET_VES_EVENTS:
                    ret = self._getVesEvents(*args)
                elif order == Orders.GET_VES_POS:
                    ret = self._getVesPos(*args)
                elif order == Orders.GET_RAFT_POS:
                    ret = self._getRaftPos(*args)
                elif order == Orders.GET_RAFT_COUNTS:
                    ret = self._getRaftCounts(*args)
                elif order == Orders.GET_RAFT_EVENTS:
                    ret = self._getRaftEvents(*args)
                elif order == Orders.GET_VERTS_V:
                    ret = self._getVertsV(*args)
                elif order is Orders.GET_DATA:
                    name, = args
                    ret = getattr(self, name)
                elif order is Orders.EXIT:
                    break
                else:
                    raise NotImplementedError(f'No such order: {order}')
                self._queue_snd.put(ret)
        utils.SetVerbosity(1)

    def _getMesh(self, hdf, hdfgroup, meshGroup):
        self._allElems = {
            Loc.VERT: {},  # {elem_ind: (x, y, z)}
            Loc.TRI: {},  # {elem_ind: (v1, v2, v3)}
            Loc.TET: {},  # {elem_ind: (v1, v2, v3, v4)}
        }
        triGrids = {}
        compSurfaces = {}
        allTetSizes = []
        bbox = [None, None]
        for groupName in meshGroup:
            group = meshGroup[groupName]
            if 'topology' in group:
                vertInds = np.array(group['vertInds'])
                elemInds = np.array(group['elemInds'])
                positions = np.array(group['XYZ'])
                topology = np.array(group['topology'])
                loc_id = group.attrs['loc_id']

                groupSurface = set()

                for vInd, pos in zip(vertInds, positions):
                    if vInd not in self._allElems[Loc.VERT]:
                        self._allElems[Loc.VERT][vInd] = pos

                tpeMap = {6: Loc.TET, 4: Loc.TRI, 1: Loc.VERT}
                triGrid = []
                for eInd, (tpe, *inds) in zip(elemInds, topology):
                    loc = tpeMap[tpe]
                    self._allElems[loc][eInd] = tuple(vertInds[inds])
                    if loc is Loc.TRI:
                        triGrid.append(tuple(vertInds[inds]))
                    if loc is Loc.TET:
                        tris = set()
                        for t in range(4):
                            tris.add(frozenset(tuple(vertInds[inds[0:t]]) + tuple(vertInds[inds[t + 1:]])))
                        groupSurface ^= tris
                        # Compute tetsize
                        vertices = np.array([self._allElems[Loc.VERT][vertInds[ind]] for ind in inds])
                        vmin = np.min(vertices, axis=0)
                        vmax = np.max(vertices, axis=0)
                        allTetSizes.append(vmax - vmin)
                        bbox[0] = np.min([bbox[0], vmin], axis=0) if bbox[0] is not None else vmin
                        bbox[1] = np.max([bbox[1], vmax], axis=0) if bbox[1] is not None else vmax
                if len(triGrid) > 0:
                    triGrids.setdefault(loc_id, [])
                    triGrids[loc_id] += triGrid
                if len(groupSurface) > 0:
                    compSurfaces[loc_id] = groupSurface

        # Add Endocytic zones
        for zoneName, zoneDct in hdfgroup.staticData.get('EndocyticZones', {}).items():
            triGrids[zoneName] = [self._allElems[Loc.TRI].get(tidx, tuple()) for tidx in zoneDct['tris']]

        # Compute tri 2 tets mapping
        self._tri2Tets = {}
        for triIdx, tri in self._allElems[Loc.TRI].items():
            self._tri2Tets.setdefault(frozenset(tri), [])
        for tetIdx, tet in self._allElems[Loc.TET].items():
            for t in range(4):
                face = frozenset(tuple(tet[0:t]) + tuple(tet[t + 1:]))
                tetLst = self._tri2Tets.get(face, None)
                if tetLst is not None:
                    tetLst.append(tetIdx)

        meshSurface = set()
        for _, surf in compSurfaces.items():
            meshSurface ^= surf

        # Remove triangles that are already part of a patch
        meshSurface -= set(frozenset(vinds) for _, vinds in self._allElems[Loc.TRI].items())

        if len(meshSurface) > 0:
            surfName = 'mesh_surface'
            while surfName in triGrids:
                surfName += '_'
            triGrids[surfName] = [tuple(vinds) for vinds in meshSurface]

        # Average tet size
        avgTetSize = np.mean(allTetSizes)

        return (self._allElems, triGrids, compSurfaces, avgTetSize, bbox)

    def _getModelObjects(self, hdf, hdfgroup, results):
        """Send all different species saved in results"""
        import steps.API_2.model as model

        self._timeSteps = None
        self._spec2Rs = {}
        self._linkspec2RS = {}
        self._ves2SpecsRs = {}
        self._spec2RaftRs = {}
        self._ves2Rs = {}
        self._raft2Rs = {}
        self._vesInds2Tpe = {}
        self._raftInds2Tpe = {}
        self._vesEventsRs = {}
        self._raftEventsRs = {}
        self._vertsPot = {}
        ps2Inds = {}
        ls2Inds = {}
        ves2IndsLocs = {}
        raft2IndsLocs = {}
        raft2Specs = {}
        str2TpeMap = {'Tri': Loc.TRI, 'Tet': Loc.TET}
        allVes = set()
        allRaft = set()
        allSpecs = set()
        maxTime = 0
        for sel in results:
            obj_type = sel.metaData.get('obj_type', None)
            obj_id = sel.metaData.get('obj_id', None)
            loc_type = sel.metaData.get('loc_type', None)
            loc_id = sel.metaData.get('loc_id', None)
            ves_type = sel.metaData.get('vesicle_type', None)
            ves_loc = sel.metaData.get('vesicle_loc', None)
            raft_type = sel.metaData.get('raft_type', None)
            pointspec_type = sel.metaData.get('species_type', None)
            linkspec_type = sel.metaData.get('linkspecies_type', None)
            props = sel.metaData.get('property', None)

            # Check that all timesteps are the same and that they all cover the same period
            if self._timeSteps is None:
                self._timeSteps = sel.time[self._rInd]
            elif len(self._timeSteps) != len(sel.time[self._rInd]):
                warnings.warn(
                    'Some result selectors have different number of timepoints.')
            elif (self._timeSteps != sel.time[self._rInd]).any():
                warnings.warn(f'Some result selectors were saved at different timepoints.')

            maxTind = len(sel.time[self._rInd]) - 1
            self._maxTind = min(self._maxTind, maxTind) if self._maxTind is not None else maxTind
            if obj_id is None and loc_type is not None and loc_id is not None and props is not None:
                # Membrane potential
                for i, (prop, etpe, eid) in enumerate(zipNone(props, loc_type, loc_id)):
                    if prop == 'V' and etpe == 'Vert':
                        self._vertsPot[eid] = (sel, i)

            # Tet and tri species
            if all(x is not None for x in [obj_type, obj_id, loc_type, loc_id]):
                for i, (tpe, oid, etpe, eid) in enumerate(zipNone(obj_type, obj_id, loc_type, loc_id)):
                    # TODO Inprovement: also load link species from there when it is saved with e.g. sim.comp.VESICLES().L1.Pos
                    if tpe == model.Species._elemStr:
                        if etpe in str2TpeMap:
                            allSpecs.add(oid)

                            # TODO Optimization: make this structure more efficient?
                            loc = str2TpeMap[etpe]
                            struct = self._spec2Rs.setdefault(loc, {})
                            elems, colInds = struct.setdefault(oid, {}).setdefault(sel, ([], []))
                            elems.append(eid)
                            colInds.append(i)
            if ves_type is not None and props is not None and linkspec_type is None:
                # Vesicle positions
                for i, (vesTpe, prop, oid, lid) in enumerate(zipNone(ves_type, props, obj_id, loc_id)):
                    if oid is None and vesTpe is not None and prop == 'Pos':
                        allVes.add(vesTpe)
                        self._ves2Rs.setdefault(vesTpe, {}).setdefault(sel, []).append(i)
                        for vesDct in sel.data[self._rInd, :, i]:
                            if isinstance(vesDct, dict):
                                for idx in vesDct.keys():
                                    ves2IndsLocs.setdefault(vesTpe, {})[idx] = lid
                                    self._vesInds2Tpe[idx] = vesTpe
                if ves_loc is not None:
                    # Vesicle specs
                    for i, (vesTpe, vesLoc, prop, objTpe, oid, psTpe) in enumerate(
                            zipNone(ves_type, ves_loc, props, obj_type, obj_id, pointspec_type)):
                        if psTpe is not None:
                            # Individual point spec positions
                            allSpecs.add(psTpe)
                            self._ves2SpecsRs.setdefault(vesTpe, {}).setdefault(
                                (prop, vesLoc, 'dct'), {}).setdefault(psTpe, {}).setdefault(sel, []).append(i)
                        elif objTpe == 'Spec' and vesLoc is not None and vesTpe is not None:
                            # Bulk point spec positions
                            allSpecs.add(oid)
                            self._ves2SpecsRs.setdefault(vesTpe, {}).setdefault(
                                (prop, vesLoc, 'lst'), {}).setdefault(oid, {}).setdefault(sel, []).append(i)
            if raft_type is not None and props is not None:
                for i, (raftTpe, prop, objTpe, oid,
                        lid) in enumerate(zipNone(raft_type, props, obj_type, obj_id, loc_id)):
                    if oid is None and raftTpe is not None and prop == 'Pos':
                        # Rafts positions
                        allRaft.add(raftTpe)
                        self._raft2Rs.setdefault(raftTpe, {}).setdefault(lid, {}).setdefault(sel,
                                                                                             []).append(i)
                        for raftDct in sel.data[self._rInd, :, i]:
                            if isinstance(raftDct, dict):
                                for idx in raftDct.keys():
                                    raft2IndsLocs.setdefault(raftTpe, {})[idx] = lid
                                    self._raftInds2Tpe[idx] = raftTpe
                    if oid is not None and objTpe == 'Spec':
                        # Raft specs
                        if prop in ['Pos', 'Count']:
                            allSpecs.add(oid)
                            raft2Specs.setdefault(raftTpe, set()).add(oid)
                            self._spec2RaftRs.setdefault(oid, {}).setdefault(sel, {}).setdefault(
                                raftTpe, {}).setdefault(lid, {}).setdefault(prop, []).append(i)
            if linkspec_type is not None and props is not None:
                for i, (lsTpe, prop) in enumerate(zipNone(linkspec_type, props)):
                    if prop in ['Pos', 'LinkedTo']:
                        self._linkspec2RS.setdefault(lsTpe, {}).setdefault(prop, {}).setdefault(sel,
                                                                                                []).append(i)
                        for veslsDct in sel.data[self._rInd, :, i]:
                            if isinstance(veslsDct, dict):
                                for _, lsDct in veslsDct.items():
                                    if isinstance(lsDct, dict):
                                        ls2Inds.setdefault(lsTpe, set())
                                        ls2Inds[lsTpe] |= set(lsDct.keys())
            # Events
            if obj_type is not None and obj_id is not None and props is not None:
                for i, (tpe, oid, prop) in enumerate(zipNone(obj_type, obj_id, props)):
                    if prop == 'Events':
                        if tpe == 'Exocytosis':
                            self._vesEventsRs.setdefault(Event.EXOCYTOSIS, []).append((sel, i))
                        elif tpe == 'Endocytosis':
                            self._vesEventsRs.setdefault(Event.ENDOCYTOSIS, []).append((sel, i))
                        elif tpe == 'RaftEndocytosis':
                            self._raftEventsRs.setdefault(Event.RAFT_ENDOCYTOSIS, []).append((sel, i))

        VesInfos = []
        RaftInfos = []
        for vesTpe in filter(self.isIncluded, allVes):
            indsLocs = ves2IndsLocs.get(vesTpe, {})
            allVesSpecs = set()
            for _, specDct in self._ves2SpecsRs.get(vesTpe, {}).items():
                allVesSpecs |= set(filter(self.isIncluded, specDct.keys()))
            if len(indsLocs) > 0:
                # Fetch ves diameter from model data
                diam = hdfgroup.staticData.get('Vesicles', {}).get(vesTpe, {}).get('Diameter', None)
                if diam is not None:
                    self._ves2Rad[vesTpe] = diam / 2
                    VesInfos.append((vesTpe, diam / 2, indsLocs, allVesSpecs))

        for raftTpe in filter(self.isIncluded, allRaft):
            indsLocs = raft2IndsLocs.get(raftTpe, {})
            allRaftSpecs = set(filter(self.isIncluded, raft2Specs.get(raftTpe, set())))
            # Only add a raft if it appears at least once
            if len(indsLocs) > 0:
                # Fetch raft diameter from model data
                diam = hdfgroup.staticData.get('Rafts', {}).get(raftTpe, {}).get('Diameter', None)
                if diam is not None:
                    RaftInfos.append((raftTpe, diam / 2, indsLocs, allRaftSpecs))

        PathInfos = hdfgroup.staticData.get('VesiclePaths', {})

        SpecInfos = sorted(filter(self.isIncluded, allSpecs))
        for loc, spec2RS in self._spec2Rs.items():
            for specTpe in spec2RS.keys():
                if self.isIncluded(specTpe) and specTpe not in SpecInfos:
                    SpecInfos.append(specTpe)

        LinkSpecInfos = [(lsTpe, ls2Inds[lsTpe]) for lsTpe in self._linkspec2RS.keys()
                         if self.isIncluded(lsTpe)]

        return (SpecInfos, LinkSpecInfos, VesInfos, RaftInfos, PathInfos)

    def _getElemSpecCounts(self, spec, loc, tind):
        loc = Loc(loc)

        allElems = []
        allCounts = []
        if loc in self._spec2Rs and spec in self._spec2Rs[loc]:
            for sel, (elems, inds) in self._spec2Rs[loc][spec].items():
                for idx, cnt in zip(elems, sel.data[self._rInd, tind, inds]):
                    allElems.append(idx)
                    allCounts.append(cnt)

        return (allElems, allCounts)

    def _getVesInSpecCounts(self, ves, tind):
        allCounts = {}

        specDct = self._ves2SpecsRs.get(ves, {}).get(('Count', 'in', 'lst'), {})
        for spec, rsDct in specDct.items():
            for sel, colInds in rsDct.items():
                for data in sel.data[self._rInd, tind, colInds]:
                    for vesIdx, cnt in data.items():
                        allCounts.setdefault(vesIdx, {})[spec] = cnt

        return allCounts

    def _getVesSurfSpecRelPos(self, ves, tind):
        allPositions = {}

        propDct = self._ves2SpecsRs.get(ves, {})
        for (prop, vesLoc, dataTpe), specDct in propDct.items():
            if vesLoc == 'surf' and prop.startswith('Pos'):
                for spec, rsDct in specDct.items():
                    for sel, colInds in rsDct.items():
                        for data in sel.data[self._rInd, tind, colInds]:
                            for vesIdx, posData in data.items():
                                if len(posData) > 0:
                                    if dataTpe == 'dct':
                                        positions = np.array(list(posData.values()))
                                        idxs = list(posData.keys())
                                    elif dataTpe == 'lst':
                                        positions = np.array(posData)
                                        idxs = list(range(len(posData)))
                                    else:
                                        raise NotImplementedError()
                                    # Convert positions to full spherical coordinates
                                    if prop == 'Pos':
                                        positions = cartesian2Spherical(positions -
                                                                        self._currVesPos[(ves, tind)][vesIdx])
                                    elif prop == 'PosSpherical':
                                        positions = np.hstack((self._ves2Rad[ves] * np.ones(
                                            (len(positions), 1)), positions))
                                    allPositions.setdefault(
                                        vesIdx, {})[spec] = {idx: pos
                                                             for idx, pos in zip(idxs, positions)}
                                else:
                                    allPositions.setdefault(vesIdx, {})[spec] = {}

        return allPositions

    def _getVesLinkSpecRelPos(self, lstpe, tind):
        LinkSpecPos = {}
        if lstpe in self._linkspec2RS:
            for sel, colInds in self._linkspec2RS[lstpe].get('Pos', {}).items():
                for dct in sel.data[self._rInd, tind, colInds]:
                    for vidx, lsDct in dct.items():
                        for lsIdx, val in lsDct.items():
                            ves = (self._vesInds2Tpe[vidx], vidx)
                            vesPos = self._currVesPos[(ves[0], tind)][vidx]
                            LinkSpecPos.setdefault(
                                lsIdx, [ves, None, None])[1] = cartesian2Spherical(np.array(val) - vesPos)
            for sel, colInds in self._linkspec2RS[lstpe].get('LinkedTo', {}).items():
                for dct in sel.data[self._rInd, tind, colInds]:
                    for vidx, lsDct in dct.items():
                        for lsIdx, val in lsDct.items():
                            ves = (self._vesInds2Tpe[vidx], vidx)
                            LinkSpecPos.setdefault(lsIdx, [ves, None, None])[2] = val

        return LinkSpecPos

    def _getVesPos(self, ves, tind):
        if (ves, tind) not in self._currVesPos:
            self._currVesPos[(ves, tind)] = {}
            if ves in self._ves2Rs:
                for sel, colInds in self._ves2Rs[ves].items():
                    for dct in sel.data[self._rInd, tind, colInds]:
                        for idx, pos in dct.items():
                            self._currVesPos[(ves, tind)][idx] = np.array(pos)

        return self._currVesPos[(ves, tind)]

    def _getRaftPos(self, raftTpe, tind):
        raftTimeKey = (raftTpe, tind)
        self._currRaftPos[raftTimeKey] = {}
        for loc, selDct in self._raft2Rs.get(raftTpe, {}).items():
            for sel, colInds in selDct.items():
                for dct in sel.data[self._rInd, tind, colInds]:
                    for idx, pos in dct.items():
                        self._currRaftPos[raftTimeKey][idx] = np.array(pos)

        return self._currRaftPos[raftTimeKey]

    def _getRaftCounts(self, raftTpe, tind):
        raftCounts = {}
        # Species on rafts
        for spec, specDct in self._spec2RaftRs.items():
            for sel, raftDct in specDct.items():
                for loc, propDct in raftDct.get(raftTpe, {}).items():
                    colInds = propDct.get('Count', None)
                    if colInds is not None:
                        for data in sel.data[self._rInd, tind, colInds]:
                            for raftIdx, specCount in data.items():
                                raftCounts.setdefault(raftIdx, {})[spec] = specCount

        return raftCounts

    def _processEvent(self, sel, tind, evTpe, event, eventDct):
        t1, t2 = sel.time[self._rInd, tind - 1:tind + 1]
        if evTpe == Event.ENDOCYTOSIS:
            time, tidx, vidx = event
            if vidx in self._vesInds2Tpe:
                vesTpe = self._vesInds2Tpe[vidx]
                triPos = np.mean([self._allElems[Loc.VERT][vidx] for vidx in self._allElems[Loc.TRI][tidx]],
                                 axis=0)
                evTime = (time - t1) / (t2 - t1)
                eventDct.setdefault(vesTpe, {}).setdefault(vidx, []).append(
                    (evTpe.value, evTime, tidx, triPos))
            else:
                warnings.warn(
                    f'Vesicle with unique identifier {vidx} was created and destroyed between two saving time steps. It will not be present in the visualization.'
                )
        elif evTpe == Event.EXOCYTOSIS:
            time, vidx, tidx, ridx = event
            if vidx in self._vesInds2Tpe:
                vesTpe = self._vesInds2Tpe[vidx]
                triPos = np.mean([self._allElems[Loc.VERT][vidx] for vidx in self._allElems[Loc.TRI][tidx]],
                                 axis=0)
                evTime = (time - t1) / (t2 - t1)
                eventDct.setdefault(vesTpe, {}).setdefault(vidx, []).append(
                    (evTpe.value, evTime, tidx, triPos, ridx))
            else:
                warnings.warn(
                    f'Vesicle with unique identifier {vidx} was created and destroyed between two saving time steps. It will not be present in the visualization.'
                )
        elif evTpe == Event.RAFT_ENDOCYTOSIS:
            time, ridx, tidx, vidx = event
            if ridx in self._raftInds2Tpe:
                raftTpe = self._raftInds2Tpe[ridx]
                triPos = np.mean([self._allElems[Loc.VERT][vidx] for vidx in self._allElems[Loc.TRI][tidx]],
                                 axis=0)
                evTime = (time - t1) / (t2 - t1)
                eventDct.setdefault(raftTpe, {}).setdefault(ridx, []).append(
                    (evTpe.value, evTime, tidx, triPos, vidx))
            else:
                warnings.warn(
                    f'Raft with unique identifier {ridx} was created and destroyed between two saving time steps. It will not be present in the visualization.'
                )
        else:
            raise NotImplementedError(evTpe, event)

    def _loadAllVesEvents(self, tind):
        if tind not in self._vesEvents:
            events = self._vesEvents.setdefault(tind, {})
            if tind > 0:
                for evTpe, selLst in self._vesEventsRs.items():
                    for sel, colInd in selLst:
                        for event in sel.data[self._rInd, tind, colInd]:
                            self._processEvent(sel, tind, evTpe, event, events)

    def _loadAllRaftEvents(self, tind):
        if tind not in self._raftEvents:
            events = self._raftEvents.setdefault(tind, {})
            if tind > 0:
                for evTpe, selLst in self._raftEventsRs.items():
                    for sel, colInd in selLst:
                        for event in sel.data[self._rInd, tind, colInd]:
                            self._processEvent(sel, tind, evTpe, event, events)

    def _getRelatedEvents(self, otherEvents, tpe, idx2Tpe, compatibleEvents):
        for otherTpe, evDct in otherEvents.items():
            for idx, evLst in evDct.items():
                for evTpe, evTime, *evInfo in evLst:
                    if Event(evTpe) in [Event.RAFT_ENDOCYTOSIS, Event.EXOCYTOSIS]:
                        tidx, triPos, oidx = evInfo
                        if evTpe in compatibleEvents and idx2Tpe.get(oidx, None) == tpe:
                            yield oidx, evTpe, evTime, tidx, triPos, idx

    def _getVesEvents(self, ves, tind):
        self._loadAllVesEvents(tind)
        self._loadAllRaftEvents(tind)
        events = self._vesEvents[tind].get(ves, {})
        for vidx, *ev in self._getRelatedEvents(self._raftEvents[tind], ves, self._vesInds2Tpe,
                                                [Event.RAFT_ENDOCYTOSIS.value]):
            events.setdefault(vidx, []).append(ev)
        return events

    def _getRaftEvents(self, raft, tind):
        self._loadAllVesEvents(tind)
        self._loadAllRaftEvents(tind)
        events = self._raftEvents[tind].get(raft, {})
        for ridx, *ev in self._getRelatedEvents(self._vesEvents[tind], raft, self._raftInds2Tpe,
                                                [Event.EXOCYTOSIS.value]):
            events.setdefault(ridx, []).append(ev)
        return events

    def _getVertsV(self, vertInds, tind):
        verts = []
        for v in vertInds:
            sel, colInd = self._vertsPot.get(v, (None, None))
            if sel is not None:
                #TODO Optimization: better data access
                verts.append(float(sel.data[self._rInd, tind, colInd]))
            else:
                return None
        return np.array(verts)


D2BQueue, B2DQueue = Queue(), Queue()


def D2BQueue_getter():
    return D2BQueue


def B2DQueue_getter():
    return B2DQueue


class QueueManager(BaseManager):
    pass


QueueManager.register('get_D2BQueue', callable=D2BQueue_getter)
QueueManager.register('get_B2DQueue', callable=B2DQueue_getter)

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Start a data loading server for Blender visualization.')

    # HDF5 File
    parser.add_argument(
        'HDFPath',
        action='store',
        help='The path prefix to a STEPS HDF5 file (see the documentation of steps.API_2.saving.HDF5Handler)')

    # Server options
    parser.add_argument('--port',
                        type=int,
                        action='store',
                        help='Port to connect to on the data loading server',
                        default=57395)
    parser.add_argument('--authkey',
                        action='store',
                        help='Authentication key to connect to the data loading server',
                        default='STEPSBlender')

    args = parser.parse_args()

    m = QueueManager(address=('', args.port), authkey=args.authkey.encode('utf-8'))
    m.start()

    dataLoader = _HDF5BlenderDataLoader(args.HDFPath, m.get_D2BQueue(), m.get_B2DQueue())
    print(f'Data loading server listening on port {args.port}')
    dataLoader.serve()
    m.shutdown()
