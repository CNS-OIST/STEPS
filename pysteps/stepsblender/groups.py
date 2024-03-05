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

import numpy as np

from . import objects
from . import state
from . import utils

from .utils import Event, Loc, spherical2Cartesian


# class MeshGroup(objects.BlenderCollection):
class MeshGroup(objects.BlenderCollection, state.State):
    _meshVertices = None
    _meshSurfaces = None
    _compSurfaces = None

    # Only used for parameter listing
    sampleMesh: objects.STEPSMeshObject = None

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)
        self._meshes = {}
        for name, tris in self._meshSurfaces.items():
            self._meshes[name] = self._getParam(
                name,
                objects.STEPSMeshObject,
                name=name,
                color=utils.GetColor(v=0.5),
                _scale=self.parent.scale,
                _verts=self._meshVertices,
                _tris=tris,
            )
        self._compMeshes = {}
        # Compute connected compartments (not separated by a patch)
        # We consider that if, among the triangles shared by two compartments, at least one
        # is part of a patch, then the compartments are fully separated.
        patchTris = set()
        for name, tris in self._meshSurfaces.items():
            patchTris |= set(frozenset(tri) for tri in tris)

        self._comp2Group = {}
        finalGroups = {}
        compSurfaces = list(self._compSurfaces.items())
        while len(compSurfaces) > 0:
            name, surf = compSurfaces[0]
            merged = False
            for i, (name2, surf2) in enumerate(compSurfaces[1:]):
                inter = surf & surf2
                if len(inter) > 0 and not any(tri in patchTris for tri in inter):
                    newName = f'{name}_{name2}'
                    self._comp2Group[name] = newName
                    self._comp2Group[name2] = newName
                    compSurfaces = [(newName, surf ^ surf2)] + compSurfaces[1:i + 1] + compSurfaces[i + 2:]
                    merged = True
                    break
            if not merged:
                finalGroups[name] = surf
                compSurfaces = compSurfaces[1:]

        for name, tris in finalGroups.items():
            # These meshes are not user parameterizable
            self._compMeshes[name] = self._getParam(
                name,
                objects.STEPSMeshObject,
                name=name,
                material=None,
                solidifySurface=False,
                _scale=self.parent.scale,
                _verts=self._meshVertices,
                _tris=tris,
            )
            obj = self._compMeshes[name].blenderObj
            obj.hide_viewport = True
            obj.hide_render = True

    def getSurfMesh(self, name):
        return self._meshes[name]

    def getCompMesh(self, name):
        while name in self._comp2Group:
            name = self._comp2Group[name]
        return self._compMeshes[name]

    def _getState(self, tind):
        state = {}
        for name, meshobj in self._meshes.items():
            potVals = self.parent._getVertsV(meshobj.mesh._vertInds, tind)
            if potVals is not None:
                state.setdefault(name, {})['V'] = potVals
        return state

    def _stateInterpolation(self, state1, state2, ratio):
        state = {}
        for name, dct in state1.items():
            for prop, vals1 in dct.items():
                vals2 = state2[name][prop]
                state.setdefault(name, {})[prop] = vals1 + (vals2 - vals1) * ratio
        return state

    def _updateDisplay(self, scene, depg, state):
        for name, dct in state.items():
            for prop, vals in dct.items():
                meshobj = self._meshes.get(name, None)
                if meshobj is not None and all(v is not None for v in vals):
                    meshobj.mesh.updateVertProp(scene, depg, prop, vals)


class VesiclePathGroup(objects.BlenderCollection):
    _pathData = None

    # Only used for parameter listing
    sampleVesiclePath: objects.STEPSVesiclePath = None

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)
        self._paths = {}
        for pathName, path in self._pathData.items():
            self._paths[pathName] = self._getParam(
                pathName,
                objects.STEPSVesiclePath,
                name=pathName,
                _data=path,
            )


class SpeciesGroup(objects.BlenderCollection, state.State):
    _specData = None

    # Only used for parameter listing
    sampleSpecies: objects.BlenderSpecies = None

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)
        self._species = {}
        for spec in self._specData:
            self._species[spec] = self._getParam(
                spec,
                objects.BlenderSpecies,
                name=spec,
                color=utils.GetColor(v=0.5),
            )

    def _getState(self, tind):
        state = {}
        for specName, spec in self._species.items():
            state.setdefault(specName, {})
            for loc in [Loc.TRI, Loc.TET]:
                state[specName][loc] = self.parent._getElemSpecCount(specName, loc, tind)
        return state

    def updateState(self, scene, depg, blenderLoader):
        super().updateState(scene, depg, blenderLoader)
        for specName, spec in self._species.items():
            cnt = 0
            for loc, specData in self.state[specName].items():
                if loc in [Loc.TET, Loc.TRI]:
                    idx, cnts = specData
                cnt += sum(cnts)
            spec._updateCount(cnt)

    @staticmethod
    def _specNbInterpolation(c1, c2, ratio):
        v = c1 + (c2 - c1) * ratio
        if v - int(v) == 0.5:
            return int(v) if c1 > c2 else int(v) + 1
        else:
            return round(v)

    def _stateInterpolation(self, state1, state2, ratio):
        state = {}
        for specName, locDct1 in state1.items():
            for loc, vals1 in locDct1.items():
                state.setdefault(specName, {})

                vals2 = state2[specName][loc]
                if loc in [Loc.TRI, Loc.TET]:
                    inds1, cnts1 = vals1
                    inds2, cnts2 = vals2
                    cnts = []
                    for idx, c1, c2 in zip(inds1, cnts1, cnts2):
                        if ratio < 0.5:
                            cnts.append(c1)
                        else:
                            cnts.append(self._specNbInterpolation(c1, c2, 2 * (ratio - 0.5)))
                    state[specName][loc] = (inds1, cnts)
        return state

    def _updateDisplay(self, scene, depg, state):
        for specName, spec in self._species.items():
            positions = []
            for tri, cnt in zip(*state[specName][Loc.TRI]):
                p0, p1, p2 = [self.parent._allElems[Loc.VERT][v] for v in self.parent._allElems[Loc.TRI][tri]]
                for i in range(int(cnt)):
                    s, t = np.random.random(2)
                    u = s**0.5
                    v = u * t
                    positions.append((1 - u) * p0 + (u - v) * p1 + v * p2)

            for tet, cnt in zip(*state[specName][Loc.TET]):
                p0, p1, p2, p3 = [
                    self.parent._allElems[Loc.VERT][v] for v in self.parent._allElems[Loc.TET][tet]
                ]
                for i in range(int(cnt)):
                    s, t, u = np.random.random(3)
                    if s + t > 1:
                        s = 1 - s
                        t = 1 - t
                    if t + u > 1:
                        tmp = u
                        u = 1 - s - t
                        t = 1 - tmp
                    elif s + t + u > 1:
                        tmp = u
                        u = s + t + u - 1
                        s = 1 - t - tmp
                    a = 1 - s - t - u
                    positions.append(a * p0 + s * p1 + t * p2 + u * p3)

            spec._setPositions(scene, depg, self.parent.scale * np.array(positions))


class LinkSpeciesGroup(objects.BlenderCollection, state.State):
    _linkSpecData = None

    # Only used for parameter listing
    sampleLinkSpecies: objects.BlenderLinkSpecies = None

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)
        self._linkSpecs = {}
        for i, (linkspec, indexes) in enumerate(self._linkSpecData):
            self._linkSpecs[linkspec] = self._getParam(
                linkspec,
                objects.BlenderLinkSpecies,
                name=linkspec,
                _indexes=indexes,
                color=utils.GetColor(v=0.5),
            )

    def _getState(self, tind):
        state = {}
        idx2lsTpe = {}
        for lspecName, lspec in self._linkSpecs.items():
            state[lspecName] = self.parent._getVesLinkSpecRelPos(lspecName, tind)
            idx2lsTpe = {**idx2lsTpe, **{idx: lspecName for idx in state[lspecName]}}
        for lspecName, poslink in state.items():
            for idx, (_, _, link) in poslink.items():
                if link is not None:
                    # Remove the link part if the link species it is linked to was not saved or loaded
                    poslink[idx][2] = (idx2lsTpe[link], link) if link in idx2lsTpe else None
        return state

    def updateState(self, scene, depg, blenderLoader):
        super().updateState(scene, depg, blenderLoader)
        for lspecName, lspec in self._linkSpecs.items():
            cnt = len(self.state[lspecName])
            lspec._updateCount(cnt)

    def _stateInterpolation(self, state1, state2, ratio):
        state = {}
        for lspecName, poslink1 in state1.items():
            poslink2 = state2[lspecName]
            poslink = state.setdefault(lspecName, {})
            for idx1, (ves1, pos1, link1) in poslink1.items():
                if idx1 in poslink2:
                    ves2, pos2, link2 = poslink2[idx1]
                    pos = utils.sphericalInterpolation(pos1, pos2, ratio)
                    if link1 != link2:
                        link = link1 if ratio <= 0.5 else link2
                    else:
                        link = link1
                    poslink[idx1] = (ves1, pos, link)
                elif ratio <= 0.5:
                    # Disappearing linkspec
                    poslink[idx1] = (ves1, pos1, link1)
            # Appearing linkspecs
            if ratio >= 0.5:
                for idx in set(poslink2.keys()) - set(poslink1.keys()):
                    poslink[idx] = poslink2[idx]
        return state

    def _updateDisplay(self, scene, depg, state):
        for lspecName, poslink in state.items():
            positions = []
            linkPositions = {}
            for idx, ((vesName1, vesIdx1), spos, link) in poslink.items():
                vpos1 = self.parent.Vesicles.getVesPos(vesName1, vesIdx1)
                p1 = self.parent.scale * (vpos1 + spherical2Cartesian(spos))
                positions.append(p1)
                if link is not None:
                    linktpe, linkidx = link
                    (vesName2, vesIdx2), spos2, _ = state[linktpe][linkidx]
                    vpos2 = self.parent.Vesicles.getVesPos(vesName2, vesIdx2)
                    p2 = self.parent.scale * (vpos2 + spherical2Cartesian(spos2))
                    linkPositions[idx] = (p1, p2)

            self._linkSpecs[lspecName]._setPositions(scene, depg, np.array(positions))
            self._linkSpecs[lspecName]._setLinkPositions(scene, depg, linkPositions)


class VesicleGroup(objects.BlenderCollection, state.State):
    _vesData = None

    # Only used for parameter listing
    sampleVesicle: objects.BlenderVesicles = None

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)
        self._vesicles = {}
        self._ves2Rad = {}
        for i, (ves, radius, indexesLocations, specs) in enumerate(self._vesData):
            meshLocations = {
                idx: self.parent.Meshes.getCompMesh(loc)
                for idx, loc in indexesLocations.items()
            }
            species = [self.parent.Species._species[specName] for specName in specs]
            self._ves2Rad[ves] = radius
            self._vesicles[ves] = self._getParam(
                ves,
                objects.BlenderVesicles,
                name=ves,
                color=utils.GetColor(v=0.1),
                _indexes=set(indexesLocations),
                _specs=species,
                _locations=meshLocations,
                _radius=self.parent.scale * radius,
            )

    def _getState(self, tind):
        state = {}
        for vesName, ves in self._vesicles.items():
            # TODO Optimization: group all calls in one?
            pos = self.parent._getVesPositions(vesName, tind)
            countIn = self.parent._getVesInSpecCounts(vesName, tind)
            posSurf = self.parent._getVesSurfSpecRelPos(vesName, tind)
            events = self.parent._getVesEvents(vesName, tind)
            vesDct, evDct = state.setdefault(vesName, ({}, events))
            for idx, p in pos.items():
                ci = countIn.get(idx, {})
                ps = posSurf.get(idx, {})
                ev = events.get(idx, [])
                vesDct[idx] = (p, ci, ps)
        return state

    def _getTriNormal(self, triIdx):
        elems = self.parent._getAllElems()
        verts = elems[Loc.TRI].get(triIdx, None)
        if verts is not None:
            p1, p2, p3 = [elems[Loc.VERT][v] for v in verts]
            norm = np.cross(p1 - p2, p3 - p2)
            return norm / np.linalg.norm(norm)
        return None

    def _animateExocytosis(self, vesName, pos1, triPos, triIdx, ratio):
        pos = None
        if ratio <= 0.5:
            pos = pos1 + 2 * ratio * (triPos - pos1)
        elif ratio <= 0.75:
            triNorm = self._getTriNormal(triIdx)
            if np.dot(triNorm, triPos - pos1) < 0:
                triNorm = -triNorm
            dest = triPos + triNorm * self._ves2Rad[vesName]
            pos = triPos + 4 * (ratio - 0.5) * (dest - triPos)
        return pos

    def _animateEndocytosis(self, vesName, pos2, triPos, triIdx, ratio):
        pos = None
        if 0.25 <= ratio <= 0.5:
            triNorm = self._getTriNormal(triIdx)
            if np.dot(triNorm, triPos - pos2) < 0:
                triNorm = -triNorm
            orig = triPos + triNorm * self._ves2Rad[vesName]
            pos = orig + 4 * (ratio - 0.25) * (triPos - orig)
        elif ratio > 0.5:
            pos = triPos + (2 * ratio - 1) * (pos2 - triPos)
        return pos

    def _stateInterpolation(self, state1, state2, ratio):
        state = {}
        for vesName, (st1, ev1) in state1.items():
            st2, ev2 = state2[vesName]
            st, _ = state.setdefault(vesName, ({}, ev2))
            for idx, (pos1, cnts1, surfPos1) in st1.items():
                st2Tpl = st2.get(idx, None)
                if st2Tpl is not None:
                    pos2, cnts2, surfPos2 = st2Tpl

                    pos = pos1 + (pos2 - pos1) * ratio
                    cnts = {
                        spec: SpeciesGroup._specNbInterpolation(cnt1, cnts2[spec], ratio)
                        for spec, cnt1 in cnts1.items()
                    }

                    surfPos = {}
                    for spec, sdata1 in surfPos1.items():
                        immob = self._vesicles[vesName].isSpecImmobile(spec)
                        sdata2 = surfPos2[spec]
                        spDct = surfPos.setdefault(spec, {})
                        for sidx1, spos1 in sdata1.items():
                            spos2 = sdata2.get(sidx1, None)
                            if spos2 is not None:
                                spDct[sidx1] = spos1 if immob else utils.sphericalInterpolation(
                                    spos1, spos2, ratio)
                            elif ratio <= 0.5 or immob:
                                # Disappearing point specs
                                spDct[sidx1] = spos1
                        if ratio > 0.5 and not immob:
                            # Appearing point specs
                            for sidx2 in set(sdata2.keys()) - set(sdata1.keys()):
                                spDct[sidx2] = sdata2[sidx2]

                    st[idx] = (pos, cnts, surfPos)
                elif idx in ev2:
                    # Exocytosis
                    for evType, evRatio, triIdx, triPos, ridx in ev2[idx]:
                        if evType == Event.EXOCYTOSIS.value:
                            pos = self._animateExocytosis(vesName, pos1, triPos, triIdx, ratio)
                            if pos is not None:
                                st[idx] = (pos, cnts1, surfPos1)
                        else:
                            raise NotImplementedError(evType)
                elif ratio <= 0.5:
                    # disapearing vesicles
                    st[idx] = (pos1, cnts1, surfPos1)
            # appearing vesicles
            for idx in set(st2.keys()) - set(st1.keys()):
                pos2, cnts2, surfPos2 = st2[idx]
                if idx in ev2:
                    # Endocytosis and Raft endocytosis
                    for evType, evRatio, *evInfo in ev2[idx]:
                        if evType == Event.RAFT_ENDOCYTOSIS.value:
                            triIdx, triPos, ridx = evInfo
                            pos = self._animateEndocytosis(vesName, pos2, triPos, triIdx, ratio)
                            if pos is not None:
                                st[idx] = (pos, cnts2, surfPos2)
                        elif evType == Event.ENDOCYTOSIS.value:
                            triIdx, triPos = evInfo
                            pos = self._animateEndocytosis(vesName, pos2, triPos, triIdx, ratio)
                            if pos is not None:
                                st[idx] = (pos, cnts2, surfPos2)
                        else:
                            raise NotImplementedError(evType)
                elif ratio >= 0.5:
                    st[idx] = st2[idx]
        return state

    def updateState(self, scene, depg, blenderLoader):
        super().updateState(scene, depg, blenderLoader)
        for vesName, ves in self._vesicles.items():
            specCounts = {}
            scaledVesPos = {}
            vesDct, events = self.state[vesName]
            for idx, (pos, counts, surfPos) in vesDct.items():
                scaledVesPos[idx] = self.parent.scale * pos
                specCounts.setdefault(Loc.VES_IN, {})[idx] = counts
                specCounts.setdefault(Loc.VES_SURF,
                                      {})[idx] = {spec: len(spos)
                                                  for spec, spos in surfPos.items()}
            # Update vesicle position here otherwise it would trigger a depgraph update because of the
            # boolean modifiers
            ves._setPositions(scene, depg, scaledVesPos)
            ves._updateSpecCounts(specCounts)
            ves._setEventStatus(scene, depg, events)

    def _updateDisplay(self, scene, depg, state):
        s = self.parent.scale
        for vesName, (st, events) in state.items():
            scaledSpecPos = {}
            for idx, (pos, _, surfPos) in st.items():
                scaledSpecPos.setdefault(Loc.VES_SURF, {})[idx] = {
                    spec: s * (spherical2Cartesian(np.array(list(spos.values()))) + pos)
                    for spec, spos in surfPos.items() if len(spos) > 0
                }
            self._vesicles[vesName]._setSpecPositions(scene, depg, scaledSpecPos)

    def getVesStateInfo(self, vesTpe):
        """Return the state information that will be used for state-dependent vesicles"""
        vesInfos, vesEvents = self.state.get(vesTpe, ({}, {}))
        for idx, (pos, specIn, specSurf) in vesInfos.items():
            yield idx, pos, specIn, specSurf

    def getVesPos(self, vesTpe, vesIdx):
        return self._state[vesTpe][0][vesIdx][0]

    def getVesRad(self, vesTpe):
        return self._ves2Rad[vesTpe]

    def isVesUnderEvent(self, vesTpe, vesIdx):
        return len(self._state[vesTpe][1].get(vesIdx, [])) > 0


class RaftGroup(objects.BlenderCollection, state.State):
    _raftData = None

    # Only used for parameter listing
    sampleRaft: objects.BlenderRafts = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)
        self._rafts = {}
        for i, (raftTpe, radius, indexesLocations, specs) in enumerate(self._raftData):
            meshLocations = {
                idx: self.parent.Meshes.getSurfMesh(loc)
                for idx, loc in indexesLocations.items()
            }
            species = [self.parent.Species._species[specName] for specName in specs]
            self._rafts[raftTpe] = self._getParam(
                raftTpe,
                objects.BlenderRafts,
                name=raftTpe,
                color=utils.GetColor(v=0.7),
                _indexes=set(indexesLocations),
                _specs=species,
                _locations=meshLocations,
                _radius=radius * self.parent.scale,
            )

    def _getState(self, tind):
        state = {}
        for raftTpe, raft in self._rafts.items():
            pos = self.parent._getRaftPositions(raftTpe, tind)
            counts = self.parent._getRaftCounts(raftTpe, tind)
            events = self.parent._getRaftEvents(raftTpe, tind)
            raftDct, evDct = state.setdefault(raftTpe, ({}, events))
            for idx, p in pos.items():
                cnt = counts.get(idx, {})
                raftDct[idx] = (p, cnt)
        return state

    def _animateRaftEndocytosis(self, pos1, triPos, ratio):
        pos = None
        if ratio <= 0.25:
            pos = pos1 + (triPos - pos1) * ratio * 4
        elif ratio <= 0.5:
            pos = triPos
        return pos

    def _animateExocytosis(self, pos2, triPos, ratio):
        pos = None
        if 0.5 <= ratio < 0.75:
            pos = triPos
        elif ratio >= 0.75:
            pos = triPos + (pos2 - triPos) * (ratio - 0.75) * 4
        return pos

    def _stateInterpolation(self, state1, state2, ratio):
        state = {}
        for raftTpe, (vals1, evs1) in state1.items():
            vals2, evs2 = state2[raftTpe]
            vals, evs = state.setdefault(raftTpe, ({}, evs2))
            for idx, v1 in vals1.items():
                p1, cnt1 = v1
                v2 = vals2.get(idx, None)
                if v2 is not None:
                    p2, cnt2 = v2
                    p = p1 + (p2 - p1) * ratio
                    cnt = {
                        spec: SpeciesGroup._specNbInterpolation(c1, cnt2.get(spec, 0), ratio)
                        for spec, c1 in cnt1.items()
                    }
                    vals[idx] = (p, cnt)
                elif idx in evs2:
                    for evType, evRatio, triIdx, triPos, vidx in evs2[idx]:
                        if evType == Event.RAFT_ENDOCYTOSIS.value:
                            pos = self._animateRaftEndocytosis(p1, triPos, ratio)
                            if pos is not None:
                                vals[idx] = (pos, cnt1)
                        else:
                            raise NotImplementedError(evType)
                elif ratio <= 0.5:
                    # disapearing rafts
                    vals[idx] = v1
            # appearing rafts
            for idx in set(vals2.keys()) - set(vals1.keys()):
                p2, cnt2 = vals2[idx]
                if idx in evs2:
                    for evType, evRatio, triIdx, triPos, vidx in evs2[idx]:
                        if evType == Event.EXOCYTOSIS.value:
                            pos = self._animateExocytosis(p2, triPos, ratio)
                            if pos is not None:
                                vals[idx] = (pos, cnt2)
                elif ratio >= 0.5:
                    vals[idx] = vals2[idx]
        return state

    def updateState(self, scene, depg, blenderLoader):
        super().updateState(scene, depg, blenderLoader)
        for raftName, raft in self._rafts.items():
            specCounts = {}
            for idx, (pos, counts) in self.state[raftName][0].items():
                specCounts.setdefault(Loc.RAFT_IN, {})[idx] = counts
            raft._updateSpecCounts(specCounts)

    def _updateDisplay(self, scene, depg, state):
        for raftTpe, (vals, evs) in state.items():
            pos = {}
            counts = {}
            for idx, (p, cnt) in vals.items():
                pos[idx] = self.parent.scale * p
                counts[idx] = cnt
            self._rafts[raftTpe]._setPositions(scene, depg, pos)

    def getRaftStateInfo(self, raftTpe):
        """Return the state information that will be used for state-dependent rafts"""
        raftInfos, raftEvents = self.state.get(raftTpe, ({}, {}))
        for idx, (pos, specIn) in raftInfos.items():
            yield idx, pos, specIn
