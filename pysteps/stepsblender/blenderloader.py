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

try:
    import bpy
except ImportError:
    pass

import atexit
from multiprocessing.managers import BaseManager
import numpy as np
import os
import sys
from typing import Annotated

from . import groups
from . import objects
from . import state
from . import utils

from .utils import Orders, Loc


class QueueManager(BaseManager):
    pass


QueueManager.register('get_D2BQueue')
QueueManager.register('get_B2DQueue')


class HDF5BlenderLoader(objects.BlenderCollection, state.State):
    server: Annotated[str, 'Address of the data loading server'] = 'localhost'
    port: Annotated[int, 'Port to connect to on the data loading server'] = 57395
    authkey: Annotated[str, 'Authentication key to connect to the data loading server'] = 'STEPSBlender'

    dbInd: Annotated[
        str,
        'Unique run group identifier, only needed if several run groups exist in the file (see the documentation of steps.API_2.sim.Simulation.toDB)'] = None
    rInd: Annotated[int, 'Run index, defaults to -1 (last run in the run group)'] = -1

    include: Annotated[
        str,
        'Comma-separated list of regular expressions that controls which STEPS objects get loaded in Blender'] = ''
    exclude: Annotated[
        str,
        'Comma-separated list of regular expressions that controls which STEPS objects are prevented from being loaded in Blender'] = ''
    render: Annotated[
        bool,
        'Whether the script should launch Blender in background mode and automatically render each frame'] = False
    renderStart: Annotated[int,
                           'Start frame for rendering. If -1, will use the start frame from the file'] = -1
    renderEnd: Annotated[
        int, 'End frame for rendering (exclusive). Iif -1, will use the end frame from the file'] = -1
    renderStep: Annotated[int,
                          'Frame step for rendering, allows to skip some frames if it is higher than 1'] = 1
    outputPath: Annotated[
        str,
        'Path to which rendered frames should be saved when using the --render option, defaults to working directory'] = '.'

    scale: Annotated[float, 'Spatial scale factor'] = None

    timeScale: Annotated[int, 'Temporal scale factor'] = 1
    timeInterpFunc = utils.InterpolationFunction
    intersectAlgo: Annotated[
        str,
        'Intersection computation algorithm, defaults to "FAST" but uses "EXACT" when rendering to prevent visual bugs, "NONE" to turn it off'] = 'FAST'

    specScaleFactor: Annotated[float, 'Size of species compared to the mean tetrahedron size'] = 0.025
    background_color: Annotated[utils.colorType, 'Color of the background'] = (0.025, 0.025, 0.025, 1)
    cycles_shadows: Annotated[bool, 'Make global light cast shadows with the cycles render engine'] = False

    Meshes: groups.MeshGroup = None
    Species: groups.SpeciesGroup = None
    LinkSpecies: groups.LinkSpeciesGroup = None
    Vesicles: groups.VesicleGroup = None
    Rafts: groups.RaftGroup = None
    VesiclePaths: groups.VesiclePathGroup = None

    def __init__(self, parameters={}, **kwargs):
        parameters = utils.HierarchicalParameters(parameters)
        super().__init__(parameters=parameters, name='STEPSObjects', **kwargs)

        self._tetSize = None
        self._bbox = None
        self._simTimeSteps = None

        self._allElems = {}
        self._meshes = {}

        m = QueueManager(address=(self.server, self.port), authkey=self.authkey.encode('utf-8'))
        m.connect()

        self._queue_snd, self._queue_rcv = m.get_B2DQueue(), m.get_D2BQueue()

        atexit.register(self._exitServer)

        envInfo = {'numpy_version': np.__version__}
        self._queue_snd.put((self.dbInd, self.rInd, self.include, self.exclude, envInfo))
        if self._queue_rcv.get() == Orders.EXIT:
            print('Data loading server exited, exiting Blender.', file=sys.stderr)
            sys.exit()

        self._load()

    def _getState(self, tind):
        """Return a state dictionary for a given STEPS time step index"""
        return {'time': self._simTimeSteps[tind]}

    def _stateInterpolation(self, state1, state2, ratio):
        """Return a state dictionary that is an interpolation between state1 and state2
        ratio == 0.5 should return a state that is midway between state1 and state2"""
        t1, t2 = state1['time'], state2['time']
        return {'time': t1 + (t2 - t1) * ratio}

    def _exitServer(self):
        self._queue_snd.put((Orders.EXIT, ))

    def _getData(self, name):
        self._queue_snd.put((Orders.GET_DATA, name))
        return self._queue_rcv.get()

    def _STEPS2BlenderTime(self, t):
        return t * self.timeScale

    def _Blender2STEPSTime(self, t):
        return t // self.timeScale

    def _getAllElems(self):
        return self._allElems

    @property
    def isFromScratch(self):
        return bpy.data.filepath == ''

    def _setupGlobalProperties(self):
        if self.isFromScratch:
            # Display the materials correctly
            for area in bpy.context.screen.areas:
                if area.type == 'VIEW_3D':
                    space = area.spaces.active
                    if space.type == 'VIEW_3D':
                        space.shading.type = 'MATERIAL'

            bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[
                0].default_value = self.background_color

            # Change light to sun
            bpy.data.lights['Light'].type = 'SUN'
            bpy.data.lights['Light'].energy = 3
            bpy.data.lights['Light'].cycles.cast_shadow = self.cycles_shadows
            bpy.data.objects['Light'].rotation_euler = (0, 0, 0)

            # Set animation start and end frames
            bpy.context.scene.frame_end = self._STEPS2BlenderTime(self._getData('_maxTind'))

            # Move camera to fit all the mesh
            for _, mesh in self.Meshes._meshes.items():
                mesh.blenderObj.select_set(True)
            bpy.ops.view3d.camera_to_view_selected()
            for _, mesh in self.Meshes._meshes.items():
                mesh.blenderObj.select_set(False)

    def _getAndAddMesh(self):
        self._queue_snd.put((Orders.GET_MESH, ))
        self._allElems, triGrids, compSurfaces, avgTetSize, bbox = self._queue_rcv.get()

        self._tetSize = avgTetSize
        self._bbox = bbox

        # Auto set scale
        if self.scale is None:
            self.scale = 10 / np.max(self._bbox[1] - self._bbox[0])

        self.Meshes = self._getParam(
            'Meshes',
            groups.MeshGroup,
            name='Meshes',
            scale=self.scale,
            _meshVertices=self._allElems[Loc.VERT],
            _meshSurfaces=triGrids,
            _compSurfaces=compSurfaces,
        )

    def _getAndSetupModel(self):
        self._queue_snd.put((Orders.GET_MODEL, ))
        allSpecs, allLinkSpecs, allVes, allRafts, allVesPaths = self._queue_rcv.get()

        specRadius = self._tetSize * self.scale * self.specScaleFactor

        self.Species = self._getParam(
            'Species',
            groups.SpeciesGroup,
            name='Species',
            radius=specRadius,
            _specData=allSpecs,
        )

        self.LinkSpecies = self._getParam(
            'LinkSpecies',
            groups.LinkSpeciesGroup,
            name='LinkSpecies',
            radius=specRadius,
            bevel_depth=0.8 * specRadius,
            _linkSpecData=allLinkSpecs,
        )

        self.Vesicles = self._getParam(
            'Vesicles',
            groups.VesicleGroup,
            name='Vesicles',
            _vesData=allVes,
        )

        self.Rafts = self._getParam(
            'Rafts',
            groups.RaftGroup,
            name='Rafts',
            _raftData=allRafts,
        )

        self.VesiclePaths = self._getParam(
            'VesiclePaths',
            groups.VesiclePathGroup,
            name='VesiclePaths',
            _pathData=allVesPaths,
            _scale=self.scale,
        )

        self._simTimeSteps = self._getData('_timeSteps')

    def _getElemSpecCount(self, spec, loc, tind):
        self._queue_snd.put((Orders.GET_ELEM_SPEC_COUNT.value, spec, loc.value, tind))
        return self._queue_rcv.get()

    def _getVesLinkSpecRelPos(self, lspec, tind):
        self._queue_snd.put((Orders.GET_VES_LINKSPEC_REL_POS.value, lspec, tind))
        return self._queue_rcv.get()

    def _getVesPositions(self, ves, tind):
        self._queue_snd.put((Orders.GET_VES_POS.value, ves, tind))
        return self._queue_rcv.get()

    def _getVesEvents(self, ves, tind):
        self._queue_snd.put((Orders.GET_VES_EVENTS.value, ves, tind))
        return self._queue_rcv.get()

    def _getVesInSpecCounts(self, ves, tind):
        self._queue_snd.put((Orders.GET_VES_IN_SPEC_COUNT.value, ves, tind))
        return self._queue_rcv.get()

    def _getVesSurfSpecRelPos(self, ves, tind):
        self._queue_snd.put((Orders.GET_VES_SURF_SPEC_REL_POS.value, ves, tind))
        return self._queue_rcv.get()

    def _getRaftPositions(self, raft, tind):
        self._queue_snd.put((Orders.GET_RAFT_POS.value, raft, tind))
        return self._queue_rcv.get()

    def _getRaftCounts(self, raft, tind):
        self._queue_snd.put((Orders.GET_RAFT_COUNTS.value, raft, tind))
        return self._queue_rcv.get()

    def _getRaftEvents(self, raft, tind):
        self._queue_snd.put((Orders.GET_RAFT_EVENTS.value, raft, tind))
        return self._queue_rcv.get()

    def _getVertsV(self, vertInds, tind):
        self._queue_snd.put((Orders.GET_VERTS_V.value, vertInds, tind))
        return self._queue_rcv.get()

    def _preFrameChange(self, scene, depg):
        if 0 <= scene.frame_current <= scene.frame_end:
            tind = scene.frame_current
            print(f'=== New frame {tind} ===')
            self.updateState(scene, depg, self)
            # Update meshes
            self.Meshes.updateState(scene, depg, self)
            # Position vesicle first so that their position can be used in setting species
            self.Vesicles.updateState(scene, depg, self)
            # Position rafts
            self.Rafts.updateState(scene, depg, self)
            # Then position species
            self.Species.updateState(scene, depg, self)
            # Then position link species
            self.LinkSpecies.updateState(scene, depg, self)

            # Update all state-dependent objects
            for obj in self._getAllChildren(state.StateUser):
                obj.updateBlenderObjects(scene, depg)

            for layer in scene.view_layers:
                layer.update()

    def _postFrameChange(self, scene, depg):
        if 0 <= scene.frame_current <= scene.frame_end:
            self.Meshes.updateDisplay(scene, depg)
            self.Vesicles.updateDisplay(scene, depg)
            self.Rafts.updateDisplay(scene, depg)
            self.Species.updateDisplay(scene, depg)
            self.LinkSpecies.updateDisplay(scene, depg)

    def _setupCallbacks(self):
        # Important, we need to lock the interface during rendering or blender can SEGFAULT
        bpy.context.scene.render.use_lock_interface = True
        bpy.app.handlers.frame_change_pre.clear()
        bpy.app.handlers.frame_change_pre.append(self._preFrameChange)
        bpy.app.handlers.frame_change_post.clear()
        bpy.app.handlers.frame_change_post.append(self._postFrameChange)

        # Add callbacks for state dependent objects
        for obj in self._getAllChildren(state.StateUser):
            obj._handleSubscriptions(self)

    def _cleanBlendFile(self):
        if self.isFromScratch and 'Cube' in bpy.data.meshes:
            bpy.data.meshes.remove(bpy.data.meshes['Cube'])

    def _load(self):
        self._cleanBlendFile()

        self._getAndAddMesh()

        self._getAndSetupModel()

        bpy.context.scene.frame_current = 0
        layer = bpy.context.view_layer
        layer.update()

        self._setupCallbacks()

        self._setupGlobalProperties()

        if self.render:
            fstart = self.renderStart if self.renderStart >= 0 else bpy.context.scene.frame_start
            fend = self.renderEnd if self.renderEnd >= 0 else (bpy.context.scene.frame_end + 1)
            nbdigits = len(str(fend))
            for f in range(fstart, fend, self.renderStep):
                bpy.context.scene.render.filepath = os.path.join(self.outputPath, f'{f:0{nbdigits}}.png')
                bpy.context.scene.frame_current = f
                bpy.ops.render.render(write_still=True)
