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

from . import objects


class State:

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._state = {}
        self._interpRatio = 0
        self._stateCache = {}
        self._users = []

    def _addSubscriber(self, user):
        self._users.append(user)

    def _getState(self, tind):
        """Return a state dictionary for a given STEPS time step index"""
        raise NotImplementedError()

    def _stateInterpolation(self, state1, state2, ratio):
        """Return a state dictionary that is an interpolation between state1 and state2
        ratio == 0.5 should return a state that is midway between state1 and state2"""
        raise NotImplementedError()

    def _updateDisplay(self, scene, depg, state):
        raise NotImplementedError()

    def updateDisplay(self, scene, depg):
        """Update the way the state is displayed
        This will be called after the Blender dependency graph is updated.
        Things like setting the position of particles should be done here.
        """
        self._updateDisplay(scene, depg, self._state)

    def _getCachedState(self, stind):
        if stind not in self._stateCache:
            self._stateCache[stind] = self._getState(stind)
        return self._stateCache[stind]

    def updateState(self, scene, depg, blenderLoader):
        """Trigger queries and update the state
        Display modifications can also be done here but only if they do not involve calls
        to an object's `evaluated_get` method.
        Display modifications that would involve dependency graph updates should be done here.
        """
        blcurr = scene.frame_current
        ststart = blenderLoader._Blender2STEPSTime(blcurr)
        blstart = blenderLoader._STEPS2BlenderTime(ststart)
        stend = ststart + 1
        blend = blenderLoader._STEPS2BlenderTime(stend)

        ratio = (blcurr - blstart) / (blend - blstart)
        self._interpRatio = blenderLoader.timeInterpFunc(ratio)
        state1 = self._getCachedState(ststart)
        if blcurr > blstart:
            state2 = self._getCachedState(stend)
            self._state = self._stateInterpolation(state1, state2, self._interpRatio)
        else:
            self._state = state1

        for user in self._users:
            user._signalUpdatedState(self, scene, depg)

    @property
    def state(self):
        return self._state


class StateUser:
    stateClasses = []

    def _handleSubscriptions(self, parent):
        for cls in self.stateClasses:
            for obj in parent._getAllChildren(cls):
                obj._addSubscriber(self)

    def _signalUpdatedState(self, caller, scene, depg):
        pass

    def updateBlenderObjects(self, scene, depg):
        pass


class StateDependentObject(objects.BlenderObject, StateUser):
    """Base class for objects whose material can change as a funciton of its state

    This class can be used in place of `BlenderObject` to customize the material of the object
    as a function of the state of the STEPS object.

    There are two possible ways to achieve this:
        - Use a small number of materials that are already declared in the Blender file;
        - Programatically customize the material of the object.

    The first way is done by setting the `setMaterialFunc` parameter to a function with the following
    signature::

        def myCustomMatSetFunc(scene, depg, state) -> str

    The first two parameters are passed from the `bpy.app.handlers.frame_change_pre` callback function
    (see https://docs.blender.org/api/current/bpy.app.handlers.html). The third parameter contains the
    state of the object (see table in :py:class:`StateDependentMesh`). The function should return a
    string that is the name of a material in the Blender file.

    The second way is done by setting the `updateMaterialFunc` parameter to a function with the
    following signature::

        def myCustomMatUpdateFunc(scene, depg, material, state) -> None

    The first two parameters are identical. The third parameter is the :py:class:`bpy.types.Material`
    object that is attributed to the object and should be modified by the function. The last parameter
    is the state of the object (see tables in :py:class:`StateDependentMesh`).
    """
    updateMaterialFunc = None
    setMaterialFunc = None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.updateMaterialFunc is not None and self.setMaterialFunc is not None:
            raise ValueError(f'Cannot supply both "updateMaterialFunc" and "setMaterialFunc" simultaneously')
        self._currState = {}

    def updateBlenderObjects(self, scene, depg):
        if self.updateMaterialFunc is not None:
            self.updateMaterialFunc(scene, depg, self.material.blenderObj, self._currState)
        elif self.setMaterialFunc is not None:
            matName = self.setMaterialFunc(scene, depg, self._currState)
            self.blenderObj.material_slots[0].material = bpy.data.materials[matName]
        self._currState = {}


class StateDependentSeparateObjects(objects.SeparateObjects, StateUser):
    """Base class for state-dependent objects like Vesicles or Rafts

    This class can be used in place of `SeparateObjects` to customize the material of each object
    as a function of the state of the STEPS object (its position, the species on its surface, etc.).

    There are two possible ways to achieve this:
        - Use a small number of materials that are already declared in the Blender file;
        - Programatically customize the material of each object.

    The first way is done by setting the `setMaterialFunc` parameter to a function with the following
    signature::

        def myCustomMatSetFunc(scene, depg, state) -> str

    The first two parameters are passed from the `bpy.app.handlers.frame_change_pre` callback function
    (see https://docs.blender.org/api/current/bpy.app.handlers.html). The third parameter contains the
    state of the object (see tables in :py:class:`StateDependentVesicle` or 
    :py:class:`StateDependentRaft`). The function should return a string that is the name of a
    material in the Blender file.

    The second way is done by setting the `updateMaterialFunc` parameter to a function with the
    following signature::

        def myCustomMatUpdateFunc(scene, depg, material, state) -> None

    The first two parameters are identical. The third parameter is the :py:class:`bpy.types.Material`
    object that is attributed to the object and should be modified by the function. The last parameter
    is the state of the object (see tables in :py:class:`StateDependentVesicle` or 
    :py:class:`StateDependentRaft`).
    """
    updateMaterialFunc = None
    setMaterialFunc = None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.updateMaterialFunc is not None and self.setMaterialFunc is not None:
            raise ValueError(f'Cannot supply both "updateMaterialFunc" and "setMaterialFunc" simultaneously')
        self._stateDict = {}
        if self.updateMaterialFunc is not None:
            self._allMaterials = {}
            mat = self.obj.material.blenderObj
            for idx, obj in self._objects.items():
                # Create as many materials as there are vesicles
                matName = f'{mat.name}_{idx}'
                if matName not in bpy.data.materials:
                    self._allMaterials[idx] = mat.copy()
                    self._allMaterials[idx].name = matName
                else:
                    self._allMaterials[idx] = bpy.data.materials[matName]
                obj.blenderObj.material_slots[0].link = 'OBJECT'
                obj.blenderObj.material_slots[0].material = self._allMaterials[idx]
        if self.setMaterialFunc is not None:
            for idx, obj in self._objects.items():
                obj.blenderObj.material_slots[0].link = 'OBJECT'

    def updateBlenderObjects(self, scene, depg):
        if self.updateMaterialFunc is not None:
            for idx, dct in self._stateDict.items():
                self.updateMaterialFunc(scene, depg, self._allMaterials[idx], dct)
        elif self.setMaterialFunc is not None:
            for idx, dct in self._stateDict.items():
                matName = self.setMaterialFunc(scene, depg, dct)
                self._objects[idx].blenderObj.material_slots[0].material = bpy.data.materials[matName]
        self._stateDict = {}
