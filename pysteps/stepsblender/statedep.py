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

from . import blenderloader
from . import groups
from . import objects
from . import state


class StateDependentVesicles(objects.BlenderVesicles, state.StateDependentSeparateObjects):
    """Subclass of BlenderVesicles that allows to dynamically change material of vesicles

    See :py:class:`StateDependentSeparateObjects` for information about material set and update
    functions.

    The `state` parameter given to both functions is a dictionary that contains the following values: 
    +-----------------+---------------------------------+-----------------------------------------+
    | Name            | Example                         | Description                             |
    +=================+=================================+=========================================+
    | idx             | 45                              | Index of the specific vesicle           |
    | position        | [x, y, z]                       | 3D position of the vesicle in cartesian |
    |                 |                                 | coordinates                             |
    | name            | 'Ves1'                          | Name of the vesicle type in STEPS       |
    | surface_species | {'S1':[                         | Dictionary of species name to list      |
    |                 |   (r1, theta1, phi1),           | of spherical coordinate positions       |
    |                 |   (r2, theta2, phi2),           | relative to the vesicle                 |
    |                 | ...], ...}                      |                                         |
    | inside_species  | {'S1': 10, ...}                 | Dictionary of species name to count     |
    | link_species    | {'L1':                          | Dictionary of link species name to      |
    |                 |   {123: {                       | dictionary of specific link species     |
    |                 |     'position': [r, theta, phi] | to 3D position in relative spherical    |
    |                 |     'link': 456},               | coordinates and corresponding link      |
    |                 |   ...}, ...}                    | species                                 |
    +-----------------+---------------------------------+-----------------------------------------+

    Example::

        def redBoundVesicles(scene, depg, material, state):
            bsdf = material.node_tree.nodes['Principled BSDF']
            if 'L1' in state['link_species']:
                bsdf.inputs['Emission Color'].default_value = (1, 0, 0, 1) # Red
            else:
                bsdf.inputs['Emission Color'].default_value = (0, 1, 0, 1) # Green

        parameters = {'Ves1': StateDependentVesicles, 'updateMaterialFunc': redBoundVesicles}

        loader = HDF5BlenderLoader(parameters=parameters, ...)

    This example will display all vesicles that have link species `'L1'` on their surface as red
    and all other vesicles as green.
    """
    stateClasses = [groups.LinkSpeciesGroup, groups.VesicleGroup]

    def _signalUpdatedState(self, caller, scene, depg):
        if isinstance(caller, groups.VesicleGroup):
            for idx, pos, specIn, specSurf in caller.getVesStateInfo(self._name):
                dct = self._stateDict.setdefault(idx, self._getDefaultDict())
                dct['idx'] = idx
                dct['position'] = pos
                dct['name'] = self._name
                dct['inside_species'] = specIn
                dct['surface_species'] = specSurf
        elif isinstance(caller, groups.LinkSpeciesGroup):
            for lspecTpe, lspecDct in caller.state.items():
                for lsIdx, ((vesTpe, vesIdx), lsPos, lsLink) in lspecDct.items():
                    if vesTpe == self._name:
                        dct = self._stateDict.setdefault(vesIdx, self._getDefaultDict()).setdefault(
                            'link_species', {})
                        lsDct = {'position': lsPos, 'link': lsLink}
                        dct.setdefault(lspecTpe, {}).setdefault(lsIdx, lsDct)

    def _getDefaultDict(self):
        return {
            'idx': None,
            'position': None,
            'name': None,
            'surface_species': {},
            'inside_species': {},
            'link_species': {}
        }


class StateDependentRafts(objects.BlenderRafts, state.StateDependentSeparateObjects):
    """Subclass of BlenderRafts that allows to dynamically change material of rafts

    See :py:class:`StateDependentSeparateObjects` for information about material set and update
    functions.

    The `state` parameter given to both functions is a dictionary that contains the following values: 
    +-----------------+------------------------+--------------------------------------+
    | Name            | Example                | Description                          |
    +=================+========================+======================================+
    | idx             | 45                     | Index of the specific raft           |
    | position        | [x, y, z]              | 3D position of the raft in cartesian |
    |                 |                        | coordinates                          |
    | name            | 'Raft1'                | Name of the raft type in STEPS       |
    | surface_species | {'S1': 0, 'S2':5, ...} | Dictionary of species name to counts |
    +-----------------+------------------------+--------------------------------------+

    Example::

        def redS1Rafts(scene, depg, material, state):
            bsdf = material.node_tree.nodes['Principled BSDF']
            if 'S1' in state['surface_species']:
                bsdf.inputs['Emission Color'].default_value = (1, 0, 0, 1) # Red
            else:
                bsdf.inputs['Emission Color'].default_value = (0, 1, 0, 1) # Green

        parameters = {'Raft1': StateDependentRafts, 'updateMaterialFunc': redS1Rafts}

        loader = HDF5BlenderLoader(parameters=parameters, ...)

    This example will display all rafts that have species `'S1'` on their surface as red
    and all other rafts as green.
    """
    stateClasses = [groups.RaftGroup]

    def _signalUpdatedState(self, caller, scene, depg):
        for raftIdx, pos, counts in caller.getRaftStateInfo(self._name):
            dct = self._stateDict.setdefault(raftIdx, self._getDefaultDict())
            dct['idx'] = raftIdx
            dct['position'] = pos
            dct['name'] = self._name
            dct['surface_species'] = counts

    def _getDefaultDict(self):
        return {'idx': None, 'position': None, 'name': None, 'surface_species': {}}


class StateDependentMesh(objects.STEPSMeshObject, state.StateDependentObject):
    """Subclass of STEPSMeshObject that allows to dynamically change the material of a mesh

    See :py:class:`StateDependentObject` for information about material set and update functions.

    The `state` parameter given to both functions is a dictionary that contains the following values:

    +------------+---------+-----------------+
    | Name       | Example | Description     |
    +============+=========+=================+
    | time       | 0.23    | Simulation time |
    +------------+---------+-----------------+

    Example::

        def meshEmission(scene, depg, material, state):
            bsdf = material.node_tree.nodes['Principled BSDF']

            if 0.001 <= state['time'] < 0.002:
            if 'S1' in state['surface_species']:
                bsdf.inputs['Emission Strength'].default_value = 5
            else:
                bsdf.inputs['Emission Strength'].default_value = 0

        parameters = {'Membrane': StateDependentMesh, 'updateMaterialFunc': meshEmission}

        loader = HDF5BlenderLoader(parameters=parameters, ...)

    This example will increase the mesh emission strength between 1 and 2ms to visually convey e.g.
    a current injection.
    """
    stateClasses = [blenderloader.HDF5BlenderLoader]

    def _signalUpdatedState(self, caller, scene, depg):
        self._currState['time'] = caller.state['time']
