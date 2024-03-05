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
    import bmesh
    import mathutils
except ImportError:
    pass

import colorsys
import numpy as np
import re
from typing import Annotated

from . import utils

from .utils import Loc

####################################################################################################

_FAR_LOCATION = (0, 0, 1e3)

####################################################################################################


class BlenderWrapper(utils.HierarchicalParamReader):
    """Base class for all STEPS Blender classes

    This class automatically checks if the Blender item is already in the Blender file
    If not, it creates it with the appropriate Blender data type.

    The ._name attribute corresponds to the name of the object in the blender file
    It can be supplied directly as a keyword parameter. If not, a unique name will
    be created by combining the name of the parent object with the name of the attribute
    that holds the current object in the parent.

    With this hierarchy of objects (a V above means this object was directly named):
       V          V
    Species --> Spec1 ---> obj ------> mesh
                      \            \-> material
                       \-> emitter --> mesh
                                   \-> material
    We get the corresponding names:
    Species --> Spec1 ---> Spec1_obj ------> Spec1_obj_mesh
                      \              \-----> Spec1_obj_material
                       \-> Spec1_emitter --> Spec1_emitter_mesh
                                         \-> Spec1_emitter_material
    """

    def __init__(self, blendContName, name=None, parameters={}, **kwargs):
        self._blendContName = blendContName
        self._name = name
        self._hidden = False

        if parameters is not None:
            super().__init__(parameters=parameters, **kwargs)

        if self._name is None:
            self._name = self._getName()

        if parameters is not None:
            # Create or link Blender object
            if self._name not in self.blenderDataCont:
                self.blenderObj = self.CreateBlenderObj(self._name)
                self.setUp(self.blenderObj, fromScratch=True)
            else:
                self.blenderObj = self.blenderDataCont[self._name]
                self.setUp(self.blenderObj, fromScratch=False)

    def _getName(self):
        if self._name is None:
            if self.parent is None or not isinstance(self.parent, BlenderWrapper):
                return ''
            else:
                return self.parent._getName() + '_' + self.nameInParent
        else:
            return self._name

    def CreateBlenderObj(self, name):
        """Create and return the actual blender object"""
        return self.blenderDataCont.new(name)

    def setUp(self, obj, fromScratch):
        pass

    def setHidden(self, hidden=True, obj=None):
        if obj is None:
            obj = self.blenderObj
        if self._hidden != hidden:
            obj.hide_viewport = hidden
            obj.hide_render = hidden
            self._hidden = hidden

    @property
    def blenderDataCont(self):
        return getattr(bpy.data, self._blendContName)

    # We do not hold pointers to blender objects because they might become invalidated by some
    # Blender operations. Instead, we querry the object every time.
    @property
    def blenderObj(self):
        if self._name is None:
            raise AttributeError()
        try:
            return getattr(bpy.data, self._blendContName)[self._name]
        except KeyError:
            raise AttributeError()

    @blenderObj.setter
    def blenderObj(self, obj):
        self._name = obj.name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, v):
        self._name = v
        self.blenderObj.name = v

    @staticmethod
    def BlenderCopy(contName, obj, name):
        newObj = obj.copy()
        newObj.name = name
        return BlenderWrapper(contName, name=name, parameters=None), newObj

    def blenderCopy(self, name, obj=None):
        if obj is None:
            obj = self.blenderObj
        return BlenderWrapper.BlenderCopy(self._blendContName, obj, name)

    # TMP
    def __repr__(self):
        return f'{self.__class__.__name__}({self._name if hasattr(self, "name") else id(self)})'


####################################################################################################


class BlenderCollection(BlenderWrapper):
    """Base class for Blender collections used to group STEPS objects

    Subclasses should implement setUp(self, coll, fromScratch)
    """

    def __init__(self, **kwargs):
        self._blenderChildren = []
        super().__init__('collections', **kwargs)

    def setUp(self, coll, fromScratch):
        """Sets up the given blender collection object

            :param coll: the Blender collection object
            :type coll: :py:class:`bpy.types.Collection`
            :param fromScratch: Whether the object was created from scratch
            :type fromScratch: bool
        """
        if fromScratch:
            if self.parent is not None:
                self.parent.addChild(self)
            else:
                bpy.context.scene.collection.children.link(coll)

        for obj in self._blenderChildren:
            self.addChild(obj)

    def addChild(self, obj):
        if hasattr(self, 'blenderObj'):
            if isinstance(obj, BlenderObject) and obj.blenderObj.name not in self.blenderObj.objects:
                self.blenderObj.objects.link(obj.blenderObj)
            elif isinstance(obj, BlenderCollection) and obj.blenderObj.name not in self.blenderObj.children:
                self.blenderObj.children.link(obj.blenderObj)
        else:
            self._blenderChildren.append(obj)

    def _getCollPath(self):
        if self.parent is not None:
            yield from self.parent._getCollPath()
        yield self._name

    def SetActiveLayerCollection(self):
        layer_coll = bpy.context.view_layer.layer_collection
        for collName in self._getCollPath():
            layer_coll = layer_coll.children[collName]
        bpy.context.view_layer.active_layer_collection = layer_coll


####################################################################################################


class ShaderNodeGroup(utils.HierarchicalParamReader):
    """Base class for custom shader node group"""

    def getNodeOutput(self, nodes, node_tree, inpt):
        """This method takes its input socket as argument and should return the output socket"""
        raise NotImplementedError()


class ShaderNodeMathRescale(ShaderNodeGroup):
    InMin: Annotated[float, 'Minimum value of input'] = 0
    InMax: Annotated[float, 'Maximum value of input'] = 1
    OutMin: Annotated[float, 'Minimum value of output'] = 0
    OutMax: Annotated[float, 'Maximum value of output'] = 1

    def getNodeOutput(self, nodes, node_tree, inpt):
        multAdd = nodes.new('ShaderNodeMath')
        multAdd.operation = 'MULTIPLY_ADD'
        dIn = self.InMax - self.InMin
        dOut = self.OutMax - self.OutMin
        multAdd.inputs[1].default_value = dOut / dIn
        multAdd.inputs[2].default_value = self.OutMin - self.InMin * dOut / dIn

        node_tree.links.new(multAdd.inputs[0], inpt)

        return multAdd.outputs[0]


class ShaderNodeColorMap(ShaderNodeGroup):
    rescaler: ShaderNodeMathRescale = ShaderNodeMathRescale
    colormap: Annotated[
        str,
        'Name of the matplotlib colormap, see https://matplotlib.org/stable/gallery/color/colormap_reference.html'] = 'viridis'
    npoints: Annotated[int, 'Number of points added to the Blender colorramp'] = 16

    def getNodeOutput(self, nodes, node_tree, inpt):
        try:
            import matplotlib
        except ImportError:
            raise ImportError(
                'The matplotlib python package needs to be installed to use ShaderNodeColorMap.')

        rescaleOutput = self.rescaler.getNodeOutput(nodes, node_tree, inpt)
        ramp = nodes.new('ShaderNodeValToRGB')
        cmap = matplotlib.colormaps[self.colormap]
        for i in range(0, self.npoints):
            pos = i / (self.npoints - 1)
            if i == len(ramp.color_ramp.elements):
                elem = ramp.color_ramp.elements.new(pos)
            else:
                elem = ramp.color_ramp.elements[i]
                elem.position = pos
            elem.color = cmap(pos)

        node_tree.links.new(ramp.inputs[0], rescaleOutput)

        return ramp.outputs['Color']


####################################################################################################


class BlenderMaterial(BlenderWrapper):
    """Base class for Blender materials used to visualize STEPS objects

    Subclasses should implement setUp(self, mat, fromScratch)
    """

    def __init__(self, **kwargs):
        super().__init__('materials', **kwargs)


BlenderMaterial.setUp.__doc__ = \
"""Sets up the given blender material object

    :param mat: the Blender material object
    :type mat: :py:class:`bpy.types.Material`
    :param fromScratch: Whether the object was created from scratch
    :type fromScratch: bool
"""


class DefaultBSDFMaterial(BlenderMaterial):
    color: utils.colorType = (0, 0, 0, 0)
    alpha: utils.alphaType = 1
    emission: utils.emissionType = 0

    def setUp(self, mat, fromScratch):
        if fromScratch:
            mat.use_nodes = True

            bsdf = mat.node_tree.nodes.get("Principled BSDF")
            bsdf.inputs['Base Color'].default_value = self.color
            em_socket_name = 'Emission Color' if 'Emission Color' in bsdf.inputs else 'Emission'
            bsdf.inputs[em_socket_name].default_value = self.color
            bsdf.inputs['Alpha'].default_value = self.alpha
            bsdf.inputs['Emission Strength'].default_value = self.emission
            if self.alpha < 1:
                mat.blend_method = 'BLEND'
            mat.shadow_method = 'NONE'


class MeshMaterial(DefaultBSDFMaterial):
    alpha = 0.5

    def setUp(self, mat, fromScratch):
        if fromScratch:
            mat.show_transparent_back = True
            mat.use_backface_culling = True
            mat.use_nodes = True
            node_tree = mat.node_tree
            nodes = node_tree.nodes

            bsdf = nodes['Principled BSDF']
            bsdf.inputs['Base Color'].default_value = self.color
            em_socket_name = 'Emission Color' if 'Emission Color' in bsdf.inputs else 'Emission'
            bsdf.inputs[em_socket_name].default_value = self.color
            bsdf.inputs['Alpha'].default_value = self.alpha
            bsdf.inputs['Emission Strength'].default_value = self.emission
            if self.alpha < 1:
                mat.blend_method = 'BLEND'
            mat.shadow_method = 'NONE'

            transp = nodes.new('ShaderNodeBsdfTransparent')
            mix = nodes.new('ShaderNodeMixShader')
            facing = nodes.new('ShaderNodeLayerWeight')
            mat_out = nodes['Material Output']

            node_tree.links.clear()

            node_tree.links.new(mat_out.inputs['Surface'], mix.outputs['Shader'])
            node_tree.links.new(mix.inputs['Fac'], facing.outputs['Facing'])
            node_tree.links.new(mix.inputs[1], transp.outputs['BSDF'])
            node_tree.links.new(mix.inputs[2], bsdf.outputs['BSDF'])


class StateDepMeshMaterial(MeshMaterial):
    attrName: Annotated[str, 'Attribute that should be used by the material'] = 'V'
    bsdfAttr: Annotated[str, 'BSDF attribute that should be modified by the state'] = 'Base Color'
    attrPipeline: ShaderNodeGroup = ShaderNodeColorMap.using(
        rescaler=ShaderNodeMathRescale.using(InMin=-0.065, InMax=0))

    def setUp(self, mat, fromScratch):
        super().setUp(mat, fromScratch)
        if fromScratch:
            node_tree = mat.node_tree
            nodes = node_tree.nodes
            bsdf = nodes['Principled BSDF']

            attr = nodes.new('ShaderNodeAttribute')
            attr.attribute_name = self.attrName

            output = self.attrPipeline.getNodeOutput(nodes, node_tree, attr.outputs['Fac'])

            node_tree.links.new(bsdf.inputs[self.bsdfAttr], output)


class SpeciesMaterial(DefaultBSDFMaterial):
    emission = 1


class VesiclePathMaterial(DefaultBSDFMaterial):
    color = (0.5, 0.5, 0.5, 1)
    emission = 1


class VesicleMaterial(BlenderMaterial):
    color: utils.colorType = None
    outline_color: utils.colorType = None
    alpha: utils.alphaType = 0.3

    fresnel_IOR: Annotated[float, 'Fresnel Index Of Refraction'] = 1.05
    fresnel_multiplier: Annotated[float, 'Strength of the outer rim'] = 2

    def setUp(self, mat, fromScratch):
        if fromScratch:
            if self.outline_color is None:
                hsvColor = colorsys.rgb_to_hsv(*self.color[:3])
                self.outline_color = colorsys.hsv_to_rgb(hsvColor[0], 0.7, 0.7) + self.color[3:]

            mat.use_nodes = True
            node_tree = mat.node_tree
            nodes = node_tree.nodes

            bsdf = nodes['Principled BSDF']
            bsdf.inputs['Base Color'].default_value = self.color
            em_socket_name = 'Emission Color' if 'Emission Color' in bsdf.inputs else 'Emission'
            bsdf.inputs[em_socket_name].default_value = self.outline_color
            bsdf.inputs['Alpha'].default_value = self.alpha
            if self.alpha < 1:
                mat.blend_method = 'BLEND'
            mat.shadow_method = 'NONE'

            fresnel = nodes.new('ShaderNodeFresnel')
            fresnel.inputs['IOR'].default_value = self.fresnel_IOR

            mult = nodes.new('ShaderNodeMath')
            mult.operation = 'MULTIPLY'
            mult.inputs[0].default_value = self.fresnel_multiplier

            node_tree.links.new(mult.inputs[1], fresnel.outputs['Fac'])
            node_tree.links.new(bsdf.inputs['Emission Strength'], mult.outputs['Value'])


class RaftMaterial(BlenderMaterial):
    color: utils.colorType = None
    outline_color: utils.colorType = None
    alpha: utils.alphaType = 0.3
    emission: utils.emissionType = 2

    outline_frac: Annotated[float, 'Fraction of the raft radius that defines the raft border'] = 0.9

    _radius: float = None

    def setUp(self, mat, fromScratch):
        if fromScratch:
            if self.outline_color is None:
                hsvColor = colorsys.rgb_to_hsv(*self.color[:3])
                self.outline_color = colorsys.hsv_to_rgb(hsvColor[0], 0.7, 0.7) + self.color[3:]

            mat.use_nodes = True
            node_tree = mat.node_tree
            nodes = node_tree.nodes

            bsdf = nodes['Principled BSDF']
            bsdf.inputs['Base Color'].default_value = self.color
            em_socket_name = 'Emission Color' if 'Emission Color' in bsdf.inputs else 'Emission'
            bsdf.inputs[em_socket_name].default_value = self.outline_color
            bsdf.inputs['Alpha'].default_value = self.alpha
            if self.alpha < 1:
                mat.blend_method = 'BLEND'
            mat.shadow_method = 'NONE'

            geom = nodes.new('ShaderNodeNewGeometry')
            obj = nodes.new('ShaderNodeObjectInfo')

            dist = nodes.new('ShaderNodeVectorMath')
            dist.operation = 'DISTANCE'
            node_tree.links.new(dist.inputs[0], obj.outputs['Location'])
            node_tree.links.new(dist.inputs[1], geom.outputs['Position'])

            gt = nodes.new('ShaderNodeMath')
            gt.operation = 'GREATER_THAN'
            gt.inputs[1].default_value = self._radius * self.outline_frac
            node_tree.links.new(gt.inputs[0], dist.outputs['Value'])

            mult = nodes.new('ShaderNodeMath')
            mult.operation = 'MULTIPLY'
            mult.inputs[0].default_value = self.emission
            node_tree.links.new(mult.inputs[1], gt.outputs['Value'])

            node_tree.links.new(bsdf.inputs['Emission Strength'], mult.outputs['Value'])


####################################################################################################


class BlenderMesh(BlenderWrapper):
    """Base class for Blender mesh that are used to visualize STEPS objects
    """

    def __init__(self, **kwargs):
        super().__init__('meshes', **kwargs)


class DefaultPointMesh(BlenderMesh):
    radius: Annotated[float, 'Radius of the sphere'] = 0.05
    subdivisions: Annotated[int, 'Number of subdivisions'] = 2

    def setUp(self, mesh, fromScratch):
        if fromScratch:
            bm = bmesh.new()
            bmesh.ops.create_icosphere(bm, subdivisions=self.subdivisions, radius=self.radius)
            bm.to_mesh(mesh)
            bm.free()


class STEPSVesicleMesh(BlenderMesh):
    subdivisions: Annotated[int, 'Number of subdivisions'] = 5

    _radius: float = None

    def setUp(self, mesh, fromScratch):
        if fromScratch:
            bm = bmesh.new()
            bmesh.ops.create_icosphere(bm, subdivisions=self.subdivisions, radius=self._radius)
            bm.to_mesh(mesh)
            bm.free()


class STEPSMesh(BlenderMesh):
    smooth_angle: Annotated[float, 'Angle for smooth shading (in radians)'] = 1.2

    _scale: float = 1
    _verts = None
    _tris = None
    _vertInds = None

    __VERT_INDS_ATTR_NAME = 'STEPSVertInds'
    __VERT_STEPS_PROPS = ['V']

    def _STEPS2BlenderMesh(self, allVerts, tris):
        v2p = {}
        vertices = []
        faces = []
        vertInds = []
        for verts in tris:
            for v in verts:
                if v not in v2p:
                    v2p[v] = len(vertices)
                    vertInds.append(v)
                    vertices.append(allVerts[v] * self._scale)
            faces.append(tuple(v2p[v] for v in verts))
        return vertices, faces, vertInds

    def setUp(self, mesh, fromScratch):
        if fromScratch:
            verts, faces, self._vertInds = self._STEPS2BlenderMesh(self._verts, self._tris)

            mesh.from_pydata(verts, [], faces)
            mesh.update()

            # Recalculate normals
            bm = bmesh.new()
            bm.from_mesh(mesh)
            bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
            bm.to_mesh(mesh)

            # Mesh attributes
            vinds = mesh.attributes.new(STEPSMesh.__VERT_INDS_ATTR_NAME, type='INT', domain='POINT')
            vinds.data.foreach_set('value', self._vertInds)
            for prop in STEPSMesh.__VERT_STEPS_PROPS:
                mesh.attributes.new(prop, type='FLOAT', domain='POINT')

            # Smooth shading
            if self.smooth_angle is not None:
                for f in mesh.polygons:
                    f.use_smooth = True
                mesh.use_auto_smooth = True
                mesh.auto_smooth_angle = self.smooth_angle
        else:
            # Retrieve vert indices from mesh
            self._vertInds = np.zeros(len(mesh.vertices), dtype=np.int32)
            mesh.attributes[STEPSMesh.__VERT_INDS_ATTR_NAME].data.foreach_get('value', self._vertInds)

    def updateVertProp(self, scene, depg, propName, values):
        self.blenderObj.attributes[propName].data[0].value = values[0]
        self.blenderObj.attributes[propName].data.foreach_set('value', values)


####################################################################################################


class BlenderCurve(BlenderWrapper):
    """Base class for Blender mesh that are used to visualize STEPS objects
    """

    def __init__(self, **kwargs):
        super().__init__('curves', **kwargs)

    def CreateBlenderObj(self, name):
        """Create and return the actual blender object"""
        return self.blenderDataCont.new(name, 'CURVE')


BlenderMesh.setUp.__doc__ = \
    """Sets up the given blender mesh object

    :param mat: the Blender mesh object
    :type mat: :py:class:`bpy.types.Mesh`
    :param fromScratch: Whether the object was created from scratch
    :type fromScratch: bool
    """


class STEPSLinkSpeciesCurve(BlenderCurve):
    bevel_depth: Annotated[float, 'Width of the link, defaults to 0.8 times the species radius'] = None

    def setUp(self, curve, fromScratch):
        if fromScratch:
            splines = curve.splines.new('POLY')
            splines.points.add(1)
            splines.points[0].co = (-1, 0, 0, 0)
            splines.points[1].co = (1, 0, 0, 0)

            curve.bevel_depth = self.bevel_depth


class STEPSVesiclePathCurve(BlenderCurve):
    path_thickness: Annotated[float, 'Thickness to the path'] = 0.01

    _data = None
    _scale: float = None

    def setUp(self, curve, fromScratch):
        if fromScratch:
            curve.dimensions = '3D'
            for idx1, (pos1, conns) in self._data.items():
                for idx2, _ in conns.items():
                    pos2, _ = self._data[idx2]
                    splines = curve.splines.new('POLY')
                    splines.points.add(1)
                    splines.points[0].co = np.array(pos1 + [0]) * self._scale
                    splines.points[1].co = np.array(pos2 + [0]) * self._scale

            curve.bevel_depth = self.path_thickness


####################################################################################################


class BlenderObject(BlenderWrapper):
    """Base class for Blender objects that are used to visualize STEPS objects
    """
    mesh: BlenderMesh = None
    material: BlenderMaterial = None

    def __init__(self, **kwargs):
        super().__init__('objects', **kwargs)

    def CreateBlenderObj(self, name):
        """Create and return the actual blender object"""
        return self.blenderDataCont.new(name, self.mesh.blenderObj if self.mesh is not None else None)

    def setUp(self, obj, fromScratch):
        """Sets up the given blender object

            :param obj: the Blender object
            :type obj: :py:class:`bpy.types.Object`
            :param fromScratch: Whether the object was created from scratch
            :type fromScratch: bool
        """
        if self.mesh is not None:
            obj.data = self.mesh.blenderObj
        if self.material is not None and self.material._name not in obj.data.materials:
            obj.data.materials.clear()
            obj.data.materials.append(self.material.blenderObj)
        if fromScratch and self.parent is not None:
            self.parent.addChild(self)


class STEPSMeshObject(BlenderObject):
    mesh: BlenderMesh = STEPSMesh
    material: BlenderMaterial = MeshMaterial
    solidifySurface: Annotated[bool, 'Whether the mesh surface should be solidified'] = True

    def setUp(self, obj, fromScratch):
        """Sets up the given blender object

            :param obj: the Blender object
            :type obj: :py:class:`bpy.types.Object`
            :param fromScratch: Whether the object was created from scratch
            :type fromScratch: bool
        """
        super().setUp(obj, fromScratch)
        if fromScratch and self.solidifySurface:
            obj.modifiers.new('solidify', type='SOLIDIFY')


class STEPSVesiclePath(BlenderObject):
    mesh: BlenderMesh = STEPSVesiclePathCurve
    material: BlenderMaterial = VesiclePathMaterial


####################################################################################################


class BlenderObjectSet(BlenderCollection):
    obj: BlenderObject = BlenderObject

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)

        self.obj.blenderObj.location = _FAR_LOCATION

    def _setPositions(self, scene, depg, positions):
        pass


class ParticleSystem(BlenderObjectSet):
    _emitter: BlenderObject = BlenderObject

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)

        self._emitter.blenderObj.location = _FAR_LOCATION

        _emitter = self._emitter.blenderObj

        if fromScratch:
            _emitter.modifiers.clear()
            _emitter.modifiers.new("particles", type='PARTICLE_SYSTEM')

        # particle system:
        particleSystem = _emitter.particle_systems[0]
        settings = particleSystem.settings

        # Reset particle system settings
        settings.emit_from = 'VERT'
        settings.physics_type = 'NEWTON'
        settings.particle_size = 1
        settings.render_type = 'OBJECT'
        settings.instance_object = self.obj.blenderObj
        settings.show_unborn = True
        settings.use_dead = True
        settings.count = 0
        settings.frame_start = -1
        settings.frame_end = -1
        settings.display_percentage = 100
        settings.normal_factor = 0
        settings.mass = 0
        settings.effector_weights.gravity = 0
        settings.effector_weights.all = 0
        self._emitter.blenderObj.show_instancer_for_render = False
        settings.use_rotations = True
        settings.use_dynamic_rotation = False
        settings.angular_velocity_mode = 'NONE'
        settings.size_random = 0
        settings.use_scale_instance = True
        settings.use_size_deflect = True

        particleSystem.point_cache.use_disk_cache = False
        particleSystem.point_cache.use_library_path = False

    def _updateCount(self, cnt):
        settings = self._emitter.blenderObj.particle_systems[0].settings
        settings.count = int(cnt)

    def _setPositions(self, scene, depg, positions):
        psys = self._emitter.blenderObj.evaluated_get(depg).particle_systems[0]
        psys.particles.foreach_set("location", positions.flatten())


class SeparateObjects(BlenderObjectSet):
    _indexes = None

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)

        self._hiddenCol = BlenderCollection(parameters=self._parameters,
                                            name=f'{self._name}_instances',
                                            parent=self)

        self._objects = {}
        hiddenColObj = self._hiddenCol.blenderObj
        mainObj = self.obj.blenderObj
        addedObjs = []
        for idx in self._indexes:
            obj_name = f'{self._name}_{idx}'
            if fromScratch:
                self._objects[idx], newObj = self.obj.blenderCopy(obj_name, obj=mainObj)
                addedObjs.append(newObj)
            else:
                self._objects[idx] = BlenderWrapper(self.obj._blendContName, name=obj_name, parameters=None)

        for obj in addedObjs:
            hiddenColObj.objects.link(obj)

    def _setPositions(self, scene, depg, positions):
        for idx, obj in self._objects.items():
            blenderObj = obj.blenderObj
            try:
                blenderObj.location = positions[idx]
                obj.setHidden(False, blenderObj)
            except KeyError:
                obj.setHidden(True, blenderObj)


class BlenderSpecies(ParticleSystem):
    obj: BlenderObject = BlenderObject.using(mesh=DefaultPointMesh, material=SpeciesMaterial)
    _emitter: BlenderObject = BlenderObject.using(mesh=DefaultPointMesh, material=None)


class BlenderLinks(SeparateObjects):
    obj: BlenderObject = BlenderObject.using(mesh=STEPSLinkSpeciesCurve, material=SpeciesMaterial)

    def _setPositions(self, scene, depg, linkPositions):
        v0 = mathutils.Vector((1, 0, 0))
        for idx, obj in self._objects.items():
            blenderObj = obj.blenderObj
            try:
                p1, p2 = linkPositions[idx]
                v1 = mathutils.Vector(p2 - p1)
                blenderObj.scale = (np.linalg.norm(p2 - p1) / 4, 1, 1)
                blenderObj.rotation_euler = v0.rotation_difference(v1).to_euler()
                blenderObj.location = p1 + (p2 - p1) / 4
                obj.setHidden(False, blenderObj)
            except KeyError:
                blenderObj.rotation_euler = (0, 0, 0)
                obj.setHidden(True, blenderObj)


class BlenderLinkSpecies(BlenderCollection):
    specs: BlenderSpecies = BlenderSpecies
    links: BlenderLinks = BlenderLinks

    def _setPositions(self, scene, depg, positions):
        self.specs._setPositions(scene, depg, positions)

    def _setLinkPositions(self, scene, depg, linkPos):
        self.links._setPositions(scene, depg, linkPos)

    def _updateCount(self, cnt):
        self.specs._updateCount(cnt)


class BlenderVesicleRafts(SeparateObjects):
    obj: BlenderObject = BlenderObject
    _locations = {}
    _specs = []

    def __init__(self, **kwargs):
        self._specSystems = {}
        # Turn the boolean modifier off for vesicles by defaults to avoid useless computation
        self._defaultBooleanModifOn = isinstance(self, BlenderRafts)
        super().__init__(**kwargs)

    def _getParticleSys(self, fromScratch, obj, specObj, psys_name, tpe='HAIR', seed=0):
        if fromScratch:
            obj.blenderObj.modifiers.new(psys_name, type='PARTICLE_SYSTEM')
        elif psys_name not in obj.blenderObj.particle_systems:
            raise Exception(f'Particle systen {psys_name} was not found in object {obj.name}')
        psys = obj.blenderObj.particle_systems[psys_name]

        settings = psys.settings

        # Reset particle system settings
        psys.seed = seed
        settings.type = tpe
        settings.render_type = 'OBJECT'
        settings.instance_object = specObj.blenderObj
        settings.particle_size = 1
        settings.count = 0
        settings.use_modifier_stack = True
        if tpe == 'HAIR':
            settings.emit_from = 'VOLUME'
            settings.hair_length = 1
            settings.distribution = 'RAND'
        if tpe == 'EMITTER':
            settings.emit_from = 'VERT'
            settings.physics_type = 'NEWTON'
            settings.show_unborn = True
            settings.use_dead = True
            settings.frame_start = -1
            settings.frame_end = -1
            settings.display_percentage = 100
            settings.normal_factor = 0
            settings.mass = 0
            settings.effector_weights.gravity = 0
            settings.effector_weights.all = 0
            settings.use_rotations = True
            settings.use_dynamic_rotation = False
            settings.angular_velocity_mode = 'NONE'
            settings.size_random = 0
            settings.use_scale_instance = True
            settings.use_size_deflect = True

            psys.point_cache.use_disk_cache = False
            psys.point_cache.use_library_path = False

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)

        if self.parent.parent.intersectAlgo != "NONE":
            for idx, obj in self._objects.items():
                blendObj = obj.blenderObj
                if fromScratch:
                    blendObj.modifiers.clear()
                    # Intersection with compartments for vesicle and patches for rafts
                    boolean = blendObj.modifiers.new('boolean', type='BOOLEAN')
                    boolean.object = self._locations[idx].blenderObj
                    boolean.operation = 'INTERSECT'
                else:
                    boolean = blendObj.modifiers['boolean']
                # Modify the algo and visibility even when loading from file
                boolean.solver = self.parent.parent.intersectAlgo
                boolean.show_viewport = self._defaultBooleanModifOn
                boolean.show_render = self._defaultBooleanModifOn
                obj._booleanModifOn = self._defaultBooleanModifOn

    def _updateSpecCounts(self, allCounts):
        for loc, counts in allCounts.items():
            for idx, cntDct in counts.items():
                for spec, cnt in cntDct.items():
                    objects, name = self._specSystems.get(loc, {}).get(spec, (None, None))
                    if name is not None:
                        ss = objects[idx].blenderObj.particle_systems[name]
                        # The seed needs to be changed for new positions to be generated
                        ss.seed = ss.seed + 1
                        ss.settings.count = cnt


class BlenderVesicles(BlenderVesicleRafts):
    obj: BlenderObject = BlenderObject.using(mesh=STEPSVesicleMesh, material=VesicleMaterial)
    innerSpecMargin: Annotated[float,
                               'Outer fraction of the vesicle radius that is free of inner species'] = 0.3
    immobileSpecs: Annotated[
        str,
        'Comma-separated list of species (without spaces) that should not be animated in between saving time points'] = ''

    def _getInnerObj(self, obj, fromScratch):
        innerObjName = f'{obj.name}_inner'
        if fromScratch:
            innerObj, blendInnerObj = obj.blenderCopy(innerObjName)
            self._hiddenCol.blenderObj.objects.link(blendInnerObj)
        else:
            if innerObjName not in self._hiddenCol.blenderObj.objects:
                raise Exception(f'{innerObjName} could not be found in {self._hiddenCol.blenderObj.name}')
            innerObj = BlenderWrapper(self.obj._blendContName, name=innerObjName, parameters=None)
            blendInnerObj = innerObj.blenderObj

        blendInnerObj.show_instancer_for_viewport = False
        blendInnerObj.show_instancer_for_render = False

        blendInnerObj.scale = (1 - self.innerSpecMargin, ) * 3

        blendInnerObj.constraints.clear()
        constr = blendInnerObj.constraints.new(type='COPY_LOCATION')
        constr.target = obj.blenderObj

        return innerObj

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)

        self._immobileSpecs = [re.compile(reg)
                               for reg in self.immobileSpecs.split(',')] if self.immobileSpecs != '' else []

        self._innerObjs = {}
        for idx, obj in self._objects.items():
            if len(self._specs) > 0:
                self._innerObjs[idx] = self._getInnerObj(obj, fromScratch)
            # particle systems:
            for i, spec in enumerate(self._specs):
                # Surface particles
                psys_name = f'{spec._name}_particles_surf'
                psys = self._getParticleSys(fromScratch,
                                            obj,
                                            spec.obj,
                                            psys_name,
                                            tpe='EMITTER',
                                            seed=3 * idx * len(self._specs) + i)
                self._specSystems.setdefault(Loc.VES_SURF, {})[spec._name] = (self._objects, psys_name)

                # Inner particles
                psys_name = f'{spec._name}_particles'
                psys = self._getParticleSys(fromScratch,
                                            self._innerObjs[idx],
                                            spec.obj,
                                            psys_name,
                                            seed=2 * idx * len(self._specs) + i)
                self._specSystems.setdefault(Loc.VES_IN, {})[spec._name] = (self._innerObjs, psys_name)

    def isSpecImmobile(self, spec):
        return any(reg.match(spec) for reg in self._immobileSpecs)

    def _setSpecPositions(self, scene, depg, positions):
        for loc, vesDct in positions.items():
            for idx, specDct in vesDct.items():
                for spec, poss in specDct.items():
                    objects, name = self._specSystems.get(loc, {}).get(spec, (None, None))
                    if name is not None:
                        eobj = objects[idx].blenderObj.evaluated_get(depg)
                        psys = eobj.particle_systems[name]
                        if self.parent.isVesUnderEvent(self._name, idx):
                            newPositions = []
                            vesPos = self.parent.getVesPos(self._name, idx) * self.parent.parent.scale
                            for pos in poss:
                                try:
                                    found, projPos, norm, fidx = eobj.closest_point_on_mesh(pos - vesPos)
                                except RuntimeError:
                                    found = False
                                if found:
                                    newPositions.append(np.array(projPos) + vesPos)
                                else:
                                    newPositions.append(_FAR_LOCATION)
                            poss = np.array(newPositions)
                        psys.particles.foreach_set("location", poss.flatten())

    def _setPositions(self, scene, depg, positions):
        super()._setPositions(scene, depg, positions)
        for idx, innerObj in self._innerObjs.items():
            innerObj.setHidden(self._objects[idx]._hidden)

    def _setEventStatus(self, scene, depg, events):
        comps = set()
        # Only turn on boolean intersection modifiers if the vesicle is undergoing some event
        for idx, obj in self._objects.items():
            if idx in events:
                if not obj._booleanModifOn:
                    for obj2 in [obj, self._innerObjs[idx]]:
                        boolean = obj2.blenderObj.modifiers['boolean']
                        boolean.show_viewport = True
                        boolean.show_render = True
                    obj._booleanModifOn = True

                    comps.add(self._locations[idx].blenderObj)
            elif obj._booleanModifOn:
                for obj2 in [obj, self._innerObjs[idx]]:
                    boolean = obj2.blenderObj.modifiers['boolean']
                    boolean.show_viewport = False
                    boolean.show_render = False
                obj._booleanModifOn = False

        # If boolean modifiers were turned on, we need to update the display of the corresponding
        # object it intersects with (the compartment in which the vesicle is). If we do not do this,
        # the boolean intersection is computed incorrectly.
        # TODO Remove this part when the above Blender issue is fixed
        for comp in comps:
            comp.hide_viewport = False
            comp.hide_render = False
        if len(comps) > 0:
            for layer in scene.view_layers:
                layer.update()
        for comp in comps:
            comp.hide_viewport = True
            comp.hide_render = True


class BlenderRafts(BlenderVesicleRafts):
    obj: BlenderObject = BlenderObject.using(mesh=STEPSVesicleMesh, material=RaftMaterial)

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)

        for idx, obj in self._objects.items():
            if fromScratch:
                # Always snap to mesh surface
                obj.blenderObj.constraints.clear()
                constr = obj.blenderObj.constraints.new(type='SHRINKWRAP')
                constr.target = self._locations[idx].blenderObj

            for i, spec in enumerate(self._specs):
                psys_name = f'{spec._name}_particles'
                psys = self._getParticleSys(fromScratch,
                                            obj,
                                            spec.obj,
                                            psys_name,
                                            seed=2 * idx * len(self._specs) + i)
                self._specSystems.setdefault(Loc.RAFT_IN, {})[spec._name] = (self._objects, psys_name)
