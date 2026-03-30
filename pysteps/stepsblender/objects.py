####################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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
import warnings

from . import utils

from .utils import Loc, progress

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

    def setAlphaBlendShadowMethod(self, mat, shadow_method='NONE'):
        if bpy.app.version < (4, 2, 0):
            if self.alpha < 1:
                mat.blend_method = 'BLEND'
            # This property was deprecated in 4.2 and removed in 4.3
            mat.shadow_method = shadow_method
        else:
            if self.alpha < 1:
                mat.surface_render_method = 'BLENDED'
            if shadow_method == 'OPAQUE':
                mat.use_transparent_shadow = False


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
            self.setAlphaBlendShadowMethod(mat)


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
            self.setAlphaBlendShadowMethod(mat)

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
    facing_blend: Annotated[float, 'Gradient of the inner color'] = 0.4

    shadow_method: Annotated[str, 'Blender shadow method (only for Blender < 4.2)'] = 'NONE'

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
            self.setAlphaBlendShadowMethod(mat, self.shadow_method)

            fresnel = nodes.new('ShaderNodeFresnel')
            fresnel.inputs['IOR'].default_value = self.fresnel_IOR

            layer_weight = nodes.new('ShaderNodeLayerWeight')
            layer_weight.inputs['Blend'].default_value = self.facing_blend

            mult = nodes.new('ShaderNodeMath')
            mult.operation = 'MULTIPLY'
            mult.inputs[0].default_value = self.fresnel_multiplier

            add = nodes.new('ShaderNodeMath')
            add.operation = 'ADD'

            node_tree.links.new(mult.inputs[1], fresnel.outputs['Fac'])
            node_tree.links.new(add.inputs[0], mult.outputs['Value'])
            node_tree.links.new(add.inputs[1], layer_weight.outputs['Facing'])
            node_tree.links.new(bsdf.inputs['Emission Strength'], add.outputs['Value'])


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
            self.setAlphaBlendShadowMethod(mat)

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
            if self.smooth_angle is not None and bpy.app.version < (4, 1, 0):
                for f in mesh.polygons:
                    f.use_smooth = True
                # These properties were removed in Blender 4.1
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


class STEPSPathLinkCurve(BlenderCurve):
    bevel_depth: Annotated[float, 'Width of the link between vesicle and vesicle path'] = 0.01

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

    cast_shadows: Annotated[str,
        'Whether the object casts shadows. Possible values: ["ON", "OFF"] (only for Blender >= 4.2)'] = "ON"

    def __init__(self, **kwargs):
        super().__init__('objects', **kwargs)

    def setShadowVisibility(self, obj):
        assert(self.cast_shadows in ['ON', 'OFF'])
        if bpy.app.version >= (4, 2, 0):
            obj.visible_shadow = self.cast_shadows == "ON"
        # Otherwise, this is handled by the material

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
        if fromScratch:
            self.setShadowVisibility(obj)
        if fromScratch and self.parent is not None:
            self.parent.addChild(self)


class STEPSMeshObject(BlenderObject):
    mesh: BlenderMesh = STEPSMesh
    material: BlenderMaterial = MeshMaterial
    surfaceThickness: Annotated[float,
        'The thickness of the mesh surface, should be greater than 0 if rafts are present on the surface'] = 0.01

    # Default cast shadows to false for meshes
    cast_shadows: Annotated[str,
        'Whether the object casts shadows. Possible values: ["ON", "OFF"] (only for Blender >= 4.2)'] = "OFF"

    def setUp(self, obj, fromScratch):
        """Sets up the given blender object

            :param obj: the Blender object
            :type obj: :py:class:`bpy.types.Object`
            :param fromScratch: Whether the object was created from scratch
            :type fromScratch: bool
        """
        super().setUp(obj, fromScratch)
        if fromScratch:
            if self.surfaceThickness > 0:
                solid = obj.modifiers.new('solidify', type='SOLIDIFY')
                solid.thickness = self.surfaceThickness
            # Smooth shading
            if bpy.app.version >= (4, 1, 0) and self.mesh.smooth_angle is not None:
                # Since Blender 4.1, the smooth shading has to be done on the object and not the mesh.
                # However, the new way to set it (bpy.ops.object.modifier_add_node_group) does not work
                # in scripts in Blender 4.1 (see https://projects.blender.org/blender/blender/issues/117399).
                # So the code below (adapted from https://blenderartists.org/t/asset-loading-is-unfinished-warning/1510459/8)
                # creates the geometry node required for smooth shading
                modifier = obj.modifiers.new("Smooth by Angle", "NODES")
                node_group = bpy.data.node_groups.new("Smooth by Angle", "GeometryNodeTree")
                node_group.interface.new_socket("Geometry", in_out="INPUT", socket_type="NodeSocketGeometry")
                input_node = node_group.nodes.new("NodeGroupInput")
                input_node.select = False
                node_group.interface.new_socket("Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry")
                output_node = node_group.nodes.new("NodeGroupOutput")
                output_node.is_active_output = True
                output_node.select = False
                node_group.is_modifier = True
                modifier.node_group = node_group

                nodes = node_group.nodes
                geom_in = nodes.get("Group Input")
                geom_out = nodes.get("Group Output")
                # Edge Angle node
                node_edge_angle = nodes.new("GeometryNodeInputMeshEdgeAngle")
                node_edge_angle.label = "Edge Angle"
                # Compare less or equal node
                node_less_than_or_equal = nodes.new("FunctionNodeCompare")
                node_less_than_or_equal.operation = "LESS_EQUAL"
                node_less_than_or_equal.inputs["B"].default_value = self.mesh.smooth_angle
                # Shade Smooth node
                node_set_smooth_edge = nodes.new("GeometryNodeSetShadeSmooth")
                node_set_smooth_edge.label = "Set Shade Smooth (Edge)"
                node_set_smooth_edge.domain = "EDGE"
                # Shade Smooth node
                node_set_smooth_face = nodes.new("GeometryNodeSetShadeSmooth")
                node_set_smooth_face.label = "Set Shade Smooth (Face)"
                node_set_smooth_face.domain = "FACE"
                node_set_smooth_face.inputs["Shade Smooth"].default_value = True
                # Linking all the nodes
                node_group.links.new(node_edge_angle.outputs["Unsigned Angle"], node_less_than_or_equal.inputs["A"])
                node_group.links.new(geom_in.outputs["Geometry"], node_set_smooth_edge.inputs["Geometry"])
                node_group.links.new(node_less_than_or_equal.outputs["Result"], node_set_smooth_edge.inputs["Shade Smooth"])
                node_group.links.new(node_set_smooth_edge.outputs["Geometry"], node_set_smooth_face.inputs["Geometry"])
                node_group.links.new(node_set_smooth_face.outputs["Geometry"], geom_out.inputs["Geometry"])


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
        for idx in progress(self._indexes, 'Add individual objects'):
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
    _linkScale: float = 0.5 # Default for link species

    def _setPositions(self, scene, depg, linkPositions):
        v0 = mathutils.Vector((1, 0, 0))
        for idx, obj in self._objects.items():
            blenderObj = obj.blenderObj
            try:
                p1, p2 = linkPositions[idx]
                v1 = mathutils.Vector(p2 - p1)
                blenderObj.scale = (np.linalg.norm(p2 - p1) / 2 * self._linkScale, 1, 1)
                blenderObj.rotation_euler = v0.rotation_difference(v1).to_euler()
                blenderObj.location = p1 + (p2 - p1) / 2 * self._linkScale
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

    def _setupParticleSys(self, fromScratch, obj, specObj, psys_name, tpe='HAIR', seed=0):
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
                intersectAlgo = self.parent.parent.intersectAlgo
                if intersectAlgo == "FAST" and bpy.app.version >= (5, 0, 0):
                    # "FAST" was renamed to "FLOAT" in Blender 5.0.0
                    # see https://projects.blender.org/blender/blender/pulls/141686
                    intersectAlgo = "FLOAT"
                try:
                    boolean.solver = intersectAlgo
                except TypeError as ex:
                    warnings.warn(
                        f'The --intersectAlgo option has an invalid value ({intersectAlgo}) '
                        f'that resulted in the following exception: {ex}'
                    )
                boolean.show_viewport = self._defaultBooleanModifOn
                boolean.show_render = self._defaultBooleanModifOn
                obj._booleanModifOn = self._defaultBooleanModifOn

    def _updateSpecCounts(self, allCounts):
        #TODO: Check if we need to do something related to immobilespecs
        totCnt = {loc: {spec._name: 0 for spec in self._specs} for loc in allCounts.keys()}
        for loc, counts in allCounts.items():
            for idx, cntDct in counts.items():
                for spec, cnt in cntDct.items():
                    totCnt[loc][spec] += cnt
        for loc, specs in totCnt.items():
            for spec, cnt in specs.items():
                objects, name = self._specSystems[loc].get(spec, (None, None))
                if name is not None:
                    ss = self.obj.blenderObj.particle_systems[name]
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

    # Only used for parameter listing
    pathLinks: BlenderLinks = None

    def setUp(self, coll, fromScratch):
        super().setUp(coll, fromScratch)

        self._immobileSpecs = [re.compile(reg)
                               for reg in self.immobileSpecs.split(',')] if self.immobileSpecs != '' else []

        defRad = self.parent.getVesRad(self._name)
        for idx, obj in progress(self._objects.items(), 'Set vesicle radii'):
            # Scale vesicles to their real radius
            rad = self.parent.getVesRad(self._name, idx)
            if rad != defRad:
                obj.blenderObj.scale = (rad / defRad,) * 3

        # particle systems:
        for i, spec in progress(enumerate(self._specs), 'Add vesicle species'):
            # Surface particles
            psys_name = f'{spec._name}_particles_surf'
            self._setupParticleSys(fromScratch,
                                   self.obj,
                                   spec.obj,
                                   psys_name,
                                   tpe='EMITTER',
                                   seed=3 * idx * len(self._specs) + i)
            # TODO: Do we still need self._objects there?
            self._specSystems.setdefault(Loc.VES_SURF, {})[spec._name] = (self._objects, psys_name)

            # Surface particles
            psys_name = f'{spec._name}_particles'
            self._setupParticleSys(fromScratch,
                                   self.obj,
                                   spec.obj,
                                   psys_name,
                                   tpe='EMITTER',
                                   seed=2 * idx * len(self._specs) + i) # TODO: Probably don't need seed anymore
            self._specSystems.setdefault(Loc.VES_IN, {})[spec._name] = (self._objects, psys_name)

        self.pathLinks = self._getParam(
            'pathLinks',
            BlenderLinks.using(
                obj=BlenderObject.using(
                    mesh=STEPSPathLinkCurve,
                    material=VesiclePathMaterial,
                ),
                _linkScale = 1,
            ),
            name=f'{self._name}_path_links',
            _indexes=self._objects.keys(),
            color=(0.5, 0.5, 0.5, 1),
        )

    def isSpecImmobile(self, spec):
        return any(reg.match(spec) for reg in self._immobileSpecs)

    def _setSpecPositions(self, scene, depg, positions):
        allPos = {}
        for loc, vesDct in positions.items():
            for idx, specDct in vesDct.items():
                if self.parent.isVesUnderEvent(self._name, idx):
                    # vesPos = self.parent.getVesPos(self._name, idx) * self.parent.parent.scale
                    eobj = self._objects[idx].blenderObj.evaluated_get(depg)
                    vesPos = np.array(eobj.location)
                    for spec, poss in specDct.items():
                        if len(poss) == 0:
                            continue
                        newPositions = []
                        for pos in poss:
                            if loc == Loc.VES_IN and utils.point_in_obj(pos - vesPos, eobj):
                                newPositions.append(pos)
                            else:
                                try:
                                    found, projPos, norm, fidx = eobj.closest_point_on_mesh(pos - vesPos)
                                except RuntimeError:
                                    found = False
                                if found:
                                    newPositions.append(np.array(projPos) + vesPos)
                                else:
                                    newPositions.append(_FAR_LOCATION)
                        allPos.setdefault(loc, {}).setdefault(spec, []).extend(newPositions)
                else:
                    for spec, poss in specDct.items():
                        if len(poss) == 0:
                            continue
                        allPos.setdefault(loc, {}).setdefault(spec, []).extend(poss)

        for loc, posDct in allPos.items():
            for spec, poss in posDct.items():
                objects, name = self._specSystems[loc].get(spec, (None, None))
                if name is not None:
                    eobj = self.obj.blenderObj.evaluated_get(depg)
                    psys = eobj.particle_systems[name]
                    psys.particles.foreach_set("location", np.array(poss).flatten())

    def _setPositions(self, scene, depg, positions):
        super()._setPositions(scene, depg, positions)

    def _setPathLinksPositions(self, scene, depg, pathLinkPos):
        self.pathLinks._setPositions(scene, depg, pathLinkPos)

    def _setEventStatus(self, scene, depg, events):
        comps = set()
        # Only turn on boolean intersection modifiers if the vesicle is undergoing some event
        for idx, obj in self._objects.items():
            if idx in events:
                if not obj._booleanModifOn:
                    boolean = obj.blenderObj.modifiers['boolean']
                    boolean.show_viewport = True
                    boolean.show_render = True
                    obj._booleanModifOn = True

                    comps.add(self._locations[idx].blenderObj)
            elif obj._booleanModifOn:
                boolean = obj.blenderObj.modifiers['boolean']
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

        for i, spec in progress(enumerate(self._specs), 'Add raft species'):
            #TODO: ?
            psys_name = f'{spec._name}_particles'
            self._setupParticleSys(fromScratch,
                                    self.obj,
                                    spec.obj,
                                    psys_name,
                                   tpe='EMITTER',
                                    )
            self._specSystems.setdefault(Loc.RAFT_IN, {})[spec._name] = (self._objects, psys_name)

        for idx, obj in self._objects.items():
            if fromScratch:
                # Always snap to mesh surface
                obj.blenderObj.constraints.clear()
                constr = obj.blenderObj.constraints.new(type='SHRINKWRAP')
                constr.target = self._locations[idx].blenderObj

    def _setPositions(self, scene, depg, positions):
        super()._setPositions(scene, depg, positions)

    def _setSpecPositions(self, scene, depg, counts):
        allPos = {}
        for idx, cnts in counts.items():
            eobj = self._objects[idx].blenderObj.evaluated_get(depg)
            rpos = np.array(eobj.location)
            for spec, cnt in cnts.items():
                poss = utils.get_points_in_sphere(rpos, self.parent.getRaftRad(self.name), cnt)
                newPositions = []
                for pos in poss:
                    try:
                        found, projPos, norm, fidx = eobj.closest_point_on_mesh(pos - rpos)
                    except RuntimeError:
                        found = False
                    if found:
                        newPositions.append(np.array(projPos) + rpos)
                    else:
                        newPositions.append(_FAR_LOCATION)
                allPos.setdefault(spec, []).extend(newPositions)

        for spec, poss in allPos.items():
            objects, name = self._specSystems[Loc.RAFT_IN].get(spec, (None, None))
            if name is not None:
                eobj = self.obj.blenderObj.evaluated_get(depg)
                psys = eobj.particle_systems[name]
                psys.particles.foreach_set("location", np.array(poss).flatten())
