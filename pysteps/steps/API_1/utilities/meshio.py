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
from __future__ import  print_function

"""
.. Note::

    This module is preliminary, means some of the functions are still under development.
    Code modification / debugging is wellcomed.
    Please email steps.dev@gmail.com if you would like to share you changes with others.

Mesh Import and Export Utilities

The MeshIO utilities are designed for importing geometry data from outside mesh-generators,
as well as creating STEPS mesh object and exporting it using STEPS' own format.

Current importing support includes:

* Abaqus data format (http://www.simulia.com)

  Abaqus format is supported by CUBIT (http://cubit.sandia.gov/)
  and Netgen (http://sourceforge.net/projects/netgen-mesher/)

* Tengen's data formats (http://tetgen.berlios.de/)

  Currently, only the .node, .ele and .face files are supported,
  other Tetgen files is not used in STEPS.

Each import function (either importTetgen or importAbaqus) will read and import data
from given data file(s), returns a Tetmesh object for STEPS simulation, and three
ElementProxy objects where geometry data and id mapping are stored.

Once a Tetmesh object has been created the data can be saved to two separate files
by the saveTetmesh() function:

* An xml annotated file containing the nodes, triangles and tetrahedra.

* A text file containing further information needed by STEPS solvers found
  in Tetmesh object constructor.

These files can later be loaded by the loadTetmesh() method and reconstruct
the Tetmesh object. This is intened to drastically reduce mesh-loading time
for large meshes (over ~100,000 voxels). By storing all data required by STEPS
internallyin these two files this infomation does not have to be found each time
by the Tetmesh object constructor.

**IMPORTANT NOTICE**

User is recommanded to save the tetmesh objects using the saveTetmesh() method
and recreate the objects from the saved files, instead of creating the objects
via the importing functions each time, if the tetmesh objects are intented to
be used in multiple simulations. Since the importing functions require a massive
amount of time to create the Tetmesh object, comparing to the loadTetmesh() method.

"""

import time
import re
import steps.API_1.geom as stetmesh
import os.path as opath
from steps.API_1.utilities.steps_shadow import *

#############################################################################################

class ElementProxy:
    """ Element Proxy Object

    The Element Proxy object contains data imported from outside files (Abaqus, Tetgen),
    as well as the mapping from the element's import id to its STEPS id, and vice versa.
    it also records block and/or group information contained in the imported file, which
    could be used for data tracing and compartment/patch creactions::

        ElementProxy(type, unitlength)

    Construct an empty ElementProxy object.

    Parameters:

    * type (string): The type of the elements stored in this object.

    * unitlength (int): The length of vector units required for each element data.
      e.g. A node with 3 coordinates requires unitlength of 3,
      a tetrahedron with 4 nodes requires unitlength of 4.
      A triangle requires unitlength of 3.

    Example::

        nodeproxy = ElementProxy('node', 3)

    """

    def __init__(self, type, unitlength):
        self.type = type
        self.unitlength = unitlength
        self.idcounter = 0
        self.data = []
        self.importid = []
        self.stepsid = {}
        self.blocks = {}
        self.groups = {}
        self.tempblockname = ''
        self.tempblockstart = -1

    def insert(self, import_id, import_data):
        """
        Insert an element to the Element Map object, a STEPS id will be assigned automatically.

        Parameters:

        * import_id (int): The id of the element given by the importing file.
        * import_data (list): a list of data belongs to the element,
          e.g. for a node, a 3 units list of its coordinates.

        Example:

        nodeproxy.insert(1, [0.1, 0.2, 0.4])

        Insert a node with import id 1 and coordinates x = 0.1, y = 0.2, z = 0.4 to nodeproxy.
        type nodeproxy.getSTEPSID(1) will return the STEPS id of this element.

        """
        self.data += import_data
        self.importid.append(import_id)
        self.stepsid[import_id] = self.idcounter
        self.idcounter += 1

    def getType(self):
        """
        Return the type of elements stored in this Element Map object.

        """
        return self.type

    def getSize(self):
        """
        Return the number of elements stored in the Element Map object.

        """
        return self.idcounter

    def getDataFromSTEPSID(self, steps_id):
        """
        Return the data of element with STEPS id steps_id.

        Parameters:

        * steps_id (int): The STEPS id of the element, this id is automatically allocated when
          inserting the element to the Element Map (see insert() method).

        """
        return self.data[steps_id * self.unitlength : (steps_id + 1) * self.unitlength]

    def getDataFromImportID(self, import_id):
        """
        Return the data of element with import id import_id.

        Parameters:

        import_id (int):
        The import id of the element, this id is given by outside generator
        and contained in the imported file.

        """
        return self.getDataFromSTEPSID(self.stepsid[import_id])

    def getAllData(self):
        """
        Return all data as an one dimentional list.
        This list can be used to construct the Tetmesh.

        Example::

            nodedata = nodeproxy.getAllData()
            tetdata = tetproxy.getAllData()
            mesh = steps.geom.Tetmesh(nodedata, tetdata)

        """
        return self.data

    def getSTEPSID(self, import_id):
        """
        Return the STEPS id of a element from its import id.

        Parameters:

        * import_id (int): The import id of the element.

        """
        return self.stepsid[import_id]

    def getImportID(self, steps_id):
        """
        Return the import id of a element from its STEPS id.

        Parameters:

        * steps_id (int): The STEPS id of the element.

        """

        return self.importid[steps_id]

    def addGroup(self, group_name, group_ids):
        """
        Add a group of elements in the Element Map object.

        An group is defined as a list of element ids. Unlike blocks (see blockBegin() method),
        the ids of a group need not be continual. For this reason, user should provide a list
        which contains all ids of elements belong to this group. Each group has a unique name
        for data access.

        No uniqueness or order is predefined for a group, this can be defined by the user for
        specific usage.

        Groups can be used to construct compartments and patches in STEPS, please refer to the
        user manual for examples about using groups and blocks.

        Parameters:

        * group_name (string): Name of the element group.
        * group_ids (list): A list of ids of elements belong to the group.

        """

        self.groups[group_name] = group_ids

    def addToGroup(self, group_name, id):
        if (group_name in self.groups):
            self.groups[group_name].append(id)
        else:
            self.groups[group_name] = [id]

    def getGroups(self):
        """
        Return all groups stored in the Element Map object.

        Example:

        The following group dictionary contains two groups::

            "Group1" contains element ids 1,2,6,4,8
            "Group2" contains element ids 3,1,2,1

            groups = {"Group1":[1,2,6,4,8], "Group2":[3,1,2,1]}

        Note that element 1 appears twice in "Group2".

        For meshes imported from Abaqus files, the keys of the groups
        are fetched from the ELSET id.

        For meshes imported from Gmsh2.2 files, the keys of the groups
        are pairs of (tag_idx, tag_value), where tag_idx is the gmsh tag
        index position of the element, and tag_value is the 
        tag value of the element. For example, if a tetrahedron is defined
        as below in Gmsh

            1 4 2 2 3 16 17 13 18

        It has two tags, tag[0] = 2 and tag[1] = 3, thus the element is 
        contained in both groups[(0, 2)] and groups[(1, 3)].
        Note that if a specfic name string is mapped to a value in the 
        physical tag (tag[0]) in the file, the name string can also be used
        to access the same group, for example, if the 
        following is defined in the file

            $PhysicalNames
            1
            3 2 "inner"
            $EndPhysicalNames

        groups[(0, 2)] can also be accessed using groups[(0, "inner")]
        """

        return self.groups

    def getGroup(self, group_key):
        """
        Return the group stored in the Element Map object with key = group_key.
        If the key is not present in the Element Map object, an empty list is
        returned.

        Example:
            g1 = getGroup("Group1")
            g2 = getGroup("Group2")
        """

        if group_key not in self.groups.keys():
            return []
        else:
            return self.groups[group_key]

    def blockBegin(self, block_name):
        """
        Notify the Element Map object that a new block with name block_name begins from
        the next element.

        A block is a special type of group whose element's STEPS ids is continual,
        i.e. can be represented by a begin id and an end id. An unique name should be
        given to each block, which can be used to access associated block information
        (i.e. a two unit list contains the block begin id and block end id).

        If another block was initialized before calling the blockEnd() method,
        this method will also end the previous block and the its data will be stored.

        A dictionary of blocks can be converted to a dictionary of group via the
        blocksToGroups() method.

        Notice:

        Although a block will be finished and stored when a new block is notified without
        calling the blockEnd() method, user is recommanded to call blockEnd() in the usual
        case.

        Parameters:

        * block_name (string): Name of the element block.

        """

        if (self.tempblockstart != -1):
            self.blocks[self.tempblockname] = [self.tempblockstart, self.idcounter - 1]
        self.tempblockname = block_name
        self.tempblockstart = self.idcounter

    def blockEnd(self):
        """
        Notify the Element Map object that the current element block should be ended.
        Block information (i.e. a dictionary element blockname:[begin_id, end_id]) will
        be added to the ElementProxy object.

        Notice:

        If the blockEnd() method is called before any block begins, or a block has been
        ended and there is no new block started, no information will be added to the block
        dictionary.

        """
        if (self.tempblockstart != -1):
            self.blocks[self.tempblockname] = [self.tempblockstart, self.idcounter - 1]
            self.tempblockstart = -1

    def getBlocks(self):
        """
        Return the block dictionary stored in the Element Map object.

        A block dictionary uses the unique block names as keys, and a list of the start and
        end ids as values.

        Example:

        The following block dictionary contains two blocks: "Block1" starts from element with
        id of 0 and end at element 100, "Block2" starts from element 101 and ends at element 110::

            blocks = { "Block1":[0, 100], "Block2":[100, 110] }

        To access individual blocks, for example "Block1", type::

            block1 = blocks["Block1"]

        User is recommanded to check Python's manual for the usage of dictionary.

        """

        return self.blocks


    def blocksToGroups(self):
        """
        Convert the block dictionary to a group dictionary.

        Return:

        A dictionary of groups.

        """

        it = iter(self.blocks)
        converted_groups = {}
        for key in it:
            blockrange = self.blocks[key]
            converted_groups[key] = list(range(blockrange[0], blockrange[1] + 1))
        return converted_groups



#############################################################################################

def importTetGen(pathroot, scale, verbose = False):
    """
    Read a TetGen-generated or refined mesh from a set of files.

    The TetGen file format for describing a mesh actually consists of
    multiple files with varying suffixes. Currently, this routine only
    reads meshes consisting of three files:

    * <input>.node: describing the tetrahedral mesh node points.
    * <input>.ele: describing tetrahedral elements, each of which
      consists of 4 pointers into the <input>.node file. (TetGen
      also supports 10-node elements; these 6 extra nodes are obviously
      not read by STEPS.)
    * <input>.face: describing triangular faces, each of which
      consists of 3 pointers into the <input>.node file. This file is optional.

    Other files are .vol (list of maximum volumes), .var (variant constraints)
    .neigh (list of neighbours), .smesh (simple PLC descriptions) and .edge
    (list of boundary edges) files. None of these seem relevant for our
    use cases, so we don't load them even when they are there. In particular,
    the .neigh file is computed by STEPS itself.

    Please refer to the TetGen manual (pages 31-40 in the last edition)
    for more information on these file formats

    tetgen.berlios.de/files/tetgen-manual.pdf

    (See the documentation for steps.geom.tetmesh for details about the mesh object.)

    PARAMETERS:

    * pathroot
      The root of the path name. E.g. mesh/torus would make this
      routine try to read files mesh/torus.node, mesh/torus.ele
      and optionally for mesh/torus.face

    * scale: LENGTH scale from the importing mesh to real geometry. e.g. a radius
      of 10 in the importing file to a radius of 1 micron in STEPS, scale is 1e-7.

    RETURNS:

    mesh, nodeproxy, tetproxy, triproxy

    mesh: The STEPS TetMesh object
    nodeproxy: Element Map for nodes
    tetproxy: Element Map for tetrahedrons
    triproxy: Element Map for triangles

   IMPORTANT NOTICE:

   User is recommanded to save the tetmesh objects using the saveTetmesh() method
   and recreate the objects from the saved files, instead of creating the objects
   via the importing functions each time, if the tetmesh objects are intented to
   be used in multiple simulations. Since the importing functions require a massive
   amount of time to create the Tetmesh object, comparing to the loadTetmesh() method.

    """
    nodefname = pathroot + '.node'
    elefname = pathroot + '.ele'
    facefname = pathroot + '.face'

    nodeproxy = ElementProxy('node', 3)
    tetproxy = ElementProxy('tet', 4)
    triproxy = ElementProxy('tri', 3)

    # Is there a .node file?
    if not opath.isfile(nodefname):
        print(nodefname)
        return None
    if not opath.isfile(elefname):
        print(elefname)
        return None
    if not opath.isfile(facefname):
        facefname = ''

    # Try to read the node file.
    nodefile = open(nodefname, 'r')
    # First line is:  <x> <y> <z> [att<# of points> <dimension (3)> <# of attributes>
    #                <boundary marker (0 or 1)>
    line = nodefile.readline()
    while (line[0] == '#' or line[0] == '\n'):
        line = nodefile.readline()
    tokens = line.split()
    assert len(tokens) == 4
    nnodes = int(tokens[0])
    assert nnodes > 0
    ndims = int(tokens[1])
    assert ndims == 3
    nattribs = int(tokens[2])
    bmarkers = int(tokens[3])
    while (nodeproxy.getSize() != nnodes):
        line = nodefile.readline()
        commentstart = line.find('#')
        if commentstart != -1:
            line = line[0:commentstart]
        # Remaing lines: <point #>ributes]
        #                [boundary marker]
        tokens = line.split()
        if len(tokens) == 0:
            continue
        nodeid = int(tokens[0])
        node = [0.0, 0.0, 0.0]
        node[0] = float(tokens[1]) * scale
        node[1] = float(tokens[2]) * scale
        node[2] = float(tokens[3]) * scale
        nodeproxy.insert(nodeid, node)
    # Close the file.
    nodefile.close()

    # Try to read the tet file.
    elefile = open(elefname, 'r')
    line = elefile.readline()
    while (line[0] == '#' or line[0] == '\n'):
        line = elefile.readline()
    tokens = line.split()
    assert len(tokens) == 3
    ntets = int(tokens[0])
    assert ntets > 0
    ndims = int(tokens[1])
    assert ndims == 4
    nattribs = int(tokens[2])
    while (tetproxy.getSize() != ntets):
        line = elefile.readline()
        commentstart = line.find('#')
        if commentstart != -1:
            line = line[0:commentstart]
        # Remaing lines: <point #>ributes]
        #                [boundary marker]
        tokens = line.split()
        if len(tokens) == 0:
            continue
        tetid = int(tokens[0])
        tet = [0, 0, 0, 0]
        tet[0] = nodeproxy.getSTEPSID(int(tokens[1]))
        tet[1] = nodeproxy.getSTEPSID(int(tokens[2]))
        tet[2] = nodeproxy.getSTEPSID(int(tokens[3]))
        tet[3] = nodeproxy.getSTEPSID(int(tokens[4]))
        tetproxy.insert(tetid, tet)
        if nattribs == 1:
            tetproxy.addToGroup(tokens[5], tetproxy.getSTEPSID(tetid))
    elefile.close()

    # Try to read the tri file.
    if (facefname != ''):
        facefile = open(facefname, 'r')
        line = facefile.readline()
        while (line[0] == '#' or line[0] == '\n'):
            line = facefile.readline()
        tokens = line.split()
        assert len(tokens) == 2
        ntris = int(tokens[0])
        assert ntris > 0
        nattribs = int(tokens[1])
        while (triproxy.getSize() != ntris):
            line = facefile.readline()
            commentstart = line.find('#')
            if commentstart != -1:
                line = line[0:commentstart]
            # Remaing lines: <point #>ributes]
            #                [boundary marker]
            tokens = line.split()
            if len(tokens) == 0:
                continue
            triid = int(tokens[0])
            tri = [0, 0, 0]
            tri[0] = nodeproxy.getSTEPSID(int(tokens[1]))
            tri[1] = nodeproxy.getSTEPSID(int(tokens[2]))
            tri[2] = nodeproxy.getSTEPSID(int(tokens[3]))
            triproxy.insert(triid, tri)
            if nattribs == 1:
                triproxy.addToGroup(tokens[4], triproxy.getSTEPSID(triid))
        # Close the file.
        facefile.close()

    if (verbose): print("Read TetGen files succesfully")

    nodedata = nodeproxy.getAllData()
    tetdata = tetproxy.getAllData()
    tridata = triproxy.getAllData()

    if (verbose): print("creating Tetmesh object in STEPS...")
    mesh = stetmesh.Tetmesh(nodedata, tetdata, tridata)
    if (verbose): print("Tetmesh object created.")
    return mesh, nodeproxy, tetproxy, triproxy


#############################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#############################################################################################

def parseAbaqusLine(line):
    line = line.replace("\n", "")
    line = line.replace("\r", "")
    type = "undefined"
    result = {}
    if line == "":
        pass
    elif (line.find('**', 0, 2) == 0):
        type = "comment"
    # keyword
    elif (line.find('*', 0, 1) == 0):
        type = "keyword"
        line = line.split(",")
        for seg in line:
            if seg[0] == "*":
                result["keyword"] = seg[1:].upper()
            elif "=" in seg:
                seg = seg.replace(" ", "").split("=")
                result[seg[0].upper()] = seg[1]
            else:
                result[seg.upper()] = None
    # data
    else:
        type = "data"
        result["data"] = line.replace(" ", "").split(",")
    return type, result

#############################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#############################################################################################

def importAbaqus(filename, scale, ebs = None, shadow_mesh = None, verbose = False):
    """
    Read a ABAQUS-formated mesh file, return the created steps.geom.Tetmesh object,
    the element mapping for nodes, tetraedrons and triangles.

    PARAMETERS:

    * filename: the Abaqus filename (or path) including any suffix.
    * scale: LENGTH scale from the importing mesh to real geometry. e.g. a radius
      of 10 in the importing file to a radius of 1 micron in STEPS, scale is 1e-7.
    * ebs: specify the names of selected element blocks which are included in the mesh.
    * shadow_mesh: name of the ShadowMesh file exported using the STEPS-CUBIT supporting toolkit, can also be the ShadowMesh object itself.

    RETURNS:

    mesh, nodeproxy, tetproxy, triproxy

    * mesh: The STEPS TetMesh object
    * nodeproxy: Element Map for nodes
    * tetproxy: Element Map for tetrahedrons
    * triproxy: Element Map for triangles

   IMPORTANT NOTICE:

   User is recommanded to save the tetmesh objects using the saveTetmesh() method
   and recreate the objects from the saved files, instead of creating the objects
   via the importing functions each time, if the tetmesh objects are intented to
   be used in multiple simulations. Since the importing functions require a massive
   amount of time to create the Tetmesh object, comparing to the loadTetmesh() method.

    """
    if (verbose): print("Reading Abaqus file...")
    btime = time.time()

    abaqusfile = open(filename, 'r')

    nodeproxy = ElementProxy('node', 3)
    tetproxy = ElementProxy('tet', 4)
    triproxy = ElementProxy('tri', 3)

    line = abaqusfile.readline()
    # check we have the right kind of CUBIT output file here
    #assert(line.find('*HEADING'))

    if ebs == None:
        include = True
    else:
        include = False
    currmap = None

    # add tri data record for remove duplicated triangles
    recorded_tris = {}

    # read line one by one until the end
    while(line):
        type, result = parseAbaqusLine(line)
        # status check
        if type == "comment":
            pass
        elif type == "keyword" and (result["keyword"] == "NODE"):
            if (currmap != None):
                currmap.blockEnd()
            if "ELSET" in result.keys():
                nodeproxy.blockBegin(result["ELSET"])
            else:
                nodeproxy.blockBegin("AllNodes")
            currmap = nodeproxy
        elif type == "keyword" and (result["keyword"] == "ELEMENT"):
            if (currmap != None):
                currmap.blockEnd()
            if "ELSET" in result.keys():
                blockname = result["ELSET"]
            else:
                blockname = "AllElements"
            if (ebs == None) or (blockname in ebs):
                include = True
            else:
                include  = False
            if (result["TYPE"] == 'C3D4') and include:
                currmap = tetproxy
                tetproxy.blockBegin(blockname)
            elif (result["TYPE"] == 'STRI3') and include:
                currmap = triproxy
                triproxy.blockBegin(blockname)
            elif currmap != None:
                currmap.blockEnd()
                currmap = None
        elif (type == "keyword") and (currmap != None):
            currmap.blockEnd()
            currmap = None
        elif type == "data" and currmap != None:
            if (currmap == nodeproxy):
                nodeid = int(result["data"][0])
                node = [0.0, 0.0, 0.0]
                # set to the right scale
                node[0] = float(result["data"][1])*scale
                node[1] = float(result["data"][2])*scale
                node[2] = float(result["data"][3])*scale
                currmap.insert(nodeid, node)
            elif (currmap == tetproxy) and include:
                nodeid = int(result["data"][0])
                node = [0,0,0,0]
                node[0] = nodeproxy.getSTEPSID(int(result["data"][1]))
                node[1] = nodeproxy.getSTEPSID(int(result["data"][2]))
                node[2] = nodeproxy.getSTEPSID(int(result["data"][3]))
                node[3] = nodeproxy.getSTEPSID(int(result["data"][4]))
                currmap.insert(nodeid, node)
            elif (currmap == triproxy) and include:
                nodeid = int(result["data"][0])
                node = [0,0,0]
                node[0] = nodeproxy.getSTEPSID(int(result["data"][1]))
                node[1] = nodeproxy.getSTEPSID(int(result["data"][2]))
                node[2] = nodeproxy.getSTEPSID(int(result["data"][3]))
                if (node[0], node[1], node[2]) in recorded_tris:
                    if (verbose): print("Triangle: ", (node[0], node[1], node[2]), " with index ", recorded_tris[(node[0], node[1], node[2])], " has been imported, ignore duplicated triangle ", nodeid)
                else:
                    currmap.insert(nodeid, node)
                    recorded_tris[(node[0], node[1], node[2])] = nodeid

        line = abaqusfile.readline()

    if (currmap != None):
        currmap.blockEnd()

    abaqusfile.close()

    nodedata = nodeproxy.getAllData()
    tetdata = tetproxy.getAllData()
    tridata = triproxy.getAllData()

    if (verbose): print("Number of nodes imported: ", nodeproxy.getSize())
    if (verbose): print("Number of tetrahedrons imported: ", tetproxy.getSize())
    if (verbose): print("Number of triangles imported: ", triproxy.getSize())

    if (verbose): print("creating Tetmesh object in STEPS...")
    mesh = stetmesh.Tetmesh(nodedata, tetdata, tridata)
    if (verbose): print("Tetmesh object created.")

    if shadow_mesh != None:
        if isinstance(shadow_mesh, str):
            shadow_mesh = ShadowMesh.importFrom(shadow_mesh)

        if (verbose): print("Importing data from shadow mesh.")
        for c in shadow_mesh.comps.values():
            if (verbose): print("Import compartment ", c.name)
            steps_indices = [tetproxy.getSTEPSID(i) for i in c.indices]
            comp = stetmesh.TmComp(c.name, mesh, steps_indices)
            for v in c.vsyss:
                comp.addVolsys(v)
        if len(shadow_mesh.patches.keys()) != 0:
            print("ImportAbaqus does not support patch importing, use ImportAbaqus2 instead.")
        for roi in shadow_mesh.rois.items():
            roi_name = roi[0]
            roi_type = roi[1]["Type"]
            roi_import_indices = roi[1]["Indices"]
            if (verbose): print("Import ROI data ", roi[0])
            if roi_type == ELEM_VERTEX:
                steps_indices = [nodeproxy.getSTEPSID(i) for i in roi_import_indices]
            elif roi_type == ELEM_TET:
                steps_indices = [tetproxy.getSTEPSID(i) for i in roi_import_indices]
            elif roi_type == ELEM_TRI:
                print("ImportAbaqus does not support triangle ROI importing, use ImportAbaqus2 instead.")
            else:
                steps_indices = roi_import_indices
            mesh.addROI(roi_name, roi_type, steps_indices)

    return mesh, nodeproxy, tetproxy,triproxy

def importAbaqus2(tetfilename, trifilename, scale, shadow_mesh = None, verbose = False):
    """
        Read two ABAQUS-formated mesh files, one with tetrahedron data and the other with triangle data, return the created steps.geom.Tetmesh object,
        the element mapping for nodes, tetraedrons and triangles.

        PARAMETERS:

        * tetfilename: the Abaqus filename for tetrahedron data including any suffix.
        * trifilename: the Abaqus filename for triangle data including any suffix.
        * scale: LENGTH scale from the importing mesh to real geometry. e.g. a radius
          of 10 in the importing file to a radius of 1 micron in STEPS, scale is 1e-7.
        * shadow_mesh: name of the ShadowMesh file exported using the STEPS-CUBIT supporting toolkit, can also be the ShadowMesh object itself.

        RETURNS:

        mesh, nodeproxy, tetproxy, triproxy

        * mesh: The STEPS TetMesh object
        * nodeproxy: Element Map for nodes
        * tetproxy: Element Map for tetrahedrons
        * triproxy: Element Map for triangles

        IMPORTANT NOTICE:

        User is recommanded to save the tetmesh objects using the saveTetmesh() method
        and recreate the objects from the saved files, instead of creating the objects
        via the importing functions each time, if the tetmesh objects are intented to
        be used in multiple simulations. Since the importing functions require a massive
        amount of time to create the Tetmesh object, comparing to the loadTetmesh() method.

    """
    if (verbose): print("Reading Abaqus file...")
    btime = time.time()

    tetfile = open(tetfilename, 'r')

    nodeproxy = ElementProxy('node', 3)
    tetproxy = ElementProxy('tet', 4)
    triproxy = ElementProxy('tri', 3)

    vert_idxs = {}

    line = tetfile.readline()
    # check we have the right kind of CUBIT output file here
    #assert(line.find('*HEADING'))

    currmap = None
    # read line one by one until the end
    while(line):
        type, result = parseAbaqusLine(line)
        # status check
        if type == "comment":
            pass
        elif type == "keyword" and (result["keyword"] == "NODE"):
            if (currmap != None):
                currmap.blockEnd()
            #if (verbose): print('Found *NODE section, start reading nodes.')
            if "ELSET" in result.keys():
                nodeproxy.blockBegin(result["ELSET"])
            else:
                nodeproxy.blockBegin("AllNodes")
            currmap = nodeproxy
        elif type == "keyword" and (result["keyword"] == "ELEMENT"):
            if (currmap != None):
                currmap.blockEnd()
            if "ELSET" in result.keys():
                blockname = result["ELSET"]
            else:
                blockname = "AllElements"
            if result["TYPE"] == 'C3D4':
                currmap = tetproxy
                tetproxy.blockBegin(blockname)
            elif result["TYPE"] == 'STRI3':
                print("This importing function is not designed for file with both tetrahedrons and triangles, try importAbaqus instead.")
                return
            elif currmap != None:
                currmap.blockEnd()
                currmap = None
        elif (type == "keyword") and (currmap != None):
            currmap.blockEnd()
            currmap = None
        elif type == "data" and currmap != None:
            if (currmap == nodeproxy):
                nodeid = int(result["data"][0])
                node = [0.0, 0.0, 0.0]
                # set to the right scale
                node[0] = float(result["data"][1])*scale
                node[1] = float(result["data"][2])*scale
                node[2] = float(result["data"][3])*scale
                currmap.insert(nodeid, node)
                vert_idxs[tuple(node)] = currmap.getSize() - 1

            elif (currmap == tetproxy):
                nodeid = int(result["data"][0])
                node = [0,0,0,0]
                # set to the right scale
                node[0] = nodeproxy.getSTEPSID(int(result["data"][1]))
                node[1] = nodeproxy.getSTEPSID(int(result["data"][2]))
                node[2] = nodeproxy.getSTEPSID(int(result["data"][3]))
                node[3] = nodeproxy.getSTEPSID(int(result["data"][4]))
                currmap.insert(nodeid, node)

        line = tetfile.readline()

    if (currmap != None):
        currmap.blockEnd()

    tetfile.close()

    trifile = open(trifilename, 'r')

    line = trifile.readline()
    # check we have the right kind of CUBIT output file here
    #assert(line.find('*HEADING'))
    # add tri data record for remove duplicated triangles
    recorded_tris = {}
    currmap = None
    # read line one by one until the end
    while(line):
        # status check
        type, result = parseAbaqusLine(line)
        if type == "comment":
            pass
        elif type == "keyword" and (result["keyword"] == "NODE"):
            if (currmap != None):
                currmap.blockEnd()
            #if (verbose): print('Found *NODE section, start reading nodes.')
            if "ELSET" in result.keys():
                nodeproxy.blockBegin(result["ELSET"])
            else:
                nodeproxy.blockBegin("AllNodes")
            currmap = nodeproxy
        elif type == "keyword" and (result["keyword"] == "ELEMENT"):
            if (currmap != None):
                currmap.blockEnd()
            if "ELSET" in result.keys():
                blockname = result["ELSET"]
            else:
                blockname = "AllElements"
            if result["TYPE"] == 'STRI3':
                currmap = triproxy
                triproxy.blockBegin(blockname)

        elif (type == "keyword") and (currmap != None):
            currmap.blockEnd()
            currmap = None

        elif type == "data":
            if (currmap == nodeproxy):
                nodeid = int(result["data"][0])
                node = [0.0, 0.0, 0.0]
                # set to the right scale
                node[0] = float(result["data"][1])*scale
                node[1] = float(result["data"][2])*scale
                node[2] = float(result["data"][3])*scale
                if tuple(node) not in vert_idxs.keys():
                    print("coordinates ", node, " is not in the tet file, abort function.")
                    raise
            elif (currmap == triproxy):
                nodeid = int(result["data"][0])
                node = [0,0,0]
                # set to the right scale
                node[0] = nodeproxy.getSTEPSID(int(result["data"][1]))
                node[1] = nodeproxy.getSTEPSID(int(result["data"][2]))
                node[2] = nodeproxy.getSTEPSID(int(result["data"][3]))
                if (node[0], node[1], node[2]) in recorded_tris:
                    if (verbose): print("Triangle: ", (node[0], node[1], node[2]), " with index ", recorded_tris[(node[0], node[1], node[2])], " has been imported, ignore duplicated triangle ", nodeid)
                else:
                    currmap.insert(nodeid, node)
                    recorded_tris[(node[0], node[1], node[2])] = nodeid

        line = trifile.readline()

    if (currmap != None):
        currmap.blockEnd()

    trifile.close()


    nodedata = nodeproxy.getAllData()
    tetdata = tetproxy.getAllData()
    tridata = triproxy.getAllData()

    if (verbose): print("Number of nodes imported: ", nodeproxy.getSize())
    if (verbose): print("Number of tetrahedrons imported: ", tetproxy.getSize())
    if (verbose): print("Number of triangles imported: ", triproxy.getSize())

    if (verbose): print("creating Tetmesh object in STEPS...")
    mesh = stetmesh.Tetmesh(nodedata, tetdata, tridata)
    if (verbose): print("Tetmesh object created.")

    if shadow_mesh != None:
        if isinstance(shadow_mesh, str):
            shadow_mesh = ShadowMesh.importFrom(shadow_mesh)

        if (verbose): print("Importing data from shadow mesh.")
        temp_comps = {}
        for c in shadow_mesh.comps.values():
            if (verbose): print("Import compartment ", c.name)
            steps_indices = [tetproxy.getSTEPSID(i) for i in c.indices]
            comp = stetmesh.TmComp(c.name, mesh, steps_indices)
            temp_comps[c] = comp
            for v in c.vsyss:
                comp.addVolsys(v)
        for p in shadow_mesh.patches.values():
            if (verbose): print("Import patch ", p.name)
            steps_indices = [triproxy.getSTEPSID(i) for i in p.indices]
            icomp = temp_comps[p.icomp]
            ocomp = None
            if p.ocomp != None:
                ocomp = temp_comps[p.ocomp]
            patch = stetmesh.TmPatch(p.name, mesh, steps_indices, icomp, ocomp)
            for s in p.ssyss:
                patch.addSurfsys(s)
        for roi in shadow_mesh.rois.items():
            roi_name = roi[0]
            roi_type = roi[1]["Type"]
            roi_import_indices = roi[1]["Indices"]
            if (verbose): print("Import ROI data ", roi[0])
            if roi_type == ELEM_VERTEX:
                steps_indices = [nodeproxy.getSTEPSID(i) for i in roi_import_indices]
            elif roi_type == ELEM_TET:
                steps_indices = [tetproxy.getSTEPSID(i) for i in roi_import_indices]
            elif roi_type == ELEM_TRI:
                steps_indices = [triproxy.getSTEPSID(i) for i in roi_import_indices]
            else:
                steps_indices = roi_import_indices
            mesh.addROI(roi_name, roi_type, steps_indices)

    return mesh, nodeproxy, tetproxy,triproxy


#############################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#############################################################################################

def importGmsh(filename, scale, verbose = False):
    """
    Read a Gmsh-formated mesh file, return the created steps.geom.Tetmesh object,
    the element mapping for nodes, tetraedrons and triangles.

    PARAMETERS:

    * filename: the Abaqus filename (or path) including any suffix.
    * scale: LENGTH scale from the importing mesh to real geometry. e.g. a radius
      of 10 in the importing file to a radius of 1 micron in STEPS, scale is 1e-7.

    RETURNS:

    mesh, nodeproxy, tetproxy, triproxy

    * mesh: The STEPS TetMesh object
    * nodeproxy: Element Map for nodes
    * tetproxy: Element Map for tetrahedrons
    * triproxy: Element Map for triangles

   IMPORTANT NOTICE:

   User is recommanded to save the tetmesh objects using the saveTetmesh() method
   and recreate the objects from the saved files, instead of creating the objects
   via the importing functions each time, if the tetmesh objects are intented to
   be used in multiple simulations. Since the importing functions require a massive
   amount of time to create the Tetmesh object, comparing to the loadTetmesh() method.

    """
    if (verbose): print("Reading Gmsh file...")
    btime = time.time()

    meshfile = open(filename, 'r')

    nodeproxy = ElementProxy('node', 3)
    tetproxy = ElementProxy('tet', 4)
    triproxy = ElementProxy('tri', 3)

    physical_name_mapping = {}

    with open(filename, 'r') as meshfile:
        line = meshfile.readline()
        while line:
            if line == '$MeshFormat\n':
                line = meshfile.readline()
                linesec = line.split()
                if linesec[0] != "2.2" or linesec[1] != "0":
                    raise Exception('This importer only supports Gmsh 2.2 ASCII format.')
                line = meshfile.readline()
                assert(line == '$EndMeshFormat\n')
            elif line == '$PhysicalNames\n':
                line = meshfile.readline()
                nphysical_names = int(line)
                for n in range(nphysical_names):
                    line = meshfile.readline()
                    linesec = line.split()
                    physical_name_mapping[(int(linesec[0]), int(linesec[1]))] = linesec[2].replace('"', '')
                line = meshfile.readline()
                assert(line == '$EndPhysicalNames\n')
            elif line == '$Nodes\n':
                line = meshfile.readline()
                n_nodes = int(line)
                # read nodes
                for i in range(n_nodes):
                    line = meshfile.readline()
                    linesec = line.split()
                    id = int(linesec[0])
                    node = [0.0,0.0,0.0]
                    node[0] = float(linesec[1])*scale
                    node[1] = float(linesec[2])*scale
                    node[2] = float(linesec[3])*scale
                    nodeproxy.insert(id, node)
                line = meshfile.readline()
                assert(line == '$EndNodes\n')
            elif line == '$Elements\n':
                line = meshfile.readline()
                n_elems = int(line)
                # read elem
                for i in range(n_elems):
                    line = meshfile.readline()
                    linesec = line.split()
                    id = int(linesec[0])
                    type = int(linesec[1])
                    ntags = int(linesec[2])
                    # triangle
                    if type == 2:
                        node = [0,0,0]
                        node[0] = nodeproxy.getSTEPSID(int(linesec[-3]))
                        node[1] = nodeproxy.getSTEPSID(int(linesec[-2]))
                        node[2] = nodeproxy.getSTEPSID(int(linesec[-1]))
                        triproxy.insert(id, node)
                        steps_id = triproxy.getSTEPSID(id)
                        for tag in range(ntags):
                            tag_id = int(linesec[3+tag])
                            triproxy.addToGroup((tag, tag_id), steps_id)
                    # tet
                    if type == 4:
                        group_id = linesec[4]
                        node = [0,0,0,0]
                        node[0] = nodeproxy.getSTEPSID(int(linesec[-4]))
                        node[1] = nodeproxy.getSTEPSID(int(linesec[-3]))
                        node[2] = nodeproxy.getSTEPSID(int(linesec[-2]))
                        node[3] = nodeproxy.getSTEPSID(int(linesec[-1]))
                        tetproxy.insert(id, node)
                        steps_id = tetproxy.getSTEPSID(id)
                        for tag in range(ntags):
                            tag_id = int(linesec[3+tag])
                            tetproxy.addToGroup((tag, tag_id), steps_id)
                line = meshfile.readline()
                assert(line == '$EndElements\n')
            line = meshfile.readline()
    meshfile.close()

    for physical_tag in physical_name_mapping:
        dimension = physical_tag[0]
        tag = physical_tag[1]
        tag_name = physical_name_mapping[(dimension, tag)]
        if dimension == 0:
            nodeproxy.addGroup((0, tag_name), nodeproxy.getGroup((0, tag)))
        elif dimension == 2:
            triproxy.addGroup((0, tag_name), triproxy.getGroup((0, tag)))
        elif dimension == 3:
            tetproxy.addGroup((0, tag_name), tetproxy.getGroup((0, tag)))
                                    
    if (verbose): print("Read Msh file succesfully")

    nodedata = nodeproxy.getAllData()
    tetdata = tetproxy.getAllData()
    tridata = triproxy.getAllData()

    if (verbose): print("Number of nodes imported: ", nodeproxy.getSize())
    if (verbose): print("Number of tetrahedrons imported: ", tetproxy.getSize())
    if (verbose): print("Number of triangles imported: ", triproxy.getSize())

    if (verbose): print("creating Tetmesh object in STEPS...")
    mesh = stetmesh.Tetmesh(nodedata, tetdata, tridata)
    if (verbose): print("Tetmesh object created.")

    return mesh, nodeproxy, tetproxy,triproxy

#############################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#############################################################################################

def importVTK(filename, scale, verbose = False):
    """
    Read a VTK-formated mesh file, return the created steps.geom.Tetmesh object,
    the element mapping for nodes, tetraedrons and triangles.

    PARAMETERS:

    * filename: the VTK filename (or path) including any suffix.
    * scale: LENGTH scale from the importing mesh to real geometry. e.g. a radius
      of 10 in the importing file to a radius of 1 micron in STEPS, scale is 1e-7.

    RETURNS:

    mesh, nodeproxy, tetproxy, triproxy

    * mesh: The STEPS TetMesh object
    * nodeproxy: Element Map for nodes
    * tetproxy: Element Map for tetrahedrons
    * triproxy: Element Map for triangles

   IMPORTANT NOTICE:

   User is recommanded to save the tetmesh objects using the saveTetmesh() method
   and recreate the objects from the saved files, instead of creating the objects
   via the importing functions each time, if the tetmesh objects are intented to
   be used in multiple simulations. Since the importing functions require a massive
   amount of time to create the Tetmesh object, comparing to the loadTetmesh() method.

    """
    if (verbose): print("Reading VTK file...")
    #btime = time.time()

    vtkfile = open(filename, 'r')

    nodeproxy = ElementProxy('node', 3)
    tetproxy = ElementProxy('tet', 4)
    triproxy = ElementProxy('tri', 3)
    elements = []

    warningprinted = False

    line = vtkfile.readline()

    file_it = iter(vtkfile)
    for line in file_it:
        if line.strip() == "":
            continue
        tokens = line.split()
        if tokens[0] == "DATASET" and tokens[1] != "UNSTRUCTURED_GRID":
            print("Error: not an unstructured grid!")
            break
        elif tokens[0] == "POINTS":
            i = 0
            n_points = int(tokens[1])
            while i < n_points:
                line = next(file_it)
                tokens = line.split()
                for j in range(len(tokens) // 3):
                    nodeproxy.insert(i, [float(tokens[j])*scale, float(tokens[j + 1])*scale, float(tokens[j + 2])*scale])
                    i += 1
        elif tokens[0] == "CELLS":
            for i in range(int(tokens[1])):
                line = next(file_it)
                tokens = line.split()
                elements.append(tokens)
        elif tokens[0] == "CELL_TYPES":
            for i in range(int(tokens[1])):
                line = next(file_it)
                element = elements[i]
                if int(line) == 5 and int(element[0]) == 3:
                    triproxy.insert(i, [int(element[1]), int(element[2]), int(element[3])])
                elif int(line) == 10 and int(element[0]) == 4:
                    tetproxy.insert(i, [int(element[1]), int(element[2]), int(element[3]), int(element[4])])
                elif not warningprinted:
                    print("Warning: at least one element is neither a triangle nor a tetrahedron")
                    warningprinted = True                 
            break
    vtkfile.close()
    nodedata = nodeproxy.getAllData()
    tetdata = tetproxy.getAllData()
    tridata = triproxy.getAllData()

    if (verbose): print("creating Tetmesh object in STEPS...")
    mesh = stetmesh.Tetmesh(nodedata, tetdata, tridata)
    if (verbose): print("Tetmesh object created.")
    return mesh, nodeproxy, tetproxy, triproxy

#############################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#############################################################################################

def saveMesh(pathname, tetmesh):

    """
    Save a STEPS Tetmesh in an XML file

   This file stores the basic information about the mesh which tends to be
   common information for any software that supports tetrahedral meshes.

   * NODES are stored by cartesian coordinates.
   * TRIANGLES are stored by the indices of their 3 nodes.
   * TETRAHEDRONS are stored by the indices of their 4 nodes.

   The XML file also stores infomation about any compartments or
   patches created in STEPS (class steps.geom.TmComp steps.geom.TmPatch
   respectively).

   * COMPARTMENT(S) are stored by:

     * their string identification.
     * a list of any volume systems added to the compartment at time of saving.
     * a list of tetrahedrons belonging to the compartment

   * PATCH(ES) are stored by:

     * their string identification.
     * a list of any surface systems added to the patch at time of saving.
     * the inner compartment id.
     * the outer compartment id (if it exists).
     * a list of trianlges belonging to this patch.

    PARAMETERS:

    * pathname:

      the root of the path to store the files.

      e.g. 'meshes/spine1' will save data in /meshes/spine1.xml

    * tetmesh:

      A valid STEPS Tetmesh object (of class steps.geom.Tetmesh).
      This mesh can be made in a variety of ways, e.g. to save a mesh loaded from tetgen::

          >>> import meshio
          >>> ### Use cubit script to create steps.geom.Tetmesh object from tetgen output file ###
          >>> mymesh = meshio.readTetgen(tetgenfilename)
          >>> ### Save this mesh in XML (and ASCII) format for quick-loading in future ###
          >>> meshio.saveMesh('/meshes/spine1', mymesh[0])

    """

    # Perform a basic test on the object itself
    if not isinstance(tetmesh, stetmesh.Tetmesh):
        raise TypeError(f"Expected a steps.API_1.geom.Tetmesh object, got {tetmesh} instead.")

    # Following will raise IOError if pathname not a valid directory
    xmlfile = open(pathname+'.xml', 'w')

    xmlfile.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
    xmlfile.write('<tetmesh>\n')


    ### Write the node information

    nverts = tetmesh.countVertices()

    xmlfile.write('\t<nodes size="'+str(nverts)+'">\n')

    for node in range(nverts):
        # Write indices and coordinates to xml file
        xmlfile.write('\t\t<node idx = "' + str(node) + '">\n')
        coords = str(tetmesh.getVertex(node))[1:-1]
        xmlfile.write('\t\t\t<coords>'+ coords + '</coords>\n')
        xmlfile.write('\t\t</node>\n')
    xmlfile.write('\t</nodes>\n')


    ### Write the triangle information

    ntris = tetmesh.countTris()
    xmlfile.write('\t<triangles size="' +str(ntris)+'">\n')

    for tri in range(ntris):
        # Write indices and nodes to xml file
        xmlfile.write('\t\t<tri idx="' + str(tri) + '">\n')
        nodes = str(tetmesh.getTri(tri))[1:-1]
        xmlfile.write('\t\t\t<nodes>' + nodes + '</nodes>\n')
        xmlfile.write('\t\t</tri>\n')
    xmlfile.write('\t</triangles>\n')


    ### Write the tetrahedron information

    ntets = tetmesh.countTets()
    xmlfile.write('\t<tetrahedrons size="' +str(ntets)+'">\n')

    for tet in range(ntets):
        # Write indices and nodes to xml file
        xmlfile.write ('\t\t<tet idx="' + str(tet) + '">\n')
        nodes = str(tetmesh.getTet(tet))[1:-1]
        xmlfile.write('\t\t\t<nodes>' + nodes + '</nodes>\n')
        xmlfile.write('\t\t</tet>\n')
        # Write the volume, barycenter, trineighbours, tetneighbours to xml file
    xmlfile.write('\t</tetrahedrons>\n')


    ### Write the comp and patch information.
    # TODO: Changes need to be made to steps code to make it easier to get the tet
    # and tri members. Currently the mesh returns base pointer (Comp or Patch)
    # which cannot be used to find the indices.

    comps = tetmesh.getAllComps()
    ncomps = comps.__len__()
    xmlfile.write('\t<compartments size="' + str(ncomps) + '">\n')

    if (ncomps > 0) :
        ids = []
        vsys=[]
        tets = []
        for c in range(ncomps):
            ids.append(comps[c].getID())
            vsys.append(comps[c].getVolsys())
            tets.append([])
        assert(tets.__len__() == ncomps)
        # Only choice right now is to loop over all tets and compare comp to tet id
        for tet in range(ntets):
            comptemp = tetmesh.getTetComp(tet)
            if not comptemp: continue
            idtemp = comptemp.getID()
            for c in range(ncomps):
                if idtemp == ids[c]:
                    tets[c].append(tet)
                    break
        # Now we have the tet members of each comp, we can write this to xml
        for c in range(ncomps):
            xmlfile.write('\t\t<comp idx="' +str(c) + '">\n')
            xmlfile.write('\t\t\t<id>' + ids[c] + '</id>\n')
            xmlfile.write('\t\t\t<volsys>')
            for v in vsys[c]:
                if (v != vsys[c][0]) : xmlfile.write(',')
                xmlfile.write(v)
            xmlfile.write('</volsys>\n')
            xmlfile.write('\t\t\t<tets>')
            for t in tets[c]:
                if (t != tets[c][0]) : xmlfile.write(',')
                xmlfile.write(str(t))
            xmlfile.write('</tets>\n')
            xmlfile.write('\t\t</comp>\n')
    xmlfile.write('\t</compartments>\n')

    patches = tetmesh.getAllPatches()
    npatches = patches.__len__()
    xmlfile.write('\t<patches size="' + str(npatches) + '">\n')

    if (npatches > 0) :
        ids = []
        ssys = []
        icomp=[]
        ocomp=[]
        tris = []
        for p in range(npatches):
            ids.append(patches[p].getID())
            ssys.append(patches[p].getSurfsys())
            icomp.append(patches[p].getIComp().getID())
            if (not patches[p].getOComp()): ocomp.append('null')
            else : ocomp.append(patches[p].getOComp().getID())
            tris.append([])
        assert(ids.__len__() == ssys.__len__() == icomp.__len__() == ocomp.__len__() == tris.__len__() == npatches)

        for tri in range(ntris):
            patchtemp = tetmesh.getTriPatch(tri)
            if not patchtemp: continue
            idtemp = patchtemp.getID()
            for p in range(npatches):
                if idtemp == ids[p]:
                    tris[p].append(tri)
                    break
        # Write all to xml
        for p in range(npatches):
            xmlfile.write('\t\t<patch idx = "'+str(p) + '">\n')
            xmlfile.write('\t\t\t<id>' + ids[p] + '</id>\n')
            xmlfile.write('\t\t\t<surfsys>')
            for s in ssys[p]:
                if (s != ssys[p][0]) : xmlfile.write(',')
                xmlfile.write(s)
            xmlfile.write('</surfsys>\n')
            xmlfile.write('\t\t\t<icomp>' + icomp[p] + '</icomp>\n')
            xmlfile.write('\t\t\t<ocomp>' + ocomp[p] + '</ocomp>\n')
            xmlfile.write('\t\t\t<tris>')
            for t in tris[p]:
                if (t != tris[p][0] ) : xmlfile.write(',')
                xmlfile.write(str(t))
            xmlfile.write('</tris>\n')
            xmlfile.write('\t\t</patch>\n')
    xmlfile.write('\t</patches>\n')

    # add write ROI data
    roi_ids = tetmesh.getAllROINames()
    xmlfile.write('\t<ROI_records size="' + str(len(roi_ids)) + '">\n')
    for n in roi_ids:
        type = tetmesh.getROIType(n)
        xmlfile.write('\t\t<ROI id="'+ n + '" type="' + str(type) + '">\n')
        indices = tetmesh.getROIData(n)
        xmlfile.write('\t\t\t<indices>')
        for i in indices:
            if (i != indices[0] ) : xmlfile.write(',')
            xmlfile.write(str(i))
        xmlfile.write('</indices>\n')
        xmlfile.write('\t\t</ROI>\n')
    xmlfile.write('\t</ROI_records>\n')
    xmlfile.write('</tetmesh>')

    xmlfile.close()

#############################################################################################

def loadMesh(pathname, scale=1, strict=False):
    """
    Load a mesh in STEPS from the XML file.

    PARAMETERS:

    * pathname: the root of the path where the file(s) are stored.
      e.g. with 'meshes/spine1' this function will look for the file /meshes/spine1.xml

    * scale: optionally rescale the mesh on loading by given factor.

    * strict: apply strict(-er) checking to the input XML.

    RETURNS: A tuple (mesh, comps, patches)

    * mesh
      The STEPS Tetmesh object (steps.geom.Tetmesh)
    * comps
      A list of any compartment objects (steps.geom.TmComp) from XML file
    * patches
      A list of any patches objects (steps.geom.TmPatch) from XML file
    """

    # More robust/relaxed attribute parsing:
    element_rx=re.compile(r'\s*<(\/?\w+)((?:\s+\w+\s*=\s*(?:\w+|"[^"]*"|\'[^\']*\'))*)\s*/?\s*>')
    attrib_rx=re.compile(r'\s*(\w+)\s*=\s*(?:(\w+)|"([^"]*)"|\'([^\']*)\')')

    def parse_xml_element(s):
        (tag,attrstr) = element_rx.match(s).groups()
        attrs = {}
        for a,u,d,q in attrib_rx.findall(attrstr):
            attrs[a] = u+d+q
        return (tag,attrs)

    # Try to open the XML file. An error will be raisen if it doesn't exist
    xmlfile = open(pathname+'.xml', 'r')

    # Perform a basic check to see if we have the expected kind of file which has not been altered.
    info = xmlfile.readline()
    if(xmlfile.readline().rstrip() != '<tetmesh>'):
        print('XML file is not a recognised STEPS mesh file')
        return

    # Collect basic node information and perform some checks on the data read from XML file
    nodeinfo = xmlfile.readline().strip()
    (tag,attrs) = parse_xml_element(nodeinfo)
    assert(tag == 'nodes')
    assert('size' in attrs)
    nnodes = int(attrs['size'])

    nodes_out = [0.0]*(nnodes*3)
    for i in range(nnodes):
        idxtemp = xmlfile.readline().strip()
        assert(int(idxtemp[13:-2]) == i)
        coordtemp = xmlfile.readline().strip()
        assert(coordtemp[:8] == '<coords>' and coordtemp[-9:] == '</coords>')
        coordtemp = coordtemp[8:-9].split(', ')
        nodes_out[i*3], nodes_out[(i*3)+1], nodes_out[(i*3)+2] = float(coordtemp[0]), float(coordtemp[1]), float(coordtemp[2])
        assert(xmlfile.readline().strip() == '</node>')

    assert(xmlfile.readline().strip() == '</nodes>')

    # Now read triangle information from xml file and text file if we have it
    triinfo = xmlfile.readline().strip()
    (tag,attrs) = parse_xml_element(triinfo)
    assert(tag == 'triangles')
    assert('size' in attrs)
    ntris = int(attrs['size'])

    tris_out = [0]*(ntris*3) # numpy.zeros(ntris*3, dtype = 'int')

    for i in range(ntris):
        idxtemp = xmlfile.readline().strip()
        if strict:
            (tag,attrs) = parse_xml_element(idxtemp)
            assert(tag == 'tri')
            assert(int(attrs['idx']) == i)

        nodetemp = xmlfile.readline().strip()
        assert (nodetemp[:7] == '<nodes>' and nodetemp[-8:] == '</nodes>')
        nodetemp = nodetemp[7:-8].split(', ')
        tris_out[i*3], tris_out[(i*3)+1], tris_out[(i*3)+2] = int(nodetemp[0]), int(nodetemp[1]), int(nodetemp[2])
        assert(xmlfile.readline().strip() == '</tri>')
        # Now read the text file, if it exists, and get the extra information

    assert(xmlfile.readline().strip() == '</triangles>')

    # Now read tet information from xml file and text file if we have it
    tetinfo = xmlfile.readline().strip()
    (tag,attrs) = parse_xml_element(tetinfo)
    assert(tag == 'tetrahedrons')
    assert('size' in attrs)
    ntets = int(attrs['size'])

    tets_out = [0]*(ntets*4)	# numpy.zeros(ntets*4, dtype = 'int')
    for i in range(ntets):
        idxtemp = xmlfile.readline().strip()
        if strict:
            (tag,attrs) = parse_xml_element(idxtemp)
            assert(tag == 'tet')
            assert(int(attrs['idx']) == i)

        nodetemp = xmlfile.readline().strip()
        assert (nodetemp[:7] == '<nodes>' and nodetemp[-8:] == '</nodes>')
        nodetemp = nodetemp[7:-8].split(', ')
        tets_out[i*4], tets_out[(i*4)+1], tets_out[(i*4)+2], tets_out[(i*4)+3] = int(nodetemp[0]), int(nodetemp[1]), int(nodetemp[2]), int(nodetemp[3])
        assert(xmlfile.readline().strip() == '</tet>')

    assert(xmlfile.readline().strip() == '</tetrahedrons>')

    # Rescale coordinates if requested.
    if scale!=1:
        for i in range(len(nodes_out)):
            nodes_out[i] *= scale

    # We have all the information now. Time to make the Tetmesh object.
    mesh = stetmesh.Tetmesh(nodes_out, tets_out, tris_out)

    # Now fetch any comp and patch information from XML file
    compinfo = xmlfile.readline().strip()
    (tag,attrs) = parse_xml_element(compinfo)
    assert(tag == 'compartments')
    assert('size' in attrs)
    ncomps = int(attrs['size'])
    comps_out = []
    for i in range(ncomps):
        idxtemp = xmlfile.readline().strip()
        if strict:
            (tag,attrs) = parse_xml_element(idxtemp)
            assert(tag == 'comp')
            assert(int(attrs['idx']) == i)

        idtemp = xmlfile.readline().strip()
        assert(idtemp[:4] == '<id>' and idtemp[-5:] == '</id>')
        idtemp = idtemp[4:-5]
        volsystemp = xmlfile.readline().strip()
        assert(volsystemp[:8] == '<volsys>' and volsystemp[-9:] == '</volsys>')
        volsystemp = volsystemp[8:-9].split(',')
        if (volsystemp[0] == '') : volsystemp = []
        tettemp = xmlfile.readline().strip()
        assert(tettemp[:6] == '<tets>' and tettemp[-7:] == '</tets>')
        tettemp = tettemp[6:-7].split(',')
        nctets = tettemp.__len__()
        ctets = [0]*nctets		# numpy.zeros(nctets, dtype = 'int')
        for ct in range(nctets): ctets[ct] = int(tettemp[ct])
        c_out = stetmesh.TmComp(idtemp, mesh, ctets)
        for v in volsystemp: c_out.addVolsys(v)
        comps_out.append(c_out)
        assert(xmlfile.readline().strip() == '</comp>')
    assert(xmlfile.readline().strip() == '</compartments>')

    # Retrieve patch info
    patchinfo = xmlfile.readline().strip()
    (tag,attrs) = parse_xml_element(patchinfo)
    assert(tag == 'patches')
    assert('size' in attrs)
    npatches = int(attrs['size'])
    patches_out = []
    for i in range(npatches):
        idxtemp = xmlfile.readline().strip()
        if strict:
            (tag,attrs) = parse_xml_element(idxtemp)
            assert(tag == 'patch')
            assert(int(attrs['idx']) == i)

        idtemp = xmlfile.readline().strip()
        assert(idtemp[:4] == '<id>' and idtemp[-5:] == '</id>')
        idtemp = idtemp[4:-5]
        surfsystemp = xmlfile.readline().strip()
        assert(surfsystemp[:9] == '<surfsys>' and surfsystemp[-10:] == '</surfsys>')
        surfsystemp = surfsystemp[9:-10].split(',')
        if (surfsystemp[0] == '') : surfsystemp = []
        icomptemp = xmlfile.readline().strip()
        assert(icomptemp[:7] == '<icomp>' and icomptemp[-8:] == '</icomp>')
        icomptemp = icomptemp[7:-8]
        ocomptemp = xmlfile.readline().strip()
        assert(ocomptemp[:7] == '<ocomp>' and ocomptemp[-8:] == '</ocomp>')
        ocomptemp = ocomptemp[7:-8]
        tritemp = xmlfile.readline().strip()
        assert(tritemp[:6] == '<tris>' and tritemp[-7:] == '</tris>')
        tritemp = tritemp[6:-7].split(',')
        nptris = tritemp.__len__()
        ptris = [0]*nptris		# numpy.zeros(nptris, dtype='int')
        for pt in range(nptris): ptris[pt] = int(tritemp[pt])
        if (ocomptemp != 'null'): p_out = stetmesh.TmPatch(idtemp, mesh, ptris, mesh.getComp(icomptemp), mesh.getComp(ocomptemp))
        else :  p_out = stetmesh.TmPatch(idtemp, mesh, ptris, mesh.getComp(icomptemp))
        for s in surfsystemp: p_out.addSurfsys(s)
        patches_out.append(p_out)
        assert(xmlfile.readline().strip() == '</patch>')
    assert(xmlfile.readline().strip() == '</patches>')

    # Retrieve ROI info
    info = xmlfile.readline().strip()
    (tag,attrs) = parse_xml_element(info)
    if tag == 'ROI_records':
        ROIinfo = info
        assert('size' in attrs)
        nROI = int(attrs['size'])

        for r in range(nROI):
            ROItemp = xmlfile.readline().strip()
            (tag,attrs) = parse_xml_element(ROItemp)
            assert(tag == 'ROI')
            assert('id' in attrs)
            assert('type' in attrs)
            id = attrs['id']
            type = int(attrs['type'])

            idxtemp = xmlfile.readline().strip()
            assert(idxtemp.find("<indices>") != -1)
            assert(idxtemp.find("</indices>") != -1)
            idxtemp = idxtemp.replace("<indices>", "").replace("</indices>", "").split(",")
            indices = [int(i) for i in idxtemp]
            assert(xmlfile.readline().strip() == '</ROI>')
            mesh.addROI(id, type, indices)
        assert(xmlfile.readline().strip() == '</ROI_records>')
        info = xmlfile.readline().strip()

    xmlfile.close()

    assert(info == '</tetmesh>')
    return (mesh,comps_out,patches_out)
