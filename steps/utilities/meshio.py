# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
# Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
# See the file AUTHORS for details.
#
# This file is part of STEPS.
#
# STEPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STEPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#  Last Changed Rev:  $Rev: 314 $
#  Last Changed Date: $Date: 2010-03-28 08:47:50 +0900 (Sun, 28 Mar 2010) $
#  Last Changed By:   $Author: wchen $

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
import steps.geom as stetmesh
import os.path as opath

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

    def getGroups(self):
        """
        Return all groups stored in the Element Map object.

        Example:
        
        The following group dictionary contains two groups::
            
            "Group1" contains element ids 1,2,6,4,8
            "Group2" contains element ids 3,1,2,1

            groups = {"Group1":[1,2,6,4,8], "Group2":[3,1,2,1]}

        Note that element 1 appears twice in "Group2".

        """

        return self.groups
        
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

def importTetGen(pathroot):
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
        node[0] = float(tokens[1])
        node[1] = float(tokens[2])
        node[2] = float(tokens[3])
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
    # Close the file.
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
        # Close the file.
        facefile.close()

    print("Read TetGen files succesfully")    
    
    nodedata = nodeproxy.getAllData()
    tetdata = tetproxy.getAllData() 
    tridata = triproxy.getAllData() 

    print("creating Tetmesh object in STEPS...")
    mesh = stetmesh.Tetmesh(nodedata, tetdata, tridata)
    print("Tetmesh object created.")
    return mesh, nodeproxy, tetproxy, triproxy
    

#############################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#############################################################################################

def importAbaqus(filename, scale):
    """ 
    Read a ABAQUS-formated mesh file, return the created steps.geom.Tetmesh object, 
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
    print("Reading Abaqus file...")
    btime = time.time()
    
    abaqusfile = open(filename, 'r')
    
    nodeproxy = ElementProxy('node', 3)
    tetproxy = ElementProxy('tet', 4)
    triproxy = ElementProxy('tri', 3)

    line = abaqusfile.readline()
    # check we have the right kind of CUBIT output file here
    assert(line == '*HEADING\n')
    
    currmap = None
    # read line one by one until the end
    while(line):
        # status check
        if (line.find('**', 0, 2) == 0):
            pass
        elif (line.find('*NODE', 0, 5) == 0):
            if (currmap != None):
                currmap.blockEnd()
            print('Found *NODE section, start reading nodes.')
            lineset = line.split('=') 
            if (len(lineset) == 1):
                nodeproxy.blockBegin('AllNodes')
            else:
                nodeproxy.blockBegin(lineset[1].strip('\n'))
            currmap = nodeproxy
        elif (line.find('*ELEMENT', 0, 8) == 0):
            if (currmap != None):
                currmap.blockEnd()
            print('Found *Element section, start reading elements.')
            lineset = line.split(',') 
            elementtype = lineset[1].split('=')[1]
            blockname = 'AllElements'
            if (len(lineset) == 3):
                blockname = lineset[2].split('=')[1].strip('\n')
            if (elementtype == 'C3D4'):
                currmap = tetproxy
                tetproxy.blockBegin(blockname)
            elif (elementtype == 'STRI3'):
                currmap = triproxy
                triproxy.blockBegin(blockname)
        elif (line.find('**', 0, 2) == -1 and line[0] == '*'):
            if (currmap != None):
                currmap.blockEnd()
                currmap = None
        else:
            if (currmap == nodeproxy):
                linesec = line.replace(',', ' ')
                linesec = linesec.split() 
                nodeid = int(linesec[0])
                node = [0.0, 0.0, 0.0]
                # set to the right scale
                node[0] = float(linesec[1])*scale
                node[1] = float(linesec[2])*scale
                node[2] = float(linesec[3])*scale 
                currmap.insert(nodeid, node)
            elif (currmap == tetproxy): 
                linesec = line.replace(',', ' ')
                linesec = linesec.split() 
                nodeid = int(linesec[0])
                node = [0,0,0,0]
                # set to the right scale
                node[0] = nodeproxy.getSTEPSID(int(linesec[1]))
                node[1] = nodeproxy.getSTEPSID(int(linesec[2]))
                node[2] = nodeproxy.getSTEPSID(int(linesec[3])) 
                node[3] = nodeproxy.getSTEPSID(int(linesec[4])) 
                currmap.insert(nodeid, node)
            elif (currmap == triproxy):
                linesec = line.replace(',', ' ')
                linesec = linesec.split() 
                nodeid = int(linesec[0])
                node = [0,0,0]
                # set to the right scale
                node[0] = nodeproxy.getSTEPSID(int(linesec[1]))
                node[1] = nodeproxy.getSTEPSID(int(linesec[2]))
                node[2] = nodeproxy.getSTEPSID(int(linesec[3])) 
                currmap.insert(nodeid, node)           

        line = abaqusfile.readline()

    if (currmap != None):
        currmap.blockEnd()

    abaqusfile.close()
    print("Read Abaqus file succesfully")    
    
    nodedata = nodeproxy.getAllData()
    tetdata = tetproxy.getAllData() 
    tridata = triproxy.getAllData() 

    print("creating Tetmesh object in STEPS...")
    mesh = stetmesh.Tetmesh(nodedata, tetdata, tridata)
    print("Tetmesh object created.")
    
    return mesh, nodeproxy, tetproxy,triproxy

#############################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#############################################################################################

def saveMesh(pathname, tetmesh):
    
    """ 
    Save a STEPS Tetmesh in two separate files 
    
    #. An XML file.
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
        
    #. An ASCII file storing information important to STEPS internally. This
       information must be found by STEPS once from the basic mesh infomation and
       is vital for simulations in STEPS. This can take a significant amount of 
       time for larger meshes, so storing this information in this way can drastically
       reduce future mesh loading times. The information stored is:
      
       * each triangle's area.
       * each triangle's normal.
       * each triangle's two (or one for surface tris) tetrahedron neighbours.
       * each tetrahedron's volume.
       * each tetrahedron's barycenter.
       * each tetrahedron's four neighbouring triangles.
       * each tetrahedron's four (or fewer for surface tets) tetrahedron neighbours.
    
    PARAMETERS:
    
    * pathname: 
      
      the root of the path to store the files. 
      
      e.g. 'meshes/spine1' will save data in /meshes/spine1.xml and /meshes/spine1.txt
      
    * tetmesh:
    
      A valid STEPS Tetmesh object (of class steps.geom.Tetmesh). 
      This mesh can be made in a variety of ways, e.g. to save a mesh loaded from tetgen::
      
          >>> import meshio
          >>> ### Use cubit script to create steps.geom.Tetmesh object from tetgen output file ###
          >>> mymesh = meshio.readTetgen(tetgenfilename)
          >>> ### Save this mesh in XML (and ASCII) format for quick-loading in future ###
          >>> meshio.saveMesh('/meshes/spine1', mymesh[0])
        
    """
    
    # Performa a basic test on the object itself
    if (tetmesh.__str__()[1:19] != 'steps.geom.Tetmesh'):
        print("2nd parameter not a valid steps.geom.Tetmesh object.")
        return 0
    
    # Following will throw IOError if pathname not a valid directory
    xmlfile = open(pathname+'.xml', 'w')
    textfile = open(pathname+'.txt', 'w')
    
    xmlfile.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
    xmlfile.write('<tetmesh>\n')
    
    
    ### Write the node information  
    
    nverts = tetmesh.countVertices()
    
    xmlfile.write('\t<nodes size = "'+str(nverts)+'">\n')
    textfile.write(str(nverts)+'\n\n')
    
    for node in range(nverts):
        # Write indices and coordinates to xml file
        xmlfile.write('\t\t<node idx = "' + str(node) + '">\n')
        coords = str(tetmesh.getVertex(node)).strip(')').strip('(')
        xmlfile.write('\t\t\t<coords>'+ coords + '</coords>\n')
        xmlfile.write('\t\t</node>\n')
    xmlfile.write('\t</nodes>\n')
    
    
    ### Write the triangle information
    
    ntris = tetmesh.countTris()
    xmlfile.write('\t<triangles size = "' +str(ntris)+'">\n')
    textfile.write(str(ntris)+'\n')
    
    for tri in range(ntris):
        # Write indices and nodes to xml file
        xmlfile.write('\t\t<tri idx = "' + str(tri) + '">\n')
        nodes = str(tetmesh.getTri(tri)).strip(')').strip('(')
        xmlfile.write('\t\t\t<nodes>' + nodes + '</nodes>\n')
        xmlfile.write('\t\t</tri>\n')
        # Write the area, normal, tetneighbours to text file
        textfile.write(str(tetmesh.getTriArea(tri)) + " ")
        norm = tetmesh.getTriNorm(tri)
        textfile.write(str(norm[0]) + " " + str(norm[1]) + " " + str(norm[2]) + " ")
        tetneighb = tetmesh.getTriTetNeighb(tri)
        textfile.write(str(tetneighb[0]) + " " + str(tetneighb[1]) + " ")
        textfile.write('\n')
    xmlfile.write('\t</triangles>\n')
    
    
    ### Write the tetrahedron information
    
    ntets = tetmesh.countTets()
    xmlfile.write('\t<tetrahedrons size = "' +str(ntets)+'">\n')
    textfile.write('\n' +str(ntets)+'\n')
    
    for tet in range(ntets):
        # Write indices and nodes to xml file
        xmlfile.write ('\t\t<tet idx = "' + str(tet) + '">\n')
        nodes = str(tetmesh.getTet(tet)).strip(')').strip('(')
        xmlfile.write('\t\t\t<nodes>' + nodes + '</nodes>\n')
        xmlfile.write('\t\t</tet>\n')
        # Write the volume, barycenter, trineighbours, tetneighbours to xml file
        textfile.write(str(tetmesh.getTetVol(tet)) + " ")
        baryc = tetmesh.getTetBarycenter(tet)
        textfile.write(str(baryc[0]) + " " + str(baryc[1]) + " " + str(baryc[2]) + " ")
        trineighb = tetmesh.getTetTriNeighb(tet)
        textfile.write(str(trineighb[0])+" "+str(trineighb[1])+" "+str(trineighb[2])+" "+str(trineighb[3])+" ")
        tetneighb = tetmesh.getTetTetNeighb(tet)
        textfile.write(str(tetneighb[0])+" "+str(tetneighb[1])+" "+str(tetneighb[2])+" "+str(tetneighb[3])+" ")
        textfile.write('\n')
    xmlfile.write('\t</tetrahedrons>\n')
    
    
    ### Write the comp and patch information. 
    # TODO: Changes need to be made to steps code to make it easier to get the tet
    # and tri members. Currently the mesh returns base pointer (Comp or Patch) 
    # which cannot be used to find the indices. 
    
    comps = tetmesh.getAllComps()
    ncomps = comps.__len__()
    xmlfile.write('\t<compartments size = "' + str(ncomps) + '">\n')
    
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
            xmlfile.write('\t\t<comp idx = "' +str(c) + '">\n')
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
    xmlfile.write('\t<patches size = "' + str(npatches) + '">\n')
    
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
    
    
    xmlfile.write('</tetmesh>')
    
    xmlfile.close()
    textfile.close()

#############################################################################################

def loadMesh(pathname):
    """ 
    Load a mesh in STEPS from the XML (and ASCII) file. This will
    work with just the XML file, but this is akin to creating the mesh in STEPS
    from scratch and really negates the use of storing the mesh infomation at all.
    For maximum benefit the XML file should be accompanied by the ASCII file, which
    contains all the internal information.
     
    PARAMETERS:
    
    * pathname: the root of the path where the file(s) are stored. 
      e.g. with 'meshes/spine1' this function will look for files /meshes/spine1.xml 
      and /meshes/spine1.txt
    
    RETURNS: A tuple (mesh, comps, patches)
    
    * mesh 
      The STEPS Tetmesh object (steps.geom.Tetmesh)
    * comps
      A list of any compartment objects (steps.geom.TmComp) from XML file 
    * patches
      A list of any patches objects (steps.geom.TmPatch) from XML file 
    """


    # Try to open the XML file. An error will be thrown if it doesn't exist
    xmlfile = open(pathname+'.xml', 'r')

    # Try to open the text file. A warning will be shown and a flag set if it doesn't exist 
    havetxt = True
    try :
        textfile = open(pathname+'.txt', 'r')
    except:
        havetxt = False
    if (havetxt == False) : print("WARNING: text file not found. Will construct mesh from information in XML file only.")
    
    # Perform a basic check to see if we have the expected kind of file which has not been altered.
    info = xmlfile.readline()
    if(xmlfile.readline().rstrip() != '<tetmesh>'):
        print('XML file is not a recognised STEPS mesh file')
        return
    
    # Collect basic node information and perform some checks on the data read from XML file
    nodeinfo = xmlfile.readline().strip()
    assert(nodeinfo.__len__() > 17)
    assert(nodeinfo[-2:] == '">')
    nnodes = int(nodeinfo[15:-2])
    if (havetxt): assert (nnodes == int(textfile.readline()))
    
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
    assert(triinfo.__len__() > 21)
    assert(triinfo[-2:] == '">')
    ntris = int(triinfo[19:-2])
    if (havetxt) : 
        textfile.readline()
        assert (ntris == int(textfile.readline()))

    tris_out = [0]*(ntris*3) # numpy.zeros(ntris*3, dtype = 'int')
    if (havetxt) :
        triareas_out = [0.0]*ntris 		# numpy.zeros(ntris)
        trinorms_out = [0.0]*(ntris*3)	# numpy.zeros(ntris*3)
        tritetns_out = [0]*(ntris*2)	# numpy.zeros(ntris*2, dtype = 'int')
    
    for i in range(ntris): 
        idxtemp = xmlfile.readline().strip()
        assert(int(idxtemp[12:-2]) == i)
        nodetemp = xmlfile.readline().strip()
        assert (nodetemp[:7] == '<nodes>' and nodetemp[-8:] == '</nodes>')
        nodetemp = nodetemp[7:-8].split(', ')
        tris_out[i*3], tris_out[(i*3)+1], tris_out[(i*3)+2] = int(nodetemp[0]), int(nodetemp[1]), int(nodetemp[2])
        assert(xmlfile.readline().strip() == '</tri>')
        # Now read the text file, if it exists, and get the extra information
        if (havetxt):
            line = textfile.readline().rstrip().split(" ")
            assert(line.__len__() == 6)
            triareas_out[i], trinorms_out[i*3], trinorms_out[(i*3)+1], trinorms_out[(i*3)+2] = float(line[0]), float(line[1]), float(line[2]), float(line[3])
            tritetns_out[i*2], tritetns_out[(i*2)+1] = int(line[4]), int(line[5])       
    assert(xmlfile.readline().strip() == '</triangles>')
    
    # Now read tet information from xml file and text file if we have it
    tetinfo = xmlfile.readline().strip()
    assert(tetinfo.__len__() > 24)
    assert(tetinfo[-2:] == '">')
    ntets = int(tetinfo[22:-2])
    if (havetxt) :
        textfile.readline()
        assert(ntets == int(textfile.readline()))   
    
    tets_out = [0]*(ntets*4)	# numpy.zeros(ntets*4, dtype = 'int')
    if (havetxt): 
        tetvols_out =  [0.0]*ntets		# numpy.zeros(ntets)
        tetbarycs_out = [0.0]*(ntets*3)	# numpy.zeros(ntets*3)
        tettrins_out = 	[0]*(ntets*4)	# numpy.zeros(ntets*4, dtype = 'int')
        tettetns_out = 	[0]*(ntets*4)	# numpy.zeros(ntets*4, dtype = 'int')
    for i in range(ntets): 
        idxtemp = xmlfile.readline().strip()
        assert(int(idxtemp[12:-2]) == i)
        nodetemp = xmlfile.readline().strip()
        assert (nodetemp[:7] == '<nodes>' and nodetemp[-8:] == '</nodes>')      
        nodetemp = nodetemp[7:-8].split(', ')
        tets_out[i*4], tets_out[(i*4)+1], tets_out[(i*4)+2], tets_out[(i*4)+3] = int(nodetemp[0]), int(nodetemp[1]), int(nodetemp[2]), int(nodetemp[3])
        assert(xmlfile.readline().strip() == '</tet>')
        # Read the text file if we have it and get further information
        if (havetxt):
            line = textfile.readline().rstrip().split(" ")
            assert (line.__len__() == 12)
            tetvols_out[i], tetbarycs_out[i*3], tetbarycs_out[(i*3)+1], tetbarycs_out[(i*3)+2] = float(line[0]), float(line[1]), float(line[2]), float(line[3])
            tettrins_out[i*4], tettrins_out[(i*4)+1], tettrins_out[(i*4)+2], tettrins_out[(i*4)+3] = int(line[4]), int(line[5]), int(line[6]), int(line[7])
            tettetns_out[i*4], tettetns_out[(i*4)+1], tettetns_out[(i*4)+2], tettetns_out[(i*4)+3] = int(line[8]), int(line[9]), int(line[10]), int(line[11])
    assert(xmlfile.readline().strip() == '</tetrahedrons>')

    # We have all the information now. Time to make the Tetmesh object. New constructor keeps the order, which is:
    # nodes, tris, tri areas, tri normals, tri tet neighbours, tets, tet volumes, tet barycenters, tet tri neighbs, tet tet neighbs.
    mesh = stetmesh.Tetmesh(nodes_out, tris_out, triareas_out, trinorms_out, tritetns_out, tets_out, tetvols_out, tetbarycs_out, tettrins_out, tettetns_out)
    
    # Now fetch any comp and patch information from XML file
    compinfo = xmlfile.readline().strip()
    assert(compinfo.__len__() > 24)
    assert(compinfo[-2:] == '">')
    ncomps = int(compinfo[22:-2])
    comps_out = []
    for i in range(ncomps):
        idxtemp = xmlfile.readline().strip()
        assert(int(idxtemp[13:-2]) == i)
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
    assert(patchinfo.__len__() > 19)
    assert(patchinfo[-2:] == '">')
    npatches = int(patchinfo[17:-2])
    patches_out = []
    for i in range(npatches):
        idxtemp = xmlfile.readline().strip()
        assert(int(idxtemp[14:-2]) == i)
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

    
    return (mesh,comps_out,patches_out)
#############################################################################################

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
