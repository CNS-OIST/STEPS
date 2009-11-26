# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2009 Stefan Wils. All rights reserved.
#
# This file is part of STEPS.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
#
# $Id: meshio.py 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

""" Mesh saving and loading tool. 
Once a Tetmesh object has been created this file allows the user to save the 
data in two files: 
	An xml annotated file containing the nodes, triangles and tetrahedra.
	A text file containing further information needed by STEPS solvers found in
	Tetmesh object constructor. 
This is intened to drastically reduce mesh-loading time for large meshes 
(over ~100,000 voxels) by storing all data required by STEPS internally so this 
infomation does not have to be found each time by the Tetmesh object constructor. 
"""

import numpy
import steps.geom as stetmesh

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def saveXML(pathname, tetmesh):
    
	""" Save a STEPS Tetmesh in two separate files: 
	1) an XML file.
	This file stores the basic information about the mesh which tends to be
	common information for any software that supports tetrahedral meshes.
	- NODES are stored by cartesian coordinates.
	- TRIANGLES are stored by the indices of their 3 nodes.
	- TETRAHEDRONS are sotred by the indices of their 4 nodes.
	
	 The XML file also stores infomation about any compartments or 
	 patches created in STEPS (class steps.geom.TmComp steps.geom.TmPatch
	 respectively). 
	 - COMPARTMENT(S) are stored by: 
		their string identification.
		a list of any volume systems added to the compartment at time of saving.
		a list of tetrahedrons belonging to the compartment
	 - PATCH(ES) are stored by:
		their string identification.
		a list of any surface systems added to the patch at time of saving.
		the inner compartment id.
		the outer compartment id (if it exists).
		a list of trianlges belonging to this patch.
		
	2) An ASCII file storing information important to STEPS internally. This
	information must be found by STEPS once from the basic mesh infomation and
	is vital for simulations in STEPS. This can take a significant amount of 
	time for larger meshes, so storing this information in this way can drastically
	reduce future mesh loading times. The information stored is:
	- each triangle's area.
	- each triangle's normal.
	- each triangle's two (or one for surface tris) tetrahedron neighbours.
	- each tetrahedron's volume.
	- each tetrahedron's barycenter.
	- each tetrahedron's four neighbouring triangles.
	- each tetrahedron's four (or fewer for surface tets) tetrahedron neighbours.
	
	PARAMETERS:
		1) pathname: the root of the path to store the files. 
			e.g. 'meshes/spine1' will save data in /meshes/spine1.xml and /meshes/spine1.txt
		2) tetmesh: A valid STEPS Tetmesh object (of class steps.geom.Tetmesh). This mesh 
				can be made in a variety of ways, e.g. to save a CUBIT mesh: 
			>>> import cubit
			>>> ### Use cubit script to create steps.geom.Tetmesh object from CUBIT output file ###
			>>> mymesh = cubit.makeMesh(cubitfilename, 1e-6)
			>>> import meshio
			>>> ### Save this mesh in XML (and ASCII) format for quick-loading in future ###
			>>> meshio.saveXML('/meshes/spine1', mymesh)
		
	"""
	
	
	if (tetmesh.__str__()[1:19] != 'steps.geom.Tetmesh'):
		print "2nd parameter not a valid steps.geom.Tetmesh object."
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
	
	for node in xrange(nverts):
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
	
	for tri in xrange(ntris):
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
	
	for tet in xrange(ntets):
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
		for c in xrange(ncomps):
			ids.append(comps[c].getID())
			vsys.append(comps[c].getVolsys())
			tets.append([])
		assert(tets.__len__() == ncomps)
		# Only choice right now is to loop over all tets and compare comp to tet id
		for tet in xrange(ntets):
			comptemp = tetmesh.getTetComp(tet)
			if not comptemp: continue
			idtemp = comptemp.getID()
			for c in xrange(ncomps):
				if idtemp == ids[c]:
					tets[c].append(tet)
					break
		# Now we have the tet members of each comp, we can write this to xml
		for c in xrange(ncomps):
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
		for p in xrange(npatches):
			ids.append(patches[p].getID())
			ssys.append(patches[p].getSurfsys())
			icomp.append(patches[p].getIComp().getID())
			if (not patches[p].getOComp()): ocomp.append('null')
			else : ocomp.append(patches[p].getOComp().getID())
			tris.append([])
		assert(ids.__len__() == ssys.__len__() == icomp.__len__() == ocomp.__len__() == tris.__len__() == npatches)
		
		for tri in xrange(ntris):
			patchtemp = tetmesh.getTriPatch(tri)
			if not patchtemp: continue
			idtemp = patchtemp.getID()
			for p in xrange(npatches):
				if idtemp == ids[p]:
					tris[p].append(tri)
					break
		# Write all to xml
		for p in xrange(npatches):
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

###############################################################

def loadXML(pathname):
	""" Load a mesh in STEPS from the XML (and ASCII) file. This will
	work with just the XML file, but this is akin to creating the mesh in STEPS
	from scratch and really negates the use of storing the mesh infomation at all.
	For maximum benefit the XML file should be accompanied by the ASCII file, which
	contains all the internal information.
	 
	PARAMETERS:
		1) pathname: the root of the path where the file(s) are stored. 
			e.g. with 'meshes/spine1' this function will look for files /meshes/spine1.xml 
			and /meshes/spine1.txt
	
	RETURNS: A tuple (mesh, comps, patches)
		mesh 
			The STEPS Tetmesh object (steps.geom.Tetmesh)
		comps
			A list of any compartment objects (steps.geom.TmComp) from XML file 
		patches
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
	if (havetxt == False) : print "WARNING: text file not found. Will construct mesh from information in XML file only."
	
	# Perform a basic check to see if we have the expected kind of file which has not been altered.
	info = xmlfile.readline()
	if(xmlfile.readline().rstrip() != '<tetmesh>'):
		print 'XML file is not a recognised STEPS mesh file'
		return
	
	# Collect basic node information and perform some checks on the data read from XML file
	nodeinfo = xmlfile.readline().strip()
	assert(nodeinfo.__len__() > 17)
	assert(nodeinfo[-2:] == '">')
	nnodes = int(nodeinfo[15:-2])
	if (havetxt): assert (nnodes == int(textfile.readline()))
	
	nodes_out = numpy.zeros(nnodes*3)
	for i in xrange(nnodes):
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

	tris_out = numpy.zeros(ntris*3, dtype = 'int')
	if (havetxt) :
		triareas_out = numpy.zeros(ntris)
		trinorms_out = numpy.zeros(ntris*3)
		tritetns_out = numpy.zeros(ntris*2, dtype = 'int')
	
	for i in xrange(ntris):	
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
	
	tets_out = numpy.zeros(ntets*4, dtype = 'int')
	if (havetxt): 
		tetvols_out = numpy.zeros(ntets)
		tetbarycs_out = numpy.zeros(ntets*3)
		tettrins_out = numpy.zeros(ntets*4, dtype = 'int')
		tettetns_out = numpy.zeros(ntets*4, dtype = 'int')
	for i in xrange(ntets):	
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
	for i in xrange(ncomps):
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
		ctets = numpy.zeros(nctets, dtype = 'int')
		for ct in xrange(nctets): ctets[ct] = int(tettemp[ct])
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
	for i in xrange(npatches):
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
		ptris = numpy.zeros(nptris, dtype='int')
		for pt in xrange(nptris): ptris[pt] = int(tritemp[pt])
		if (ocomptemp != 'null'): p_out = stetmesh.TmPatch(idtemp, mesh, ptris, mesh.getComp(icomptemp), mesh.getComp(ocomptemp))
		else :  p_out = stetmesh.TmPatch(idtemp, mesh, ptris, mesh.getComp(icomptemp))
		for s in surfsystemp: p_out.addSurfsys(s)
		patches_out.append(p_out)
		assert(xmlfile.readline().strip() == '</patch>')
	assert(xmlfile.readline().strip() == '</patches>')

	
	return (mesh,comps_out,patches_out)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
