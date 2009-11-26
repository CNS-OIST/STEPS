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

""" Supoprt for CUBIT (http://cubit.sandia.gov/), a powerful geometry and 
mesh-generation toolkit (license required) which supports generation of
tetrahedral meshes that this file imports to STEPS. Output from CUBIT should be
saved in ABAQUS format as this is the only format this script currently 
supports. """
 

import time
import steps.geom as stetmesh					

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def makeMesh(filename, scale):

	""" Read a CUBIT-generated tetrahedral mesh saved in ABAQUS format and 
	return the created steps.geom.Tetmesh object.
	
	PARAMETERS:
		1) filename: the CUBIT exported filename (or path) including any suffix.
			Must be in abacus format.
		2) scale: LENGTH scale from the cubit mesh to real geometry. e.g. a radius 
		   of 10 in CUBIT to a radius of 1 micron in STEPS, scale is 1e-7.
	
	INDEXING:
	Within STEPS, indexing of nodes, triangles and tetrahedrons must start at 0,
	whereas indexing in CUBIT begins at 1. So the resulting indices in STEPS 
	are equal to the indices in CUBIT minus 1. Therefore information about the 
	mesh elements can be passed between STEPS and CUBIT (e.g. a compartment in 
	STEPS could be highlighted and viewed in CUBIT) but the user should take 
	care to allow for this difference of 1.
	note: The triangle indices will not be related to each other because this 
	information is not exported by CUBIT.
	
	For large meshes (over ~100,000 voxels) creation of STEPS Tetmesh can be time consuming.
	Therefore it is advised to only use this method once, and use meshio.py to save mesh in
	XML and ASCII format for quick loading in future.
	Example use:
		>>> import cubit
		>>> cubitfilename = '/cubit_meshes/spine.inp'
		>>> ### Use this function to create steps.geom.Tetmesh object from CUBIT output file ###
		>>> mymesh = cubit.makeMesh(cubitfilename, 1e-6)
		>>> import meshio
		>>> ### Save this mesh in XML (and ASCII) format for quick-loading in future ###
		>>> meshio.saveXML('/meshes/spine1', mymesh)
		>>> ### Mesh saved in /meshes/spine1.xml and /meshes/spine1.txt ###
	
	"""
			 
	print "Reading CUBIT file..."
	btime = time.time()
	
	# Try to open the CUBIT file. An error will be thrown if it doesn't exist
	cubfile = open(filename, 'r')
	
	# 1st line is: *HEADING
	# 2nd line is: cubit(<filename>): <date>: <time>
	# 3rd line is: *NODE
	line = cubfile.readline()
	# check we have the right kind of CUBIT output file here
	assert(line == '*HEADING\n')
	line = cubfile.readline() 
	line = cubfile.readline()
	assert(line == '*NODE\n')
	
	# OK, we have the right kind of file here, lets make the data structures
	# Problem is we don't know how big these are going to be at this point
	# For the numpy problem in STEPS 0.5. keep these structures as simple lists
	pnts = []
	tets = []
	# TODO: Look at exporting surface triangles in CUBIT
	
	line = cubfile.readline()
	line= line.replace(',', ' ')
	line = line.split()
	while(line[0] != '*ELEMENT'):
		assert(len(line) == 4)
		nodidx = int(line[0])
		# Create the node to add to the list of nodes
		node = [0.0, 0.0, 0.0]
		# set to the right scale
		node[0] = float(line[1])*scale
		node[1] = float(line[2])*scale
		node[2] = float(line[3])*scale 
		# And add to the list of points. That's all there is to it
		pnts+= node
		# Fetch the next line
		line = cubfile.readline()
		line= line.replace(',', ' ')
		line = line.split()
	assert(pnts.__len__()%3 == 0)
	
	# Now we are on the tets. We can just read to the end of the file now
	line = cubfile.readline()
	while(line):
		line = line.replace(',', ' ')
		line = line.split()
		assert(len(line) == 5)
		tetidx = int(line[0])
		# create the element to add to the list of tets 
		element = [0, 0, 0, 0]
		# Must subtract 1 to get proper indexing (starting at 0, not 1 as in CUBIT output)
		element[0] = int(line[1])-1		
		element[1] = int(line[2])-1
		element[2] = int(line[3])-1
		element[3] = int(line[4])-1
		tets+=element
		line = cubfile.readline()
	assert(tets.__len__()%4 == 0)
	
	print "Read CUBIT file succesfully"
	
	# Output the number of nodes and tets
	print "Number of nodes: ", pnts.__len__()/3
	print "Number of tets: ", tets.__len__()/4
	
	print "creating Tetmesh object in STEPS..."
	mesh = stetmesh.Tetmesh(pnts, tets)
	
	# Check mesh was created properly from this minimal information
	assert(pnts.__len__()/3 == mesh.nverts)
	assert(tets.__len__()/4 == mesh.ntets)
	
	return mesh


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
