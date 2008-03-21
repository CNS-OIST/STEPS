# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
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
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""Classes and routines for handling (static) tetrahedral meshes.

Construction of a tetrahedral mesh using the classes in this module should
be performed in the following order:
    
    1/ Create a TetMesh object and provide its initialization method with
       the raw geometric information, which includes:
       
            - an N_pnt * 3 floating point array of corner points
            
            - an N_tet * 4 integer array of indices. Each row of this
              array corresponds to 1 tetrahedron, given by 4 indices into 
              the list of corner points.
              
            - an N_tri * 3 integer array of indices. Each row of this
              array corresponds to one triangle, given by 3 indices into
              the list of corner points. It should be noted that not every
              'possible' triangle should be given, only those triangles
              that are of special importance, i.e. triangles belonging 
              to a patch that separates two compartments.
    
       The object will copy these arrays (rather than storing references to
       them) and will pre-compute some auxiliary tables which are useful
       for subsequent manipulation and visualization of the tetrahedral 
       mesh.
    
    2/ Annotate the tetrahedrons; this means creating Comp objects
       and specifying which tetrahedrons belong to a given compartment.
    
    3/ Annotate the triangles; this means creating Patch objects and 
       specifying which triangles belong to a given patch.

TODO: make iterators for running over tetrahedrons and triangles, that
give access to immediate properties such as volumes, neighbours, quality
etcetera.
"""


import math
import numpy
import scipy.io.mio as mio

import steps.error as serr
import steps.math.tetrahedron as stet
import steps.math.triangle as stri
import wm


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class _tet_iterator(object):
    
    """Auxiliary class for iterating over tetrahedrons in a mesh.
    """ 
    
    def __init__(self, tetmesh, tetlist):
        self._mesh = tetmesh
        self._tets = tetlist
        self._ntets = len(self._tets)
        # The iterator index (not tetrahedron index!)
        self._iidx = 0
    
    def __iter__(self):
        return self

    def next(self):
        if self._iidx == self._ntets:
            raise StopIteration
        t = Tet(self._mesh, self._tets[self._iidx])
        self._iidx = self._iidx + 1
        return t
    
    def reset(self):
        self._iidx = 0


class Tet(object):
    
    
    def __init__(self, tetmesh, tetidx):
        self._mesh = tetmesh
        self._tidx = tetidx
        self._nodes = numpy.empty((1, 4),dtype=int)
        self._nodes[0,:] = self._mesh._tets[self._tidx,:]
        self._nodes2 = self._mesh._tets[self._tidx,:]
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getIdx(self):
        return self._tidx
    
    idx = property(getIdx)
    
    
    def getNode0(self):
        return self._nodes2[0]
    
    node0 = property(getNode0)
    
    
    def getNode1(self):
        return self._nodes2[1]
    
    node1 = property(getNode1)
    
    
    def getNode2(self):
        return self._nodes2[2]
    
    node2 = property(getNode2)
    
    
    def getNode3(self):
        return self._nodes2[3]
    
    node3 = property(getNode3)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getQualityAR(self):
        """Compute the quality of the tetrahedron. 
        
        This method uses the aspect ratio (AR) metric for tetrahedron
        quality, given by dividing the length of the longest edge with
        the smallest altitude. The smaller this value, the more regular 
        the tetrahedron.
         
        TODO: implement getQualityAR
        """
        return 0.0
    
    qar = property(getQualityAR)
    
    
    def getQualityRER(self):
        """Compute the quality of the tetrahedron. 
        
        This method uses the radius-edge ratio (RER) metric for tetrahedron
        quality, given by dividing the radius of the tetrahedron's 
        circumsphere with the length of the shortest edge. 
        
        The smaller this value, the more regular the tetrahedron. The 
        lowest possible value of this metric is given by computing the 
        RER for a fully regular tetrahedron:
        
            Q = sqrt(6)/4 ~ 0.612
         
        This is a slightly weaker metric than getQualityAR, because 
        certain slivers (degenerate tetrahedrons) can still have a fairly 
        small value. 
        """
        return stet.circumradius(self._mesh._pnts, self._nodes) / \
            stet.shortestEdge(self._mesh._pnts, self._nodes)
    
    qrer = property(getQualityRER)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getBarycenter(self):
        return stet.barycenter(self._mesh._pnts, self._nodes)
    
    barycenter = property(getBarycenter)
    
    
    def getTriBarycenter(self, i):
        tri = { 0: [0,1,2], 1: [0,1,3], 2:[0,2,3], 3:[1,2,3] } 
        return stri.barycenter(self._mesh._pnts, self._nodes[:,tri[i]])
    
    def getTri0Barycenter(self):
        return stri.barycenter(self._mesh._pnts, self._nodes[:,[0,1,2]])
    
    def getTri1Barycenter(self):
        return stri.barycenter(self._mesh._pnts, self._nodes[:,[0,1,3]])
    
    def getTri2Barycenter(self):
        return stri.barycenter(self._mesh._pnts, self._nodes[:,[0,2,3]])
    
    def getTri3Barycenter(self):
        return stri.barycenter(self._mesh._pnts, self._nodes[:,[1,2,3]])
    
    tri0barycenter = property(getTri0Barycenter)
    tri1barycenter = property(getTri1Barycenter)
    tri2barycenter = property(getTri2Barycenter)
    tri3barycenter = property(getTri3Barycenter)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getVol(self):
        return numpy.asscalar(stet.vol(self._mesh._pnts, self._nodes))
    
    vol = property(getVol)
    
    
    def getArea(self, i):
        tri = { 0: [0,1,2], 1: [0,1,3], 2:[0,2,3], 3:[1,2,3] }
        n = numpy.asscalar(stri.area(self._mesh._pnts, self._nodes[:,tri[i]]))
        return n
    
    def getArea0(self):
        n = numpy.asscalar(stri.area(self._mesh._pnts, self._nodes[:,[0,1,2]]))
        return n

    def getArea1(self):
        n = numpy.asscalar(stri.area(self._mesh._pnts, self._nodes[:,[0,1,3]]))
        return n
    
    def getArea2(self):
        n = numpy.asscalar(stri.area(self._mesh._pnts, self._nodes[:,[0,2,3]]))
        return n
    
    def getArea3(self):
        n = numpy.asscalar(stri.area(self._mesh._pnts, self._nodes[:,[1,2,3]]))
        return n
    
    area0 = property(getArea0)
    area1 = property(getArea1)
    area2 = property(getArea2)
    area3 = property(getArea3)
    
    
    def getTriDist(self, i):
        tmp = self.getTriBarycenter(i) - self.barycenter
        return math.sqrt((tmp * tmp).sum(axis = 1))
    
    def getTri0Dist(self):
        tmp = self.tri0barycenter - self.barycenter
        return math.sqrt((tmp * tmp).sum(axis = 1))
    
    def getTri1Dist(self):
        tmp = self.tri1barycenter - self.barycenter
        return math.sqrt((tmp * tmp).sum(axis = 1))
    
    def getTri2Dist(self):
        tmp = self.tri2barycenter - self.barycenter
        return math.sqrt((tmp * tmp).sum(axis = 1))
    
    def getTri3Dist(self):
        tmp = self.tri3barycenter - self.barycenter
        return math.sqrt((tmp * tmp).sum(axis = 1))
    
    tri0dist = property(getTri0Dist)
    tri1dist = property(getTri1Dist)
    tri2dist = property(getTri2Dist)
    tri3dist = property(getTri3Dist)
    
    
    def getTetDist(self, i):
        ntet = self.getNextTet(i)
        if ntet == None: return 0.0
        tmp = self.barycenter - ntet.barycenter
        return math.sqrt((tmp * tmp).sum(axis = 1))
    
    def getTet0Dist(self):
        ntet = self.ntet0
        if ntet == None: return 0.0
        tmp = self.barycenter - ntet.barycenter
        return math.sqrt((tmp * tmp).sum(axis = 1))
    
    def getTet1Dist(self):
        ntet = self.ntet1
        if ntet == None: return 0.0
        tmp = self.barycenter - ntet.barycenter
        return math.sqrt((tmp * tmp).sum(axis = 1))
    
    def getTet2Dist(self):
        ntet = self.ntet2
        if ntet == None: return 0.0
        tmp = self.barycenter - ntet.barycenter
        return math.sqrt((tmp * tmp).sum(axis = 1))
    
    def getTet3Dist(self):
        ntet = self.ntet3
        if ntet == None: return 0.0
        tmp = self.barycenter - ntet.barycenter
        return math.sqrt((tmp * tmp).sum(axis = 1))
    
    tet0dist = property(getTet0Dist)
    tet1dist = property(getTet1Dist)
    tet2dist = property(getTet2Dist)
    tet3dist = property(getTet3Dist)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def isInside(self, p):
        return stet.inside(self._mesh._pnts[self._nodes2], numpy.array(p))
    
    
    def getRanPnt(self, rng, n = 1):
        return stet.ranpnt(rng, self._mesh._pnts[self._nodes2], n)
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getComp(self):
        return self._mesh._tet_comps[self._tidx]

    comp = property(getComp)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getNextTetIdx(self, i):
        return self._mesh._tet_tet_neighbours[self._tidx,i]
    
    def getNextTet0Idx(self):
        return self._mesh._tet_tet_neighbours[self._tidx,0]
    
    def getNextTet1Idx(self):
        return self._mesh._tet_tet_neighbours[self._tidx,1]
    
    def getNextTet2Idx(self):
        return self._mesh._tet_tet_neighbours[self._tidx,2]
    
    def getNextTet3Idx(self):
        return self._mesh._tet_tet_neighbours[self._tidx,3]
    
    ntet0idx = property(getNextTet0Idx)
    ntet1idx = property(getNextTet1Idx)
    ntet2idx = property(getNextTet2Idx)
    ntet3idx = property(getNextTet3Idx)
    
    
    def getNextTet(self, i):
        nidx = self.getNextTetIdx(i)
        if nidx == -1: return None
        return Tet(self._mesh, nidx)
    
    def getNextTet0(self):
        nidx = self.ntet0idx
        if nidx == -1: return None
        return Tet(self._mesh, nidx)
    
    def getNextTet1(self):
        nidx = self.ntet1idx
        if nidx == -1: return None
        return Tet(self._mesh, nidx)
    
    def getNextTet2(self):
        nidx = self.ntet2idx
        if nidx == -1: return None
        return Tet(self._mesh, nidx)
    
    def getNextTet3(self):
        nidx = self.ntet3idx
        if nidx == -1: return None
        return Tet(self._mesh, nidx)
    
    ntet0 = property(getNextTet0)
    ntet1 = property(getNextTet1)
    ntet2 = property(getNextTet2)
    ntet3 = property(getNextTet3)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getNextTriIdx(self, i):
        return self._mesh._tet_tri_neighbours[self._tidx,i]
    
    def getNextTri0Idx(self):
        return self._mesh._tet_tri_neighbours[self._tidx,0]
    
    def getNextTri1Idx(self):
        return self._mesh._tet_tri_neighbours[self._tidx,1]
    
    def getNextTri2Idx(self):
        return self._mesh._tet_tri_neighbours[self._tidx,2]
    
    def getNextTri3Idx(self):
        return self._mesh._tet_tri_neighbours[self._tidx,3]
    
    ntri0idx = property(getNextTri0Idx)
    ntri1idx = property(getNextTri1Idx)
    ntri2idx = property(getNextTri2Idx)
    ntri3idx = property(getNextTri3Idx)
    
    
    def getNextTri(self, i):
        nidx = self.getNextTriIdx(i)
        if nidx == -1: return None
        return Tet(self._mesh, nidx)
    
    def getNextTri0(self):
        nidx = self.ntri0idx
        if nidx == -1: return None
        return Tet(self._mesh, nidx)
    
    def getNextTri1(self):
        nidx = self.ntri1idx
        if nidx == -1: return None
        return Tet(self._mesh, nidx)
    
    def getNextTri2(self):
        nidx = self.ntri2idx
        if nidx == -1: return None
        return Tet(self._mesh, nidx)
    
    def getNextTri3(self):
        nidx = self.ntri3idx
        if nidx == -1: return None 
        return Tet(self._mesh, nidx)
    
    ntri0 = property(getNextTri0)
    ntri1 = property(getNextTri1)
    ntri2 = property(getNextTri2)
    ntri3 = property(getNextTri3)
        


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class _tri_iterator(object):
    
    """Auxiliary class for iterating over triangles in a mesh.
    """ 
    
    def __init__(self, tetmesh, trilist):
        self._mesh = tetmesh
        self._tris = trilist
        self._ntris = len(self._tris)
        # The iterator index (not triangle index!)
        self._iidx = 0
    
    def __iter__(self):
        return self
    
    def next(self):
        if self._iidx == self._ntris:
            raise StopIteration
        t = Tri(self._mesh, self._tris[self._iidx])
        self._iidx = self._iidx + 1
        return t
    
    def reset(self):
        self._iidx = 0


class Tri(object):
    
    def __init__(self, tetmesh, triidx):
        self._mesh = tetmesh
        self._tidx = triidx
        self._nodes = numpy.empty((1, 3),dtype=int)
        self._nodes[0,:] = self._mesh._tris[self._tidx,:]
        self._nodes2 = self._mesh._tris[self._tidx,:]
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getIdx(self):
        return self._tidx
    
    idx = property(getIdx)
    
    
    def getNode0(self):
        return self._nodes2[0]
    
    node0 = property(getNode0)
    
    
    def getNode1(self):
        return self._nodes2[1]
    
    node1 = property(getNode1)
    
    
    def getNode2(self):
        return self._nodes2[2]
    
    node2 = property(getNode2)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getArea(self):
        return numpy.asscalar(stri.area(self._mesh._pnts, self._nodes))
    
    area = property(getArea)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getPatch(self):
        return self._mesh._tri_patches[self._tidx]
    
    patch = property(getPatch)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getInnerTetIdx(self):
        return self._mesh._tri_tet_neighbours[self._tidx,0]
    
    def getOuterTetIdx(self):
        return self._mesh._tri_tet_neighbours[self._tidx,1]
    
    itetidx = property(getInnerTetIdx)
    otetidx = property(getOuterTetIdx)
    
    
    def getInnerTet(self):
        iidx = self.itetidx
        if iidx == -1: return None
        return Tet(self._mesh, iidx)
    
    def getOuterTet(self):
        oidx = self.otetidx
        if oidx == -1: return None
        return Tet(self._mesh, oidx)
    
    itet = property(getInnerTet)
    otet = property(getOuterTet)
    
    
    def getBarycenter(self):
        return stri.barycenter(self._mesh._pnts, self._nodes)
    
    barycenter = property(getBarycenter)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class TetMesh(wm.Geom):


    """The main container class for static tetrahedral meshes.
    
    This class stores the node points, tetrahedrons and boundary triangles
    that comprise the tetrahedral mesh. In addition, it also precomputes
    some auxiliary data for the mesh as a whole:
    
        - Rectangular, axis-aligned bounding box.
        - Overall volume
    
    Auxiliary data is also stored for the tetrahedrons:
    
        - Volume of each tetrahedron.
        - For each tetrahedron, the indices of its 4 neighbouring tets. 
          If there is no neighbour (i.e. if the tetrahedron lies on the 
          border), this index will be -1. The sequence of neighbours is
          determined by the following common boundary triangles: (0,1,2);
          (0,1,3); (0,2,3); (1,2,3).
        - For each tetrahedron, the indices of its 4 neighbouring 
          boundary triangles. If their is no neighbouring triangle (i.e.
          if the triangle has not been explicitly added as a boundary 
          triangle) this index is -1. The sequence of neighbours is also
          determined by (0,1,2); (0,1,3); (0,2,3); (1,2,3).
        - Tetrahedron annotation: this is basically the compartment (Comp
          object) that a tetrahedron belongs to.
    
    And finally, also for the triangles:
    
        - Area of each triangle.
        - Normal vector for each triangle, normalized to length 1.0.
        - For each triangle, the indices of its inside and outside 
          neighbouring tetrahedron. If this tetrahedron does not exist 
          (because the triangle lies on the outer boundary), this index 
          will be -1.
        - Triangle annotation: this is basically the Patch object that a 
          triangle belongs to.
    
    NOTES: 
        - Keep in mind that these auxiliary tables take up extra memory!
        - Adding/deleting/moving node points after initiation is 
          currently not implemented.
        - Adding/deleting of tetrahedrons after initiation is currently 
          not implemented.
        - Deleting triangles is currently not implemented.
    
    SEE ALSO:
    
        The docstring of module steps.geom.tetmesh, for a larger scale 
        overview of the role this class plays in representing tetrahedral
        meshes.
    """
    
    
    def __init__(self, pnts, tets, tris = None, **params):
        """Initialize a TetMesh object.
        
        PARAMETERS:
            pnts
                An X*3 array of floats, storing the coordinates of node points.
            tets
                A Y * 4 array of integers, storing the corner nodes for each 
                tetrahedron as an index into array nodes.
            tris
                An (optional) Z * 3 array of integers, storing the corner 
                nodes for triangles (i.e. special boundaries between or 
                around tetrahedrons.
        
        RETURNS:
            ---
        
        RAISES:
            ---
        """
        # Call the initializer for the parent class.
        wm.Geom.__init__(self)
        
        # SET UP: CORNER POINTS
        assert pnts.shape[1] == 3
        assert pnts.dtype == float
        #npnts = pnts.shape[0]
        self._pnts = pnts.copy()

        # Find minimal and maximal boundary values.
        self.__minx = self._pnts[:,0].min()
        self.__miny = self._pnts[:,1].min()
        self.__minz = self._pnts[:,2].min()
        self.__maxx = self._pnts[:,0].max()
        self.__maxy = self._pnts[:,1].max()
        self.__maxz = self._pnts[:,2].max()
        
        # SET UP: TETRAHEDRONS
        assert tets.shape[1] == 4
        assert tets.dtype == int
        ntets = tets.shape[0]
        self._tets = tets.copy()
        # Check range; check duplicates
        self._tet_vols = stet.vol(self._pnts, self._tets)
        self._tet_tri_neighbours = numpy.ones((ntets, 4), dtype = int) * -1
        self._tet_tet_neighbours = numpy.empty((ntets, 4), dtype = int)
        self.__set_tet_tet_neighbours()
        self._tet_comps = [ None ] * ntets
        
        # SETUP: TRIANGLES
        self._tris = numpy.empty((0, 3), dtype=int)
        self._tri_areas = numpy.empty((0), dtype=float)
        self._tri_norms = numpy.empty((0,3), dtype=float)
        self._tri_tet_neighbours = numpy.empty((0,2), dtype=int)
        self._tri_patches = [ ]
        if tris != None: 
            self.addTris(tris)
    
    
    def __set_tet_tet_neighbours(self):
        """Find neighbouring tetrahedrons.
        
        NOTES:
            - Internal function.
            - Might benefit from reimplementation in C++.
        """
        #
        ntets = self._tets.shape[0]
        tmp1 = numpy.asarray(xrange(0, ntets))
        tmp2 = numpy.ones(ntets, dtype=int)
        
        # Create array for triangle formed by edges (0,1,2).
        t012 = self._tets[:,(0,1,2)]
        t012.sort()
        t012 = numpy.c_[ t012, tmp1, tmp2 * 0 ]
        
        # Create array for triangle formed by edges (0,1,3).
        t013 = self._tets[:,(0,1,3)]
        t013.sort()
        t013 = numpy.c_[ t013, tmp1, tmp2 ]
        
        # Create array for triangle formed by edges (0,2,3).
        t023 = self._tets[:,(0,2,3)]
        t023.sort()
        t023 = numpy.c_[ t023, tmp1, tmp2 * 2 ]
        
        # Create array for triangle formed by edges (1,2,3).
        t123 = self._tets[:,(1,2,3)]
        t123.sort()
        t123 = numpy.c_[ t123, tmp1, tmp2 * 3 ]
        
        # Concatenate arrays and sort by first three columns only.
        # The fourth column is the tetrahedron it belongs to, the 
        # fifth column denotes which boundary.
        tn = numpy.r_[ t012, t013, t023, t123 ]
        tn = tn[numpy.lexsort((tn[:,2], tn[:,1], tn[:,0])),:]
        
        # Loop over this matrix to detect duplicate entries. These duplicate
        # entries signal a shared triangle; single entries denote a boundary
        # triangle.
        ix1 = 0
        ix2 = 1
        nx = tn.shape[0]
        while ix1 < nx:
            # Check if we've reached the last one.
            if ix1 == nx - 1:
                # Special case: the last triangle has been reached; it should 
                # be a boundary triangle, else it couldn't be the last one.
                self._tet_tet_neighbours[tn[ix1,3], tn[ix1,4]] = -1
                break
            # Shared or unique boundary?
            if (tn[ix1, 0:3] == tn[ix2, 0:3]).all():
                # We're dealing with a shared triangle.
                tet1 = tn[ix1, 3]
                tet2 = tn[ix2, 3]
                self._tet_tet_neighbours[tet1, tn[ix1, 4]] = tet2
                self._tet_tet_neighbours[tet2, tn[ix2, 4]] = tet1
                ix1 = ix1 + 2
            else:
                # We're dealing with a boundary triangle.
                self._tet_tet_neighbours[tn[ix1,3], tn[ix1,4]] = -1
                ix1 = ix1 + 1
            # Update other counter.
            ix2 = ix1 + 1
    
    
    def __set_tet_tri_neighbours(self, tri_idcs):
        """Resolve neighbouring tetrahedrons for a set of triangles.
        
        PARAMETERS:
            tri_idcs
                Indices of triangles whose neighbours should be resolved.
        """
        # First, we make a dictionary of the triangles and their indices.
        tris = self._tris[tri_idcs, :]
        tris.sort()
        ntris = len(tri_idcs)
        tri_dict = { }
        for ii in xrange(0, ntris):
            tri_dict[(tris[ii,0], tris[ii,1], tris[ii,2])] = tri_idcs[ii]
        # tri_dict now can be used to find the original ('real') index of
        # a triangle defined by its corner points.
        
        #
        ntets = self._tets.shape[0]
        
        # Find all (0,1,2)
        unr012 = [ x for x in xrange(0, ntets) \
            if self._tet_tri_neighbours[x,0] == -1 ]
        # DEBUG 09-Oct-2007:
        # tri012 = self._tets[unr012,(0,1,2)
        # doesn't work in Numpy, so we have to turn it into two statements.
        # The same changes made in the cases for (0,1,3), (0, 2, 3) and
        # (1, 2, 3).
        tri012 = self._tets[unr012,:]
        tri012 = tri012[:,(0,1,2)]
        tri012.sort()
        n012 = len(unr012)
        for ii in xrange(0, n012):
            teidx = unr012[ii]
            tridx = tri_dict.get((tri012[teidx, 0], \
                                  tri012[teidx, 1], \
                                  tri012[teidx, 2]), -1)
            if tridx == -1:
                continue
            if self._tri_tet_neighbours[tridx,0] == -1:
                self._tri_tet_neighbours[tridx,0] = teidx
            elif self._tri_tet_neighbours[tridx,1] == -1:
                self._tri_tet_neighbours[tridx,1] = teidx
            else:
                assert False
            self._tet_tri_neighbours[teidx,0] = tridx
        
        # Find all (0,1,3)
        unr013 = [ x for x in xrange(0, ntets) \
            if self._tet_tri_neighbours[x,1] == -1 ]
        tri013 = self._tets[unr013,:]
        tri013 = tri013[:,(0,1,3)]
        tri013.sort()
        n013 = len(unr013)
        for ii in xrange(0, n013):
            teidx = unr013[ii]
            tridx = tri_dict.get((tri013[teidx, 0], \
                                  tri013[teidx, 1], \
                                  tri013[teidx, 2]), -1)
            if tridx == -1:
                continue
            if self._tri_tet_neighbours[tridx,0] == -1:
                self._tri_tet_neighbours[tridx,0] = teidx
            elif self._tri_tet_neighbours[tridx,1] == -1:
                self._tri_tet_neighbours[tridx,1] = teidx
            else:
                assert False
            self._tet_tri_neighbours[teidx,1] = tridx
        
        # Find all (0,2,3)
        unr023 = [ x for x in xrange(0, ntets) \
            if self._tet_tri_neighbours[x,2] == -1 ]
        tri023 = self._tets[unr023,:]
        tri023 = tri023[:,(0,2,3)]
        tri023.sort()
        n023 = len(unr023)
        for ii in xrange(0, n023):
            teidx = unr023[ii]
            tridx = tri_dict.get((tri023[teidx, 0], \
                                  tri023[teidx, 1], \
                                  tri023[teidx, 2]), -1)
            if tridx == -1:
                continue
            if self._tri_tet_neighbours[tridx,0] == -1:
                self._tri_tet_neighbours[tridx,0] = teidx
            elif self._tri_tet_neighbours[tridx,1] == -1:
                self._tri_tet_neighbours[tridx,1] = teidx
            else:
                assert False
            self._tet_tri_neighbours[teidx,2] = tridx
        
        # Find all (1,2,3)
        unr123 = [ x for x in xrange(0, ntets) \
            if self._tet_tri_neighbours[x,3] == -1 ]
        tri123 = self._tets[unr123,:]
        tri123 = tri123[:,(1,2,3)]
        tri123.sort()
        n123 = len(unr123)
        for ii in xrange(0, n123):
            teidx = unr123[ii]
            tridx = tri_dict.get((tri123[teidx, 0], \
                                  tri123[teidx, 1], \
                                  tri123[teidx, 2]), -1)
            if tridx == -1:
                continue
            if self._tri_tet_neighbours[tridx,0] == -1:
                self._tri_tet_neighbours[tridx,0] = teidx
            elif self._tri_tet_neighbours[tridx,1] == -1:
                self._tri_tet_neighbours[tridx,1] = teidx
            else:
                assert False
            self._tet_tri_neighbours[teidx,3] = tridx
        
            
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def checkConsistency(self):
        """Perform a consistency check on the entire mesh.
        """
        pass
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getBoundMin(self):
        """Get the minimal coordinate of the rectangular bounding box
        around the entire mesh.
        
        RETURNS:
            A tuple (x,y,z).
        """
        return ( self.__minx, self.__miny, self.__minz )
    
    
    def getBoundMax(self):
        """Get the maximal coordinate of the rectangular bounding box
        around the entire mesh.
        
        RETURNS:
            A tuple (x,y,z).
        """
        return ( self.__maxx, self.__maxy, self.__maxz )
    
    
    def countPnts(self):
        """Return the number of node points in the mesh.
        """
        return self._pnts.shape[0]
    
    npnts = property(countPnts)
    
        
    def getPnts(self):
        return self._pnts
    
    pnts = property(getPnts)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
            
    
    def __find_unique_tris(self, tris, incl_mesh = True):
        """Find all triangles in an array of triangle indices that are unique.
        
        PARAMETERS:
            tris
                An N * 3 array of integer (index) values.
            incl_mesh
                Flag to control whether the triangles that have already 
                been added to the mesh should be included (== True) or 
                not (== False). This will most often be true of course.
        
        RETURNS:
            An array of indices pointing to the unique triangles in the
            input.
        
        NOTES:
            Internal function.
        """
        # The dictionary that will collect triangles.
        etris = { }
        
        # Do we also want to check vis-a-vis existing triangles in the mesh?
        notris = self._tris.shape[0]
        if (incl_mesh == True) and (notris > 0):
            # Create a sorted array of old triangles.
            otris = self._tris.copy()
            otris.sort()
            otris = otris[numpy.lexsort((otris[:,2],otris[:,1],otris[:,0])),:]
            # Put them in the dictionary and mark them as 'already existing',
            # which means marking them with -1.
            for ii in notris:
                etris.setdefault((otris[ii,0], otris[ii,1], otris[ii,2]), -1)
        
        # Create a sorted array of new triangles.
        nntris = tris.shape[0]
        ntris = tris.copy()
        ntris.sort()
        ntris = ntris[numpy.lexsort((ntris[:,2],ntris[:,1],ntris[:,0])),:]
        for ii in xrange(0, nntris):
            etris.setdefault((ntris[ii,0], ntris[ii,1], ntris[ii,2]), ii)
        
        # Build a set of unique triangles.
        # SOMETHING WRONG HERE
        utris = [ x for x in etris.itervalues() if x >= 0 ]
        utris.sort()
        return utris
        
    
    def addTris(self, tris):
        """Add a number of boundary triangles to the mesh.
        """
        # Check basic assumptions about input.
        assert tris.shape[1] == 3
        assert tris.dtype == int
        # Filter out duplicates.
        utris = tris[self.__find_unique_tris(tris), :]
        # Check trivial case: return.
        ntris = utris.shape[0]
        if ntris == 0: 
            return
        
        # SIMPLE STUFF.
        # Add the triangles themselves.
        self._tris = numpy.r_[ self._tris, utris ]
        # Add their areas.
        self._tri_areas = \
            numpy.r_[ self._tri_areas, stri.area(self._pnts, utris) ]
        # Add their normal vectors.
        self._tri_norms = \
            numpy.r_[ self._tri_norms, stri.normal(self._pnts, utris) ]
        
        # MORE DIFFICULT STUFF.
        # Find the neighbouring tetrahedrons, depending on the current
        # triangle order.
        self._tri_tet_neighbours = \
            numpy.r_[ self._tri_tet_neighbours, \
                       numpy.ones((ntris, 2), dtype=int) * -1 ]
        xrstop = self._tris.shape[0]
        xrstart = xrstop - ntris
        self.__set_tet_tri_neighbours(xrange(xrstart,xrstop))
        
        # OTHER STUFF.
        # Set their annotation to the default value: None.
        self._tri_patches += [None] * ntris
        
    
    def countTris(self):
        """Return the number of boundary triangles in the mesh.
        """
        return self._tris.shape[0]
    
    ntris = property(countTris)
    
    
    def getTris(self):
        return _tri_iterator(self, xrange(0,self.ntris))
    
    tris = property(getTris)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def countTets(self):
        """Return the number of tetrahedrons in the mesh.
        """
        return self._tets.shape[0]
    
    ntets = property(countTets)
    
    
    #def getTets(self):
    #    return self._tets
    #
    #tets = property(getTets)
    
    
    def getTets(self):
        return _tet_iterator(self, xrange(0,self.ntets))
    
    tets = property(getTets)
    
    
    def findTets(self, p):
        r = [ ]
        for i in self.tets:
            if i.isInside(p): r.append(i.idx)
        return r
    
    def findTet(self, p):
        for i in self.tets:
            if i.isInside(p): return i.idx
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #


    def findCompByPoint(self, p, inclusive = True):
        """Return the compartment(s) to which one or more 3D point(s) belong.
        
        PARAMETERS:
            p
                The 3-dimensional point(s) that should be tested. This 
                parameter can take a lot of forms. The simplest is a 
                single tuple of real values (x, y, z). A list of such tuples
                can also be given. Alternatively, the method accept a
                N*3 Numpy array object. Internally, the input parameter
                gets temporarily converted into a Numpy array anyway.
                
            inclusive
                A boolean value. If True (as it is by default), a point  
                lying on the edge of a compartment (i.e. lying on a patch) 
                is considered to be on the inside. Patches are often
                shared by two compartments, so the method simply returns 
                the last compartment that tested positively. When False, 
                the method returns None for such edge points.
                
        RETURNS:
            A list of references to compartments (or None values for points 
            not inside any compartment). This is always a list, even if the
            parameter was a single tuple.
        
        RAISES:
            ---
        """
        pass


    def findCompByTet(self, t):
        """Return the compartment(s) to which one or more tetrahedrons belong.
        
        PARAMETERS:
            t
                The tetrahedron(s) that should be tested. This parameter can
                take multiple forms. The simplest form is a single integer
                index value, indicating a single tetrahedron to be tested.
                A sequence of such indices can also be given.
        
        RETURNS:
            A list of references to compartments (or None for voxel indices
            that are not inside any compartments). This is always a list,
            even if the parameter was a single integer index value.
        
        RAISES:
            ---
        """
        pass


    def findTetByPoint(self, p, inclusive = True):
        """Return the voxel(s) to which one or more 3D point(s) belong.
        
        PARAMETERS:
            p
                The 3-dimensional point(s) that should be tested. This 
                parameter can take a lot of forms. The simplest is a 
                single tuple of real values (x, y, z). A list of such tuples
                can also be given. Alternatively, the method accept a
                N*3 Numpy array object. Internally, the input parameter
                gets temporarily converted into a Numpy array anyway.

            inclusive
                A boolean value. If True (as it is by default), a point  
                lying on the edge of a voxel is considered to be on the 
                inside. As triangles, edges or corner points are usually
                shared by two or more voxels, the method simply returns 
                the last voxel that tested positively. When False, the 
                method returns -1 for such points.

        RETURNS:
            An array of tetrahedral voxel indices (or -1 for points not 
            inside any tetrahedron).
        
        RAISES:
            ---
        """
        pass
        
        
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def allBorderTris(self):
        pass
    
    
    def allBorderTriIndices(self, reorient = True):
        pass
        
        
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Comp(wm.Comp):


    """Provides annotation for a group of tetrahedrons in a TetMesh.
    """
    
    
    def __init__(self, id, container, tets):
        """Initialize a tetrahedral mesh compartment.
        
        PARAMETERS:
            id
                The identifier of the compartment.
            container
                A TetMesh object.
            tets
                A sequence of indices to tetrahedrons in the tetmesh.
        
        RETURNS:
            ---
        
        RAISES:
            steps.error.ArgumentError
                When specified tetrahedrons do not exist; 
                when the identifier is not valid.
        
        See also:
            steps.geom.wm.Comp.__init__()
        """        
        # Store all tetrahedron indices as a set.
        tetset = set(tets)
        # Check index range.
        if (min(tetset) < 0) or (max(tetset) >= container.ntets):
            raise serr.ArgumentError, 'Invalid tetrahedron index specified.'
        # Make a list.
        self._tet_indices = list(tetset)
        
        # Check if specified tetrahedrons belong to another comp already.
        for ii in self._tet_indices:
            if container._tet_comps[ii] != None:
                raise serr.ArgumentError, \
                    'Cannot add tetrahedron to compartment.'

        # Call parent object initializer.
        wm.Comp.__init__(self, id, container)
        
        # Compute volume.
        self._vol = self.container._tet_vols[self._tet_indices].sum()
        
        # Perform annotation of tetrahedrons.
        for i in self._tet_indices:
            self.container._tet_comps[i] = self
        
        # Compute bounds.
        _tets = self.container._tets[self._tet_indices,:]
        pts = set(_tets[:,0]) | set(_tets[:,1]) \
            | set(_tets[:,2]) | set(_tets[:,3])
        pts = list(pts)
        xs = self.container.pnts[pts,0]
        self.__minx = xs.min()
        self.__maxx = xs.max()
        ys = self.container.pnts[pts,1]
        self.__miny = ys.min()
        self.__maxy = ys.max()
        zs = self.container.pnts[pts,2]
        self.__minz = zs.min()
        self.__maxz = zs.max()
        
    
    def setVol(self, vol):
        """Method disabled. Volume is computed automatically.
        
        See also:
            steps.geom.wm.Comp.setVol()
        """
        pass
    
    
    def checkConsistency(self):
        """Perform a consistency check on the compartment.
        """
        pass
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def _addIPatch(self, patch):
        """Internal function.
        
        REWRITE INTERNAL COMMENTS
        """
        # Call the parent method.
        super(Comp, self)._addIPatch(patch)
        # Nothing else needs to be done (in the current implementation).
        
    
    def _delIPatch(self, patch):
        """Internal function.
        """
        # Call original method.
        super(Comp, self)._delIPatch(patch)
        # Nothing else needs to be done (in the current implementation).
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def _addOPatch(self, patch):
        """Internal function.
        """
        # Call original method.
        super(Comp, self)._addOPatch(patch)
        # Nothing else needs to be done (in the current implementation).
    
    
    def _delOPatch(self, patch):
        """Internal function.
        """
        # Call original method.
        super(Comp, self)._delOPatch(patch)
        # Nothing else needs to be done (in the current implementation).
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getBoundMin(self):
        """Get the minimal coordinate of the rectangular bounding box
        around the compartment.
        
        RETURNS:
            A tuple (x, y, z).
        """
        return ( self.__minx, self.__miny, self.__minz )
    
    
    def getBoundMax(self):
        """Get the maximal coordinate of the rectangular bounding box
        around the compartment.
        
        RETURNS:
            A tuple (x, y, z).
        """
        return ( self.__maxx, self.__maxy, self.__maxz )
        

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def isTetInside(self, t):
        """Test whether a tetrahedron(s) (specified by index) is inside
        this compartment.
        
        PARAMETERS:
            t
                A tetrahedral index or sequence of tetrahedral indices.
                
        RETURNS:
            A boolean value, or a list of boolean values.
        
        RAISES:
            ---
        """
        pass
    
    
    def findTetByPoint(self, p, inclusive = True):
        """Return the voxel(s) to which one or more 3D point(s) belong.
        
        When a point is not inside any tetrahedron, or when it lies right on 
        the edge between two or more neighbouring tetrahedrons, the method 
        returns -1 for that point.
        
        PARAMETERS:
            p
                The 3-dimensional point(s) that should be tested. This 
                parameter can take a lot of forms. The simplest is a 
                single tuple of real values (x, y, z). A list of such tuples
                can also be given. Alternatively, the method accept a
                N*3 Numpy array object. Internally, the input parameter
                gets temporarily converted into a Numpy array anyway.

        RETURNS:
            An array of tetrahedral voxel indices (or -1 for points not 
            inside any tetrahedron).
        
        RAISES:
            ---
        """
        pass
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    def getTets(self):
        return _tet_iterator(self.container, self._tet_indices)
    
    tets = property(getTets)
    
    def allPatches(self):
        """Return a list of all patches surrounding this compartment.
        """
        pass
    
    
    def allTets(self):
        """Return an array of all tetrahedrons in this compartments.
        """
        pass
    
    
    def allTetIndices(self):
        """Return an array of indices to all tetrahedrons in this compartment.
        """
        pass
    
    
    def allBorderTris(self, reorient = True):
        """Return an array of all border triangles surrounding this
        compartment.
        """
        pass
    
    
    def allBorderTriIndices(self):
        """Return an array of indices to all triangles surrounding this
        compartment.
        """
        pass


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Patch(wm.Patch):


    """
    """
    
    
    def __init__(self, id, container, tris, icomp, ocomp = None):
        """Initialize a tetrahedral mesh patch.
        
        Currently, if some of the specified triangles have already been 
        assigned to another patch, the initializer will raise an error.
        It also assumes that the tetrahedral neighbour(s) of each patch 
        triangle belong to the inner (and outer) compartments specified
        in parameters icomp (and ocomp).
        
        PARAMETERS:
            id
                The identifier of the patch.
            container
                A TetMesh object.
            tris
                A sequence of indices to triangles in the tetmesh.
            icomp
                A pointer to a valid steps.geom.tetmesh.Comp object.
            ocomp
                A pointer to a valid steps.geom.tetmesh.Comp object; this
                parameter can be None, if the patch lies on the outermost
                boundary of the mesh.
        
        RETURNS:
            ---
        
        RAISES:
            steps.error.ArgumentError
                When specified triangles do not exist; when specified
                triangles have already been assigned to another patch;
                when the neighbouring tetrahedrons do not match the
                specified compartments; when the identifier is not valid.
        
        See also:
            steps.geom.wm.Patch.__init__()
        """
        # Store all triangle indices as a set.
        triset = set(tris)
        # Check index range.
        if (min(triset) < 0) or (max(triset) >= container.ntris):
            raise serr.ArgumentError, 'Invalid triangle index specified.'
        # Make a list.
        self._tri_indices = list(triset)
        
        # Check if specified triangles belong to another patch already.
        for ii in self._tri_indices:
            if container._tri_patches[ii] != None:
                raise serr.ArgumentError, 'Cannot add triangle to patch.'
        
        # Check if the triangle neighbours belong to the correct compartments.
        # (Taking into account that their relative position might need to be
        # flipped later on... we already build a flip list at this point that
        # we will use later in this initializer.)
        fliplist = [ ]
        for ii in self._tri_indices:
            icmp = container._tri_tet_neighbours[ii,0]
            if icmp != -1:
                icmp = container._tet_comps[icmp]
            else:
                icmp = None
            ocmp = container._tri_tet_neighbours[ii,1]
            if ocmp != -1:
                ocmp = container._tet_comps[ocmp]
            else:
                ocmp = None
            if (icmp, ocmp) == (icomp, ocomp): continue
            if (icmp, ocmp) == (ocomp, icomp): 
                fliplist.append(ii)
                continue
            raise serr.ArgumentError, 'Triangle cannot belong to this patch.'
        
        # Okay, everything should work out smoothly from here on, so let's
        # start changing data structures.
        #
        # First, we call parent object initializer.
        wm.Patch.__init__(self, id, container, icomp, ocomp)
        
        # Next, we compute area.
        self._area = self.container._tri_areas[self._tri_indices].sum()
        
        # Annotate the triangles as belonging to this patch.
        for i in self._tri_indices:
            self.container._tri_patches[i] = self
        
        # Switch internal and external tetrahedrons for triangles where this
        # is necessary.
        if len(fliplist) > 0:
            tmp = self.container._tri_tet_neighbours[fliplist,0]
            self.container._tri_tet_neighbours[fliplist,0] = \
                self.container._tri_tet_neighbours[fliplist,1]
            self.container._tri_tet_neighbours[fliplist,1] = tmp
        
        # Check the triangle orientation, and flip vertices if needed.
        # First, for each triangle, compute the vector from the inner
        # neighbouring tetrahedron barycenter to the triangle barycenter.
        # Then, check if the dot product of this vector with the triangle
        # normal is positive. If not, the triangle vertices need to be
        # flipped and the normal has to be negated.
        tris = self.container._tris[self._tri_indices,:]
        tri_barycenters = stri.barycenter(self.container.pnts, tris)
        tets = self.container._tri_tet_neighbours[self._tri_indices,0]
        tets = self.container._tets[tets,:]
        tet_barycenters = stet.barycenter(self.container.pnts, tets)
        vecs1 = tet_barycenters - tri_barycenters
        vecs2 = self.container._tri_norms[self._tri_indices,:]
        fliplist2 = numpy.array(self._tri_indices)
        fliplist2 = fliplist2[(vecs1 * vecs2).sum(axis=1) < 0.0]
        tmp = self.container._tris[fliplist2,0]
        self.container._tris[fliplist2,0] = self.container._tris[fliplist2,1]
        self.container._tris[fliplist2,1] = tmp
        self.container._tri_norms[fliplist2,:] = \
            - self.container._tri_norms[fliplist2,:]
    
    
    def setArea(self, area):
        """Method disabled. Area is computed automatically.
        
        See also:
            steps.geom.wm.Patch.setArea()
        """
        pass


    def checkConsistency(self):
        """Perform a consistency check on the patch.
        """
        pass


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def fromCDF(cdf):
    """Construct a mesh from NetCDF.
    
    TODO:
        Implement!
    """
    pass
    

def toCDF(tetmesh):
    """Save the mesh in NetCDF.
    
    TODO:
        Implement!
    """
    pass


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def toMatlab(filename, tetmesh, prefix=''):
    """Save the mesh in a Matlab file.
    
    Currently, only the following information is saved:
        1/ the node point coordinates
        2/ the boundary triangles
        3/ the tetrahedrons themselves
    
    PARAMETERS:
        filename
            Path and name of the destination Matlab file.
        tetmesh
            A steps.geom.tetmesh.TetMesh object.
        prefix
            Optional parameter to specify whether the variable names 
            should get a prefix (to avoid name clashes in Matlab).
    
    RETURNS:
        ---
    
    RAISES:
        ---
    """
    # If prefix has been defined, check whether we should add 
    # an extra underscore.
    prefix2 = prefix
    if len(prefix2) > 0:
        if prefix2[-1] != "_": 
            prefix2 = prefix2 + "_"
    # Create a dictionary of the stuff we're going to save.
    contents = { }
    contents[prefix2 + "nodes"] = tetmesh.pnts
    contents[prefix2 + "tris"] = tetmesh.tris
    contents[prefix2 + "tets"] = tetmesh.tets
    # Save to a Matlab file.
    mio.savemat(filename, contents)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
