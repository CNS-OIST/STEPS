# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# STEPS - STochastic Engine for Pathway Simulation
# Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
#
# $Id$
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


"""Classes and routines for handling (static) tetrahedral meshes.
"""


import core
import numpy
import tetrahedron as stet


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class TetMesh(core.Container):


    """The main container class for tetrahedral meshes.
    
    To construct a tetrahedral mesh, first create a TetMesh object and
    provide it will all the points, tetrahedrons and triangles. Next,
    these tetrahedrons and triangles are ordered and marked in compartments.
    """
    
    
    def __init__(self, pnts, tets, tris, **params):
        """Initialize a TetMesh object.
        
        Parameters:
            pnts
                An X*3 array of floats, storing the coordinates of node points.
            tets
                A Y * 4 array of integers, storing the corner nodes for each 
                tetrahedron as an index into array nodes.
            tris
                A Z * 3 array of integers, storing the corner nodes for 
                triangles (i.e. special boundaries between or around
                tetrahedrons.
        
        Returns:
            ---
        
        Raises:
            ---
        """
        # Call the initializer for the parent class.
        super(self).__init__()
        
        # Copy the 3 crucial tables.
        assert pnts.shape[1] == 3
        assert pnts.dtype == float
        self._pnts = pnts.copy()
        assert tets.shape[1] == 4
        assert tets.dtype == int
        self._tets = tets.copy()
        assert tris.shape[1] == 3
        assert tris.dtype == int
        self._tris = tris.copy()
        
        # Find minimal and maximal boundary values.
        self.__minx = self._pnts[:,0].min()
        self.__miny = self._pnts[:,1].min()
        self.__minz = self._pnts[:,2].min()
        self.__maxx = self._pnts[:,0].max()
        self.__maxy = self._pnts[:,1].max()
        self.__maxz = self._pnts[:,2].max()
    
    
    def checkConsistency(self):
        """Perform a consistency check on the entire mesh.
        """
        pass
    
    
    def getPnts(self):
        return self._pnts
    
    pnts = property(getPnts)
    
    
    def getTets(self):
        return self._tets
    
    tets = property(getTets)
    
    
    def getTris(self):
        return self._tris
    
    tris = property(getTris)
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getBoundMin(self):
        """Get the minimal coordinate of the rectangular bounding box
        around the entire mesh.
        
        Returns:
            A tuple (x,y,z).
        """
        return ( self.__minx, self.__miny, self.__minz )
    
    
    def getBoundMax(self):
        """Get the maximal coordinate of the rectangular bounding box
        around the entire mesh.
        
        Returns:
            A tuple (x,y,z).
        """
        return ( self.__maxx, self.__maxy, self.__maxz )
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #


    def findCompByPoint(self, p, inclusive = True):
        """Return the compartment(s) to which one or more 3D point(s) belong.
        
        Parameters:
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
                
        Returns:
            A list of references to compartments (or None values for points 
            not inside any compartment). This is always a list, even if the
            parameter was a single tuple.
        
        Raises:
            ---
        """
        pass


    def findCompByTet(self, t):
        """Return the compartment(s) to which one or more tetrahedrons belong.
        
        Parameters:
            t
                The tetrahedron(s) that should be tested. This parameter can
                take multiple forms. The simplest form is a single integer
                index value, indicating a single tetrahedron to be tested.
                A sequence of such indices can also be given.
        
        Returns:
            A list of references to compartments (or None for voxel indices
            that are not inside any compartments). This is always a list,
            even if the parameter was a single integer index value.
        
        Raises:
            ---
        """
        pass


    def findTetByPoint(self, p, inclusive = True):
        """Return the voxel(s) to which one or more 3D point(s) belong.
        
        Parameters:
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

        Returns:
            An array of tetrahedral voxel indices (or -1 for points not 
            inside any tetrahedron).
        
        Raises:
            ---
        """
        pass
        
        
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def allBorderTris(self):
        pass
    
    
    def allBorderTriIndices(self, reorient = True):
        pass
        
        
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class Comp(core.Comp):


    """
    """
    
    
    def __init__(self, id, container, tets):
        """Initialize a tetrahedral mesh compartment.
        """
        # Call parent object initializer.
        super(self).__init__(id, container)
        
        # Store all tetrahedron indices as a set.
        self._tets_idx = set(tets)
        
        # Next, copy the tetrahedrons
        
        # Compute volume.
        vols = stet.vol(self.container.pnts, self.container.tets[self._tets,:])
        self._vol = vols.sum()
        
        # Compute bounds.
        pts = set(self.container.tets[self._tets,0]) \
            | set(self.container.tets[self._tets,1]) \
            | set(self.container.tets[self._tets,2]) \
            | set(self.container.tets[self._tets,3])
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
        """
        pass
    
    
    def checkConsistency(self):
        """Perform a consistency check on the compartment.
        """
        pass
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def getBoundMin(self):
        """Get the minimal coordinate of the rectangular bounding box
        around the compartment.
        """
        pass
    
    
    def getBoundMax(self):
        """Get the maximal coordinate of the rectangular bounding box
        around the compartment.
        """
        pass 
        

    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
    def isTetInside(self, t):
        """Test whether a tetrahedron(s) (specified by index) is inside
        this compartment.
        
        Parameters:
            t
                A tetrahedral index or sequence of tetrahedral indices.
                
        Returns:
            A list of boolean values.
        
        Raises:
            ---
        """
        pass
    
    
    def findTetByPoint(self, p, inclusive):
        """Return the voxel(s) to which one or more 3D point(s) belong.
        
        When a point is not inside any tetrahedron, or when it lies right on 
        the edge between two or more neighbouring tetrahedrons, the method 
        returns -1 for that point.
        
        Parameters:
            p
                The 3-dimensional point(s) that should be tested. This 
                parameter can take a lot of forms. The simplest is a 
                single tuple of real values (x, y, z). A list of such tuples
                can also be given. Alternatively, the method accept a
                N*3 Numpy array object. Internally, the input parameter
                gets temporarily converted into a Numpy array anyway.

        Returns:
            An array of tetrahedral voxel indices (or -1 for points not 
            inside any tetrahedron).
        
        Raises:
            ---
        """
        pass
    
    
    #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
    
    
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


class Patch(core.Patch):


    """
    """
    
    
    def __init__(self, id, container):
        super(self).__init__(id, container)


    def checkConsistency(self):
        """Perform a consistency check on the patch.
        """
        pass


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def fromCDF(cdf):
    pass
    

def toCDF(geom):
    pass


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
