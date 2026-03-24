# cython:language_level=3str
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

"""
This file is the user-interface file for all objects related to STEPS 4. All
objects are directly derived from the corresponding Cython objects.
"""

# The following fixes compatibility issues with mpi4py >= 4.0.0 and Open MPI
# versions that do not support MPI-4 yet.
# See https://github.com/mpi4py/mpi4py/issues/525
# This could be removed once common Open MPI packages support MPI-4
cdef extern from *:
    """
    #include <mpi.h>
    
    #if (MPI_VERSION < 3) && !defined(PyMPI_HAVE_MPI_Message)
    typedef void *PyMPI_MPI_Message;
    #define MPI_Message PyMPI_MPI_Message
    #endif
    
    #if (MPI_VERSION < 4) && !defined(PyMPI_HAVE_MPI_Session)
    typedef void *PyMPI_MPI_Session;
    #define MPI_Session PyMPI_MPI_Session
    #endif
    """

from cython.operator cimport dereference as deref
cimport mpi4py.MPI as MPI

include "cysteps_mpi.pyx"
from steps_dist_tetmesh cimport *
from steps_dist_solver cimport *

from typing import NamedTuple

_py_MembraneResistivity = NamedTuple('MembraneResistivity', [('resistivity', float), ('reversal_potential', float)])

# ======================================================================================================================
# Python bindings to namespace steps::dist
# ======================================================================================================================

cdef class _py_Library(_py__base):
    "Python wrapper class for MPI environment"
    cdef Library *ptrx(self):
        return <Library*> self._ptr

    def __init__(self, MPI.Comm comm):
        """Construct a Library

         Args:
             comm: MPI communicator. Default is COMM_WORLD
         """

        import sys
        cdef int argc = len(sys.argv)
        cdef std.vector[char*] argv
        argv.reserve(argc)
        for arg in sys.argv:
          argv.push_back(arg.encode("utf-8"))
        cdef char** argv_ptr = argv.data()
        self._ptr = new Library(&argc, &argv_ptr, comm.ob_mpi)

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_DistMesh(_py_Geom):
    "Python wrapper class for distributed mesh"
# ----------------------------------------------------------------------------------------------------------------------

    cdef DistMesh *ptrx(self):
        return <DistMesh*> self._ptr

    def __init__(self, _py_Library library, str path, float scale=0):
        """Construct a DistMesh

        Args:
            library: MPI environment
            path: mesh location on the filesystem
            scale: LENGTH scale from the importing mesh to real geometry. e.g. a radius of 10 in the importing file to a radius of 1 micron in STEPS, scale is 1e-7.
        """
        self._ptr = new DistMesh(deref(library.ptrx()), to_std_string(path), scale)

    @staticmethod
    cdef void _setPtr(_py_DistMesh self, _py_DistMesh other):
        self._ptr = other._ptr

    @staticmethod
    def _use_gmsh():
        """Return whether omega_h uses gmsh

        For internal use only

        Syntax::

            _use_gmsh()

        Arguments:
        None

        Return:
        bool

        """
        return DistMesh.use_gmsh()

    def addDiffusionBoundary(self, str name, str comp1, str comp2, triangles=None):
        """
        Add a diffusion boundary between comp1 and comp2 to the mesh.
        Optionally supply a list of triangle ids. Otherwise all the shared triangles are used

        Syntax::

            addDiffusionBoundary(name, comp1, comp2, triangles)

        Arguments:
        str name: name of the boundary
        str comp1: first compartment
        str comp2: second compartment
        List[GO] triangles: triangles on which the diffusion boundary is applied. They must be on the boundary between the two compartments. If this is null, all the shared triangles are taken

        Return:
        None

        """
        cdef std.set[triangle_global_id_t] tris
        if triangles is not None:
            for tri in triangles:
                tris.insert(triangle_global_id_t(tri))
            self.ptrx().addDiffusionBoundary(
                diffusion_boundary_name(to_std_string(name)),
                compartment_id(to_std_string(comp1)),
                compartment_id(to_std_string(comp2)),
                tris
            )
        else:
            self.ptrx().addDiffusionBoundary(
                diffusion_boundary_name(to_std_string(name)),
                compartment_id(to_std_string(comp1)),
                compartment_id(to_std_string(comp2))
            )


    def countTets(self, bool local=False):
        """
        Returns the total number of tetrahedrons in the mesh.

        Syntax::

            countTets()

        Arguments:
        bool local

        Return:
        int
        """
        if local:
            return self.ptrx().num_elems()
        else:
            return self.ptrx().total_num_elems()

    def countTris(self, bool local=False):
        """
        Returns the total number of triangles in the mesh.

        Syntax:

        countTris()

        Arguments:
        bool local

        Return:
        int
        """
        if local:
            return self.ptrx().num_bounds()
        else:
            return self.ptrx().total_num_bounds()

    def countBars(self, bool local=False):
        """
        Returns the total number of bars in the mesh.

        Syntax:

        countBars()

        Arguments:
        bool local

        Return:
        int
        """
        if local:
            return self.ptrx().num_bars()
        else:
            return self.ptrx().total_num_bars()

    def countVertices(self, bool local=False):
        """
        Returns the total number of vertices in the mesh.

        Syntax:

        countVertices()

        Arguments:
        bool local

        Returns:
        uint
        """
        if local:
            return self.ptrx().num_verts()
        else:
            return self.ptrx().total_num_verts()

    def redistributed(self):
        return self.ptrx().redistributed()

    @property
    def num_elems(self):
        """
        The number of elements owned by this process
        """
        return self.ptrx().num_elems()

    @property
    def total_num_elems(self):
        """
        The number of elements in the entire mesh
        """
        return self.ptrx().total_num_elems()

    @property
    def num_bounds(self):
        """
        The number of boundaries owned by this process
        """
        return self.ptrx().num_bounds()

    @property
    def total_num_bounds(self):
        """
        The number of boundaries in the entire mesh
        """
        return self.ptrx().total_num_bounds()

    @property
    def num_verts(self):
        """
        The number of vertices owned by this process
        """
        return self.ptrx().num_verts()

    @property
    def total_num_verts(self):
        """
        The number of vertices in the entire mesh
        """
        return self.ptrx().total_num_verts()

    def getTetComp(self, GO idx, bool local=False):
        """
        Returns a reference to a steps.geom.DistComp object: the compartment which
        tetrahedron with index idx belongs to. Returns None if tetrahedron not
        assigned to a compartment.

        Syntax::

            getTetComp(idx)

        Arguments:
        GO idx
        bool local

        Return:
        steps.geom.DistComp

        """
        cdef DistComp* compPtr;
        if local:
            compPtr = self.ptrx().getTetComp(tetrahedron_local_id_t(<LO>(idx)) )
        else:
            compPtr = self.ptrx().getTetComp(tetrahedron_global_id_t(idx))
        return _py_DistComp.from_ptr(compPtr) if compPtr != NULL else None

    def getTriPatch(self, GO idx, bool local=False):
        """
        Returns a reference to a steps.geom.DistPatch object: the patch which
        triangle with index idx belongs to. Returns None if the triangle is not
        assigned to a compartment.

        Syntax::

            getTriPatch(idx)

        Arguments:
        GO idx
        bool local

        Return:
        steps.geom.DistPatch

        """
        cdef DistPatch* patchPtr;
        if local:
            patchPtr = self.ptrx().getTriPatch(triangle_local_id_t(<LO>(idx)))
        else:
            patchPtr = self.ptrx().getTriPatch(triangle_global_id_t(idx))
        return _py_DistPatch.from_ptr(patchPtr) if patchPtr != NULL else None

    def getTetVol(self, GO idx, bool local=False):
        """
        Returns the volume of the tetrahedron with index idx.

        Syntax::

            getTetVol(idx)

        Arguments:
        GO idx
        bool local

        Return:
        float

        """
        if local:
            return self.ptrx().getTetVol(<tetrahedron_local_id_t>(<LO>(idx)))
        else:
            return self.ptrx().getTetVol(<tetrahedron_global_id_t>(idx))

    def getTriArea(self, GO idx, bool local=False):
        """
        Returns the area of the triangle with index idx.

        Syntax::

            getTriArea(idx)

        Arguments:
        GO idx
        bool local

        Return:
        float

        """
        if local:
            return self.ptrx().getTriArea(<triangle_local_id_t>(<LO>(idx)))
        else:
            return self.ptrx().getTriArea(<triangle_global_id_t>(idx))

    def getSurfTris(self, bool local=False, bool owned=True):
        """
        Returns a list of triangles that form the mesh boundary.
        Support function for steps.utilities.visual.

        Syntax::

            getSurfTris()

        Arguments:
        bool local
        bool owned

        Return:
        list<GO>

        """
        if local:
            return [tl.get() for tl in self.ptrx().getSurfLocalTris(owned)]
        else:
            return [tg.get() for tg in self.ptrx().getSurfTris()]

    def getTet(self, GO idx, bool local=False):
        """
        Returns the tetrahedron with index idx in the container by its four vertex indices.

        Syntax::
            getTet(idx)

        Arguments:
        index_t idx
        bool local

        Return:
        list<index_t, length = 4>

        """
        if local:
            return [vl.get() for vl in self.ptrx().getTet_(<tetrahedron_local_id_t>(<LO>(idx)))]
        else:
            return [vg.get() for vg in self.ptrx().getTet_(<tetrahedron_global_id_t>(idx))]

    def getTri(self, GO idx, bool local=False):
        """
        Returns the triangle with index idx in the container by its three vertex indices.

        Syntax::

            getTri(idx)

        Arguments:
        GO idx
        bool local

        Return:
        list<GO, length = 3>

        """
        if local:
            return [vl.get() for vl in self.ptrx().getTri_(<triangle_local_id_t>(<LO>(idx)))]
        else:
            return [vg.get() for vg in self.ptrx().getTri_(<triangle_global_id_t>(idx))]

    def getBar(self, GO idx, bool local=False):
        """
        Returns the bar with index idx in the container by its two vertex indices.

        Syntax::

            getBar(idx)

        Arguments:
        GO idx
        bool local

        Return:
        list<GO, length = 3>

        """
        if local:
            return [vl.get() for vl in self.ptrx().getBarVertNeighb(<bar_local_id_t>(<LO>(idx)))]
        else:
            return [vg.get() for vg in self.ptrx().getBarVertNeighb(<bar_global_id_t>(idx))]

    def getVertex(self, GO idx, bool local=False):
        """
        Returns the coordinates of vertex with index idx in the container.

        Syntax::

            getVertex(idx)

        Arguments:
        GO idx
        bool local

        Return:
        list<float, length = 3>

        """
        if local:
            return self.ptrx().getVertex(<vertex_local_id_t>(<LO>(idx)))
        else:
            return self.ptrx().getVertex(<vertex_global_id_t>(idx))

    def getTetTetNeighb(self, GO idx, bool local=False, bool owned=True):
        """
        Returns the indices of the four neighbouring tetrahedrons of tetrahedron with index idx.
        An index of UNKNOWN_TET indicates no neighbour (tetrahedron is on the mesh border).

        Syntax::

            getTetTetNeighb(idx)

        Arguments:
        GO idx
        bool local
        bool owned

        Return:
        list<GO, length = 4>

        """
        if local:
            return [tl.get() if tl.valid() else UNKNOWN_TET for tl in self.ptrx().getTetTetNeighb(<tetrahedron_local_id_t>(<LO>(idx)), owned)]
        else:
            return [tg.get() if tg.valid() else UNKNOWN_TET for tg in self.ptrx().getTetTetNeighb(<tetrahedron_global_id_t>(idx))]

    def getTetTriNeighb(self, GO idx, bool local=False):
        """
        Returns the indices of the four neighbouring triangles of tetrahedron with index idx.

        Syntax::

            getTetTriNeighb(idx)

        Arguments:
        GO idx
        bool local

        Return:
        list<GO, length = 4>

        """
        if local:
            return [tl.get() for tl in self.ptrx().getTetTriNeighb(<tetrahedron_local_id_t>(<LO>(idx)))]
        else:
            return [tg.get() for tg in self.ptrx().getTetTriNeighb(<tetrahedron_global_id_t>(idx))]

    def getTriBars(self, GO idx, bool local=False):
        """
        Returns the indices of the bars that comprise the triangle.

        Syntax::

            getTriBars(tidx)

        Arguments:
        GO idx
        bool local

        Return:
        list<GO, length = 3>

        """
        if local:
            return [tl.get() for tl in self.ptrx().getTriBarNeighb(<triangle_local_id_t>(<LO>(idx)))]
        else:
            return [tg.get() for tg in self.ptrx().getTriBarNeighb(<triangle_global_id_t>(idx))]

    def getTriTetNeighb(self, GO idx, bool local=False, bool owned=True):
        """
        Returns the indices of the neighbouring tetrahedrons of triangle with
        index idx. If the triangle is on the mesh boundary, only one tetrahedron is returned

        Syntax::

            getTriTetNeighb(idx)

        Arguments:
        GO idx
        bool local
        bool owned

        Return:
        list<GO>

        """
        if local:
            return [tl.get() for tl in self.ptrx().getTriTetNeighb(<triangle_local_id_t>(<LO>(idx)), owned)]
        else:
            return [tg.get() for tg in self.ptrx().getTriTetNeighb(<triangle_global_id_t>(idx))]

    def getTriTriNeighbs(self, GO idx, bool local=False, bool owned=True):
        """
        Returns the indices of the neighbouring triangles of triangle with
        index idx.

        Syntax::

            getTriTriNeighbs(idx)

        Arguments:
        GO idx
        bool local
        bool owned

        Return:
        list<GO>

        """
        if local:
            return [tl.get() for tl in self.ptrx().getTriTriNeighbs(<triangle_local_id_t>(<LO>(idx)), owned)]
        else:
            return [tg.get() for tg in self.ptrx().getTriTriNeighbs(<triangle_global_id_t>(idx))]

    def getTriTriNeighb(self, GO idx, _py_DistPatch patch, bool local=False, bool owned=True):
        """
        Returns the indices of the neighbouring triangles of triangle with
        index idx within a given patch.

        Syntax::

            getTriTriNeighb(idx, patch)

        Arguments:
        GO idx
        DistPatch patch
        bool local
        bool owned

        Return:
        list<GO>

        """
        if local:
            return [tl.get() for tl in self.ptrx().getTriTriNeighbs(<triangle_local_id_t>(<LO>(idx)), deref(patch.ptrx()), owned)]
        else:
            return [tg.get() for tg in self.ptrx().getTriTriNeighbs(<triangle_global_id_t>(idx), deref(patch.ptrx()))]


    def getBoundMin(self, bool local=False):
        """
        Returns the minimal Cartesian coordinate of the rectangular bounding box
        of the mesh. 
        
        A bounding box is the smallest box that encompasses all the vertices. It is defined by 2 extreme vertices: 
        min and max. 

        Syntax::

            getBoundMin()

        Arguments:
        bool local

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMin(local)

    def getBoundMax(self, bool local=False):
        """
        Returns the maximal Cartesian coordinate of the rectangular bounding box
        of the mesh.
        
        A bounding box is the smallest box that encompasses all the vertices. It is defined by 2 extreme vertices: 
        min and max. 

        Syntax::

            getBoundMax()

        Arguments:
        bool local

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMax(local)

    def getTetBarycenter(self, GO idx, bool local=False):
        """
        Returns the barycenter of the tetrahedron with index idx.

        Syntax::

            getTetBarycenter(idx)

        Arguments:
        GO idx
        bool local

        Return:
        list<float, length = 3>

        """
        if local:
            return self.ptrx().getTetBarycenter(<tetrahedron_local_id_t>(<LO>(idx)))
        else:
            return self.ptrx().getTetBarycenter(<tetrahedron_global_id_t>(idx))

    def getTriBarycenter(self, GO idx, bool local=False):
        """
        Returns the Cartesian coordinates of the barycenter of triangle with index idx.

        Syntax::

            getTriBarycenter(idx)

        Arguments:
        GO idx
        bool local

        Return:
        list<float, length = 3>

        """
        if local:
            return self.ptrx().getTriBarycenter(<triangle_local_id_t>(<LO>(idx)))
        else:
            return self.ptrx().getTriBarycenter(<triangle_global_id_t>(idx))

    def findTetByPoint(self, std.vector[double] p, bool local=False):
        """
        Returns the index of the tetrahedron which encompasses a given point
        p (given in Cartesian coordinates x,y,z). Returns UNKNOWN_TET if p is a position
        outside the mesh.

        Syntax::

            findTetByPoint(p)

        Arguments:
        list<float, length = 3> p

        Return:
        GO

        """
        cdef tetrahedron_global_id_t globalInd
        cdef tetrahedron_local_id_t localInd
        if local:
            localInd = self.ptrx().findLocalTetByPoint(p)
            return localInd.get() if localInd.valid() else UNKNOWN_TET
        else:
            globalInd = self.ptrx().findTetByPoint(p)
            return globalInd.get() if globalInd.valid() else UNKNOWN_TET

    def findLocalTetByPointLinear(self, std.vector[double] p):
        """
        Returns the local index of the tetrahedron which encompasses a given point
        p (given in Cartesian coordinates x,y,z). Returns UNKNOWN_TET if p is a position
        outside the mesh. Linear search.

        Syntax::

            findLocalTetByPointLinear(p)

        Arguments:
        list<double, length = 3> p

        Return:
        LO

        """
        cdef tetrahedron_local_id_t localInd
        localInd = self.ptrx().findLocalTetByPointLinear(p)
        return localInd.get() if localInd.valid() else UNKNOWN_TET

    def findLocalTetByPointWalk(self, std.vector[double] p):
        """
        Returns the local index of the tetrahedron which encompasses a given point
        p (given in Cartesian coordinates x,y,z). Returns UNKNOWN_TET if p is a position
        outside the mesh. A* search with random initial seeding and restarts.

        Syntax::

            findLocalTetByPointWalk(p)

        Arguments:
        list<double, length = 3> p

        Return:
        LO

        """
        cdef tetrahedron_local_id_t localInd
        localInd = self.ptrx().findLocalTetByPointWalk(p)
        return localInd.get() if localInd.valid() else UNKNOWN_TET

    def isPointInTet(self, std.vector[double] p, GO tidx, bool local=False):
        """
        Check if point belongs to the tetrahedron or not

        Syntax::

            isPointInTet(p, tidx)

        Arguments:
        list<float, length = 3> p
        int tetrahedron tidx
        bool local

        Return:
        bool
        """
        if local:
            return self.ptrx().isPointInTet(p, tetrahedron_local_id_t(<LO>(tidx)))
        else:
            return self.ptrx().isPointInTet(p, tetrahedron_global_id_t(tidx))

    def getTaggedTetrahedrons(self, str tag, bool local=False, bool owned=True):
        """
        Returns the global indexes of all tetrahedrons corresponding to a tag

        Syntax::

            getTaggedTetrahedrons(tag)

        Arguments:
        str tag
        bool local
        bool owned

        Return:
        List[int]
        """
        if local:
            return [tl.get() for tl in self.ptrx().getTaggedLocalTetrahedrons(compartment_id(to_std_string(tag)), owned) if tl.valid()]
        else:
            return [tg.get() for tg in self.ptrx().getTaggedTetrahedrons(compartment_id(to_std_string(tag)))]

    def getTaggedTriangles(self, str tag, bool local=False, bool owned=True):
        """
        Returns the global indexes of all triangles corresponding to a tag

        Syntax::

            getTaggedTriangles(tag)

        Arguments:
        str tag
        bool local
        bool owned

        Return:
        List[int]
        """
        if local:
            return [tl.get() for tl in self.ptrx().getTaggedLocalTriangles(patch_id(to_std_string(tag)), owned) if tl.valid()]
        else:
            return [tg.get() for tg in self.ptrx().getTaggedTriangles(patch_id(to_std_string(tag)))]

    def getTaggedVertices(self, str tag, bool local=False, bool owned=True):
        """
        Returns the global indexes of all vertices corresponding to a tag

        Syntax::

            getTaggedVertices(tag)

        Arguments:
        str tag
        bool local
        bool owned

        Return:
        List[int]
        """
        if local:
            return [tl.get() for tl in self.ptrx().getTaggedLocalVertices(vertgroup_id(to_std_string(tag)), owned) if tl.valid()]
        else:
            return [tg.get() for tg in self.ptrx().getTaggedVertices(vertgroup_id(to_std_string(tag)))]

    def getTags(self, int dim):
        """Get the list of physical groups with dimension dim

        Syntax::

            getTags(dim)

        Arguments:
        int dim

        Return:
        List[str]
        """
        return [from_std_string(s) for s in self.ptrx().getTags(dim)]

    def getTetLocalIndex(self, GO idx, bool owned=True):
        """Return the local index of tetrahedron with global index idx

        Return None if the tetrahedron does not exist locally.

        Syntax::

            getTetLocalIndex(idx)

        Arguments:
        int idx
        bool owned

        Return:
        int
        """
        cdef tetrahedron_local_id_t ind = self.ptrx().getLocalIndex(tetrahedron_global_id_t(idx), owned)
        return ind.get() if ind.valid() else None

    def getTriLocalIndex(self, GO idx, bool owned=True):
        """Return the local index of triangle with global index idx

        Return None if the triangle does not exist locally.

        Syntax::

            getTriLocalIndex(idx)

        Arguments:
        int idx
        bool owned

        Return:
        int
        """
        cdef triangle_local_id_t ind = self.ptrx().getLocalIndex(triangle_global_id_t(idx), owned)
        return ind.get() if ind.valid() else None

    def getBarLocalIndex(self, GO idx, bool owned=True):
        """Return the local index of bar with global index idx

        Return None if the bar does not exist locally.

        Syntax::

            getBarLocalIndex(idx)

        Arguments:
        int idx
        bool owned

        Return:
        int
        """
        cdef bar_local_id_t ind = self.ptrx().getLocalIndex(bar_global_id_t(idx), owned)
        return ind.get() if ind.valid() else None

    def getVertLocalIndex(self, GO idx, bool owned=True):
        """Return the local index of vertex with global index idx

        Return None if the vertex does not exist locally.

        Syntax::

            getVertLocalIndex(idx)

        Arguments:
        int idx
        bool owned

        Return:
        int
        """
        cdef vertex_local_id_t ind = self.ptrx().getLocalIndex(vertex_global_id_t(idx), owned)
        return ind.get() if ind.valid() else None

    def getTetGlobalIndex(self, LO idx):
        """
        Return the global index of tetrahedron with local index idx

        Syntax::

            getTetGlobalIndex(idx)

        Arguments:
        int idx

        Return:
        int
        """
        return self.ptrx().getGlobalIndex(tetrahedron_local_id_t(idx)).get()

    def getTriGlobalIndex(self, LO idx):
        """
        Return the global index of triangle with local index idx

        Syntax::

            getTriGlobalIndex(idx)

        Arguments:
        int idx

        Return:
        int
        """
        return self.ptrx().getGlobalIndex(triangle_local_id_t(idx)).get()

    def getBarGlobalIndex(self, LO idx):
        """
        Return the global index of bar with local index idx

        Syntax::

            getBarGlobalIndex(idx)

        Arguments:
        int idx

        Return:
        int
        """
        return self.ptrx().getGlobalIndex(bar_local_id_t(idx)).get()

    def getVertGlobalIndex(self, LO idx):
        """
        Return the global index of vertex with local index idx

        Syntax::

            getTetGlobalIndex(idx)

        Arguments:
        int idx

        Return:
        int
        """
        return self.ptrx().getGlobalIndex(vertex_local_id_t(idx)).get()

    def getAllTetIndices(self, bool local=False, bool owned=True):
        """
        Returns a list of all tetrahedrons in the mesh.

        Syntax::

            getAllTetIndices()

        Arguments:
        bool local
        bool owned

        Return:
        list<int>
        """
        if local:
            return [tl.get() for tl in self.ptrx().getLocalTetIndices(owned) if tl.valid()]
        else:
            return [tg.get() for tg in self.ptrx().getAllTetIndices()]

    def getAllTriIndices(self, bool local=False, bool owned=True):
        """
        Returns a list of all triangles in the mesh.

        Syntax::

            getAllTriIndices()

        Arguments:
        bool local
        bool owned

        Return:
        list<int>
        """
        if local:
            return [tl.get() for tl in self.ptrx().getLocalTriIndices(owned) if tl.valid()]
        else:
            return [tg.get() for tg in self.ptrx().getAllTriIndices()]

    def getAllBarIndices(self, bool local=False, bool owned=True):
        """
        Returns a list of all bars in the mesh.

        Syntax::

            getAllBarIndices()

        Arguments:
        bool local
        bool owned

        Return:
        list<int>
        """
        if local:
            return [tl.get() for tl in self.ptrx().getLocalBarIndices(owned) if tl.valid()]
        else:
            return [tg.get() for tg in self.ptrx().getAllBarIndices()]

    def getAllVertIndices(self, bool local=False, bool owned=True):
        """
        Returns a list of all vertices in the mesh.

        Syntax::

            getAllVertIndices()

        Arguments:
        bool local
        bool owned

        Return:
        list<int>
        """
        if local:
            return [vl.get() for vl in self.ptrx().getLocalVertIndices(owned) if vl.valid()]
        else:
            return [vg.get() for vg in self.ptrx().getAllVertIndices()]

    def getMeshVolume(self, bool local=False):
        """
        Returns the total volume of the mesh.

        Syntax::

            getMeshVolume()

        Arguments:
        bool local

        Return:
        float

        """
        if local:
            return self.ptrx().local_measure(compartment_id(to_std_string('__MESH__')))
        else:
            return self.ptrx().total_measure(compartment_id(to_std_string('__MESH__')))

    def intersect(self, double[:, :] points, int sampling=-1, bool local=True):
        """
        Computes the intersection of line segment(s) given the vertices with the current mesh

        Args:
            points: A 2-D NumPy array (/memview) of points in the 3D space, 
                    where each element contains the 3 point coordinates
            int sampling: any value --> deterministic method (montecarlo not implemented for STEPS4)
            bool local

        Returns:
            A list of lists of tuples representing the intersected tetrahedrons, one element per line segment.
            Each tuple is made of 2 elements, a tetrahedron local identifier, and its respective intersection ratio.
        """
        if (points.strides[0] != 24 or points.strides[1] != 8):
            raise Exception("Wrong memory layout for points, np array should be [pts,3] and row major")
        if local:
            return [[(t.first.get(), t.second) for t in row] for row in self.ptrx().localIntersect(&points[0][0], points.shape[0], sampling)]
        else:
            return [[(t.first.get(), t.second) for t in row] for row in self.ptrx().intersect(&points[0][0], points.shape[0], sampling)]
        

    def intersectIndependentSegments(self, double[:, :] points, int sampling=-1, bool local=True):
        """
        Similar to the intersect method but here we deal with independent segments, i.e.
        every two points we have a segment not related to previous or following ones.
        E.g. seg0 = (points[0], points[1]), seg1 = (points[2], points[3]), etc.

        Args:
            points: A 2-D NumPy array (/memview) of points in the 3D space, 
                    where each element contains the 3 point coordinates
            int sampling: any value --> deterministic method (montecarlo not implemented for STEPS4)
            bool local

        Returns:
            A list of lists of tuples representing the intersected tetrahedrons, one element per line segment.
            Each tuple is made of 2 elements, a tetrahedron local identifier, and its respective intersection ratio.
        """
        if (points.strides[0] != 24 or points.strides[1] != 8):
            raise Exception("Wrong memory layout for points, np array should be [pts,3] and row major")
        if local:
            return [[(t.first.get(), t.second) for t in row] for row in self.ptrx().localIntersectIndependentSegments(&points[0][0], points.shape[0], sampling)]
        else:
            return [[(t.first.get(), t.second) for t in row] for row in self.ptrx().intersectIndependentSegments(&points[0][0], points.shape[0], sampling)]

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_DistComp(_py_Comp):
    "Python wrapper class for distributed compartment"
# ----------------------------------------------------------------------------------------------------------------------
    cdef DistComp *ptrx(self):
        return <DistComp*> self._ptr

    def __init__(self, str id, _py_DistMesh mesh, tets=None, physical_tag=None, conductivity=0, local=False):
        """Construct a DistComp

        The compartment will be composed of all tetrahedrons that are tagged with the physical tag associated to the
        compartment name.

        Optionally provide a volume conductivity for the compartment.

        Args:
            string id
            steps.geom.DistMesh mesh
            List[int] tets
            int physical_tag
            float conductivity
        """
        cdef std.vector[tetrahedron_global_id_t] tet_indices
        cdef std.vector[tetrahedron_local_id_t] local_tet_indices
        if physical_tag is None:
            if tets is None:
                self._ptr = new DistComp(
                    compartment_name(to_std_string(id)),
                    mesh.ptrx()[0],
                    conductivity
                )
            elif local:
                local_tet_indices.reserve(len(tets))
                for ind in tets:
                    local_tet_indices.push_back(tetrahedron_local_id_t(ind))
                self._ptr = new DistComp(
                    compartment_name(to_std_string(id)),
                    mesh.ptrx()[0],
                    <std.vector[tetrahedron_local_id_t]> local_tet_indices,
                    <double> conductivity
                )
            else:
                tet_indices.reserve(len(tets))
                for ind in tets:
                    tet_indices.push_back(tetrahedron_global_id_t(ind))
                self._ptr = new DistComp(
                    compartment_name(to_std_string(id)),
                    mesh.ptrx()[0],
                    <std.vector[tetrahedron_global_id_t]> tet_indices,
                    <double> conductivity
                )
        else:
            self._ptr = new DistComp(
                compartment_name(to_std_string(id)),
                mesh.ptrx()[0],
                <compartment_physical_tag> compartment_physical_tag(physical_tag),
                <double> conductivity
            )

    def getAllTetIndices(self, bool local=False, bool owned=True):
        """
        Returns a list of all tetrahedrons assigned to the compartment.

        Syntax::

            getAllTetIndices()

        Arguments:
        bool local
        bool owned

        Return:
        list<int>

        """
        if local:
            return [tl.get() for tl in self.ptrx().getLocalTetIndices(owned)]
        else:
            return [tg.get() for tg in self.ptrx().getAllTetIndices()]

    def getSurfTris(self, bool local=False):
        """
        Returns a list of triangles that form the compartment boundary.

        Syntax::

            getSurfTris()

        Arguments:
        bool local

        Return:
        list<GO>

        """
        if local:
            return [tl.get() for tl in self.ptrx().getSurfLocalTris()]
        else:
            return [tg.get() for tg in self.ptrx().getSurfTris()]

    def getConductivity(self):
        """
        Return the conductivity of the compartment

        Syntax::

            getConductivity()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getConductivity()

    def setConductivity(self, conductivity):
        """
        Set the conductivity of the compartment

        Syntax::

            getConductivity(conductivity)

        Arguments:
        float conductivity

        Return:
        None

        """
        return self.ptrx().setConductivity(conductivity)

    def getVol(self, bool local=False):
        """
        Get the volume of the compartment (in m^3).

        Syntax::

            getVol()

        Arguments:
        bool local

        Return:
        float

        """
        if local:
            return self.ptrx().getOwnedVol()
        else:
            return self.ptrx().getTotalVol()

    def getBoundMin(self, bool local=False):
        """
        Returns the minimal Cartesian coordinate of the rectangular bounding box
        of the compartment. 
        
        A bounding box is the smallest box that encompasses all the vertices. It is defined by 2 extreme vertices: 
        min and max. 

        Syntax::

            getBoundMin()

        Arguments:
        bool local

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMin(local)

    def getBoundMax(self, bool local=False):
        """
        Returns the maximal Cartesian coordinate of the rectangular bounding box
        of the compartment.
        
        A bounding box is the smallest box that encompasses all the vertices. It is defined by 2 extreme vertices: 
        min and max. 

        Syntax::

            getBoundMax()

        Arguments:
        bool local

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMax(local)


    @staticmethod
    cdef _py_DistComp from_ptr(DistComp *ptr):
        if (ptr == NULL):
            return None
        cdef _py_DistComp obj = _py_DistComp.__new__(_py_DistComp)
        obj._ptr = ptr
        return obj

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_DistPatch(_py_Patch):
    "Python wrapper class for distributed patch"
# ----------------------------------------------------------------------------------------------------------------------
    cdef DistPatch *ptrx(self):
        return <DistPatch*> self._ptr

    def __init__(self, str id, _py_DistMesh mesh, tris=None, _py_DistComp icomp=None, _py_DistComp ocomp=None, physical_tag=None, local=False):
        """Construct a DistPatch

        The patch will be composed of all triangles that are tagged with the physical tag associated to the
        patch name.

        Args:
            string id
            steps.geom.DistMesh container
            List[int] tris
            steps.geom.DistComp icomp
            steps.geom.DistComp ocomp
            int physical_tag
        """
        assert icomp is not None
        cdef std.vector[triangle_global_id_t] tri_indices
        cdef std.vector[triangle_local_id_t] local_tri_indices
        if physical_tag is None:
            if tris is None:
                self._ptr = new DistPatch(
                    patch_name(to_std_string(id)),
                    mesh.ptrx()[0],
                    deref(icomp.ptrx()),
                    ocomp.ptrx() if ocomp is not None else NULL
                )
            elif local:
                local_tri_indices.reserve(len(tris))
                for ind in tris:
                    local_tri_indices.push_back(triangle_local_id_t(ind))
                self._ptr = new DistPatch(
                    patch_name(to_std_string(id)),
                    mesh.ptrx()[0],
                    local_tri_indices,
                    deref(icomp.ptrx()),
                    ocomp.ptrx() if ocomp is not None else NULL
                )
            else:
                tri_indices.reserve(len(tris))
                for ind in tris:
                    tri_indices.push_back(triangle_global_id_t(ind))
                self._ptr = new DistPatch(
                    patch_name(to_std_string(id)),
                    mesh.ptrx()[0],
                    tri_indices,
                    deref(icomp.ptrx()),
                    ocomp.ptrx() if ocomp is not None else NULL
                )
        else:
            self._ptr = new DistPatch(
                patch_name(to_std_string(id)),
                mesh.ptrx()[0],
                patch_physical_tag(physical_tag),
                deref(icomp.ptrx()),
                ocomp.ptrx() if ocomp is not None else NULL
            )

    def getAllTriIndices(self, bool local=False, bool owned=True):
        """
        Returns a list of all triangles assigned to the compartment.

        Syntax::

            getAllTriIndices()

        Arguments:
        bool local
        bool owned

        Return:
        list<int>

        """
        if local:
            return [tl.get() for tl in self.ptrx().getLocalTriIndices(owned)]
        else:
            return [tg.get() for tg in self.ptrx().getAllTriIndices()]

    def getArea(self, bool local=False):
        """
        Get the area of the patch (in m^2).

        Syntax::

            getArea()

        Arguments:
        bool local

        Return:
        float

        """
        if local:
            return self.ptrx().getOwnedArea()
        else:
            return self.ptrx().getTotalArea()

    def getBoundMin(self, bool local=False):
        """
        Returns the minimal Cartesian coordinate of the rectangular bounding box
        of the patch. 
        
        A bounding box is the smallest box that encompasses all the vertices. It is defined by 2 extreme vertices: 
        min and max. 

        Syntax::

            getBoundMin()

        Arguments:
        bool local

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMin(local)

    def getBoundMax(self, bool local=False):
        """
        Returns the maximal Cartesian coordinate of the rectangular bounding box
        of the patch.
        
        A bounding box is the smallest box that encompasses all the vertices. It is defined by 2 extreme vertices: 
        min and max. 

        Syntax::

            getBoundMax()

        Arguments:
        bool local

        Return:
        list<float, length = 3>

        """
        return self.ptrx().getBoundMax(local)

    @staticmethod
    cdef _py_DistPatch from_ptr(DistPatch *ptr):
        if (ptr == NULL):
            return None
        cdef _py_DistPatch obj = _py_DistPatch.__new__(_py_DistPatch)
        obj._ptr = ptr
        return obj

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_DistMemb(_py_Memb):
    "Python wrapper class for distributed membrane"
# ----------------------------------------------------------------------------------------------------------------------
    cdef DistMemb *ptrx(self):
        return <DistMemb*> self._ptr

    def __init__(self, str id, _py_DistMesh mesh, patches, double capacitance=0):
        """Construct a DistMemb

        Args:
            string id
            steps.geom.DistMesh mesh
            List[DistPatch] patches
            float capacitance
        """
        cdef std.set[patch_id] ptchs
        for p in patches:
            ptchs.insert(patch_id(to_std_string(p.getID())))
        self._ptr = new DistMemb(membrane_id(to_std_string(id)), mesh.ptrx()[0], ptchs, capacitance)

    def getCapacitance(self):
        """
        Return the capacitance of the compartment

        Syntax::

            getCapacitance()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getCapacitance()

    def setCapacitance(self, capacitance):
        """
        Set the capacitance of the compartment

        Syntax::

            setCapacitance(capacitance)

        Arguments:
        float capacitance

        Return:
        None

        """
        return self.ptrx().setCapacitance(capacitance)

cimport steps_dist_solver
from steps_dist_solver cimport Simulation

cdef class _py_SSAMethod:
    SSA = 0
    RSSA = 1
    RLEAPING = 2

cdef class _py_SearchMethod:
    DIRECT = 0
    GIBSON_BRUCK = 1
    RLEAPING = 2

cdef class _py_DiffMethod:
    CONSTANT_DT = 0
    TAU_LEAPING_DT = 1

cdef class _py_DistributionMethod:
    UNIFORM = 0
    MULTINOMIAL = 1

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_DistTetOpSplitP(_py__base):
    """Bindings for MPI DistTetOpSplitP"""
# ----------------------------------------------------------------------------------------------------------------------

    cdef unique_ptr[Simulation] _uniqueptr

    cdef Simulation *ptrx(self):
        return <Simulation*> self._ptr

    def __init__(self, _py_Model model, _py_DistMesh mesh, _py_RNG rng, SSAMethod=_py_SSAMethod.SSA,
            searchMethod=_py_SearchMethod.GIBSON_BRUCK, diffMethod=_py_DiffMethod.CONSTANT_DT, bool indepKProcs=False, bool isEfield=True):
        """
        Construction::

            sim = steps.solver.DistTetOpSplit(model, mesh, rng)

        Create a distributed spatial stochastic solver based on operator splitting, that is that reaction events are
        partitioned and diffusion is approximated. Keyword parameters SSAMethod and searchMethod respectively set the
        SSA method (SSA or RSSA) and the next event search method (DIRECT or GIBSON_BRUCK).
        Keyword parameter diffMethod sets the diffusion method.

        Arguments:
        steps.model.Model model
        steps.geom.DistMesh mesh
        steps.rng.RNG rng
        steps.sim.SSAMethod SSAMethod
        steps.sim.NextEventSearchMethod searchMethod
        steps.sim.DiffusionMethod diffMethod
        bool indepKProcs
        bool isEfield

        """
        if model == None:
            raise TypeError('The Model object is empty.')
        if mesh == None:
            raise TypeError('The Mesh object is empty.')
        if rng == None:
            raise TypeError('The RNG object is empty.')

        self._uniqueptr = GetSimulation(deref(model.ptr()), deref(mesh.ptrx()), rng.ptr(), SSAMethod, searchMethod, diffMethod, indepKProcs, isEfield)
        self._ptr = self._uniqueptr.get()

    def getSolverName(self):
        """
        Returns a string of the solver's name.

        Syntax::

            getSolverName()

        Arguments:
        None

        Return:
        string

        """
        return from_std_string(self.ptrx().getSolverName())

    def getReacExtent(self, bool local=False):
        """
        Return the number of reaction events that have happened in the simulation.

        if all processes call this function, it will return the accumulated
        result across all processes. It can also be called in individual process with
        the local argument set to true, in which case it returns the local result of this process.

        By default it is called globally and return the accumulated result.

        Syntax::

            getReacExtent(local)

        Arguments:
        bool local (default = False)

        Return:
        index_t
        """
        return self.ptrx().getReacExtent(local)

    def getDiffExtent(self, bool local=False):
        """
        Return the number of diffusion events that have happened in the simulation.

        if all processes call this function, it will return the accumulated
        result accross all processes. It can also be called in individual process with
        the local argument set to true, in which case it returns the local result of this process.

        By default it is called globally and return the accumlated result.

        Syntax::

            getDiffExtent(local)

        Arguments:
        bool local (default = False)

        Return:
        index_t
        """
        return self.ptrx().getDiffExtent(local)

    def getEFieldTime(self, ):
        """
        Return the accumulated EField run time of the process.

        Syntax::

            getEFieldTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getEFieldTime()

    def getRDTime(self, ):
        """
        Return the accumulated reaction-diffusion run time of the process.

        Syntax::

            getRDTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getRDTime()

    def getDiffusionTime(self, ):
        """
        Return the accumulated diffusion run time of the process.

        Syntax::

            getDiffusionTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getDiffusionTime()

    def getReactionTime(self, ):
        """
        Return the accumulated reaction run time of the process.

        Syntax::

            getReactionTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getReactionTime()

    def getReactionDebugInfo(self, bool local=False):
        """
        Return debug information for the reaction operator

        If all processes call this function, it will return the accumulated
        result across all processes. It can also be called in individual process with
        the local argument set to true, in which case it returns the local result of this process.

        By default it is called globally and returns the accumulated result.

        Syntax::

            getReactionDebugInfo(local)

        Arguments:
        bool local (default = False)

        Return:
        Dict[str, float]
        """
        return {from_std_string(pair.first): pair.second for pair in self.ptrx().getReactionDebugInfo(local)}

    def getDiffusionDebugInfo(self, bool local=False):
        """
        Return debug information for the diffusion operator

        If all processes call this function, it will return the accumulated
        result across all processes. It can also be called in individual process with
        the local argument set to true, in which case it returns the local result of this process.

        By default it is called globally and returns the accumulated result.

        Syntax::

            getDiffusionDebugInfo(local)

        Arguments:
        bool local (default = False)

        Return:
        Dict[str, float]
        """
        return {from_std_string(pair.first): pair.second for pair in self.ptrx().getDiffusionDebugInfo(local)}

    def getCompSpecCount(self, str comp, str spec):
        """
        Returns the number of molecules of a species with identifier string spec 
        in compartment with identifier string comp.

        In a mesh-based simulation this is the combined count from 
        all tetrahedral elements in the compartment.

        Syntax::
            
            getCompSpecCount(comp, spec)
            
        Arguments:
        string comp
        string spec

        Return:
        float

        """
        return self.ptrx().getCompSpecCount(compartment_id(to_std_string(comp)), species_name(to_std_string(spec)))

    def getCompSpecConc(self, str comp, str spec):
        """
        Returns the concentration (in Molar units) of species with identifier string spec 
        in compartment with identifier string comp.

        Note: in a mesh-based simulation this is calculated from the combined 
        number of molecules from all tetrahedral elements in the compartment and the total 
        volume of the tetrahedrons.

        Syntax::
            
            getCompSpecConc(comp, spec)
            
        Arguments:
        string comp
        string spec

        Return:
        float

        """
        return self.ptrx().getCompSpecConc(compartment_id(to_std_string(comp)), species_name(to_std_string(spec)))

    def setCompSpecCount(self, str comp, str spec, double n, distributionMethod=_py_DistributionMethod.UNIFORM):
        """
        Set the number of molecules of a species with identifier string spec 
        in compartment with identifier string comp.

        In a mesh-based simulation this is the combined count from 
        all tetrahedral elements in the compartment.

        The distributing is weighted with the volume fraction V_tet/V_tot: bigger elements get
        a higher amount of molecules.

        distributionMethod=UNIFORM the distribution is deterministic (apart from roundings) and the number of
        molecules per element is n*V_tet/V_tot.

        distributionMethod=MULTINOMIAL the distribution is multinomial and the probability
        of putting an element in a tet is V_tet/V_tot

        Syntax::
            
            setCompSpecCount(comp, spec, n, distributionMethod)
            
        Arguments:
        string comp
        string spec
        int n
        DistributionMethod distributionMethod

        Return:
        None

        """
        self.ptrx().setCompSpecCount(compartment_id(to_std_string(comp)), species_name(to_std_string(spec)), n, distributionMethod)

    def setCompSpecConc(self, str comp, str spec, double conc, distributionMethod=_py_DistributionMethod.UNIFORM):
        """
        Sets the concentration (in Molar units) of species with identifier string spec 
        in compartment with identifier string comp to conc. In a discrete solver the 
        continuous concentration is converted to a discrete number of 
        molecules.

        In a mesh-based simulation this is the combined count from 
        all tetrahedral elements in the compartment.

        The distributing is weighted with the volume fraction V_tet/V_tot: bigger elements get
        a higher amount of molecules.

        distributionMethod=UNIFORM the distribution is deterministic (apart from roundings) and the number of
        molecules per element is n*V_tet/V_tot.

        distributionMethod=MULTINOMIAL the distribution is multinomial and the probability
        of putting an element in a tet is V_tet/V_tot

        Syntax::

            setCompSpecConc(comp, spec, conc, distributionMethod)
            
        Arguments:
        string comp
        string spec
        float conc
        DistributionMethod distributionMethod

        Return:
        None

        """
        self.ptrx().setCompSpecConc(compartment_id(to_std_string(comp)), species_name(to_std_string(spec)), conc, distributionMethod)

    def getCompSpecClamped(self, str comp, str spec):
        """
        Returns whether species with identifier string spec is clamped
        in the compartment with identifier string comp.

        Syntax::

            getCompSpecClamped(comp, spec)

        Arguments:
        string comp
        string spec

        Return:
        bool

        """
        return self.ptrx().getCompSpecClamped(compartment_id(to_std_string(comp)), species_name(to_std_string(spec)))

    def setCompSpecClamped(self, str comp, str spec, bool clamped):
        """
        Sets whether species with identifier string spec is clamped
        in the compartment with identifier string comp.

        Syntax::

            setCompSpecClamped(comp, spec, clamped)

        Arguments:
        string comp
        string spec
        bool clamped

        Return:
        None

        """
        self.ptrx().setCompSpecClamped(compartment_id(to_std_string(comp)), species_name(to_std_string(spec)), clamped)

    def getCompComplexCount(self, str c, str complex, f):
        """
        Returns the number of molecules of a complex with identifier string complex
        matching filter filt in compartment with identifier string comp.

        In a mesh-based simulation this is the combined count from 
        all tetrahedral elements in the compartment.

        Syntax::
            
            getCompComplexCount(comp, complex, filt)
            
        Arguments:
        string comp
        string complex
        list[uint[:]] filt

        Return:
        float

        """
        return self.ptrx().getCompComplexCount(compartment_id(to_std_string(c)), complex_name(to_std_string(complex)), _get_filters(f))

    def setCompComplexCount(self, str c, str complex, i, double n, distributionMethod=_py_DistributionMethod.UNIFORM):
        """
        Set the number of molecules of a complex with identifier string complex
        in state init in compartment with identifier string comp.

        In a mesh-based simulation this is the combined count from 
        all tetrahedral elements in the compartment.

        Syntax::
            
            setCompComplexCount(comp, complex, init, nspec)
            
        Arguments:
        string comp
        string complex
        uint[:] init
        int nspec

        Return:
        None

        """
        self.ptrx().setCompComplexCount(compartment_id(to_std_string(c)), complex_name(to_std_string(complex)), _get_filters(i), n, distributionMethod)

    def getCompComplexSUSCount(self, str c, str complex, f, uint m):
        """
        Returns the number of subunits in state m of a complex with identifier string complex
        matching filter filt in compartment with identifier string comp.

        In a mesh-based simulation this is the combined count from 
        all tetrahedral elements in the compartment.

        Syntax::
            
            getCompComplexSUSCount(comp, complex, filt, m)
            
        Arguments:
        string comp
        string complex
        list[uint[:]] filt
        uint m

        Return:
        float

        """
        return self.ptrx().getCompComplexSUSCount(compartment_id(to_std_string(c)), complex_name(to_std_string(complex)), _get_filters(f), complex_substate_id(m))

    def getCompReacK(self, str comp, str reac):
        """
        Gets the macroscopic reaction constant of volume reaction with identifier 
        string reac in compartment with identifier string comp to kf. The unit of the reaction 
        constant depends on the order of the reaction. 

        Note: This method gets the currently set value for the patch,
        individual triangles in the patch might have different values.

        Syntax::

            getCompReacK(comp, reac)

        Arguments:
        string comp
        string reac

        Return:
        float

        """
        return self.ptrx().getCompReacK(compartment_id(to_std_string(comp)), reaction_id(to_std_string(reac)))

    def setCompReacK(self, str comp, str reac, double kf):
        """
        Sets the macroscopic reaction constant of volume reaction with identifier 
        string reac in compartment with identifier string comp to kf. The unit of the reaction 
        constant depends on the order of the reaction. 

        Note: This method sets the surface reaction constant in all triangular elements of the patch to kf.

        Note: The default value still comes from the model description, so calling 
        reset() will return the surface reaction constant to that value.

        Syntax::

            setCompReacK(comp, reac, kf)

        Arguments:
        string comp
        string reac
        float kf

        Return:
        None

        """
        self.ptrx().setCompReacK(compartment_id(to_std_string(comp)), reaction_id(to_std_string(reac)), kf)

    def getCompReacExtent(self, str comp, str reac):
        """
        Gets the extent of volume reaction with identifier string reac in compartment with identifier string comp.

        Syntax::

            getCompReacExtent(comp, reac)

        Arguments:
        string comp
        string reac

        Return:
        uint

        """
        return self.ptrx().getCompReacExtent(compartment_id(to_std_string(comp)), reaction_id(to_std_string(reac)))

    def getCompComplexReacExtent(self, str comp, str reac):
        """
        Gets the extent of complex volume reaction with identifier string reac in compartment with identifier string comp.

        Syntax::

            getCompComplexReacExtent(comp, reac)

        Arguments:
        string comp
        string reac

        Return:
        uint

        """
        return self.ptrx().getCompComplexReacExtent(compartment_id(to_std_string(comp)), complex_reaction_id(to_std_string(reac)))

    def getCompDiffD(self, str c, str d):
        """
        Returns the diffusion constant of diffusion rule with identifier string diff
        in compartment with identifier string comp. This constant is in units m^2/s.

        The value for the compartment is
        returned, although individual or groups of tetrahedral elements may have different
        values (set with setTetDiffD).

        Syntax::

            getCompDiffD(comp, diff)

        Arguments:
        string comp
        string diff

        Return:
        float

        """
        return self.ptrx().getCompDiffD(compartment_id(to_std_string(c)), diffusion_id(to_std_string(d)))

    def setCompDiffD(self, str c, str d, double dcst):
        """
        Sets the diffusion constant of diffusion rule with identifier string diff
        in compartment with identifier string comp to dcst (in m^2/s).

        Note: This method will set the diffusion constant in all tetrahedral elements
        in the compartment.

        Note: The default value still comes from the steps.model description,
        so calling reset() will return the diffusion constants to that value.

        Syntax::

            setCompDiffD(comp, diff, dcst)

        Arguments:
        string comp
        string diff
        float dcst

        Return:
            None

        """
        self.ptrx().setCompDiffD(compartment_id(to_std_string(c)), diffusion_id(to_std_string(d)), dcst)

    def getPatchComplexCount(self, str c, str complex, f):
        """
        Returns the number of molecules of a complex with identifier string complex
        matching filter filt in patch with identifier string patch.

        In a mesh-based simulation this is the combined count from 
        all triangles in the patch.

        Syntax::
            
            getPatchComplexCount(patch, complex, filt)
            
        Arguments:
        string patch
        string complex
        list[uint[:]] filt

        Return:
        float

        """
        return self.ptrx().getPatchComplexCount(patch_id(to_std_string(c)), complex_name(to_std_string(complex)), _get_filters(f))

    def setPatchComplexCount(self, str c, str complex, i, double n, distributionMethod=_py_DistributionMethod.UNIFORM):
        """
        Set the number of molecules of a complex with identifier string complex
        in state init in patch with identifier string patch.

        In a mesh-based simulation this is the combined count from 
        all triangles in the patch.

        Syntax::
            
            setPatchComplexCount(patch, complex, init, nspec)
            
        Arguments:
        string patch
        string complex
        uint[:] init
        int nspec

        Return:
        None

        """
        self.ptrx().setPatchComplexCount(patch_id(to_std_string(c)), complex_name(to_std_string(complex)), _get_filters(i), n, distributionMethod)

    def getPatchComplexSUSCount(self, str c, str complex, f, uint m):
        """
        Returns the number of subunits in state m of a complex with identifier string complex
        matching filter filt in patch with identifier string patch.

        In a mesh-based simulation this is the combined count from 
        all triangles in the patch.

        Syntax::
            
            getPatchComplexSUSCount(patch, complex, filt, m)
            
        Arguments:
        string patch
        string complex
        list[uint[:]] filt
        uint m

        Return:
        float

        """
        return self.ptrx().getPatchComplexSUSCount(patch_id(to_std_string(c)), complex_name(to_std_string(complex)), _get_filters(f), complex_substate_id(m))


    def getPatchSpecCount(self, str patch, str spec):
        """
        Returns the number of molecules of species with identifier string spec in patch 
        with identifier string pat.Note: in a mesh-based simulation this 
        is the combined count from all triangular elements in the patch. 

        Syntax::
            
            getPatchSpecCount(patch, spec)
            
        Arguments:
        string patch
        string spec

        Return:
        float

        """
        return self.ptrx().getPatchSpecCount(patch_id(to_std_string(patch)), species_name(to_std_string(spec)))

    def setPatchSpecCount(self, str patch, str spec, double n, distributionMethod=_py_DistributionMethod.UNIFORM):
        """
        Sets the number of molecules of species with identifier string spec in patch 
        with identifier string pat to n.


        Note: In case of a mesh-based simulation the molecules, molecules are divided among triangles.

        distributionMethod=UNIFORM the distribution is deterministic (apart from roundings) and the number of
        molecules per element is n*V_tet/V_tot.

        distributionMethod=DIST_MULTINOMIAL the distribution is multinomial and the probability
        of putting an element in a tet is V_tet/V_tot

        Syntax::

            setPatcSpechCount(patch, spec, n, distributionMethod)
            
        Arguments:
        string patch
        string spec
        int n
        DistributionMethod distributionMethod

        Return:
        float

        """
        self.ptrx().setPatchSpecCount(patch_id(to_std_string(patch)), species_name(to_std_string(spec)), n, distributionMethod)

    def getPatchSpecClamped(self, str patch, str spec):
        """
        Returns whether species with identifier string spec is clamped
        in the patch with identifier string patch.

        Syntax::

            getPatchSpecClamped(patch, spec)

        Arguments:
        string patch
        string spec

        Return:
        bool

        """
        return self.ptrx().getPatchSpecClamped(patch_id(to_std_string(patch)), species_name(to_std_string(spec)))

    def setPatchSpecClamped(self, str patch, str spec, bool clamped):
        """
        Sets whether species with identifier string spec is clamped
        in the patch with identifier string patch.

        Syntax::

            setPatchSpecClamped(patch, spec, clamped)

        Arguments:
        string patch
        string spec
        bool clamped

        Return:
        None

        """
        self.ptrx().setPatchSpecClamped(patch_id(to_std_string(patch)), species_name(to_std_string(spec)), clamped)

    if USE_PETSC:

        def setMembPotential(self, str memb, double v):
            """
            Sets the potential (in volts) of membrane with string identifier memb.
            NOTE: This method will set the potential of all nodes in the volume conductor
            to the same value.

            Syntax::
                        
                setMembPotential(memb, v)
                    
            Arguments:
            string memb
            float v

            Return:
            None

            """
            self.ptrx().setMembPotential(membrane_id(to_std_string(memb)), v)

        def setMembVolRes(self, str memb, double ro):
            """
            Set the bulk electrical resistivity of the section of the mesh
            representing the volume conductor for the membrane with string identifier memb.

            Syntax::

                setMembVolRes(memb, ro)

            Arguments:
            string memb
            float ro

            Return:
            None

            """
            self.ptrx().setMembVolRes(membrane_id(to_std_string(memb)), ro)

        def setMembCapac(self, str memb, double capac):
            """
            Sets the surface capacitance (in F.m^-2) of the membrane with string identifier memb.

            Syntax::

                setMembCapac(memb, capac)

            Arguments:
            string memb
            float capac

            Return:
            None

            """
            self.ptrx().setMembCapac(membrane_id(to_std_string(memb)), capac)

        def setMembRes(self, str memb, double ro, double vrev):
            """
            Sets the surface electrical resistivity ro (in ohm.m^2) of the membrane with string identifier memb. Reversal potential vrev is required in Volts.
            
            Syntax::
                        
                setMembRes(memb, ro, vrev)
                    
            Arguments:
            string memb
            float ro
            float vrev

            Return:
            None

            """
            self.ptrx().setMembRes(membrane_id(to_std_string(memb)), ro, vrev)

        def getMembRes(self, str membrane):
            """
            Gets the resistivity and the reversal potential of the membrane 
            with string identifier membrane.

            Syntax::
                        
                getMembRes(membrane)
                    
            Arguments:
            string membrane

            Return:
            pair: (double, double)

            """
            cdef MembraneResistivity val = self.ptrx().getMembRes(membrane_id(to_std_string(membrane)))
            return _py_MembraneResistivity(val.resistivity, val.reversal_potential)

        def setMembIClamp(self, str memb, float current):
            """
            Set a current clamp on a membrane

            Syntax::

                setMembIClamp(memb, current)

            Arguments:
            str memb
            float current

            Return:
            None

            """
            self.ptrx().setMembIClamp(membrane_id(to_std_string(memb)), current)

        def getTriCapac(self, GO idx, bool local=False):
            """
            Returns the specific membrane capacitance (in farad / m^2) of triangle with index idx.

            Syntax::

                getTriCapac(idx)

            Arguments:
            index_t idx
            bool local

            Return:
            float

            """
            return self.ptrx().getTriCapac(idx, local)

        def setTriCapac(self, GO idx, double cm, bool local=False):
            """
            Sets the specific membrane capacitance (in farad / m^2) of triangle with index idx.

            Syntax::

                setTriCapac(idx, cm)

            Arguments:
            index_t idx
            float cm
            bool local

            Return:
            None

            """
            self.ptrx().setTriCapac(idx, cm, local)

        def getTriRes(self, GO idx, bool local=False):
            """
            Returns the membrane resistivity (in ohm m^2) and reversal potential (in V) for the leak current in triangle with index idx.

            Syntax::

                getTriRes(idx)

            Arguments:
            index_t idx
            bool local

            Return:
            Tuple[float, float]

            """
            cdef MembraneResistivity val = self.ptrx().getTriRes(idx, local)
            return _py_MembraneResistivity(val.resistivity, val.reversal_potential)

        def setTriRes(self, GO idx, double res, double erev, bool local=False):
            """
            Sets the membrane resistivity (in ohm m^2) and reversal potential (in V) for the leak current in triangle with index idx.

            Syntax::

                setTriRes(idx, res, erev)

            Arguments:
            index_t idx
            float res
            float erev
            bool local

            Return:
            None

            """
            self.ptrx().setTriRes(idx, res, erev, local)

        def getVertIClamp(self, GO idx, bool local=False):
            """
            Returns the current clamp on the vertex with index idx, in ampere.
            NOTE: Convention is maintained that a positive current clamp is depolarizing, a negative current clamp is hyperpolarizing.

            Syntax::

                getVertIClamp(idx)

            Arguments:
            GO idx
            bool local

            Return:
            float

            """
            return self.ptrx().getVertIClamp(idx, local)

        def setVertIClamp(self, GO idx, double current, bool local=False):
            """
            Set the current clamp on the vertex with index idx, in ampere.
            NOTE: Convention is maintained that a positive current clamp is depolarizing, a negative current clamp is hyperpolarizing.

            Syntax::

                setVertIClamp(idx, current)

            Arguments:
            GO idx
            float current
            bool local

            Return:
            None

            """
            self.ptrx().setVertIClamp(idx, current, local)

        def getVertVClamped(self, GO idx, bool local=False):
            """
            Gets voltage clamp in vertex.

            Syntax::

                getVertVClamped(idx)

            Arguments:
            GO idx
            bool local

            Return:
            bool
            """
            return self.ptrx().getVertVClamped(idx, local)

        def setVertVClamped(self, GO idx, bool clamped, bool local=False):
            """
            Sets voltage clamp in vertex.

            Syntax::

                setVertVClamped(idx, clamped)

            Arguments:
            GO idx
            bool clamped
            bool local

            Return:
            None
            """
            self.ptrx().setVertVClamped(idx, clamped, local)

        def getTriOhmicErev(self, GO idx, str ohmic_current, bool local=False):
            """
            Gets the ohmic current reversal potential of triangle in volts.

            Arguments:
                idx: Index of the triangle
                ohmic_current: name of the ohmic current
                local: whether the triangle index is local to the process or global to the mesh

            Return:
                double
            """
            return self.ptrx().getTriOhmicErev(idx, ohmic_current_id(to_std_string(ohmic_current)), local)

        def getTriVClamped(self, GO idx, bool local=False):
            """
            Gets voltage clamp in triangle.

            Syntax::

                getTriVClamped(idx)

            Arguments:
            GO idx
            bool local

            Return:
            bool
            """
            return self.ptrx().getTriVClamped(idx, local)

        def setTriVClamped(self, GO idx, bool clamped, bool local=False):
            """
            Sets voltage clamp in triangle.

            Syntax::

                setTriVClamped(idx, clamped)

            Arguments:
            GO idx
            bool clamped
            bool local

            Return:
            None
            """
            self.ptrx().setTriVClamped(idx, clamped, local)

        def getTriIClamp(self, GO idx, bool local=False):
            """
            Returns the current clamp on the triangle with index idx, in ampere.
            NOTE: Convention is maintained that a positive current clamp is depolarizing, a negative current clamp is hyperpolarizing.

            Syntax::

                getTriIClamp(idx)

            Arguments:
            GO idx
            bool local

            Return:
            float

            """
            return self.ptrx().getTriIClamp(idx, local)

        def setTriIClamp(self, GO idx, double current, bool local=False):
            """
            Set the current clamp on the triangle with index idx, in ampere.
            NOTE: Convention is maintained that a positive current clamp is depolarizing, a negative current clamp is hyperpolarizing.

            Syntax::

                setTriIClamp(idx, current)

            Arguments:
            GO idx
            float current
            bool local

            Return:
            None

            """
            self.ptrx().setTriIClamp(idx, current, local)

        def setTriOhmicErev(self, GO idx, str ohmic_current, double reversal_potential, bool local=False):
            """
            Sets the ohmic current reversal potential of triangle in volts.

            Arguments:
                idx: Index of the triangle
                ohmic_current: name of the ohmic current
                reversal_potential: value in volts to assign
                local: whether the triangle index is local to the process or global to the mesh
            """
            return self.ptrx().setTriOhmicErev(idx, ohmic_current_id(to_std_string(ohmic_current)), reversal_potential, local)

        def getTriComplexOhmicErev(self, GO idx, str ohmic_current, bool local=False):
            """
            Gets the complex ohmic current reversal potential of triangle in volts.

            Arguments:
                idx: Index of the triangle
                ohmic_current: name of the complex ohmic current
                local: whether the triangle index is local to the process or global to the mesh

            Return:
                double
            """
            return self.ptrx().getTriComplexOhmicErev(idx, complex_ohmic_current_id(to_std_string(ohmic_current)), local)

        def setTriComplexOhmicErev(self, GO idx, str ohmic_current, double reversal_potential, bool local=False):
            """
            Sets the complex ohmic current reversal potential of triangle in volts.

            Arguments:
                idx: Index of the triangle
                ohmic_current: name of the complex ohmic current
                reversal_potential: value in volts to assign
                local: whether the triangle index is local to the process or global to the mesh
            """
            return self.ptrx().setTriComplexOhmicErev(idx, complex_ohmic_current_id(to_std_string(ohmic_current)), reversal_potential, local)

        def getVertV(self, GO idx, bool local=False):
            """
            Returns the potential (in volts) of vertex element with index idx.

            Syntax::

                getVertV(idx)

            Arguments:
            GO idx
            bool local

            Return:
            float

            """
            return self.ptrx().getVertV(idx, local)

        def getTriV(self, GO idx, bool local=False):
            """
            Returns the potential (in volts) of triangle element with index idx.

            Syntax::

                getTetV(idx)

            Arguments:
            GO idx
            bool local

            Return:
            float

            """
            return self.ptrx().getTriV(idx, local)

        def getTetV(self, GO idx, bool local=False):
            """
            Returns the potential (in volts) of tetrahdron element with index idx.

            Syntax::

                getTetV(idx)

            Arguments:
            GO idx
            bool local

            Return:
            float

            """
            return self.ptrx().getTetV(idx, local)

        def setVertV(self, GO idx, double v, bool local=False):
            """
            Sets the potential (in volts) of vertex element with index idx.

            Syntax::

                setVertV(idx, v)

            Arguments:
            GO idx
            float v
            bool local
            """
            return self.ptrx().setVertV(idx, v, local)

        def setTriV(self, GO idx, double v, bool local=False):
            """
            Sets the potential (in volts) of triangle element with index idx.

            Syntax::

                setTetV(idx)

            Arguments:
            GO idx
            float v
            bool local
            """
            return self.ptrx().setTriV(idx, v, local)

        def setTetV(self, GO idx, double v, bool local=False):
            """
            Sets the potential (in volts) of tetrahdron element with index idx.

            Syntax::

                setTetV(idx)

            Arguments:
            GO idx
            float v
            bool local
            """
            return self.ptrx().setTetV(idx, v, local)

        def getTriOhmicI(self, GO idx, str oc, bool local=False):
            """
            Returns the ohmic current of triangle element with index idx, in amps.

            Syntax::

                getTriOhmicI(idx, oc)

            Arguments:
            GO idx
            string oc
            bool local

            Return:
            float

            """
            return self.ptrx().getTriOhmicI(idx, ohmic_current_id(to_std_string(oc)), local)

        def getTriComplexOhmicI(self, GO idx, str oc, bool local=False):
            """
            Returns the complex ohmic current of triangle element with index idx, in amps.

            Syntax::

                getTriComplexOhmicI(idx, oc)

            Arguments:
            GO idx
            string oc
            bool local

            Return:
            float

            """
            return self.ptrx().getTriComplexOhmicI(idx, complex_ohmic_current_id(to_std_string(oc)), local)

        def getTriGHKI(self, GO idx, str ghk, bool local=False):
            """
            Returns the GHK current of triangle element with index idx, in amps.

            Syntax::

                getTriGHKI(idx, ghk)

            Arguments:
            GO idx
            string ghk
            bool local

            Return:
            float

            """
            return self.ptrx().getTriGHKI(idx, ghk_current_id(to_std_string(ghk)), local)

        def getTriComplexGHKI(self, GO idx, str ghk, bool local=False):
            """
            Returns the complex GHK current of triangle element with index idx, in amps.

            Syntax::

                getTriComplexGHKI(idx, ghk)

            Arguments:
            GO idx
            string ghk
            bool local

            Return:
            float

            """
            return self.ptrx().getTriComplexGHKI(idx, complex_ghk_current_id(to_std_string(ghk)), local)

        def getTriSReacI(self, GO idx, str reac, bool local=False):
            """
            Syntax::

                getTriSReacI(tri, reac)

            Arguments:
            GO idx
            string reac
            bool local

            Return:
            float

            """
            return self.ptrx().getTriSReacI(idx, surface_reaction_id(to_std_string(reac)), local)

        def getTriVDepSReacI(self, GO idx, str reac, bool local=False):
            """
            Syntax::

                getTriVDepSReacI(tri, reac)

            Arguments:
            GO idx
            string reac
            bool local

            Return:
            float

            """
            return self.ptrx().getTriVDepSReacI(idx, vdep_surface_reaction_id(to_std_string(reac)), local)

        def getTriComplexSReacI(self, GO idx, str reac, bool local=False):
            """
            Syntax::

                getTriComplexSReacI(tri, reac)

            Arguments:
            GO idx
            string reac
            bool local

            Return:
            float

            """
            return self.ptrx().getTriComplexSReacI(idx, complex_surface_reaction_id(to_std_string(reac)), local)

        def getTriVDepComplexSReacI(self, GO idx, str reac, bool local=False):
            """
            Syntax::

                getTriVDepComplexSReacI(tri, reac)

            Arguments:
            GO idx
            string reac
            bool local

            Return:
            float

            """
            return self.ptrx().getTriVDepComplexSReacI(idx, vdep_complex_surface_reaction_id(to_std_string(reac)), local)

        def getTriI(self, GO idx, bool local=False):
            """
            Returns the current of triangle element with index idx, in amps.

            Syntax::

                getTriI(idx)

            Arguments:
            GO idx
            bool local

            Return:
            float

            """
            return self.ptrx().getTriI(idx, local)

        def getTetV(self, GO idx, bool local=False):
            """
            Returns the potential (in volts) of tetrahdron element with index idx.

            Syntax::

                getTetV(idx)

            Arguments:
            GO idx
            bool local

            Return:
            float

            """
            return self.ptrx().getTetV(idx, local)

        def getTetVClamped(self, GO idx, bool local=False):
            """
            Gets voltage clamp in tetrahedron.

            Syntax::

                getTetVClamped(idx)

            Arguments:
            GO idx
            bool local

            Return:
            bool
            """
            return self.ptrx().getTetVClamped(idx, local)

        def setTetVClamped(self, GO idx, bool clamped, bool local=False):
            """
            Sets voltage clamp in tetrahedron.

            Syntax::

                setTetVClamped(idx, clamped)

            Arguments:
            GO idx
            bool clamped
            bool local

            Return:
            None
            """
            self.ptrx().setTetVClamped(idx, clamped, local)

    def reset(self):
        """
        Reset the simulation to the state the solver was initialised to.

        Syntax::

            reset()

        Arguments:
        None

        Return:
        None
        """
        self.ptrx().reset()

    def run(self, double endtime):
        """
        Advance the simulation until endtime (given in seconds) is reached.
        The endtime must be larger or equal to the current simulation time.

        Syntax::

            run(endtime)

        Arguments:
        float endtime

        Return:
        None
        """
        self.ptrx().run(endtime)

    def getTime(self):
        """
        Returns the current simulation time in seconds.

        Syntax::

            getTime()

        Arguments:
        None

        Return:
        float
        """
        return self.ptrx().getTime()

    def getTetSpecCount(self, GO idx, str spec, bool local=False):
        """
        Returns the number of molecules of species with identifier string spec 
        in the tetrahedral element with index idx.

        Syntax::
            
            getTetSpecCount(idx, spec)
            
        Arguments:
        GO idx
        string spec
        bool local

        Return:
        int

        """
        return self.ptrx().getTetSpecCount(idx, species_name(to_std_string(spec)), local)

    def getTetSpecConc(self, GO idx, str spec, bool local=False):
        """
        Returns the concentration (in Molar units) of species with identifier 
        string spec in a tetrahedral element with index idx.

        Syntax::
            
            getTetSpecConc(idx, spec)
            
        Arguments:
        GO idx
        string spec
        bool local

        Return:
        float

        """
        return self.ptrx().getTetSpecConc(idx, species_name(to_std_string(spec)), local)

    def setTetSpecCount(self, GO idx, str spec, double n, bool local=False):
        """
        Sets the number of molecules of species with identifier string spec in 
        tetrahedral element with index idx to n.

        Syntax::
            
            setTetSpecCount(idx, spec, n)
            
        Arguments:
        GO idx
        string spec
        int n
        bool local

        Return:
        None

        """
        self.ptrx().setTetSpecCount(idx, species_name(to_std_string(spec)), n, local)

    def setTetSpecConc(self, GO idx, str spec, double c, bool local=False):
        """
        Sets the concentration (in Molar units) of species with identifier string spec 
        in a tetrahedral element with index idx to conc.This continuous value must be 
        converted internally to a discrete number of molecules. 

        Due to the small volumes of tetrahedral elements the difference between 'rounding 
        up' and 'rounding down' can be a large difference in concentration.

        Syntax::
            
            setTetSpecConc(idx, spec, c)
            
        Arguments:
        GO idx
        string spec
        float c
        bool local

        Return:
        None

        """
        self.ptrx().setTetSpecConc(idx, species_name(to_std_string(spec)), c, local)

    def getTetSpecClamped(self, GO idx, str spec, bool local=False):
        """
        Returns whether species with identifier string spec is clamped
        in the tetrahedral element with index idx.

        Syntax::

            getTetSpecClamped(idx, spec)

        Arguments:
        GO idx
        string spec
        bool local

        Return:
        bool

        """
        return self.ptrx().getTetSpecClamped(idx, species_name(to_std_string(spec)), local)

    def setTetSpecClamped(self, GO idx, str spec, bool clamped, bool local=False):
        """
        Sets whether species with identifier string spec is clamped
        in the tetrahedral element with index idx.

        Syntax::
 
            setTetSpecClamped(idx, spec, clamped)

        Arguments:
        GO idx
        string spec
        bool clamped
        bool local

        Return:
        None

        """
        self.ptrx().setTetSpecClamped(idx, species_name(to_std_string(spec)), clamped, local)

    def getTriSpecCount(self, GO idx, str spec, bool local=False):
        """
        Returns the number of molecules of species with identifier string spec 
        in the triangular element with index idx.

        Syntax::
            
            getTriSpecCount(idx, spec)
            
        Arguments:
        GO idx
        string spec
        bool local

        Return:
        float

        """
        return self.ptrx().getTriSpecCount(idx, species_name(to_std_string(spec)), local)

    def setTriSpecCount(self, GO idx, str spec, double n, bool local=False):
        """
        Sets the number of molecules of species with identifier string spec in 
        triangular element with index idx to n. 

        Syntax::
            
            setTriSpecCount(idx, spec, n)
            
        Arguments:
        GO idx
        string spec
        int n
        bool local

        Return:
        None

        """
        self.ptrx().setTriSpecCount(idx, species_name(to_std_string(spec)), n, local)

    def getTriSpecClamped(self, GO idx, str spec, bool local=False):
        """
        Returns whether species with identifier string spec is clamped
        in the triangular element with index idx.

        Syntax::

            getTriSpecClamped(idx, spec)

        Arguments:
        GO idx
        string spec
        bool local

        Return:
        bool

        """
        return self.ptrx().getTriSpecClamped(idx, species_name(to_std_string(spec)), local)

    def setTriSpecClamped(self, GO idx, str spec, bool clamped, bool local=False):
        """
        Sets whether species with identifier string spec is clamped
        in the triangular element with index idx.

        Syntax::
 
            setTriSpecClamped(idx, spec, clamped)
 
        Arguments:
        GO idx
        string spec
        bool clamped
        bool local

        Return:
        None

        """
        self.ptrx().setTriSpecClamped(idx, species_name(to_std_string(spec)), clamped, local)

    def getTriSReacK(self, GO idx, str reac, bool local=False):
        """
        Syntax::

            setTriSReacK(tri, reac)

        Arguments:
        GO idx
        string reac
        bool local

        Return:
        float

        """
        return self.ptrx().getTriSReacK(idx, surface_reaction_id(to_std_string(reac)), local)

    def setTriSReacK(self, GO idx, str reac, double kf, bool local=False):
        """
        Syntax::
            
            setTriSReacK(tri, reac, kf)
            
        Arguments:
        GO idx
        string reac
        float kf
        bool local

        Return:
        None

        """
        self.ptrx().setTriSReacK(idx, surface_reaction_id(to_std_string(reac)), kf, local)

    def getTriComplexSReacK(self, GO idx, str reac, bool local=False):
        """
        Syntax::

            setTriComplexSReacK(tri, reac)

        Arguments:
        GO idx
        string reac
        bool local

        Return:
        float

        """
        return self.ptrx().getTriComplexSReacK(idx, complex_surface_reaction_id(to_std_string(reac)), local)

    def setTriComplexSReacK(self, GO idx, str reac, double kf, bool local=False):
        """
        Syntax::
            
            setTriComplexSReacK(tri, reac, kf)
            
        Arguments:
        GO idx
        string reac
        float kf
        bool local

        Return:
        None

        """
        self.ptrx().setTriComplexSReacK(idx, complex_surface_reaction_id(to_std_string(reac)), kf, local)

    def getPatchSReacK(self, str patch, str reac):
        """
        Gets the macroscopic reaction constant of surface reaction with identifier 
        string sreac in patch with identifier string pat to kf. The unit of the reaction 
        constant depends on the order of the reaction. 

        Note: In a mesh-based simulation this method gets the currently set value for the patch,
        individual triangles in the patch might have different values.

        Syntax::

            getPatchSReacK(patch, reac)

        Arguments:
        string patch
        string reac

        Return:
        float

        """
        return self.ptrx().getPatchSReacK(patch_id(to_std_string(patch)), surface_reaction_id(to_std_string(reac)))

    def setPatchSReacK(self, str patch, str reac, double kf):
        """
        Sets the macroscopic reaction constant of surface reaction with identifier 
        string sreac in patch with identifier string pat to kf. The unit of the reaction 
        constant depends on the order of the reaction. 

        Note: In a mesh-based simulation this method sets the surface 
        reaction constant in all triangular elements of the patch to kf.

        Note: The default value still comes from the steps.model description, so calling 
        reset() will return the surface reaction constant to that value.

        Syntax::
            
            setPatchSReacK(patch, reac, kf)
            
        Arguments:
        string patch
        string reac
        float kf

        Return:
        None

        """
        self.ptrx().setPatchSReacK(patch_id(to_std_string(patch)), surface_reaction_id(to_std_string(reac)), kf)

    def getPatchSReacExtent(self, str patch, str reac):
        """
        Gets the extent of surface reaction with identifier string reac in patch with identifier string patch.

        Syntax::

            getPatchSReacExtent(patch, reac)

        Arguments:
        string patch
        string reac

        Return:
        uint

        """
        return self.ptrx().getPatchSReacExtent(patch_id(to_std_string(patch)), surface_reaction_id(to_std_string(reac)))

    def getPatchComplexSReacExtent(self, str patch, str reac):
        """
        Gets the extent of complex surface reaction with identifier string reac in patch with identifier string patch.

        Syntax::

            getPatchComplexSReacExtent(patch, reac)

        Arguments:
        string patch
        string reac

        Return:
        uint

        """
        return self.ptrx().getPatchComplexSReacExtent(patch_id(to_std_string(patch)), complex_surface_reaction_id(to_std_string(reac)))

    def getPatchVDepSReacExtent(self, str patch, str reac):
        """
        Gets the extent of voltage-dependent surface reaction with identifier string reac in patch with identifier string patch.

        Syntax::

            getPatchVDepSReacExtent(patch, reac)

        Arguments:
        string patch
        string reac

        Return:
        uint

        """
        return self.ptrx().getPatchVDepSReacExtent(patch_id(to_std_string(patch)), vdep_surface_reaction_id(to_std_string(reac)))

    def getPatchVDepComplexSReacExtent(self, str patch, str reac):
        """
        Gets the extent of voltage-dependent complex surface reaction with identifier string reac in patch with identifier string patch.

        Syntax::

            getPatchVDepComplexSReacExtent(patch, reac)

        Arguments:
        string patch
        string reac

        Return:
        uint

        """
        return self.ptrx().getPatchVDepComplexSReacExtent(patch_id(to_std_string(patch)), vdep_complex_surface_reaction_id(to_std_string(reac)))

    # # ---------------------------------------------------------------------------------
    # # NUMPY section - we accept numpy arrays and generically typed memory-views
    # # ---------------------------------------------------------------------------------
    def getBatchTetSpecCountsNP(self, GO[:] indices, str spec, double[:] counts, bool local=False):
        """
        Get the counts of a species s in a list of tetrahedrons.

        Syntax::
            getBatchTetSpecCountsNP(indices, spec, counts)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> counts
        bool local

        Return:
        None

        """
        self.ptrx().getBatchTetSpecCountsNP(&indices[0], indices.shape[0], species_name(to_std_string(spec)), &counts[0], counts.shape[0], local)

    def setBatchTetSpecCountsNP(self, GO[:] indices, str spec, double[:] counts, bool local=False):
        """
        Set the counts of a species s in a list of tetrahedrons.

        Syntax::
            setBatchTetSpecCountsNP(indices, spec, counts)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> counts
        bool local

        Return:
        None

        """
        self.ptrx().setBatchTetSpecCountsNP(&indices[0], indices.shape[0], species_name(to_std_string(spec)), &counts[0],
                counts.shape[0], local)

    def getBatchTetSpecConcsNP(self, GO[:] indices, str spec, double[:] concs, bool local=False):
        """
        Get the individual concentration of a species s in a list of tetrahedrons.

        Syntax::
            getBatchTetSpecCountsNP(indices, spec, concs)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> concs
        bool local

        Return:
        None

        """
        self.ptrx().getBatchTetSpecConcsNP(&indices[0], indices.shape[0], species_name(to_std_string(spec)), &concs[0], concs.shape[0], local)

    def setBatchTetSpecConcsNP(self, GO[:] indices, str spec, double[:] concs, bool local=False):
        """
        Set the concetration of a species s in a list of tetrahedrons.

        Syntax::
            setBatchTetSpecConcsNP(indices, spec, concs)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> concs
        bool local

        Return:
        None

        """
        self.ptrx().setBatchTetSpecConcsNP(&indices[0], indices.shape[0], species_name(to_std_string(spec)), &concs[0],
                concs.shape[0], local)

    def getBatchTriSpecCountsNP(self, GO[:] indices, str spec, double[:] counts, bool local=False):
        """
        Get the counts of a species s in a list of triangles.

        Syntax::
            getBatchTriSpecCountsNP(indices, spec, counts)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> counts
        bool local

        Return:
            None

        """
        self.ptrx().getBatchTriSpecCountsNP(&indices[0], indices.shape[0], species_name(to_std_string(spec)), &counts[0], counts.shape[0], local)

    def setBatchTriSpecCountsNP(self, GO[:] indices, str spec, double[:] counts, bool local=False):
        """
        Set the counts of a species s in a list of triangles.

        Syntax::
            getBatchTriSpecCountsNP(indices, spec, counts)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> counts
        bool local

        Return:
            None

        """
        self.ptrx().setBatchTriSpecCountsNP(&indices[0], indices.shape[0], species_name(to_std_string(spec)), &counts[0],
                counts.shape[0], local)

    if USE_PETSC:

        def getBatchVertVsNP(self, GO[:] indices, double[:] voltages, bool local=False):
            """
            Get the potential in a list of vertices.

            Syntax::
                getBatchVertVsNP(indices, voltages)

            Arguments:
            numpy.array<GO> indices
            numpy.array<double, length = len(indices)> voltages
            bool local

            Return:
                None

            """
            self.ptrx().getBatchVertVsNP(&indices[0], indices.shape[0], &voltages[0], voltages.shape[0], local)

        def getBatchTriVsNP(self, GO[:] indices, double[:] voltages, bool local=False):
            """
            Get the potential in a list of triangles.

            Syntax::
                getBatchTriVsNP(indices, voltages)

            Arguments:
            numpy.array<GO> indices
            numpy.array<double, length = len(indices)> voltages
            bool local

            Return:
                None

            """
            self.ptrx().getBatchTriVsNP(&indices[0], indices.shape[0], &voltages[0], voltages.shape[0], local)

        def getBatchTetVsNP(self, GO[:] indices, double[:] voltages, bool local=False):
            """
            Get the potential in a list of tetrahedra.

            Syntax::
                getBatchTetVsNP(indices, voltages)

            Arguments:
            numpy.array<GO> indices
            numpy.array<double, length = len(indices)> voltages
            bool local

            Return:
                None

            """
            self.ptrx().getBatchTetVsNP(&indices[0], indices.shape[0], &voltages[0], voltages.shape[0], local)

        def setBatchVertVsNP(self, GO[:] indices, double[:] voltages, bool local=False):
            """
            set the potential in a list of vertices.

            Syntax::
                setBatchVertVsNP(indices, voltages)

            Arguments:
            numpy.array<GO> indices
            numpy.array<double, length = len(indices)> voltages
            bool local

            Return:
                None

            """
            self.ptrx().setBatchVertVsNP(&indices[0], indices.shape[0], &voltages[0], voltages.shape[0], local)

        def setBatchTriVsNP(self, GO[:] indices, double[:] voltages, bool local=False):
            """
            set the potential in a list of triangles.

            Syntax::
                setBatchTetVsNP(indices, voltages)

            Arguments:
            numpy.array<GO> indices
            numpy.array<double, length = len(indices)> voltages
            bool local

            Return:
                None

            """
            self.ptrx().setBatchTriVsNP(&indices[0], indices.shape[0], &voltages[0], voltages.shape[0], local)

        def setBatchTetVsNP(self, GO[:] indices, double[:] voltages, bool local=False):
            """
            set the potential in a list of tetrahedra.

            Syntax::
                setBatchTetVsNP(indices, voltages)

            Arguments:
            numpy.array<GO> indices
            numpy.array<double, length = len(indices)> voltages
            bool local

            Return:
                None

            """
            self.ptrx().setBatchTetVsNP(&indices[0], indices.shape[0], &voltages[0], voltages.shape[0], local)

        def getBatchTriOhmicIsNP(self, GO[:] indices, str oc, double[:] currents, bool local=False):
            """
            Get the Ohmic currents in a list of triangles.

            Syntax::
                getBatchTriOhmicIsNP(indices, oc, currents)

            Arguments:
            numpy.array<GO> indices
            string oc
            numpy.array<double, length = len(indices)> currents
            bool local

            Return:
                None

            """
            self.ptrx().getBatchTriOhmicIsNP(&indices[0], indices.shape[0], ohmic_current_id(to_std_string(oc)), &currents[0], currents.shape[0], local)

        def getBatchTriComplexOhmicIsNP(self, GO[:] indices, str oc, double[:] currents, bool local=False):
            """
            Get the complex Ohmic currents in a list of triangles.

            Syntax::
                getBatchTriComplexOhmicIsNP(indices, oc, currents)

            Arguments:
            numpy.array<GO> indices
            string oc
            numpy.array<double, length = len(indices)> currents
            bool local

            Return:
                None

            """
            self.ptrx().getBatchTriComplexOhmicIsNP(&indices[0], indices.shape[0], complex_ohmic_current_id(to_std_string(oc)), &currents[0], currents.shape[0], local)

        def getBatchTriGHKIsNP(self, GO[:] indices, str ghk, double[:] currents, bool local=False):
            """
            Get the GHK currents in a list of triangles.

            Syntax::
                getBatchTriGHKIsNP(indices, ghk, currents)

            Arguments:
            numpy.array<GO> indices
            string ghk
            numpy.array<double, length = len(indices)> currents
            bool local

            Return:
                None

            """
            self.ptrx().getBatchTriGHKIsNP(&indices[0], indices.shape[0], ghk_current_id(to_std_string(ghk)), &currents[0], currents.shape[0], local)

        def getBatchTriComplexGHKIsNP(self, GO[:] indices, str ghk, double[:] currents, bool local=False):
            """
            Get the complex GHK currents in a list of triangles.

            Syntax::
                getBatchTriComplexGHKIsNP(indices, ghk, currents)

            Arguments:
            numpy.array<GO> indices
            string ghk
            numpy.array<double, length = len(indices)> currents
            bool local

            Return:
                None

            """
            self.ptrx().getBatchTriComplexGHKIsNP(&indices[0], indices.shape[0], complex_ghk_current_id(to_std_string(ghk)), &currents[0], currents.shape[0], local)

        def getBatchTriIsNP(self, GO[:] indices, double[:] currents, bool local=False):
            """
            Get the currents in a list of triangles.

            Syntax::
                getBatchTriIsNP(indices, currents)

            Arguments:
            numpy.array<GO> indices
            numpy.array<double, length = len(indices)> currents
            bool local

            Return:
                None

            """
            self.ptrx().getBatchTriIsNP(&indices[0], indices.shape[0], &currents[0], currents.shape[0], local)

        def getBatchTriOhmicErevsNP(self, GO[:] triangles, str ohmic_current, double[:] rv, bool local=False):
            """
            """
            self.ptrx().getBatchTriOhmicErevsNP(&triangles[0], triangles.shape[0], ohmic_current_id(to_std_string(ohmic_current)), &rv[0], rv.shape[0], local)

        def getBatchTriComplexOhmicErevsNP(self, GO[:] triangles, str ohmic_current, double[:] rv, bool local=False):
            """
            """
            self.ptrx().getBatchTriComplexOhmicErevsNP(&triangles[0], triangles.shape[0], complex_ohmic_current_id(to_std_string(ohmic_current)), &rv[0], rv.shape[0], local)

    def getTetDiffD(self, GO idx, str diff, GO direction_tet=tetrahedron_global_id_t.unknown_value(), bool local=False):
        """
        Gets the diffusion constant of diffusion rule with identifier string diff in
        tetrahedral element with index idx to dcst (in m^2/s). Specify direction_tet to get the constant only towards a given tetrahedron direction.
        Syntax::

            getTetDiffD(idx, diff, dcst, direction_tet)

        Arguments:
        GO idx
        string diff
        GO direction_tet
        bool local

        Return:
        float

        """
        return self.ptrx().getTetDiffD(idx, diffusion_id(to_std_string(diff)), direction_tet, local)

    def setTetDiffD(self, GO idx, str diff, double dcst, GO direction_tet=tetrahedron_global_id_t.unknown_value(), bool local=False):
        """
        Sets the diffusion constant of diffusion rule with identifier string diff in
        tetrahedral element with index idx to dcst (in m^2/s). Specify direction_tet to set the constant only towards a given tetrahedron direction.
        Syntax::

            setTetDiffD(idx, diff, dcst, direction_tet)

        Arguments:
        GO idx
        string diff
        float dcst
        GO direction_tet
        bool local

        Return:
        None

        """
        self.ptrx().setTetDiffD(idx, diffusion_id(to_std_string(diff)), dcst, direction_tet, local)

    def getTetReacK(self, GO idx, str reac, bool local=False):
        """
        Syntax::

            getTetReacK(idx, reac)

        Arguments:
        GO idx
        string reac
        bool local

        Return:
        float

        """
        return self.ptrx().getTetReacK(idx, reaction_id(to_std_string(reac)), local)

    def setTetReacK(self, GO idx, str reac, double kf, bool local=False):
        """
        Syntax::

            setTetReacK(idx, reac, kf)

        Arguments:
        GO idx
        string reac
        float kf
        bool local

        Return:
        None

        """
        self.ptrx().setTetReacK(idx, reaction_id(to_std_string(reac)), kf, local)

    def getTetComplexReacK(self, GO idx, str reac, bool local=False):
        """
        Syntax::

            getTetComplexReacK(idx, reac)

        Arguments:
        GO idx
        string reac
        bool local

        Return:
        float

        """
        return self.ptrx().getTetComplexReacK(idx, complex_reaction_id(to_std_string(reac)), local)

    def setTetComplexReacK(self, GO idx, str reac, double kf, bool local=False):
        """
        Syntax::

            setTetComplexReacK(idx, reac, kf)

        Arguments:
        GO idx
        string reac
        float kf
        bool local

        Return:
        None

        """
        self.ptrx().setTetComplexReacK(idx, complex_reaction_id(to_std_string(reac)), kf, local)

    def setDiffBoundarySpecDiffusionActive(self, str diffb, str spec, bool act):
        """
        Activates or inactivates diffusion across a diffusion boundary for a species.
                     
        Syntax::
                     
            setDiffBoundaryDiffusionActive(diffb, spec, act)
                     
        Arguments:
        string diffb
        string spec
        bool act
                     
        Return:
        None

        """
        self.ptrx().setDiffBoundarySpecDiffusionActive(diffusion_boundary_name(to_std_string(diffb)), species_name(to_std_string(spec)), act)

    def getDiffBoundarySpecDiffusionActive(self, str diffb, str spec):
        """
        Returns whether diffusion is active across a diffusion boundary for a species.
                     
        Syntax::
                     
            getDiffBoundaryDiffusionActive(diffb, spec)
                     
        Arguments:
        string diffb
        string spec
                     
        Return:
        bool

        """
        return self.ptrx().getDiffBoundarySpecDiffusionActive(diffusion_boundary_name(to_std_string(diffb)), species_name(to_std_string(spec)))

    def setDiffBoundarySpecDcst(self, str diffb, str spec, double dcst):
        """
        Set the diffusion constant for diffusion across a diffusion boundary for a species.

        Syntax::

            setDiffBoundaryDcst(diffb, spec, dcst)

        Arguments:
        string diffb
        string spec
        float dcst

        Return:
        None

        """
        self.ptrx().setDiffBoundarySpecDcst(diffusion_boundary_name(to_std_string(diffb)), species_name(to_std_string(spec)), dcst)

    def setDiffApplyThreshold(self, int threshold):
        """
        Set the threshold for using binomial distribution for molecule diffusion instead of
        single molecule diffusion.

        If the number of molecules in a tetrahedron await for diffusion is higher than this
        threshold, the solver will use binomial function to distribute these molecules to
        each neighboring tetrahedron. Otherwise the molecules will diffuse one by one.

        The default threshold is 10.

        Syntax::

            setDiffApplyThreshold(threshold)

        Arguments:
        int threshold

        Return:
        None
        """
        if threshold < 0:
            raise ValueError(f'The threshold cannot be negative.')
        self.ptrx().setDiffApplyThreshold(threshold)

    def setTemp(self, double t):
        """
        Set the simulation temperature. Currently, this will only
        influence the GHK flux rate, so will only influence simulations
        including membrane potential calculation.

        Syntax::

            setTemp(temp)

        Arguments:
        float temp

        Return:
        None

        """
        self.ptrx().setTemp(t)

    def getTemp(self, ):
        """
        Return the simulation temperature.

        Syntax::

            getTemp()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getTemp()

    def getDiffusionTolerance(self):
        """
        Get the diffusion tolerance if the simulation is using the TAU_LEAPING_DT diffusion method

        Syntax::

            getDiffusionTolerance()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getDiffusionTolerance()

    def setDiffusionTolerance(self, double tol):
        """
        Set the diffusion tolerance if the simulation is using the TAU_LEAPING_DT diffusion method

        Syntax::

            setDiffusionTolerance(tol)

        Arguments:
        float tol

        Return:
        None

        """
        self.ptrx().setDiffusionTolerance(tol)

    def getDiffusionNormalApproximationThreshold(self):
        """
        Get the threshold for normal approximation to skellam distribution if the simulation is using the TAU_LEAPING_DT diffusion method

        Syntax::

            getDiffusionNormalApproximationThreshold()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getDiffusionNormalApproximationThreshold()

    def setDiffusionNormalApproximationThreshold(self, double thresh):
        """
        Set the threshold for normal approximation to skellam distribution if the simulation is using the TAU_LEAPING_DT diffusion method

        Syntax::

            setDiffusionNormalApproximationThreshold(thresh)

        Arguments:
        float thresh

        Return:
        None

        """
        self.ptrx().setDiffusionNormalApproximationThreshold(thresh)

    def getDiffusionCrankNicolsonThreshold(self):
        """
        Get the threshold for using Crank-Nicolson scheme if the simulation is using the TAU_LEAPING_DT diffusion method

        Syntax::

            getDiffusionCrankNicolsonThreshold()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getDiffusionCrankNicolsonThreshold()

    def setDiffusionCrankNicolsonThreshold(self, double thresh):
        """
        Set the threshold for using Crank-Nicolson scheme if the simulation is using the TAU_LEAPING_DT diffusion method

        Syntax::

            setDiffusionCrankNicolsonThreshold(thresh)

        Arguments:
        float thresh

        Return:
        None

        """
        self.ptrx().setDiffusionCrankNicolsonThreshold(thresh)

    def getDiffusionLeapThreshold(self):
        """
        Get the minimum number of species for leaping with TAU_LEAPING_DT diffusion

        Syntax::

            getDiffusionLeapThreshold()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().getDiffusionLeapThreshold()

    def setDiffusionLeapThreshold(self, int leap_thresh):
        """
        Set the minimum number of species for leaping with TAU_LEAPING_DT diffusion

        Syntax::

            setDiffusionLeapThreshold(leap_thresh)

        Arguments:
        int leap_thresh

        Return:
        None

        """
        self.ptrx().setDiffusionLeapThreshold(leap_thresh)

    def getDiffusionMaxDtSkips(self):
        """
        Get the maximum number of dt skips if the simulation is using the TAU_LEAPING_DT diffusion method

        Syntax::

            getDiffusionMaxDtSkips()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().getDiffusionMaxDtSkips()

    def setDiffusionMaxDtSkips(self, int skips):
        """
        Set the maximum number of dt skips if the simulation is using the TAU_LEAPING_DT diffusion method

        Syntax::

            setDiffusionMaxDtSkips(skips)

        Arguments:
        int skips

        Return:
        None

        """
        self.ptrx().setDiffusionMaxDtSkips(skips)

    def getDiffusionMinDtFactor(self):
        """
        Get the factor for computing the minimum diffusion dt if the simulation is using the TAU_LEAPING_DT diffusion method

        Syntax::

            getDiffusionMinDtFactor()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getDiffusionMinDtFactor()

    def setDiffusionMinDtFactor(self, float factor):
        """
        Set the factor for computing the minimum diffusion dt if the simulation is using the TAU_LEAPING_DT diffusion method

        Syntax::

            setDiffusionMinDtFactor(factor)

        Arguments:
        float factor

        Return:
        None

        """
        self.ptrx().setDiffusionMinDtFactor(factor)

    def getReactionSSAThreshold(self):
        """
        Get the minimum leap size for using R-leaping if the simulation is using the RLEAPING reaction operator

        Syntax::

            getReactionSSAThreshold()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().getReactionSSAThreshold()

    def setReactionSSAThreshold(self, int thresh):
        """
        Set the minimum leap size for using R-leaping if the simulation is using the RLEAPING reaction operator

        Syntax::

            setReactionSSAThreshold(thresh)

        Arguments:
        int thresh

        Return:
        None

        """
        self.ptrx().setReactionSSAThreshold(thresh)

    def getReactionSSASteps(self):
        """
        Get the number of standard SSA steps to run in a row when the reaction leaps are below threshold, only available if using the RLEAPING reaction operator.

        Syntax::

            getReactionSSASteps()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().getReactionSSASteps()

    def setReactionSSASteps(self, int steps):
        """
        Set the number of standard SSA steps to run in a row when the reaction leaps are below threshold, only available if using the RLEAPING reaction operator.

        Syntax::

            setReactionSSASteps(steps)

        Arguments:
        int steps

        Return:
        None

        """
        self.ptrx().setReactionSSASteps(steps)

    def getReactionLComputePeriod(self):
        """
        Get the period at which L is computed, only available if using the RLEAPING reaction operator.

        Syntax::

            getReactionLComputePeriod()

        Arguments:
        None

        Return:
        int

        """
        return self.ptrx().getReactionLComputePeriod()

    def setReactionLComputePeriod(self, int period):
        """
        Set the period at which L is computed, only available if using the RLEAPING reaction operator.

        Syntax::

            setReactionLComputePeriod(period)

        Arguments:
        int period

        Return:
        None

        """
        self.ptrx().setReactionLComputePeriod(period)

    def getReactionTolerance(self):
        """
        Get the reaction tolerance if the simulation is using the RLEAPING reaction operator

        Syntax::

            getReactionTolerance()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getReactionTolerance()

    def setReactionTolerance(self, double tol):
        """
        Set the reaction tolerance if the simulation is using the RLEAPING reaction operator

        Syntax::

            setReactionTolerance(tol)

        Arguments:
        float tol

        Return:
        None

        """
        self.ptrx().setReactionTolerance(tol)

    def getReactionTheta(self):
        """
        Get the theta parameter if the simulation is using the RLEAPING reaction operator

        Syntax::

            getReactionTheta()

        Arguments:
        None

        Return:
        float

        """
        return self.ptrx().getReactionTheta()

    def setReactionTheta(self, double theta):
        """
        Set the theta parameter if the simulation is using the RLEAPING reaction operator

        Syntax::

            setReactionTheta(theta)

        Arguments:
        float theta

        Return:
        None

        """
        self.ptrx().setReactionTheta(theta)

    if USE_PETSC:

        def setEfieldDT(self, double efdt):
            """
            Set the stepsize for membrane potential solver (default 1us).
            This is the time for each voltage calculation step. The SSA will
            run until passing this stepsize, so in fact each membrane potential
            time step will vary slightly around the dt so as to be aligned with the SSA.

            Syntax::

                setEFieldDT(efdt)

            Arguments:
            float efdt

            Return:
            None

            """
            self.ptrx().setEfieldDT(efdt)

        def getEfieldDT(self):
            """
            Get the stepsize for the membrane potential solver.

            Syntax::

                getEFieldDT()

            Arguments:
            None

            Return:
            float

            """
            return self.ptrx().getEfieldDT()

        def setPetscOptions(self, str options):
            """
            set PETSc options

            - for the ksp: https://petsc.org/release/docs/manualpages/KSP/KSPSetFromOptions.html
            - for the pc: https://petsc.org/release/docs/manualpages/PC/PCSetFromOptions.html

            is_unpreconditioned_norm=True is on the unpreconditioned norm, False is on the preconditioned norm

            Syntax::

                setPetscOptions("-option1 val1 -option2 val2")

            Arguments:
            string options: the string with the options in the form: "-option1 val1 -option2 val2"

            Return:
            None

            """
            self.ptrx().setPetscOptions(to_std_string(options))

    def dumpDepGraphToFile(self, str path):
        """
        Dump the kproc dependency graph in a file specified by path

        Syntax::

            dumpDepGraphToFile(path)

        Arguments:
        str path

        Return:
        None

        """
        self.ptrx().dumpDepGraphToFile(to_std_string(path))

