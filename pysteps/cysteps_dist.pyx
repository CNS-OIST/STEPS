###___license_placeholder___###

"""
This file is the user-interface file for all objects related to STEPS 4. All
objects are directly derived from the corresponding Cython objects.
"""

from cython.operator cimport dereference as deref
cimport mpi4py.MPI as MPI

include "cysteps_mpi.pyx"
from steps_dist_tetmesh cimport *
from steps_dist_solver cimport *

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

    def getTetTetNeighb(self, GO idx, bool local=False):
        """
        Returns the indices of the four neighbouring tetrahedrons of tetrahedron with index idx.
        An index of UNKNOWN_TET indicates no neighbour (tetrahedron is on the mesh border).

        Syntax::

            getTetTetNeighb(idx)

        Arguments:
        GO idx
        bool local

        Return:
        list<GO, length = 4>

        """
        if local:
            return [tl.get() if tl.valid() else UNKNOWN_TET for tl in self.ptrx().getTetTetNeighb(<tetrahedron_local_id_t>(<LO>(idx)))]
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

    def getTriTetNeighb(self, GO idx, bool local=False):
        """
        Returns the indices of the neighbouring tetrahedrons of triangle with
        index idx. If the triangle is on the mesh boundary, only one tetrahedron is returned

        Syntax::

            getTriTetNeighb(idx)

        Arguments:
        GO idx
        bool local

        Return:
        list<GO>

        """
        if local:
            return [tl.get() for tl in self.ptrx().getTriTetNeighb(<triangle_local_id_t>(<LO>(idx)))]
        else:
            return [tg.get() for tg in self.ptrx().getTriTetNeighb(<triangle_global_id_t>(idx))]

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
        return self.ptrx().getGlobalIndex(tetrahedron_local_id_t(idx))

    def getTriGlobalIndex(self, LO idx):
        """
        Return the global index of triangle with local index idx

        Syntax::

            getTetGlobalIndex(idx)

        Arguments:
        int idx

        Return:
        int
        """
        return self.ptrx().getGlobalIndex(triangle_local_id_t(idx))

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
        return self.ptrx().getGlobalIndex(vertex_local_id_t(idx))

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

    def intersect(self, double[:, :] points, int sampling=-1):
        """
        Computes the intersection of line segment(s) given the vertices with the current mesh

        Args:
            points: A 2-D NumPy array (/memview) of points in the 3D space, 
                    where each element contains the 3 point coordinates
            int sampling: not specified or sampling < 1 --> use deterministic method
                          sampling > 0 --> use montecarlo method with sampling points

        Returns:
            A list of lists of tuples representing the intersected tetrahedrons, one element per line segment. 
            Each tuple is made of 2 elements, a tetrahedron global identifier, and its respective intersection ratio. 
        """
        if (points.strides[0] != 24 or points.strides[1] != 8):
            raise Exception("Wrong memory layout for points, np array should be [pts,3] and row major")
        data = self.ptrx().intersect(&points[0][0], points.shape[0], sampling)
        return [[(t.first.get(), t.second) for t in row] for row in data]

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
        cdef std.vector[triangle_global_id_t] tri_indices
        cdef std.vector[triangle_local_id_t] local_tri_indices
        if physical_tag is None:
            if tris is None:
                self._ptr = new DistPatch(
                    patch_name(to_std_string(id)),
                    mesh.ptrx()[0],
                    icomp.ptrx() if icomp is not None else NULL,
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
                    icomp.ptrx() if icomp is not None else NULL,
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
                    icomp.ptrx() if icomp is not None else NULL,
                    ocomp.ptrx() if ocomp is not None else NULL
                )
        else:
            self._ptr = new DistPatch(
                patch_name(to_std_string(id)),
                mesh.ptrx()[0],
                patch_physical_tag(physical_tag),
                icomp.ptrx() if icomp is not None else NULL,
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
from steps_dist_solver cimport TetOpSplitBase, TetOpSplit

IF USE_PETSC:
    from steps_dist_solver cimport KSPNormType

    cdef class _py_KSPNormType:
        KSP_NORM_DEFAULT = steps_dist_solver.KSP_NORM_DEFAULT
        KSP_NORM_NONE = steps_dist_solver.KSP_NORM_NONE
        KSP_NORM_PRECONDITIONED = steps_dist_solver.KSP_NORM_PRECONDITIONED
        KSP_NORM_UNPRECONDITIONED = steps_dist_solver.KSP_NORM_UNPRECONDITIONED
        KSP_NORM_NATURAL = steps_dist_solver.KSP_NORM_NATURAL

cdef class _py_SSAMethod:
    SSA = 0
    RSSA = 1

cdef class _py_SearchMethod:
    DIRECT = 0
    GIBSON_BRUCK = 1

cdef class _py_DistributionMethod:
    UNIFORM = 0
    MULTINOMIAL = 1

# ----------------------------------------------------------------------------------------------------------------------
cdef class _py_DistTetOpSplitP(_py__base):
    """Bindings for MPI DistTetOpSplitP"""
# ----------------------------------------------------------------------------------------------------------------------

    cdef TetOpSplitBase *ptrx(self):
        return <TetOpSplitBase*> self._ptr

    def __init__(self, _py_Model model, _py_DistMesh mesh, _py_RNG rng, SSAMethod=_py_SSAMethod.SSA,
            searchMethod=_py_SearchMethod.DIRECT, bool indepKProcs=False):
        """
        Construction::

            sim = steps.solver.DistTetOpSplit(model, mesh, rng)

        Create a distributed spatial stochastic solver based on operator splitting, that is that reaction events are
        partitioned and diffusion is approximated. Keyword parameters SSAMethod and searchMethod respectively set the
        SSA method (SSA or RSSA) and the next event search method (DIRECT or GIBSON_BRUCK).

        Arguments:
        steps.model.Model model
        steps.geom.DistMesh mesh
        steps.rng.RNG rng
        steps.sim.SSAMethod SSAMethod
        steps.sim.NextEventSearchMethod searchMethod
        bool indepKProcs

        """
        if model == None:
            raise TypeError('The Model object is empty.')
        if mesh == None:
            raise TypeError('The Mesh object is empty.')
        if rng == None:
            raise TypeError('The RNG object is empty.')

        if SSAMethod == _py_SSAMethod.SSA:
            if searchMethod == _py_SearchMethod.DIRECT:
                self._ptr = new TetOpSplit[steps_dist_solver.SSAMethod_SSA, steps_dist_solver.NextEventSearchMethod_Direct](
                    deref(model.ptr()), deref(mesh.ptrx()), rng.ptr(), indepKProcs
                )
            elif searchMethod == _py_SearchMethod.GIBSON_BRUCK:
                self._ptr = new TetOpSplit[steps_dist_solver.SSAMethod_SSA, steps_dist_solver.NextEventSearchMethod_GibsonBruck](
                    deref(model.ptr()), deref(mesh.ptrx()), rng.ptr(), indepKProcs
                )
            else:
                raise ValueError(f'Unknown next event search method: {searchMethod}')
        elif SSAMethod == _py_SSAMethod.RSSA:
            self._ptr = new TetOpSplit[steps_dist_solver.SSAMethod_RSSA, steps_dist_solver.NextEventSearchMethod_Direct](
                deref(model.ptr()), deref(mesh.ptrx()), rng.ptr(), indepKProcs
            )
        else:
            raise ValueError(f'Unknown SSA method: {SSAMethod}')

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

    def getCompCount(self, str comp, str spec):
        """
        Returns the number of molecules of a species with identifier string spec 
        in compartment with identifier string comp.

        In a mesh-based simulation this is the combined count from 
        all tetrahedral elements in the compartment.

        Syntax::
            
            getCompCount(comp, spec)
            
        Arguments:
        string comp
        string spec

        Return:
        float

        """
        return self.ptrx().getCompCount(to_std_string(comp), to_std_string(spec))

    def getCompConc(self, str comp, str spec):
        """
        Returns the concentration (in Molar units) of species with identifier string spec 
        in compartment with identifier string comp.

        Note: in a mesh-based simulation this is calculated from the combined 
        number of molecules from all tetrahedral elements in the compartment and the total 
        volume of the tetrahedrons.

        Syntax::
            
            getCompConc(comp, spec)
            
        Arguments:
        string comp
        string spec

        Return:
        float

        """
        return self.ptrx().getCompConc(to_std_string(comp), to_std_string(spec))

    def setCompCount(self, str comp, str spec, double n, distributionMethod=_py_DistributionMethod.UNIFORM):
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
            
            setCompCount(comp, spec, n, distributionMethod)
            
        Arguments:
        string comp
        string spec
        int n
        DistributionMethod distributionMethod

        Return:
        None

        """
        self.ptrx().setCompCount(to_std_string(comp), to_std_string(spec), n, distributionMethod)

    def setCompConc(self, str comp, str spec, double conc, distributionMethod=_py_DistributionMethod.UNIFORM):
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

            setCompConc(comp, spec, conc, distributionMethod)
            
        Arguments:
        string comp
        string spec
        float conc
        DistributionMethod distributionMethod

        Return:
        None

        """
        self.ptrx().setCompConc(to_std_string(comp), to_std_string(spec), conc, distributionMethod)

    def getPatchCount(self, str patch, str spec):
        """
        Returns the number of molecules of species with identifier string spec in patch 
        with identifier string pat.Note: in a mesh-based simulation this 
        is the combined count from all triangular elements in the patch. 

        Syntax::
            
            getPatchCount(patch, spec)
            
        Arguments:
        string patch
        string spec

        Return:
        float

        """
        return self.ptrx().getPatchCount(to_std_string(patch), to_std_string(spec))

    def setPatchCount(self, str patch, str spec, double n, distributionMethod=_py_DistributionMethod.UNIFORM):
        """
        Sets the number of molecules of species with identifier string spec in patch 
        with identifier string pat to n.


        Note: In case of a mesh-based simulation the molecules, molecules are divided among triangles.

        distributionMethod=UNIFORM the distribution is deterministic (apart from roundings) and the number of
        molecules per element is n*V_tet/V_tot.

        distributionMethod=DIST_MULTINOMIAL the distribution is multinomial and the probability
        of putting an element in a tet is V_tet/V_tot

        Syntax::

            setPatchCount(patch, spec, n, distributionMethod)
            
        Arguments:
        string patch
        string spec
        int n
        DistributionMethod distributionMethod

        Return:
        float

        """
        self.ptrx().setPatchCount(to_std_string(patch), to_std_string(spec), n, distributionMethod)

    def getPatchMaxV(self, str patch):
        """
        Returns the maximum potential across a patch

        Syntax::
            
            getPatchMaxV(patch)
            
        Arguments:
        string patch

        Return:
        float

        """
        return self.ptrx().getPatchMaxV(to_std_string(patch))

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
        self.ptrx().setMembPotential(to_std_string(memb), v)

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

    def getSolverDesc(self):
        """
        Returns a string giving a short description of the solver.

        Syntax::

            getSolverDesc()

        Arguments:
        None

        Return:
        string
        """
        return from_std_string(self.ptrx().getSolverDesc())

    def getSolverAuthors(self):
        """
        Returns a string of the solver authors names.

        Syntax::

            getSolverAuthors()

        Arguments:
        None

        Return:
        string
        """
        return from_std_string(self.ptrx().getSolverAuthors())

    def getSolverEmail(self):
        """
        Returns a string giving the author's email address.

        Syntax::

            getSolverEmail()

        Arguments:
        None

        Return:
        string

        """
        return from_std_string(self.ptrx().getSolverEmail())

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

    def getTetCount(self, GO idx, str spec, bool local=False):
        """
        Returns the number of molecules of species with identifier string spec 
        in the tetrahedral element with index idx.

        Syntax::
            
            getTetCount(idx, spec)
            
        Arguments:
        GO idx
        string spec
        bool local

        Return:
        int

        """
        return self.ptrx().getTetCount(idx, to_std_string(spec), local)

    def getTetConc(self, GO idx, str spec, bool local=False):
        """
        Returns the concentration (in Molar units) of species with identifier 
        string spec in a tetrahedral element with index idx.

        Syntax::
            
            getTetConc(idx, spec)
            
        Arguments:
        GO idx
        string spec
        bool local

        Return:
        float

        """
        return self.ptrx().getTetConc(idx, to_std_string(spec), local)

    def setTetCount(self, GO idx, str spec, double n, bool local=False):
        """
        Sets the number of molecules of species with identifier string spec in 
        tetrahedral element with index idx to n.

        Syntax::
            
            setTetCount(idx, spec, n)
            
        Arguments:
        GO idx
        string spec
        int n
        bool local

        Return:
        None

        """
        self.ptrx().setTetCount(idx, to_std_string(spec), n, local)

    def setTetConc(self, GO idx, str spec, double c, bool local=False):
        """
        Sets the concentration (in Molar units) of species with identifier string spec 
        in a tetrahedral element with index idx to conc.This continuous value must be 
        converted internally to a discrete number of molecules. 

        Due to the small volumes of tetrahedral elements the difference between 'rounding 
        up' and 'rounding down' can be a large difference in concentration.

        Syntax::
            
            setTetConc(idx, spec, c)
            
        Arguments:
        GO idx
        string spec
        float c
        bool local

        Return:
        None

        """
        self.ptrx().setTetConc(idx, to_std_string(spec), c, local)

    def getTriCount(self, GO idx, str spec, bool local=False):
        """
        Returns the number of molecules of species with identifier string spec 
        in the triangular element with index idx.

        Syntax::
            
            getTriCount(idx, spec)
            
        Arguments:
        GO idx
        string spec
        bool local

        Return:
        float

        """
        return self.ptrx().getTriCount(idx, to_std_string(spec), local)

    def setTriCount(self, GO idx, str spec, double n, bool local=False):
        """
        Sets the number of molecules of species with identifier string spec in 
        triangular element with index idx to n. 

        Syntax::
            
            setTriCount(idx, spec, n)
            
        Arguments:
        GO idx
        string spec
        int n
        bool local

        Return:
        None

        """
        self.ptrx().setTriCount(idx, to_std_string(spec), n, local)

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
        return self.ptrx().getTriOhmicI(idx, to_std_string(oc), local)

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
        return self.ptrx().getTriGHKI(idx, to_std_string(ghk), local)

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
        self.ptrx().setPatchSReacK(to_std_string(patch), to_std_string(reac), kf)

    def getBatchTetCounts(self, std.vector[GO] tets, str spec, bool local=False):
        """
        Get the counts of a species s in a list of tetrahedrons.

        Syntax::

            getBatchTetCounts(tets, spec)

        Arguments:
        list<GO> tets
        string spec
        bool local

        Return:
        list<double>

        """
        return self.ptrx().getBatchTetCounts(tets, to_std_string(spec), local)

    def setBatchTetCounts(self, std.vector[GO] tets, str spec, std.vector[double] counts, bool local=False):
        """
        Set the counts of a species s in a list of tetrahedrons individually.

        Syntax::

            setBatchTetCounts(tets, spec, counts)

        Arguments:
        list<GO> tets
        string spec
        list<double> counts
        bool local

        Return:
        None

        """
        self.ptrx().setBatchTetCounts(tets, to_std_string(spec), counts, local)

    def getBatchTriCounts(self, std.vector[GO] tris, str spec, bool local=False):
        """
        Get the counts of a species s in a list of triangles.

        Syntax::

            getBatchTriCounts(tris, spec)

        Arguments:
        list<GO> tris
        string spec
        bool local

        Return:
        list<double>

        """
        return self.ptrx().getBatchTriCounts(tris, to_std_string(spec), local)

    def setBatchTriCounts(self, std.vector[GO] tris, str spec, std.vector[double] counts, bool local=False):
        """
        Set the counts of a species s in a list of triangles.

        Syntax::

            setBatchTriCounts(tris, spec, counts)

        Arguments:
        list<GO> tris
        string spec
        list<double> counts
        bool local

        Return:
        None

        """
        return self.ptrx().setBatchTriCounts(tris, to_std_string(spec), counts, local)

    def setBatchTetConcs(self, std.vector[GO] tets, str spec, std.vector[double] concs, bool local=False):
        """
        Set the concentration of a species s in a list of tetrahedrons individually.

        Syntax::

            setBatchTetConcs(tets, spec, concs)

        Arguments:
        list<GO> tets
        string spec
        list<double> concs
        bool local

        Return:
        None

        """
        self.ptrx().setBatchTetConcs(tets, to_std_string(spec), concs, local)

    def getBatchTetConcs(self, std.vector[GO] tets, str spec, bool local=False):
        """
        Get the individual concentration of a species s in a list of tetrahedrons.

        Syntax::

            getBatchTetConcs(tets, spec)

        Arguments:
        list<GO> tets
        string spec
        bool local

        Return:
        list<double>

        """
        return self.ptrx().getBatchTetConcs(tets, to_std_string(spec), local)

    # # ---------------------------------------------------------------------------------
    # # NUMPY section - we accept numpy arrays and generically typed memory-views
    # # ---------------------------------------------------------------------------------
    def getBatchTetCountsNP(self, GO[:] indices, str spec, double[:] counts, bool local=False):
        """
        Get the counts of a species s in a list of tetrahedrons.

        Syntax::
            getBatchTetCountsNP(indices, spec, counts)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> counts
        bool local

        Return:
        None

        """
        self.ptrx().getBatchTetCountsNP(&indices[0], indices.shape[0], to_std_string(spec), &counts[0], counts.shape[0], local)

    def setBatchTetCountsNP(self, GO[:] indices, str spec, double[:] counts, bool local=False):
        """
        Set the counts of a species s in a list of tetrahedrons.

        Syntax::
            setBatchTetCountsNP(indices, spec, counts)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> counts
        bool local

        Return:
        None

        """
        self.ptrx().setBatchTetCountsNP(&indices[0], indices.shape[0], to_std_string(spec), &counts[0],
                counts.shape[0], local)

    def getBatchTetConcsNP(self, GO[:] indices, str spec, double[:] concs, bool local=False):
        """
        Get the individual concentration of a species s in a list of tetrahedrons.

        Syntax::
            getBatchTetCountsNP(indices, spec, concs)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> concs
        bool local

        Return:
        None

        """
        self.ptrx().getBatchTetConcsNP(&indices[0], indices.shape[0], to_std_string(spec), &concs[0], concs.shape[0], local)

    def setBatchTetConcsNP(self, GO[:] indices, str spec, double[:] concs, bool local=False):
        """
        Set the concetration of a species s in a list of tetrahedrons.

        Syntax::
            setBatchTetConcsNP(indices, spec, concs)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> concs
        bool local

        Return:
        None

        """
        self.ptrx().setBatchTetConcsNP(&indices[0], indices.shape[0], to_std_string(spec), &concs[0],
                concs.shape[0], local)

    def getBatchTriCountsNP(self, GO[:] indices, str spec, double[:] counts, bool local=False):
        """
        Get the counts of a species s in a list of triangles.

        Syntax::
            getBatchTriCountsNP(indices, spec, counts)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> counts
        bool local

        Return:
            None

        """
        self.ptrx().getBatchTriCountsNP(&indices[0], indices.shape[0], to_std_string(spec), &counts[0], counts.shape[0], local)

    def setBatchTriCountsNP(self, GO[:] indices, str spec, double[:] counts, bool local=False):
        """
        Set the counts of a species s in a list of triangles.

        Syntax::
            getBatchTriCountsNP(indices, spec, counts)

        Arguments:
        numpy.array<GO> indices
        string spec
        numpy.array<double, length = len(indices)> counts
        bool local

        Return:
            None

        """
        self.ptrx().setBatchTriCountsNP(&indices[0], indices.shape[0], to_std_string(spec), &counts[0],
                counts.shape[0], local)

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
            getBatchTetVsNP(indices, voltages)

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
        self.ptrx().getBatchTriOhmicIsNP(&indices[0], indices.shape[0], to_std_string(oc), &currents[0], currents.shape[0], local)

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
        self.ptrx().getBatchTriGHKIsNP(&indices[0], indices.shape[0], to_std_string(ghk), &currents[0], currents.shape[0], local)

    def setDiffBoundaryDiffusionActive(self, str diffb, str spec, bool act):
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
        self.ptrx().setDiffBoundaryDiffusionActive(to_std_string(diffb), to_std_string(spec), act)

    def getDiffBoundaryDiffusionActive(self, str diffb, str spec):
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
        return self.ptrx().getDiffBoundaryDiffusionActive(to_std_string(diffb), to_std_string(spec))

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

    IF USE_PETSC:
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

        def setEfieldTolerances(self, double atol=1e-50, double rtol=1e-8, KSPNormType norm_type = steps_dist_solver.KSP_NORM_UNPRECONDITIONED):
            """
            Set the absolute and relative tolerances for the E-field solver and on which norm they should be applied.
            See https://petsc.org/release/docs/manual/ksp/#convergence-tests

            is_unpreconditioned_norm=True is on the unpreconditioned norm, False is on the preconditioned norm

            Syntax::

                setEFielTolerances(atol, rtol)

            Arguments:
            float atol
            float rtol
            bool is_unpreconditioned_norm

            Return:
            None

            """
            self.ptrx().setEfieldTolerances(atol, rtol, norm_type)

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
        self.ptrx().setMembIClamp(to_std_string(memb), current)

