# -*- coding: utf-8 -*-

from cython.operator cimport dereference as deref
from libcpp cimport bool
from steps_common cimport *
cimport std
cimport steps_wm
cimport steps
cimport mpi4py.libmpi as libmpi
from steps_common cimport *
from steps_dist cimport *

# ======================================================================================================================
cdef extern from "util/vocabulary.hpp" namespace "steps::dist::model":
# ----------------------------------------------------------------------------------------------------------------------

    cdef cppclass vertgroup_id:
        vertgroup_id(std.string)
        std.string get()

    cdef cppclass patch_id:
        patch_id(std.string)
        std.string get()

    cdef cppclass compartment_id:
        compartment_id(std.string)
        std.string get()

    cdef cppclass membrane_id:
        membrane_id(std.string)
        std.string get()

# ======================================================================================================================
cdef extern from "geom/dist/distmesh.hpp" namespace "Omega_h":
# ----------------------------------------------------------------------------------------------------------------------
    ###### Cybinding for Distmesh ######
    cdef cppclass Library:
        Library(int*, char***, libmpi.MPI_Comm)

# ======================================================================================================================
cdef extern from "geom/dist/measure.hpp" namespace "steps::dist":
# ----------------------------------------------------------------------------------------------------------------------
    ###### Cybinding for Distmesh ######
    cdef cppclass Measure:
        double mesh_measure()
        double rank_measure()

# ======================================================================================================================
cdef extern from "geom/dist/distmesh.hpp" namespace "steps::dist":
# ----------------------------------------------------------------------------------------------------------------------
    ###### Cybinding for Distmesh ######
    cdef cppclass DistMesh:
        # All but a few functions can throw excepts- implement for all but the countXXXs functions
        DistMesh(Library, std.string, double) except +
        @staticmethod
        bool use_gmsh()
        LO num_elems()
        GO total_num_elems()
        LO num_bounds()
        GO total_num_bounds()
        LO num_verts()
        GO total_num_verts()
        void addDiffusionBoundary(diffusion_boundary_name &, compartment_id&, compartment_id&, std.set[triangle_global_id_t]) except +
        void addDiffusionBoundary(diffusion_boundary_name &, compartment_id&, compartment_id&) except +
        DistComp* getTetComp(tetrahedron_global_id_t) except +
        DistComp* getTetComp(tetrahedron_local_id_t) except +
        DistPatch* getTriPatch(triangle_global_id_t) except +
        DistPatch* getTriPatch(triangle_local_id_t) except +
        double getTetVol(tetrahedron_global_id_t) except +
        double getTetVol(tetrahedron_local_id_t) except +
        double getTriArea(triangle_global_id_t) except +
        double getTriArea(triangle_local_id_t) except +
        std.vector[triangle_global_id_t] getSurfTris() except +
        std.vector[triangle_local_id_t] getSurfLocalTris() except +
        std.vector[vertex_global_id_t] getTet_(tetrahedron_global_id_t) except +
        std.vector[vertex_local_id_t] getTet_(tetrahedron_local_id_t) except +
        std.vector[vertex_global_id_t] getTri_(triangle_global_id_t) except +
        std.vector[vertex_local_id_t] getTri_(triangle_local_id_t) except +
        std.vector[double] getVertex(vertex_global_id_t) except +
        std.vector[double] getVertex(vertex_local_id_t) except +
        std.vector[tetrahedron_global_id_t] getTetTetNeighb(tetrahedron_global_id_t) except +
        std.vector[tetrahedron_local_id_t] getTetTetNeighb(tetrahedron_local_id_t) except +
        std.vector[triangle_global_id_t] getTetTriNeighb(tetrahedron_global_id_t) except +
        std.vector[triangle_local_id_t] getTetTriNeighb(tetrahedron_local_id_t) except +
        std.vector[tetrahedron_global_id_t] getTriTetNeighb(triangle_global_id_t) except +
        std.vector[tetrahedron_local_id_t] getTriTetNeighb(triangle_local_id_t) except +
        std.vector[double] getBoundMin(bool) except +
        std.vector[double] getBoundMax(bool) except +
        std.vector[double] getTetBarycenter(tetrahedron_global_id_t) except +
        std.vector[double] getTetBarycenter(tetrahedron_local_id_t) except +
        std.vector[double] getTriBarycenter(triangle_global_id_t) except +
        std.vector[double] getTriBarycenter(triangle_local_id_t) except +
        tetrahedron_global_id_t findTetByPoint(std.vector[double], bool) except +
        std.vector[tetrahedron_global_id_t] getTaggedTetrahedrons(compartment_id&) except +
        std.vector[tetrahedron_local_id_t] getTaggedLocalTetrahedrons(compartment_id&, bool) except +
        std.vector[triangle_global_id_t] getTaggedTriangles(patch_id&) except +
        std.vector[triangle_local_id_t] getTaggedLocalTriangles(patch_id&, bool) except +
        std.vector[vertex_global_id_t] getTaggedVertices(vertgroup_id&) except +
        std.vector[vertex_local_id_t] getTaggedLocalVertices(vertgroup_id&, bool) except +
        tetrahedron_local_id_t getLocalIndex(tetrahedron_global_id_t, bool) except +
        triangle_local_id_t getLocalIndex(triangle_global_id_t, bool) except +
        vertex_local_id_t getLocalIndex(vertex_global_id_t, bool) except +
        GO getGlobalIndex(tetrahedron_local_id_t) except +
        GO getGlobalIndex(triangle_local_id_t) except +
        GO getGlobalIndex(vertex_local_id_t) except +
        std.vector[tetrahedron_global_id_t] getAllTetIndices() except +
        std.vector[tetrahedron_local_id_t] getLocalTetIndices(bool) except +
        std.vector[triangle_global_id_t] getAllTriIndices() except +
        std.vector[triangle_local_id_t] getLocalTriIndices(bool) except +
        std.vector[vertex_global_id_t] getAllVertIndices() except +
        std.vector[vertex_local_id_t] getLocalVertIndices(bool) except +
        double total_measure(compartment_id) except +
        double local_measure(compartment_id) except +
        std.vector[std.string] getTags(int) except +


# ======================================================================================================================
cdef extern from "geom/dist/distcomp.hpp" namespace "steps::dist":
# ----------------------------------------------------------------------------------------------------------------------
    ###### Cybinding for DistComp ######
    cdef cppclass DistComp:
        DistComp(compartment_name, DistMesh, double) except +
        DistComp(compartment_name, DistMesh, compartment_physical_tag, double) except +
        DistComp(compartment_name, DistMesh, std.vector[tetrahedron_global_id_t], double) except +
        std.vector[tetrahedron_global_id_t] getAllTetIndices() except +
        std.vector[tetrahedron_local_id_t] getLocalTetIndices(bool) except +
        double getConductivity()
        void setConductivity(double)
        double getVol() except+
        std.vector[double] getBoundMin(bool) except +
        std.vector[double] getBoundMax(bool) except +

# ======================================================================================================================
cdef extern from "geom/dist/distpatch.hpp" namespace "steps::dist":
# ----------------------------------------------------------------------------------------------------------------------
    ###### Cybinding for DistPatch ######
    cdef cppclass DistPatch:
        DistPatch(patch_name, DistMesh, DistComp*, DistComp*) except +
        DistPatch(patch_name, DistMesh, patch_physical_tag, DistComp*, DistComp*) except +
        DistPatch(patch_name, DistMesh, std.vector[triangle_global_id_t], DistComp*, DistComp*) except +
        std.vector[triangle_global_id_t] getAllTriIndices() except +
        std.vector[triangle_local_id_t] getLocalTriIndices(bool) except +
        double getArea() except+
        std.vector[double] getBoundMin(bool) except +
        std.vector[double] getBoundMax(bool) except +

# ======================================================================================================================
cdef extern from "geom/dist/distmemb.hpp" namespace "steps::dist":
# ----------------------------------------------------------------------------------------------------------------------
    ###### Cybinding for DistMemb ######
    cdef cppclass DistMemb:
        DistMemb(membrane_id, DistMesh, std.set[patch_id], double) except +
        DistMesh* getContainer()
        std.string getID()
        double getCapacitance()
        void setCapacitance(double)
