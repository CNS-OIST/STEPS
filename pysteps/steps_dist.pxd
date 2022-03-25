cimport std
from libcpp cimport bool

cdef extern from "Omega_h_defines.hpp" namespace "Omega_h":
    ctypedef long long int GO;
    ctypedef int LO;


# ======================================================================================================================
cdef extern from "util/vocabulary.hpp" namespace "steps::dist::mesh":
# ----------------------------------------------------------------------------------------------------------------------

    cdef cppclass compartment_name:
        compartment_name(std.string)
        std.string get()

    cdef cppclass tetrahedron_global_id_t:
        tetrahedron_global_id_t()
        tetrahedron_global_id_t(int)
        int get()
        bool valid()

    cdef cppclass tetrahedron_local_id_t:
        tetrahedron_local_id_t()
        tetrahedron_local_id_t(int)
        int get()
        bool valid()

    cdef cppclass compartment_physical_tag:
        compartment_physical_tag(int)
        int get()

    cdef cppclass patch_name:
        patch_name(std.string)
        std.string get()

    cdef cppclass triangle_global_id_t:
        triangle_global_id_t()
        triangle_global_id_t(int)
        int get()
        bool valid()

    cdef cppclass triangle_local_id_t:
        triangle_local_id_t()
        triangle_local_id_t(int)
        int get()
        bool valid()

    cdef cppclass vertex_global_id_t:
        vertex_global_id_t()
        vertex_global_id_t(int)
        int get()
        bool valid()

    cdef cppclass vertex_local_id_t:
        vertex_local_id_t()
        vertex_local_id_t(int)
        int get()
        bool valid()

    cdef cppclass patch_physical_tag:
        patch_physical_tag(int)
        int get()

    cdef cppclass diffusion_boundary_name:
        diffusion_boundary_name(std.string)
        std.string get()
