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
        @staticmethod
        int unknown_value()

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

    cdef cppclass bar_global_id_t:
        bar_global_id_t()
        bar_global_id_t(int)
        int get()
        bool valid()

    cdef cppclass bar_local_id_t:
        bar_local_id_t()
        bar_local_id_t(int)
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

    cdef cppclass species_name:
        species_name(std.string)
        std.string get()

    cdef cppclass complex_name:
        complex_name(std.string)
        std.string get()

    cdef cppclass complex_substate_id:
        complex_substate_id(int)
        int get()

    cdef cppclass ohmic_current_id:
        ohmic_current_id(std.string)
        std.string get()

    cdef cppclass complex_ohmic_current_id:
        complex_ohmic_current_id(std.string)
        std.string get()

    cdef cppclass ghk_current_id:
        ghk_current_id(std.string)
        std.string get()

    cdef cppclass complex_ghk_current_id:
        complex_ghk_current_id(std.string)
        std.string get()

    cdef cppclass reaction_id:
        reaction_id(std.string)
        std.string get()

    cdef cppclass complex_reaction_id:
        complex_reaction_id(std.string)
        std.string get()

    cdef cppclass surface_reaction_id:
        surface_reaction_id(std.string)
        std.string get()

    cdef cppclass complex_surface_reaction_id:
        complex_surface_reaction_id(std.string)
        std.string get()

    cdef cppclass vdep_surface_reaction_id:
        vdep_surface_reaction_id(std.string)
        std.string get()

    cdef cppclass vdep_complex_surface_reaction_id:
        vdep_complex_surface_reaction_id(std.string)
        std.string get()

    cdef cppclass diffusion_id:
        diffusion_id(std.string)
        std.string get()

