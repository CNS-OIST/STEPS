#pragma once

#include <iosfwd>
#include <limits>
#include <type_traits>
#include <variant>

#include "common.hpp"
#include "strong_id.hpp"
#if STEPS_USE_DIST_MESH
#include "strong_ids.hpp"
#endif  // STEPS_USE_DIST_MESH
#include "strong_string.hpp"

namespace steps {

struct tetrahedron_id_trait {};
struct tetrahedron_global_id_trait {};
struct tetrahedron_local_id_trait {};
struct triangle_id_trait {};
struct triangle_global_id_trait {};
struct triangle_local_id_trait {};
struct vertex_id_trait {};
struct bar_id_trait {};  // edge
using tetrahedron_id_t = util::strong_id<index_t, tetrahedron_id_trait>;
using tetrahedron_global_id = util::strong_id<index_t, tetrahedron_global_id_trait>;
using tetrahedron_local_id = util::strong_id<index_t, tetrahedron_local_id_trait>;
using triangle_id_t = util::strong_id<index_t, triangle_id_trait>;
using triangle_global_id = util::strong_id<index_t, triangle_global_id_trait>;
using triangle_local_id = util::strong_id<index_t, triangle_local_id_trait>;
using vertex_id_t = util::strong_id<index_t, vertex_id_trait>;
using bar_id_t = util::strong_id<index_t, bar_id_trait>;  // edge

namespace tetmesh {

// Forward declarations.
class Tetmesh;
class Memb;
class TmPatch;

}  // namespace tetmesh


#ifdef STEPS_USE_DIST_MESH

namespace dist {

/// user API i.e Python bindings
namespace model {
struct tag_species_id {};
struct tag_species_name {};
struct tag_patch_id {};

/// specie name provided by the user
using species_name = util::strong_string<tag_species_name>;

/// internal id generated for the species
using species_id = util::strong_id<osh::I32, tag_species_id>;

/// compartment label used in the mesh
struct tag_compartment_label {};
using compartment_label = util::strong_id<Omega_h::I64, tag_compartment_label>;

/// user compartment identifier
struct tag_compartment_id {};
using compartment_id = util::strong_string<tag_compartment_id>;

/// patch label used in the mesh
using patch_label = util::strong_id<Omega_h::I64, struct tag_patch_label>;

/// patch id used in the mesh
using patch_id = util::strong_string<tag_patch_id>;

// Necessary to access groups of tagged vertices
struct tag_vertgroup_id {};
using vertgroup_id = util::strong_string<tag_vertgroup_id>;

/// membrane identifier given by the user
using membrane_id = util::strong_string<struct tag_membrane_id>;

/// channel identifier given by the user
using channel_id = util::strong_string<struct tag_channel_id>;

/// ohmic current identifier given by the user
using ohmic_current_id = util::strong_string<struct tag_ohmic_current_id>;

/// ghk current identifier given by the user
using ghk_current_id = util::strong_string<struct tag_ghk_current_id>;

/// surface reaction identifier given by the user
using surface_reaction_id = util::strong_string<struct tag_surface_reaction_id>;

using region_id = std::variant<patch_id, compartment_id>;

}  // namespace model

namespace mesh {

using tetrahedron_global_id_t = util::strong_id<osh::I64, tetrahedron_id_trait>;
using triangle_global_id_t = util::strong_id<osh::I64, triangle_id_trait>;
using vertex_global_id_t = util::strong_id<osh::I64, vertex_id_trait>;
using bar_global_id_t = util::strong_id<osh::I64, struct bar_id_trait>;  // edge

/// TODO TCL add host_id_t primitive type

using tetrahedron_local_id_t = util::strong_id<osh::LO, struct tetrahedron_id_trait>;
using triangle_local_id_t = util::strong_id<osh::LO, struct triangle_id_trait>;
using vertex_local_id_t = util::strong_id<osh::LO, struct vertex_id_trait>;
using bar_local_id_t = util::strong_id<osh::LO, struct bar_id_trait>;  // edge

/// TODO TCL add host_id_t primitive type

using tetrahedron_id_t = tetrahedron_local_id_t;
using triangle_id_t = triangle_local_id_t;
using vertex_id_t = vertex_local_id_t;
using bar_id_t = bar_local_id_t;  // edge

struct tag_compartment_name {};
struct tag_compartment_physical_tag {};
struct tag_compartment_id {};
struct tag_patch_name {};
struct tag_patch_physical_tag {};
struct tag_patch_id {};
struct tag_diffusion_boundary_name {};

/// Compartment name given by the user
using compartment_name = util::strong_string<tag_compartment_name>;
using diffusion_boundary_name = util::strong_string<tag_diffusion_boundary_name>;

/// compartment physical tag defined in the mesh
using compartment_physical_tag = util::strong_id<osh::I32, tag_compartment_physical_tag>;

/// internal compartment identifier
using compartment_id = util::strong_id<osh::I32, tag_compartment_id>;

/// patch name given by the user
using patch_name = util::strong_string<tag_patch_name>;

/// patch physical tag defined in the mesh
using patch_physical_tag = util::strong_id<osh::I32, tag_patch_physical_tag>;

/// internal patch identifier
using patch_id = util::strong_id<osh::I32, tag_patch_id>;

using tetrahedron_ids = util::strong_ids<tetrahedron_id_t>;
using triangle_ids = util::strong_ids<triangle_id_t>;
using vertex_ids = util::strong_ids<vertex_id_t>;

using tetrahedron_global_ids = util::strong_ids<tetrahedron_global_id_t>;
using triangle_global_ids = util::strong_ids<triangle_global_id_t>;
using vertex_global_ids = util::strong_ids<vertex_global_id_t>;

/// unused: but I think they are better than something ambiguous
/// I suggest:
/// TODO replace tetrahedron_ids -> tetrahedron_local_ids
using tetrahedron_local_ids = util::strong_ids<tetrahedron_local_id_t>;
using triangle_local_ids = util::strong_ids<triangle_local_id_t>;
using vertex_local_ids = util::strong_ids<vertex_local_id_t>;
/// unused

}  // namespace mesh

/// internal simulation description
namespace container {

struct compartment_id_tag {};
struct species_id_tag {};
struct reaction_id_tag {};
struct diffusion_id_tag {};
struct kproc_id_tag {};
struct patch_id_tag {};

/// internal compartment identifier
using compartment_id = util::strong_id<osh::I32, compartment_id_tag>;

/// internal patch identifier
using patch_id = util::strong_id<osh::I32, patch_id_tag>;

/// internal specie identifier
using species_id = util::strong_id<osh::I32, species_id_tag>;

/// internal reaction identifier
using reaction_id = util::strong_id<osh::I64, reaction_id_tag>;

/// internal diffusion identifier
using diffusion_id = util::strong_id<osh::I64, diffusion_id_tag>;

/// internal surface reaction identifier
using surface_reaction_id = util::strong_id<osh::I64, struct surface_reaction_id_tag>;

/// internal kinetic process identifier
using kproc_id = util::strong_id<osh::I64, kproc_id_tag>;

}  // namespace container

using molecules_t = Omega_h::GO;

/**
 * Get the dimension of an entity according to the dimension of the mesh
 * \tparam Dim the mesh dimension
 * \tparam Entity vocabulary type i.e a strong_id or a strong_string
 */
template <Omega_h::Int Dim, class Entity>
struct entity_dimension {};

template <Omega_h::Int Dim>
struct entity_dimension<Dim, model::vertgroup_id> {
    static const Omega_h::Int value = Dim - 3;
};

template <Omega_h::Int Dim>
struct entity_dimension<Dim, model::patch_id> {
    static const Omega_h::Int value = Dim - 1;
};

template <Omega_h::Int Dim>
struct entity_dimension<Dim, model::compartment_id> {
    static const Omega_h::Int value = Dim;
};

}  // namespace dist

#endif  // !STEPS_USE_DIST_MESH

}  // namespace steps
