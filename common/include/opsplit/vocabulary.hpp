#pragma once

#include <Omega_h_defines.hpp>
#include <boost/variant.hpp>
#include <petscsys.h>

#include "strong_id.hpp"
#include "strong_ids.hpp"
#include "strong_string.hpp"


namespace zee {

/// user API i.e Python bindings
namespace model {

/// compartment label used in the mesh
using compartment_label = strong_id<PetscInt, struct tag_compartment_label>;

/// user compartment identifier
using compartment_id = strong_string<struct tag_compartment_id>;

/// specie identifier given to the user
using specie_id = strong_id<PetscInt, struct tag_specie_id>;

/// specie name given by the user
using specie_name = strong_string<struct tag_specie_name>;

/// patch label used in the mesh
using patch_label = strong_id<PetscInt, struct tag_patch_label>;

/// patch id used in the mesh
using patch_id = strong_string<struct tag_patch_id>;

/// region id either a patch id or a compartment id
using region_id = boost::variant<model::patch_id, model::compartment_id>;

}  // namespace model

/// internal simulation description
namespace container {

/// internal patch identifier
using patch_id = strong_id<PetscInt, struct patch_id_tag>;

/// internal compartment identifier
using compartment_id = strong_id<PetscInt, struct compartment_id_tag>;

/// internal specie identifier
using specie_id = strong_id<int, struct specie_id_tag>;

/// internal reaction identifier
using reaction_id = strong_id<PetscInt, struct reaction_id_tag>;

/// internal surface reaction identifier
using surface_reaction_id = strong_id<PetscInt, struct surface_reaction_id_tag>;

/// internal diffusion identifier
using diffusion_id = strong_id<PetscInt, struct diffusion_id_tag>;

/// internal kinetic process identifier
using kproc_id = strong_id<PetscInt, struct kproc_id_tag>;

}  // namespace container

namespace mesh {

/// internal compartment identifier
using compartment_id = strong_id<int, struct compartment_id_tag>;

using element_id = strong_id<Omega_h::LO, struct element_id_tag>;
using boundary_id = strong_id<Omega_h::LO, struct boundary_id_tag>;

using element_ids = strong_ids<element_id>;
using boundary_ids = strong_ids<boundary_id>;

}  // namespace mesh

}  // namespace zee
