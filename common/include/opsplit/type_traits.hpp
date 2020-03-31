#pragma once

#include "vocabulary.hpp"

namespace zee {

/**
 * Get the dimension of an entity according to the dimension of the mesh
 * \tparam Dim the mesh dimension
 * \tparam Entity vocabulary type i.e a strong_id or a strong_string
 */
template <Omega_h::Int Dim, class Entity>
struct entity_dimension {};

template <Omega_h::Int Dim>
struct entity_dimension<Dim, model::patch_id> {
    static const Omega_h::Int value = Dim - 1;
};

template <Omega_h::Int Dim>
struct entity_dimension<Dim, model::compartment_id> {
    static const Omega_h::Int value = Dim;
};

}  // namespace zee
