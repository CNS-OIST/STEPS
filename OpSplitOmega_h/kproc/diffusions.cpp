#include "diffusions.hpp"

#include "../simulation_data.hpp"

namespace zee {

template <osh::Int Dim>
DiffusionDiscretizedRates<Dim>::DiffusionDiscretizedRates(const osh::LOs& neighbors_per_elements,
                                                          const osh::LOs& species_per_elements)
    : super_type(species_per_elements) {
    osh::LO triangle_idx{};
    osh::LO rate_idx{};
    osh::Write<osh::LO> elements_to_triangles(ab2c().size() + 1);
    for (auto element = 0; element < neighbors_per_elements.size(); ++element) {
        const auto num_neighbors = neighbors_per_elements[element];
        for (auto specie = 0; specie < species_per_elements[element]; ++specie) {
            elements_to_triangles[triangle_idx++] = rate_idx;
            rate_idx += num_neighbors;
        }
    }
    elements_to_triangles[triangle_idx] = rate_idx;
    elements_to_triangles_ = osh::LOs(elements_to_triangles);
    diffusion_rates_ = osh::Write<osh::Real>(rate_idx);
}

template <osh::Int Dim>
PetscScalar DiffusionDiscretizedRates<Dim>::rates_max_sum() const {
    auto max_sum = osh::get_max(osh::Read<osh::Real>(ab2c()));
    auto global_max = 0.0;
    auto ierr = MPI_Allreduce(&max_sum, &global_max, 1, MPIU_SCALAR, MPI_MAX, MPI_COMM_WORLD);
    if (ierr != 0) {
        MPI_Abort(MPI_COMM_WORLD, ierr);
    }
    return global_max;
}

template <osh::LO Dim, class RNG>
Diffusions<Dim, RNG>::Diffusions(OmegaHMesh<Dim>& t_mesh,
                                 SimulationInput<RNG>& t_input,
                                 bool molecules_pools_force_dist_for_variable_sized)
    : rates_(t_mesh.neighbors_per_element(), t_input.species_per_element)
    , dv_(t_input.pools.species_per_elements())
    , total_leaving_(t_input.molecules_leaving)
    , leaving_molecules_(t_mesh.getMesh(),
                         t_input.species_per_element,
                         molecules_pools_force_dist_for_variable_sized) {}

// explicit instantiation definitions
template class DiffusionDiscretizedRates<2>;
template class DiffusionDiscretizedRates<3>;
template class Diffusions<2, std::mt19937>;
template class Diffusions<3, std::mt19937>;

}  // namespace zee
