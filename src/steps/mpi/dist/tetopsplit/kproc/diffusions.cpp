#include "diffusions.hpp"

#include "../simulation_data.hpp"
#include "geom/dist/distmesh.hpp"

namespace steps::dist::kproc {

DiffusionDiscretizedRates::DiffusionDiscretizedRates(const osh::LOs& neighbors_per_elements,
                                                     const osh::LOs& species_per_elements)
    : super_type(species_per_elements) {
    osh::LO triangle_idx{};
    osh::LO rate_idx{};
    osh::Write<osh::LO> elements_to_triangles(ab2c().size() + 1);
    for (auto element = 0; element < neighbors_per_elements.size(); ++element) {
        const auto num_neighbors = neighbors_per_elements[element];
        for (auto species = 0; species < species_per_elements[element]; ++species) {
            elements_to_triangles[triangle_idx++] = rate_idx;
            rate_idx += num_neighbors;
        }
    }
    elements_to_triangles[triangle_idx] = rate_idx;
    elements_to_triangles_ = osh::LOs(elements_to_triangles);
    diffusion_rates_ = osh::Write<osh::Real>(rate_idx, 0.0);
}

osh::Real DiffusionDiscretizedRates::rates_max_sum() const {
    return osh::get_max(osh::Reals(ab2c()));
}

Diffusions::Diffusions(DistMesh& t_mesh, SimulationInput& t_input)
    : comm_(t_mesh.comm_impl())
    , rates_(t_mesh.neighbors_per_element(),
             t_input.pools.species_per_elements() /*elements owned strictly*/)
    , total_leaving_(t_input.molecules_leaving)
    , leaving_molecules_(t_mesh, t_input.species_per_element /*elements owned or not*/) {}

osh::Real Diffusions::global_rates_max_sum() const {
    osh::Real global_max_sum;
    const osh::Real max_sum = this->rates().rates_max_sum();
    auto err = MPI_Allreduce(&max_sum, &global_max_sum, 1, MPI_DOUBLE, MPI_MAX, this->comm_);
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm_, err);
    }
    return global_max_sum;
}

}  // namespace steps::dist::kproc
