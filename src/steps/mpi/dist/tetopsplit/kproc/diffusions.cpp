#include "diffusions.hpp"

#include "../simulation_data.hpp"
#include "geom/dist/distmesh.hpp"
#include "model/diff.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "util/vocabulary.hpp"

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

osh::Real DiffusionDiscretizedRates::rates_sum_sum() const {
    return osh::get_sum(osh::Reals(ab2c()));
}

osh::Real DiffusionDiscretizedRates::rates_min_sum() const {
    return osh::get_min(osh::Reals(ab2c()));
}

Diffusions::Diffusions(DistMesh& t_mesh, const Statedef& statedef, SimulationInput& t_input)
    : comm_(t_mesh.comm_impl())
    , mesh_(t_mesh)
    , statedef_(statedef)
    , rates_(t_mesh.neighbors_per_element(), t_input.species_per_element)
    , total_leaving_(t_input.molecules_leaving)
    , leaving_molecules_(t_mesh, t_input.species_per_element /*elements owned or not*/) {}

osh::Real Diffusions::get_tet_dcst(mesh::tetrahedron_local_id_t element,
                                   container::diffusion_id diff,
                                   int i) const {
    // Try to get diffusion constant value from most specific to least specific

    // Check diffusion constant in tetrahedrons / face pairs
    auto it = dcsts_.find({element, diff, i});
    if (it != dcsts_.end()) {
        return it->second;
    }
    const auto compartment_id = mesh_.getCompartment(element);
    const auto& compartment = statedef_.getCompdef(compartment_id);
    const auto& diffusion = compartment.diffdefs().at(diff.get());

    if (i >= 0) {
        // Check diffusion constant from diffusion boundaries
        auto compartment_mid = mesh_.getRegionMeshID(element);
        auto d = mesh_.tet_neighbors_int_data()(element.get(), i);
        if (mesh_.getRegionMeshID(mesh::tetrahedron_id_t(d[0])) != compartment_mid) {
            auto diffbDcst = this->mesh_.getDiffusionBoundaryDcst(mesh::triangle_id_t(d[2]),
                                                                  compartment_mid,
                                                                  diffusion->getSpecContainerIdx());
            if (diffbDcst >= 0) {
                return diffbDcst;
            }
        }

        // Check diffusion constant in tetrahedron (without a specific face idx)
        auto it2 = dcsts_.find({element, diff, -1});
        if (it2 != dcsts_.end()) {
            return it2->second;
        }
    }

    // Check diffusion constant in compartment
    return diffusion->getDcst();
}

void Diffusions::clear_tet_dcst(mesh::tetrahedron_local_id_t element,
                                container::diffusion_id diff,
                                int i) {
    dcsts_.erase({element, diff, i});
}

void Diffusions::initialize_discretized_rates() {
    const auto& measure_info = mesh_.getMeasure();
    leaving_molecules_reset();
    osh::parallel_for(
        mesh_.num_local_elems(),
        [&measure_info, this](osh::LO e) __attribute__((always_inline, flatten)) {
            const mesh::tetrahedron_id_t elem(e);
            const auto elem_measure = measure_info.element_measure(elem);
            const auto compartment_id = this->mesh_.getCompartment(elem);
            const auto compartment_mid = this->mesh_.getRegionMeshID(elem);
            const auto& compartment = this->statedef_.getCompdef(compartment_id);
            for (const auto& diffusion: compartment.diffdefs()) {
                const auto species = diffusion->getSpecContainerIdx();
                this->rates_sum(elem, species) = 0;
                const auto num_neighbors = this->mesh_.tet_neighbors_int_data().size(elem.get());
                for (auto face_idx = 0; face_idx < num_neighbors; ++face_idx) {
                    auto d = this->mesh_.tet_neighbors_int_data()(elem.get(), face_idx);
                    auto compartment_mid2 = this->mesh_.getRegionMeshID(
                        mesh::tetrahedron_id_t(d[0]));

                    if (compartment_mid2 == compartment_mid ||
                        this->mesh_.getDiffusionBoundaryDcst(mesh::triangle_id_t(d[2]),
                                                             compartment_mid,
                                                             species) != 0.0) {
                        const auto neighbor_boundary_distance =
                            this->mesh_.tet_neighbors_real_data()(elem.get(), face_idx)[0];
                        const auto neighbor_boundary_measure =
                            this->mesh_.tet_neighbors_real_data()(elem.get(), face_idx)[1];
                        const auto dcst =
                            get_tet_dcst(elem, diffusion->getDiffContainerIdx(), face_idx);
                        const auto propensity = dcst * neighbor_boundary_measure / elem_measure /
                                                neighbor_boundary_distance;
                        this->ith_rate(elem, species, face_idx) = propensity;
                        this->rates_sum(elem, species) += propensity;
                    }
                }
            }
        });
}

osh::Real Diffusions::global_rates_max_sum() const {
    osh::Real global_max_sum;
    const osh::Real max_sum = this->rates().rates_max_sum();
    auto err = MPI_Allreduce(&max_sum, &global_max_sum, 1, MPI_DOUBLE, MPI_MAX, this->comm_);
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm_, err);
    }
    return global_max_sum;
}

osh::Real Diffusions::global_rates_mean_sum() const {
    osh::Real global_sum_sum;
    const osh::Real sum_sum = this->rates().rates_sum_sum();
    auto err = MPI_Allreduce(&sum_sum, &global_sum_sum, 1, MPI_DOUBLE, MPI_SUM, this->comm_);
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm_, err);
    }
    return global_sum_sum / static_cast<osh::Real>(mesh_.total_num_elems());
}

osh::Real Diffusions::global_rates_min_sum() const {
    osh::Real global_min_sum;
    const osh::Real min_sum = this->rates().rates_min_sum();
    auto err = MPI_Allreduce(&min_sum, &global_min_sum, 1, MPI_DOUBLE, MPI_MIN, this->comm_);
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm_, err);
    }
    return global_min_sum;
}

}  // namespace steps::dist::kproc
