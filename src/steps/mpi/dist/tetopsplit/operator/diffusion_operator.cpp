#include "diffusion_operator.hpp"

#include <random>

#include <Omega_h_for.hpp>

#include "geom/dist/distmesh.hpp"
#include "math/distributions.hpp"
#include "rng/rng.hpp"
#include "util/profile/profiler_interface.hpp"

namespace steps::dist {

DiffusionOperator::DiffusionOperator(DistMesh& t_mesh,
                                     rng::RNG& t_rng,
                                     MolState& t_pools,
                                     kproc::Diffusions& t_diffusions,
                                     kproc::KProcState& t_kproc_state)
    : mesh(t_mesh)
    , rng(t_rng)
    , pools(t_pools)
    , diffusions_(t_diffusions)
    , kproc_state_(t_kproc_state)
    , ur_distribution(0, 1) {}

void DiffusionOperator::operator()(const osh::Real opsplit_period, const osh::Real state_time) {
    Instrumentor::phase p("DiffusionOperator::operator()");

    diffusions_.reset();
    species_leaving_elements(opsplit_period, state_time);

    num_diffusions_ += diffusions_.leaving_molecules().num_diffusions();

    Instrumentor::phase_begin("sync_delta_pools()");
    diffusions_.leaving_molecules().sync_delta_pools();
    Instrumentor::phase_end("sync_delta_pools()");

    species_entering_elements();
}

void DiffusionOperator::species_leaving_elements(const osh::Real opsplit_period,
                                                 const osh::Real state_time) {
    Instrumentor::phase p("DiffusionOperator::species_leaving_elements()");

    const auto lambda = [this, opsplit_period, state_time](osh::LO elemIdx)
        __attribute__((always_inline, flatten)) {
        const auto element = mesh.owned_elems()[elemIdx];
        for (auto species: pools.species(element)) {
            const auto num_molecules = pools(element, species);
            if (num_molecules > 0) {
                species_leaving_element(
                    element, species, num_molecules, opsplit_period, state_time);
            }
        }
    };
    osh::parallel_for(mesh.owned_elems().size(), lambda, "species_leaving_elements");
}

void DiffusionOperator::species_leaving_element(mesh::tetrahedron_id_t element,
                                                container::species_id species,
                                                molecules_t num_molecules,
                                                const osh::Real opsplit_period,
                                                const osh::Real state_time) {
    const osh::Real rates_sum = diffusions_.rates_sum(element, species);
    const auto delta_pool_total = this->get_leaving_molecules(
        element, species, num_molecules, rates_sum, opsplit_period, state_time);
    if (delta_pool_total > 0) {
        // no need to update occupancy here because channels cannot diffuse and rd occupancy should
        // not track diffusion changes
        pools.add(element, species, -delta_pool_total);
        const auto& dep_map = kproc_state_.get_dependency_map_elems();
        auto ab = pools.moleculesOnElements().ab(element, species);
        const auto& kprocs = dep_map[ab];
        pools.outdated_kprocs().insert(pools.outdated_kprocs().end(), kprocs.begin(), kprocs.end());
        if (delta_pool_total <= diffusion_threshold_) {
            species_leaving_element_standard(element, species, delta_pool_total, rates_sum);
        } else {
            species_leaving_element_binomial(element, species, delta_pool_total, rates_sum);
        }
    }
}


// Possible improvement:
// the STL function "std::discrete_distribution" generates a categorical random variables.
// This function is more efficient when the number of categories (i.e. the size of rates vector) is
// big. For a small vector the current implementation is slightly faster.

void DiffusionOperator::species_leaving_element_standard(mesh::tetrahedron_id_t element,
                                                         container::species_id species,
                                                         molecules_t delta_pool_total,
                                                         osh::Real scaled_dcst) {
    while (delta_pool_total > 0) {
        const auto selector = ur_distribution(rng) * scaled_dcst;
        osh::Real partial_sum_scaled_dcst{0};
        const auto& rates = diffusions_.rates().rates(element, species);
        const auto num_rates = static_cast<osh::LO>(rates.size());
        for (auto direction = 0; direction < num_rates; ++direction) {
            auto current_scaled_dcst = rates[static_cast<size_t>(direction)];
            partial_sum_scaled_dcst += current_scaled_dcst;
            if (selector < partial_sum_scaled_dcst) {
                diffusions_.leaving_molecules().increment_ith_delta_pool(element,
                                                                         species,
                                                                         direction,
                                                                         1);
                delta_pool_total -= 1;
                break;
            }
        }
    }
}

void DiffusionOperator::species_leaving_element_binomial(mesh::tetrahedron_id_t element,
                                                         container::species_id species,
                                                         molecules_t delta_pool_total,
                                                         osh::Real scaled_dcst) {
    const auto& rates = diffusions_.rates().rates(element, species);
    for (auto e = 0; e < static_cast<osh::LO>(rates.size()); ++e) {  // loop over boundary/faces
        const auto probability_e = rates[static_cast<size_t>(e)] / scaled_dcst;
        auto delta_pool_e = delta_pool_total;
        if (probability_e < 1) {
            delta_pool_e = static_cast<molecules_t>(
                this->rng.getBinom(static_cast<uint>(delta_pool_total),
                                   static_cast<double>(probability_e)));
            delta_pool_total -= delta_pool_e;
            scaled_dcst -= rates[static_cast<size_t>(e)];
        }
        diffusions_.leaving_molecules().increment_ith_delta_pool(element, species, e, delta_pool_e);
    }
}

/**
 * Functor to update molecules pools according to diffusion increment vector
 */

struct SpeciesEnteringElement {
    SpeciesEnteringElement(const DistMesh& t_mesh,
                           MolState& t_pools,
                           kproc::PoolsIncrements& t_leaving_molecules,
                           kproc::KProcState& t_kproc_state)
        : mesh(t_mesh)
        , pools(t_pools)
        , leaving_molecules(t_leaving_molecules)
        , kproc_state_(t_kproc_state) {}

    inline void operator()(mesh::tetrahedron_id_t elem) const {
        const auto& dep_map = kproc_state_.get_dependency_map_elems();
        auto num_elements = mesh.tet_neighbors_int_data().size(elem.get());
        auto elem_comp_id = mesh.getCompartmentMeshID(elem);
        for (auto face_id = 0; face_id < num_elements; face_id++) {
            const auto& neighbor = mesh.tet_neighbors_int_data()(elem.get(), face_id);
            const auto neighbor_id = mesh::tetrahedron_id_t(neighbor[0]);
            auto comp_id = mesh.getCompartmentMeshID(neighbor_id);
            const auto neighbour_face_id = neighbor[1];
            const auto triangle_id = mesh::triangle_id_t(neighbor[2]);
            if (comp_id == elem_comp_id) {
                // Careful, pools.species is only defined for an owned element.
                for (osh::LO sp = 0; sp < leaving_molecules.size(neighbor_id.get()); ++sp) {
                    container::species_id species(sp);
                    auto num_mols = leaving_molecules.ith_delta_pool(
                        mesh::tetrahedron_id_t(neighbor_id), species, neighbour_face_id);
                    // no need to update occupancy here because channels cannot diffuse and rd
                    // occupancy should not track diffusion changes
                    pools.add(elem, species, num_mols);
                    if (num_mols > 0) {
                        auto ab = pools.moleculesOnElements().ab(elem, species);
                        const auto& kprocs = dep_map[ab];
                        pools.outdated_kprocs().insert(pools.outdated_kprocs().end(),
                                                       kprocs.begin(),
                                                       kprocs.end());
                    }
                }
            } else {
                for (osh::LO sp = 0; sp < leaving_molecules.size(neighbor_id.get()); ++sp) {
                    container::species_id species(sp);
                    auto num_mols = leaving_molecules.ith_delta_pool(
                        mesh::tetrahedron_id_t(neighbor_id), species, neighbour_face_id);
                    if (num_mols > 0) {
                        if (mesh.isActiveDiffusionBoundary(triangle_id, comp_id, species)) {
                            auto sID = mesh.convertSpeciesID(triangle_id, comp_id, species);
                            // no need to update occupancy here because channels cannot diffuse and
                            // rd occupancy should not track diffusion changes
                            pools.add(elem, sID, num_mols);
                            auto ab = pools.moleculesOnElements().ab(elem, sID);
                            const auto& kprocs = dep_map[ab];
                            pools.outdated_kprocs().insert(pools.outdated_kprocs().end(),
                                                           kprocs.begin(),
                                                           kprocs.end());
                        } else {
                            throw std::logic_error("Something is wrong here.");
                        }
                    }
                }
            }
        }
    }

  private:
    const DistMesh& mesh;
    MolState& pools;
    const kproc::PoolsIncrements& leaving_molecules;
    kproc::KProcState& kproc_state_;
};

void DiffusionOperator::species_entering_elements() {
    Instrumentor::phase p("DiffusionOperator::species_entering_elements()");

    const SpeciesEnteringElement func_per_element(mesh,
                                                  pools,
                                                  diffusions_.leaving_molecules(),
                                                  kproc_state_);

    const auto lambda = [this, &func_per_element](osh::LO ownedElemIdx)
        __attribute__((always_inline, flatten)) {
        const auto elem = mesh.owned_elems()[ownedElemIdx];
        func_per_element(mesh::tetrahedron_id_t(elem));
    };

    osh::parallel_for(mesh.owned_elems().size(), lambda, "species_entering_elements");
}

molecules_t DiffusionOperator::get_leaving_molecules(mesh::tetrahedron_id_t elem,
                                                     container::species_id species,
                                                     molecules_t num_molecules,
                                                     osh::Real sum_rates,
                                                     const osh::Real opsplit_period,
                                                     const osh::Real state_time) {
    // if probability is 0 we do not diffuse
    if (sum_rates == 0) {
        return 0;
    }

    // diffusions happen at the end of the rd_dt time step. Event_time = state_time + opsplit_period
    const auto occupancy = pools.get_occupancy_rd(elem, species, state_time + opsplit_period);

    const auto mean_population = rng.stochastic_round<molecules_t>(occupancy, num_molecules);
    return diffusions_.total_leaving()(mean_population, sum_rates, opsplit_period);
}

}  // namespace steps::dist
