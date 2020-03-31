#include "diffusion_operator.hpp"

#include <random>

#include <Omega_h_for.hpp>
#include <Omega_h_shape.hpp>

#include "../kproc/diffusions.hpp"
#include "../mesh.hpp"
#include "../mol_state.hpp"

#include "math.hpp"

namespace zee {
template <osh::LO Dim, typename RNG>
DiffusionOperator<Dim, RNG>::DiffusionOperator(OmegaHMesh<Dim>& t_mesh,
                                               RNG& t_rng,
                                               MolState& t_pools,
                                               Diffusions<Dim, RNG>& t_diffusions,
                                               osh::Real& t_time_delta)
    : mesh(t_mesh)
    , rng(t_rng)
    , pools(t_pools)
    , diffusions_(t_diffusions)
    , time_delta(t_time_delta)
    , elems2boundary_(t_mesh.getMesh().ask_down(Dim, Dim - 1).ab2b)
    , boundary2elems_(t_mesh.getMesh().ask_up(Dim - 1, Dim))
    , ur_distribution(0, 1) {}

template <osh::LO Dim, typename RNG>
void DiffusionOperator<Dim, RNG>::operator()() {
    species_leaving_elements();

    num_diffusions_ += diffusions_.leaving_molecules().num_diffusions();
    diffusions_.leaving_molecules().sync_delta_pools();

    species_entering_elements();
}

template <osh::LO Dim, typename RNG>
void DiffusionOperator<Dim, RNG>::species_leaving_elements() {
    const auto lambda = [this](osh::LO elemIdx) __attribute__((always_inline, flatten)) {
        const auto element = mesh.getOwnedElems()[elemIdx];
        const auto num_species = pools.numSpecies(element);
        for (container::specie_id specie(0); specie < num_species; ++specie) {
            const auto num_molecules = pools(element, specie);
            if (num_molecules > 0) {
                specie_leaving_element(element, specie, num_molecules);
            }
        }
    };
    osh::parallel_for(mesh.numOwnedElems(), lambda);
}

template <osh::LO Dim, typename RNG>
void DiffusionOperator<Dim, RNG>::specie_leaving_element(mesh::element_id element,
                                                         container::specie_id specie,
                                                         osh::LO num_molecules) {
    const osh::Real rates_sum = diffusions_.rates_sum(element, specie);
    const auto delta_pool_total =
        this->get_leaving_molecules(element, specie, num_molecules, rates_sum);
    if (delta_pool_total > 0) {
        pools(element, specie) -= delta_pool_total;
        if (delta_pool_total <= diffusion_threshold_) {
            specie_leaving_element_standard(element, specie, delta_pool_total, rates_sum);
        } else {
            specie_leaving_element_binomial(element, specie, delta_pool_total, rates_sum);
        }
    }
}

template <osh::LO Dim, typename RNG>
void DiffusionOperator<Dim, RNG>::specie_leaving_element_standard(mesh::element_id element,
                                                                  container::specie_id specie,
                                                                  osh::LO delta_pool_total,
                                                                  osh::Real scaled_dcst) {
    while (delta_pool_total > 0) {
        const auto selector = ur_distribution(rng) * scaled_dcst;
        PetscScalar partial_sum_scaled_dcst{0};
        const auto& rates = diffusions_.rates().rates(element, specie);
        for (auto direction = 0; direction < static_cast<osh::LO>(rates.size()); ++direction) {
            auto current_scaled_dcst = rates[static_cast<size_t>(direction)];
            partial_sum_scaled_dcst += current_scaled_dcst;
            if (selector < partial_sum_scaled_dcst) {
                diffusions_.leaving_molecules().increment_ith_delta_pool(element,
                                                                         specie,
                                                                         direction,
                                                                         1);
                delta_pool_total -= 1;
                break;
            }
        }
    }
}

template <osh::LO Dim, typename RNG>
void DiffusionOperator<Dim, RNG>::specie_leaving_element_binomial(mesh::element_id element,
                                                                  container::specie_id specie,
                                                                  osh::LO delta_pool_total,
                                                                  osh::Real scaled_dcst) {
    const auto& rates = diffusions_.rates().rates(element, specie);
    for (auto e = 0; e < static_cast<osh::LO>(rates.size()); ++e) {  // loop over boundary/faces
        const auto probability_e = rates[static_cast<size_t>(e)] / scaled_dcst;
        int delta_pool_e = delta_pool_total;
        if (probability_e < 1) {
            std::binomial_distribution<> leaving_through_e(delta_pool_total, probability_e);
            delta_pool_e = leaving_through_e(this->rng);
            delta_pool_total -= delta_pool_e;
            scaled_dcst -= rates[static_cast<size_t>(e)];
        }
        diffusions_.leaving_molecules().increment_ith_delta_pool(element, specie, e, delta_pool_e);
    }
}

/**
 * Functor to update molecules pools according to diffusion increment vector
 */
template <osh::LO Dim>
struct SpeciesEnteringElement {
    SpeciesEnteringElement(const OmegaHMesh<Dim>& t_mesh,
                           MolState& t_pools,
                           PoolsIncrements<Dim>& t_leaving_molecules,
                           const osh::LOs& t_elems2boundary,
                           const osh::Adj& t_boundary2elems)
        : mesh(t_mesh)
        , pools(t_pools)
        , leaving_molecules(t_leaving_molecules)
        , elems2boundary(t_elems2boundary)
        , boundary2elems(t_boundary2elems) {}

    void operator()(mesh::element_id elem) const {
        const auto num_neighbors = mesh.neighbors_per_element()[elem.get()];
        for (auto elem_neighbor_idx = 0; elem_neighbor_idx < num_neighbors; ++elem_neighbor_idx) {
            const auto neighbor_elem = mesh.tet_neighbors_int_data()(elem.get(),
                                                                     elem_neighbor_idx)[0];
            for (auto neighbor_neighbord_idx = 0;
                 neighbor_neighbord_idx < mesh.neighbors_per_element()[neighbor_elem];
                 ++neighbor_neighbord_idx) {
                if (mesh.tet_neighbors_int_data()(neighbor_elem, neighbor_neighbord_idx)[0] ==
                    elem) {
                    for (container::specie_id specie(0); specie < pools.numSpecies(elem);
                         ++specie) {
                        const auto face = mesh.tet_neighbors_int_data()(neighbor_elem,
                                                                        neighbor_neighbord_idx)[1];
                        const auto mols =
                            leaving_molecules.ith_delta_pool(neighbor_elem, specie, face);
                        pools(elem, specie) += mols;
                    }
                    break;
                }
            }
        }
    }

    const OmegaHMesh<Dim>& mesh;
    MolState& pools;
    const PoolsIncrements<Dim>& leaving_molecules;
    const osh::LOs& elems2boundary;
    const osh::Adj& boundary2elems;
};

template <>
struct SpeciesEnteringElement<3> {
    SpeciesEnteringElement(const OmegaHMesh<3>& t_mesh,
                           MolState& t_pools,
                           PoolsIncrements<3>& t_leaving_molecules,
                           const osh::LOs& /* t_elems2boundary */,
                           const osh::Adj& /* t_boundary2elems */)
        : mesh(t_mesh)
        , pools(t_pools)
        , leaving_molecules(t_leaving_molecules) {}

    inline void operator()(mesh::element_id elem) const {
        for (const auto& neighbor: mesh.tet_neighbors_int_data()[elem.get()]) {
            for (container::specie_id specie(0); specie < pools.numSpecies(elem); ++specie) {
                const auto neighbor_id = neighbor[0];
                const auto face_id = neighbor[1];
                pools(elem,
                      specie) += leaving_molecules.ith_delta_pool(neighbor_id, specie, face_id);
            }
        }
    }

  private:
    const OmegaHMesh<3>& mesh;
    MolState& pools;
    const PoolsIncrements<3>& leaving_molecules;
};


template <osh::LO Dim, typename RNG>
void DiffusionOperator<Dim, RNG>::species_entering_elements() {
    const SpeciesEnteringElement<Dim> func_per_element(
        mesh, pools, diffusions_.leaving_molecules(), elems2boundary_, boundary2elems_);

    const auto lambda = [this, &func_per_element](osh::LO ownedElemIdx)
        __attribute__((always_inline, flatten)) {
        const auto elem = mesh.getOwnedElems()[ownedElemIdx];
        func_per_element(mesh::element_id(elem));
    };

    osh::parallel_for(mesh.numOwnedElems(), lambda);
}

template <osh::LO Dim, typename RNG>
int DiffusionOperator<Dim, RNG>::get_leaving_molecules(mesh::element_id elem,
                                                       container::specie_id specie,
                                                       osh::LO num_molecules,
                                                       osh::Real sum_rates) {
    const auto mean_population = stochastic_round<osh::LO>(
        diffusions_.dv().occupancy(elem, specie) / time_delta, this->rng, num_molecules);
    return diffusions_.total_leaving()(mean_population, sum_rates, time_delta);
}


// explicit template instantiation definitions
template class DiffusionOperator<2, std::mt19937>;
template class DiffusionOperator<3, std::mt19937>;

}  // namespace zee
