#pragma once

#include <random>

#include <Omega_h_array.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_dist.hpp>
#include <Omega_h_mesh.hpp>

#include "../common.hpp"
#include "flat_multimap.hpp"
#include "opsplit/vocabulary.hpp"

namespace zee {

/**
 * Wrapper for
 *   - sum of diffusion rates per element and specie, handled by the inherited data-structure
 *   - diffusion rates of every element faces per specie
 *     A diffusion rate provides the probability for each specie specie to cross a specific
 *     face of an element.
 *     Only the diffusing faces are taken into account (those at the mesh boundaries are excluded).
 *     Plus, in a multi-comparments context, elements can have different number of species.
 *     3 different vectors are used to store such rate:
 *       - \a super_type::a2ab: element -> start index of its species in elements_to_triangles_
 *       - \a elements_to_triangles_: (element x specie) -> start index of their diffusing faces
 *       in \a diffusion_rates_
 *       - \a diffusing_rates_: contiguous array of rates
 *
 *
 */
template <osh::Int Dim>
class DiffusionDiscretizedRates: public zee::flat_multimap<osh::Real, 1, zee::OSH> {
  public:
    using super_type = zee::flat_multimap<osh::Real, 1, zee::OSH>;
    /**
     * \param neighbors_per_elements vector providing number of neighbors per element
     * \param species_per_elements vector providing number of species per element
     */
    DiffusionDiscretizedRates(const osh::LOs& neighbors_per_elements,
                              const osh::LOs& species_per_elements);

    inline osh::Real& rates_sum(mesh::element_id element, container::specie_id specie) noexcept {
        return this->operator()(element.get(), specie.get());
    }

    inline osh::Real rates_sum(mesh::element_id element, container::specie_id specie) const
        noexcept {
        return this->operator()(element.get(), specie.get());
    }

    inline osh::Real ith_rate(mesh::element_id element, container::specie_id specie, int i) const
        noexcept {
        const auto element_specie_index = ab(element.get(), specie.get());
        const auto rate_index = elements_to_triangles_[element_specie_index] + i;
        return diffusion_rates_[rate_index];
    }
    inline osh::Real& ith_rate(mesh::element_id element,
                               container::specie_id specie,
                               int i) noexcept {
        const auto element_specie_index = ab(element.get(), specie.get());
        const auto rate_index = elements_to_triangles_[element_specie_index] + i;
        return diffusion_rates_[rate_index];
    }

    inline const gsl::span<osh::Real> rates(mesh::element_id element,
                                            container::specie_id specie) const noexcept {
        const auto lhs_index = elements_to_triangles_[ab(element.get(), specie.get())];
        const auto rhs_index = elements_to_triangles_[ab(element.get(), specie.get() + 1)];
        return {diffusion_rates_.data() + lhs_index, diffusion_rates_.data() + rhs_index};
    }

    PetscScalar rates_max_sum() const;

  private:
    osh::LOs elements_to_triangles_;
    osh::Write<osh::Real> diffusion_rates_;
};

/**
 * Wrapper for number of molecules transferred to neighbours
 */
template <osh::Int dim>
class PoolsIncrements: public zee::flat_multimap<osh::LO, dim + 1, zee::OSH> {
  public:
    using super_type = zee::flat_multimap<osh::LO, dim + 1, zee::OSH>;

    PoolsIncrements(osh::Mesh& mesh,
                    const osh::LOs& elem2num_species,
                    bool force_dist_for_variable_sized)
        : super_type(elem2num_species)
        , synced_delta_pools_(this->ab2c().size())
        , mesh_(mesh) {
        const auto num_species_all_equal = std::adjacent_find(elem2num_species.begin(),
                                                              elem2num_species.end(),
                                                              std::not_equal_to<>()) ==
                                           elem2num_species.end();
        if (force_dist_for_variable_sized || !num_species_all_equal) {
            dist_ = osh::create_dist_for_variable_sized(mesh.ask_dist(dim),
                                                        sizes2offsets(elem2num_species))
                        .invert();
        }
    }


    PoolsIncrements(const PoolsIncrements&) = delete;

    inline void increment_ith_delta_pool(mesh::element_id element,
                                         container::specie_id specie,
                                         int element_face,
                                         osh::LO num_molecules) noexcept {
        this->operator()(element.get(),
                         specie.get())[static_cast<size_t>(element_face)] += num_molecules;
    }

    inline void reset() {
        this->assign(0);
    }

    inline void sync_delta_pools() {
        if (dist_.comm()) {
            // multiple compartment
            synced_delta_pools_ = dist_.exch(osh::LOs(this->ab2c()), 1 /* unused width */);
        } else {
            // same number of species per element
            const auto element_num_values = this->a2ab()[1] - this->a2ab()[0];
            synced_delta_pools_ = mesh_.sync_array(dim, osh::LOs(this->ab2c()), element_num_values);
        }
    }

    inline osh::LO ith_delta_pool(osh::LO element_id, container::specie_id specie, int i) const
        noexcept {
        const auto index = this->ab(element_id, specie.get());
        return synced_delta_pools_[index + i];
    }

    /**
     * \return Total number of diffusions
     */
    inline osh::LO num_diffusions() const {
        return osh::get_sum(osh::LOs(this->ab2c()));
    }

    constexpr osh::Int dims() const {
        return dim;
    }

  private:
    osh::LOs synced_delta_pools_;
    osh::Mesh& mesh_;
    osh::Dist dist_;

    static osh::LOs sizes2offsets(osh::LOs sizes) {
        osh::Write<osh::LO> offsets(sizes.size() + 1);
        int accu = 0;
        std::transform(sizes.begin(), sizes.end(), offsets.begin(), [&accu](const auto size) {
            auto res = accu;
            accu += (dim + 1) * size;
            return res;
        });
        offsets[offsets.size() - 1] = accu;
        return offsets;
    }
};

/**
 * Functor to compute number of leaving molecules of a certain specie from a given
 * triangle/tetrahedron
 */
template <typename RNG>
class LeavingMolecules {
  public:
    LeavingMolecules(RNG& t_rng)
        : rng_(t_rng) {}

    int operator()(osh::LO num_molecules,
                   osh::Real diffusion_propensity_sum,
                   osh::Real time_delta) const {
        std::binomial_distribution<> total_leaving(num_molecules,
                                                   time_delta * diffusion_propensity_sum);
        return total_leaving(rng_);
    }

  private:
    RNG& rng_;
};

/**
 * Wrapper for number of molecules transferred to neighbours
 */
class DiffusionVariables {
  public:
    DiffusionVariables(const osh::LOs& species_per_elements)
        : variables_(species_per_elements) {}

    inline osh::Real occupancy(mesh::element_id element, container::specie_id specie) const
        noexcept {
        return variables_(element.get(), specie.get())[0];
    }

    inline osh::Real& occupancy(mesh::element_id element, container::specie_id specie) noexcept {
        return variables_(element.get(), specie.get())[0];
    }

    inline osh::Real last_update_time(mesh::element_id element, container::specie_id specie) const
        noexcept {
        return variables_(element.get(), specie.get())[1];
    }

    inline osh::Real& last_update_time(mesh::element_id element,
                                       container::specie_id specie) noexcept {
        return variables_(element.get(), specie.get())[1];
    }

    inline void reset() {
        variables_.assign(0);
    }

  private:
    zee::flat_multimap<osh::Real, 2, zee::OSH> variables_;
};

/**
 * Container for all data structures required to simulate the diffusion of molecules
 * \tparam Dim Mesh dimension
 * \tparam RNG Random number generator functor
 */
template <osh::LO Dim, typename RNG>
class Diffusions {
  public:
    Diffusions(OmegaHMesh<Dim>& t_mesh,
               SimulationInput<RNG>& t_input,
               bool molecules_pools_force_dist_for_variable_sized);

    inline const DiffusionDiscretizedRates<Dim>& rates() const noexcept {
        return rates_;
    }

    inline DiffusionDiscretizedRates<Dim>& rates() noexcept {
        return rates_;
    }

    inline osh::Real& rates_sum(mesh::element_id element, container::specie_id specie) noexcept {
        return rates_.rates_sum(element, specie);
    }

    inline osh::Real rates_sum(mesh::element_id element, container::specie_id specie) const
        noexcept {
        return rates_.rates_sum(element, specie);
    }

    inline osh::Real ith_rate(mesh::element_id element, container::specie_id specie, int i) const {
        return rates_.ith_rate(element, specie, i);
    }

    inline osh::Real& ith_rate(mesh::element_id element,
                               container::specie_id specie,
                               int i) noexcept {
        return rates_.ith_rate(element, specie, i);
    }

    inline DiffusionVariables& dv() noexcept {
        return dv_;
    }

    inline const LeavingMolecules<RNG>& total_leaving() const noexcept {
        return total_leaving_;
    }

    inline const PoolsIncrements<Dim>& leaving_molecules() const noexcept {
        return leaving_molecules_;
    }

    inline PoolsIncrements<Dim>& leaving_molecules() noexcept {
        return leaving_molecules_;
    }

    inline void reset() {
        leaving_molecules_.reset();
        dv_.reset();
    }

  private:
    DiffusionDiscretizedRates<Dim> rates_;
    DiffusionVariables dv_;
    const LeavingMolecules<RNG>& total_leaving_;
    /// number of molecules leaving from every boundary/faces of every triangles/tetrahedrons
    PoolsIncrements<Dim> leaving_molecules_;
};

// explicit instantiation declarations
extern template class DiffusionDiscretizedRates<2>;
extern template class DiffusionDiscretizedRates<3>;
extern template class Diffusions<2, std::mt19937>;
extern template class Diffusions<3, std::mt19937>;

}  // namespace zee
