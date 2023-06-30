#pragma once

#include <algorithm>

#include <Omega_h_array.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_dist.hpp>
#include <Omega_h_mesh.hpp>

#include "geom/dist/distmesh.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"
#include "rng/rng.hpp"
#include "util/collections.hpp"
#include "util/flat_multimap.hpp"
#include "util/vocabulary.hpp"


namespace steps::dist::kproc {

/**
 * Wrapper for
 *   - sum of diffusion rates per element and species, handled by the inherited
 * data-structure
 *   - diffusion rates of every element faces per specie
 *     A diffusion rate provides the probability for each species to
 * cross a specific face of an element. Only the diffusing faces are taken into
 * account (those at the mesh boundaries are excluded). Plus, in a
 * multi-comparments context, elements can have different number of species. 3
 * different vectors are used to store such rate:
 *       - \a super_type::a2ab: element -> start index of its species in
 * elements_to_triangles_
 *       - \a elements_to_triangles_: (element x species) -> start index of
 * their diffusing faces in \a diffusion_rates_
 *       - \a diffusing_rates_: contiguous array of rates
 *
 *
 */
class DiffusionDiscretizedRates: public util::flat_multimap<osh::Real, 1> {
  public:
    using super_type = util::flat_multimap<osh::Real, 1>;
    /**
     * \param neighbors_per_elements vector providing number of neighbors per
     * element
     * \param species_per_elements vector providing number of species per
     * element
     */
    DiffusionDiscretizedRates(const osh::LOs& neighbors_per_elements,
                              const osh::LOs& species_per_elements);

    inline osh::Real& rates_sum(mesh::tetrahedron_id_t element,
                                container::species_id species) noexcept {
        return this->operator()(element.get(), species.get());
    }

    inline osh::Real rates_sum(mesh::tetrahedron_id_t element,
                               container::species_id species) const noexcept {
        return this->operator()(element.get(), species.get());
    }

    inline osh::Real ith_rate(mesh::tetrahedron_id_t element,
                              container::species_id species,
                              int i) const noexcept {
        const auto element_species_index = ab(element.get(), species.get());
        const auto rate_index = elements_to_triangles_[element_species_index] + i;
        return diffusion_rates_[rate_index];
    }
    inline osh::Real& ith_rate(mesh::tetrahedron_id_t element,
                               container::species_id species,
                               int i) noexcept {
        const auto element_species_index = ab(element.get(), species.get());
        const auto rate_index = elements_to_triangles_[element_species_index] + i;
        return diffusion_rates_[rate_index];
    }

    inline const gsl::span<osh::Real> rates(mesh::tetrahedron_id_t element,
                                            container::species_id species) const noexcept {
        const auto lhs_index = elements_to_triangles_[ab(element.get(), species.get())];
        const auto rhs_index = elements_to_triangles_[ab(element.get(), species.get() + 1)];
        return {diffusion_rates_.data() + lhs_index, diffusion_rates_.data() + rhs_index};
    }

    osh::Real rates_max_sum() const;

    inline void reset() {
        assign(0.0);
        osh::fill(diffusion_rates_, 0.0);
    }

  private:
    osh::LOs elements_to_triangles_;
    osh::Write<osh::Real> diffusion_rates_;
};


/**
 * Wrapper for number of molecules transferred to neighbors through boundaries
 */
class PoolsIncrements: public util::flat_multimap<molecules_t, DistMesh::dim() + 1> {
  public:
    using super_type = util::flat_multimap<molecules_t, DistMesh::dim() + 1>;

    PoolsIncrements(DistMesh& mesh, const osh::LOs& elem2num_species)
        : super_type(elem2num_species)
        , synced_delta_pools_(this->ab2c().size())
        , mesh_(mesh) {
        // find the index of the element that is not equal to the previous one
        const bool num_species_all_equal = std::adjacent_find(elem2num_species.begin(),
                                                              elem2num_species.end(),
                                                              std::not_equal_to<>()) ==
                                           elem2num_species.end();
        // 0 is the invalid value (we should not have processes without elements)
        int nspecies = 0;
        if (num_species_all_equal && elem2num_species.size() != 0) {
            nspecies = elem2num_species[0];
        }

        // MPI technique to check with just std MPI functions if nspecies is the same among all the
        // ranks. Quite clever
        int nspecies_MPI[2] = {nspecies, -nspecies};
        int err = MPI_Allreduce(MPI_IN_PLACE, nspecies_MPI, 2, MPI_INT, MPI_MIN, mesh.comm_impl());
        if (err != MPI_SUCCESS) {
            MPI_Abort(mesh_.comm_impl(), err);
        }

        if (nspecies == 0 || nspecies_MPI[0] != -nspecies_MPI[1]) {
            // Warning! sync_delta_pools relies on the fact that dist_ is
            // initialized (and its pointer to the communicator not null) only if
            // this if is true
            dist_ = mesh_.create_dist_for_variable_sized(dims(), sizes2offsets(elem2num_species))
                        .invert();
        }
    }

    PoolsIncrements(const PoolsIncrements&) = delete;

    inline void increment_ith_delta_pool(mesh::tetrahedron_id_t element,
                                         container::species_id species,
                                         int element_face,
                                         molecules_t num_molecules) noexcept {
        this->operator()(element.get(),
                         species.get())[static_cast<size_t>(element_face)] += num_molecules;
    }

    inline void reset() {
        this->assign(0);
    }

    /** Sync the molecule counts among elements (and ghost elements)
     * There are 2 ways to sync species in ghost elements:
     *
     * - sync_array: it is the fastest one and does not do alltoallv. It works
     * only if every element has the same number of species
     * - exch: it uses alltoallv under the hood. This works in any case
     *
     * dist_ is initialized only if we have variable numbers of species.
     * Therefore, checking the pointer to the communicator is a cheap trick to
     * check which function we should use
     */
    inline void sync_delta_pools() {
        if (dist_.comm()) {
            // multiple compartment
            synced_delta_pools_ = dist_.exch(util::createRead(this->ab2c()), 1 /* unused width */);
        } else {
            // same number of species per element
            const auto element_num_values = this->a2ab()[1] - this->a2ab()[0];
            synced_delta_pools_ =
                mesh_.sync_array(dims(), util::createRead(this->ab2c()), element_num_values);
        }
    }

    inline molecules_t ith_delta_pool(mesh::tetrahedron_id_t element,
                                      container::species_id species,
                                      int face) const noexcept {
        const auto index = this->ab(element.get(), species.get());
        return synced_delta_pools_[index + face];
    }

    /**
     * \return Total number of diffusions
     */
    inline molecules_t num_diffusions() const {
        return osh::get_sum(util::createRead(this->ab2c()));
    }

    static constexpr osh::Int dims() noexcept {
        return DistMesh::dim();
    }

  private:
    osh::Read<molecules_t> synced_delta_pools_;
    DistMesh& mesh_;
    osh::Dist dist_;

    static osh::LOs sizes2offsets(osh::LOs sizes) {
        osh::Write<osh::LO> offsets(sizes.size() + 1);
        osh::LO accu = 0;
        std::transform(sizes.begin(), sizes.end(), offsets.begin(), [&accu](const auto size) {
            auto res = accu;
            accu += (dims() + 1) * size;
            return res;
        });
        offsets[offsets.size() - 1] = accu;
        return offsets;
    }
};

/**
 * Functor to compute number of leaving molecules of a certain species from a
 * given triangle/tetrahedron
 */
class LeavingMolecules {
  public:
    explicit LeavingMolecules(rng::RNG& t_rng)
        : rng_(t_rng) {}

    molecules_t operator()(molecules_t num_molecules,
                           osh::Real diffusion_propensity_sum,
                           osh::Real time_delta) const {
        const auto p = std::clamp(time_delta * diffusion_propensity_sum, 0., 1.);
        return molecules_t(rng_.getBinom(static_cast<uint>(num_molecules), p));
    }

  private:
    rng::RNG& rng_;
};

/**
 * Container for all data structures required to simulate the diffusion of
 * molecules
 */
class Diffusions {
  public:
    Diffusions(DistMesh& t_mesh, SimulationInput& t_input);

    inline const DiffusionDiscretizedRates& rates() const noexcept {
        return rates_;
    }

    inline DiffusionDiscretizedRates& rates() noexcept {
        return rates_;
    }

    inline osh::Real& rates_sum(mesh::tetrahedron_id_t element,
                                container::species_id species) noexcept {
        return rates_.rates_sum(element, species);
    }

    inline osh::Real rates_sum(mesh::tetrahedron_id_t element,
                               container::species_id species) const noexcept {
        return rates_.rates_sum(element, species);
    }

    inline osh::Real ith_rate(mesh::tetrahedron_id_t element,
                              container::species_id species,
                              int i) const {
        return rates_.ith_rate(element, species, i);
    }

    inline osh::Real& ith_rate(mesh::tetrahedron_id_t element,
                               container::species_id species,
                               int i) noexcept {
        return rates_.ith_rate(element, species, i);
    }

    inline const LeavingMolecules& total_leaving() const noexcept {
        return total_leaving_;
    }

    inline const PoolsIncrements& leaving_molecules() const noexcept {
        return leaving_molecules_;
    }

    inline PoolsIncrements& leaving_molecules() noexcept {
        return leaving_molecules_;
    }

    inline void reset() {
        leaving_molecules_.reset();
    }

    /// \return maximum sum of rates on all processes
    osh::Real global_rates_max_sum() const;

  private:
    const MPI_Comm comm_;
    DiffusionDiscretizedRates rates_;
    const LeavingMolecules& total_leaving_;
    /// number of molecules leaving from every boundary/faces of every
    /// triangles/tetrahedrons
    PoolsIncrements leaving_molecules_;
};

}  // namespace steps::dist::kproc
