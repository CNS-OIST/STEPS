#pragma once

#include <random>

#include <Omega_h_array.hpp>
#include <Omega_h_shape.hpp>

#include "../common.hpp"
#include "../kproc/diffusions.hpp"

namespace zee {

template <osh::Int Dim, typename RNG>
class DiffusionOperator {
  public:
    DiffusionOperator(OmegaHMesh<Dim>& mesh,
                      RNG& t_rng,
                      MolState& t_pools,
                      Diffusions<Dim, RNG>& t_diffusions,
                      osh::Real& t_time_delta);
    void operator()();

    inline PetscInt64 getExtent() const noexcept {
        return num_diffusions_;
    }

    inline void setBinomialThreshold(osh::GO threshold) noexcept {
        diffusion_threshold_ = threshold;
    }

  private:
    /**
     * Compute leaving species on all owned elements
     */
    void species_leaving_elements();

    /**
     * Update state of owned elements to take in to account entering species
     */
    void species_entering_elements();

    /**
     * Compute specie leaving a given element (triangle/tetrahedron)
     * @param element an element identifier
     * @param specie the molecule specie
     */
    void specie_leaving_element(mesh::element_id element,
                                container::specie_id specie,
                                osh::LO num_molecules);

    void specie_leaving_element_standard(mesh::element_id element,
                                         container::specie_id specie,
                                         osh::LO delta_pool_total,
                                         osh::Real scaled_dcst);

    void specie_leaving_element_binomial(mesh::element_id element,
                                         container::specie_id specie,
                                         osh::LO delta_pool_total,
                                         osh::Real scaled_dcst);

    /**
     * Compute number of molecules leaving an element (triangle/tetrahedron)
     * \param elem an element identifier
     * \param specie the molecule specie
     * \return number of molecules leaving
     */
    int get_leaving_molecules(mesh::element_id elem,
                              container::specie_id specie,
                              osh::LO num_molecules,
                              osh::Real sum_rates);


    const OmegaHMesh<Dim>& mesh;
    RNG& rng;
    MolState& pools;
    Diffusions<Dim, RNG>& diffusions_;
    osh::Real& time_delta;
    PetscInt64 num_diffusions_{};

    osh::GO diffusion_threshold_{10};
    const osh::LOs elems2boundary_;
    const osh::Adj boundary2elems_;
    std::uniform_real_distribution<double> ur_distribution;
};

// explicit template instantiation declarations
extern template class DiffusionOperator<2, std::mt19937>;
extern template class DiffusionOperator<3, std::mt19937>;

}  // namespace zee
