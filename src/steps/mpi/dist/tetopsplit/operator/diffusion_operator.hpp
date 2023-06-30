#pragma once

#include <random>

#include <Omega_h_array.hpp>
#include <Omega_h_shape.hpp>

#include "../kproc/diffusions.hpp"
#include "../kproc/kproc_id.hpp"
#include "../kproc/kproc_state.hpp"
#include "../mol_state.hpp"
#include "geom/dist/fwd.hpp"
#include "rng/rng.hpp"

namespace steps::dist {

class DiffusionOperator {
  public:
    DiffusionOperator(DistMesh& mesh,
                      rng::RNG& t_rng,
                      MolState& t_pools,
                      kproc::Diffusions& t_diffusions,
                      kproc::KProcState& t_kproc_state);
    void operator()(osh::Real opsplit_period, osh::Real state_time);

    inline osh::I64 getExtent() const noexcept {
        return num_diffusions_;
    }

    inline void setBinomialThreshold(osh::I64 threshold) noexcept {
        diffusion_threshold_ = threshold;
    }

  private:
    /** Compute leaving species on all owned elements
     *
     * @param opsplit_period
     * @param state_time
     */
    void species_leaving_elements(osh::Real opsplit_period, osh::Real state_time);

    /**
     * Update state of owned elements to take into account entering species
     *
     * no need for opsplit_period apparently
     */
    void species_entering_elements();

    /** Compute species leaving a given element (triangle/tetrahedron)
     *
     * @param element
     * @param species
     * @param num_molecules
     * @param opsplit_period
     * @param state_time
     */
    void species_leaving_element(mesh::tetrahedron_id_t element,
                                 container::species_id species,
                                 molecules_t num_molecules,
                                 osh::Real opsplit_period,
                                 osh::Real state_time);

    void species_leaving_element_standard(mesh::tetrahedron_id_t element,
                                          container::species_id species,
                                          molecules_t delta_pool_total,
                                          osh::Real scaled_dcst);

    void species_leaving_element_binomial(mesh::tetrahedron_id_t element,
                                          container::species_id species,
                                          molecules_t delta_pool_total,
                                          osh::Real scaled_dcst);

    /** Compute number of molecules leaving an element (triangle/tetrahedron)
     *
     * \param elem an element identifier
     * \param species the molecule specie
     * @param num_molecules
     * @param sum_rates
     * @param opsplit_period
     * @param state_time
     * \return number of molecules leaving
     */
    molecules_t get_leaving_molecules(mesh::tetrahedron_id_t elem,
                                      container::species_id species,
                                      molecules_t num_molecules,
                                      osh::Real sum_rates,
                                      osh::Real opsplit_period,
                                      osh::Real state_time);

    const DistMesh& mesh;
    rng::RNG& rng;
    MolState& pools;
    kproc::Diffusions& diffusions_;
    // Needed to access the dependency_map and relate (elem_id,species) with KProcIDs (optimization
    // on updating the propensities).
    kproc::KProcState& kproc_state_;
    osh::I64 num_diffusions_{};

    osh::I64 diffusion_threshold_{10};
    std::uniform_real_distribution<double> ur_distribution;
};

}  // namespace steps::dist
