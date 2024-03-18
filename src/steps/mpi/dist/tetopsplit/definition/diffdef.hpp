#pragma once

#include "fwd.hpp"

#include "util/vocabulary.hpp"

namespace steps::dist {

/**
 * \brief State definition of a diffusion.
 *
 * The Diffdef class defines the sub biochemical container of
 * a diffusion in a compartment.
 *
 * This class corresponds to the solver::Diffdef class in STEPS.
 */

class Diffdef {
  public:
    Diffdef(const Compdef& compdef,
            container::kproc_id kproc,
            container::diffusion_id t_diffusion,
            container::species_id species,
            osh::Real t_dcst);

    inline container::kproc_id getKProcContainerIdx() const noexcept {
        return kprocContainerIdx;
    }

    inline container::diffusion_id getDiffContainerIdx() const noexcept {
        return diffusion;
    }
    inline container::species_id getSpecContainerIdx() const noexcept {
        return specContainerIdx;
    }
    inline osh::Real getDcst() const noexcept {
        return dcst;
    }

    inline const Compdef& compdef() const noexcept {
        return pCompdef;
    }

    /**
     * \param species chemical species identifier
     * \return true if the given species depends on this diffusion
     */
    inline bool depSpec(container::species_id species) const noexcept {
        return species == specContainerIdx;
    }

    void report(std::ostream& ostr) const;

  private:
    const Compdef& pCompdef;
    container::kproc_id kprocContainerIdx;
    container::diffusion_id diffusion;
    container::species_id specContainerIdx;
    osh::Real dcst;
};

}  // namespace steps::dist
