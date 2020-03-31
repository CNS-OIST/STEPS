#pragma once

#include <petscsys.h>

#include "fwd.hpp"
#include "vocabulary.hpp"

namespace zee {

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
            container::specie_id specie,
            PetscScalar t_dcst);

    inline container::kproc_id getKProcContainerIdx() const noexcept {
        return kprocContainerIdx;
    }

    inline container::diffusion_id getDiffContainerIdx() const noexcept {
        return diffusion;
    }
    inline container::specie_id getSpecContainerIdx() const noexcept {
        return specContainerIdx;
    }
    inline PetscScalar getDcst() const noexcept {
        return dcst;
    }

    inline const Compdef& compdef() const noexcept {
        return pCompdef;
    }

    /**
     * \param specie chemical specie identifier
     * \return true if the given specie depends on this diffusion
     */
    inline bool depSpec(container::specie_id specie) const noexcept {
        return specie == specContainerIdx;
    }

    void report(std::ostream& ostr) const;

  private:
    const Compdef& pCompdef;
    container::kproc_id kprocContainerIdx;
    container::diffusion_id diffusion;
    container::specie_id specContainerIdx;
    PetscScalar dcst;
};

}  // namespace zee
