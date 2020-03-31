#pragma once

#include <vector>

#include <petscsys.h>

#include "fwd.hpp"
#include "opsplit/vocabulary.hpp"

namespace zee {

/**
 * \brief State definition of a reaction.
 *
 * The Reacdef class defines the sub biochemical container of
 * a reaction in a compartment.
 *
 * A reaction is defined by
 * - lhs Left hand side species list in the reaction equation
 * - rhs Right hand side species list in the reaction equation
 * - kcst Reaction constant
 *
 * This class corresponds to the solver::Reacdef class in STEPS.
 */

class Reacdef {
  public:
    enum class PoolChangeArrayType { LHS, RHS, UPD };

    Reacdef(const Compdef& compdef,
            container::kproc_id kproc,
            container::reaction_id reaction,
            const std::vector<container::specie_id>& reactants,
            const std::vector<container::specie_id>& products,
            PetscScalar kcst);

    inline container::kproc_id getKProcContainerIdx() const noexcept {
        return kproc_id;
    }

    inline PetscScalar getKcst() const noexcept {
        return kcst;
    }

    inline PetscInt getOrder() const noexcept {
        return order;
    }

    inline const Compdef& compdef() const noexcept {
        return pCompdef;
    }

    inline bool depSpec(container::specie_id specie) const noexcept {
        return poolChangeLHS[static_cast<size_t>(specie.get())] != 0;
    }

    /**
     *
     * \return dependent species of the reaction
     */
    inline const std::vector<model::specie_id>& getUpdSpecModelIdxs() const noexcept {
        return updSpecModelIdxs;
    }

    const std::vector<PetscInt>& getPoolChangeArray(PoolChangeArrayType type) const;

    void report(std::ostream& ostr) const;

  private:
    const Compdef& pCompdef;
    const container::kproc_id kproc_id;
    const container::reaction_id reaction_id;
    const PetscScalar kcst;
    const PetscInt order;

    std::vector<PetscInt> poolChangeLHS;
    std::vector<PetscInt> poolChangeRHS;
    std::vector<PetscInt> poolChangeUPD;

    std::vector<model::specie_id> updSpecModelIdxs;
};

}  // namespace zee
