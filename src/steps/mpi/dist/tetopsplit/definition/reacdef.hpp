#pragma once

#include <vector>

#include "fwd.hpp"

#include "util/vocabulary.hpp"

namespace steps::dist {

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
 * - order sum of the reactant types (ex: 2 Ca + Na = XXX -> Order: 2)
 * - upd: update. Basically rhs + lhs
 * - dep: dependency. Basically element-wise conversion to bool of lhs
 *
 * This class corresponds to the solver::Reacdef class in STEPS.
 */

class Reacdef {
  public:
    enum class PoolChangeArrayType { LHS, RHS, UPD };
    using pool_change_t = std::vector<osh::I64>;

    Reacdef(const Compdef& compdef,
            container::kproc_id kproc,
            container::reaction_id reaction,
            const std::vector<container::species_id>& reactants,
            const std::vector<container::species_id>& products,
            osh::Real kcst);

    inline container::kproc_id getKProcContainerIdx() const noexcept {
        return kproc_id;
    }

    inline osh::Real getKcst() const noexcept {
        return kcst;
    }

    inline osh::I64 getOrder() const noexcept {
        return order;
    }

    inline const Compdef& compdef() const noexcept {
        return pCompdef;
    }

    inline bool depSpec(container::species_id species) const noexcept {
        return poolChangeLHS[static_cast<size_t>(species.get())] != 0;
    }

    /**
     *
     * \return dependent species of the reaction
     */
    inline const std::vector<model::species_id>& getUpdSpecModelIdxs() const noexcept {
        return updSpecModelIdxs;
    }

    const pool_change_t& getPoolChangeArray(PoolChangeArrayType type) const noexcept;

    const pool_change_t& getPoolChangeLHS() const noexcept {
        return poolChangeLHS;
    }

    const pool_change_t& getPoolChangeRHS() const noexcept {
        return poolChangeRHS;
    }

    const pool_change_t& getPoolChangeUPD() const noexcept {
        return poolChangeUPD;
    }

    void report(std::ostream& ostr, const mesh::tetrahedron_id_t tet_id) const;

  private:
    const Compdef& pCompdef;
    const container::kproc_id kproc_id;
    const container::reaction_id reaction_id;
    const osh::Real kcst;
    const osh::I64 order;

    /** poolChangeLHS[i] lhs of the reaction for the molecule i
     *
     * If involved in the reaction it is negative because reactants disappear in the reaction. The
     * value is the stochiometry ex: 2 Ca has -2
     */
    pool_change_t poolChangeLHS;
    /** poolChangeRHS[i] lhs of the reaction for the molecule i
     *
     * If involved in the reaction it is positive because products appear in the reaction. The value
     * is the stochiometry ex: 2 Ca has 2
     */
    pool_change_t poolChangeRHS;
    /** update
     *
     * Is is basically poolChangeLHS + poolChangeRHS
     */
    pool_change_t poolChangeUPD;

    std::vector<model::species_id> updSpecModelIdxs;
};

}  // namespace steps::dist
