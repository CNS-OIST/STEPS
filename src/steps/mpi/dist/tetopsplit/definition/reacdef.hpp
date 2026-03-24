#pragma once

#include <Omega_h_defines.hpp>
#include <vector>

#include "compdef.hpp"
#include "complexeventsdef.hpp"
#include "fwd.hpp"
#include "model/complexreac.hpp"
#include "model/fwd.hpp"
#include "statedef.hpp"

#include "util/strong_ra.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

void report_molecule(std::stringstream& s,
                     const model::species_name& name,
                     const osh::I64 stochiometry,
                     const mesh::tetrahedron_id_t tet_id);

/**
 * \brief State definition of a reaction.
 *
 * The ReacdefBase class defines the sub biochemical container of
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

template <typename ReacT>
class ReacdefBase {
  public:
    enum class PoolChangeArrayType { LHS, RHS, UPD };
    using pool_change_t = util::strongid_vector<container::species_id, osh::I64>;
    using ContReacID = typename ModID2ContID<typename Model2Def<ReacT>::model_id_type>::type;

    ReacdefBase(const Compdef& compdef,
                container::kproc_id kproc,
                ContReacID reaction,
                const ReacT& reac)
        : reac_(reac)
        , pCompdef(compdef)
        , kproc_id(kproc)
        , reaction_id(reaction)
        , kcst(reac.getKcst())
        , order(static_cast<osh::I64>(reac.getOrder()))
        , poolChangeLHS(compdef.getNSpecs(), 0)
        , poolChangeRHS(compdef.getNSpecs(), 0)
        , poolChangeUPD(compdef.getNSpecs(), 0) {
        for (const auto* spec: reac.getLHS()) {
            auto specId = compdef.getSpecContainerIdx(
                compdef.statedef().getSpecModelIdx(model::species_name(spec->getID())));
            poolChangeLHS[specId] -= 1;
            poolChangeUPD[specId] -= 1;
        }
        for (const auto* spec: reac.getRHS()) {
            auto specId = compdef.getSpecContainerIdx(
                compdef.statedef().getSpecModelIdx(model::species_name(spec->getID())));
            poolChangeRHS[specId] += 1;
            poolChangeUPD[specId] += 1;
        }

        for (auto species: poolChangeUPD.range()) {
            if (poolChangeUPD[species] != 0) {
                updSpecModelIdxs.push_back(compdef.getSpecModelIdx(species));
            }
        }
    }

    /// Reset the def-object values to model defaults
    void reset() {
        kcst = reac_.getKcst();
        extent = 0;
    }

    inline container::kproc_id getKProcContainerIdx() const noexcept {
        return kproc_id;
    }

    inline ContReacID reactionID() const noexcept {
        return reaction_id;
    }

    inline osh::Real getKcst() const noexcept {
        return kcst;
    }

    inline void setKcst(osh::Real k) noexcept {
        kcst = k;
    }

    inline osh::I64 getExtent() const noexcept {
        return extent;
    }

    inline void add_extent(int nb) noexcept {
        extent += nb;
    }

    inline osh::I64 getOrder() const noexcept {
        return order;
    }

    inline const Compdef& compdef() const noexcept {
        return pCompdef;
    }

    inline bool depSpec(container::species_id species) const noexcept {
        return poolChangeLHS[species] != 0;
    }

    /**
     *
     * \return dependent species of the reaction
     */
    inline const std::vector<model::species_id>& getUpdSpecModelIdxs() const noexcept {
        return updSpecModelIdxs;
    }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
    const pool_change_t& getPoolChangeArray(PoolChangeArrayType type) const noexcept {
        switch (type) {
        case PoolChangeArrayType::LHS:
            return poolChangeLHS;
        case PoolChangeArrayType::RHS:
            return poolChangeRHS;
        case PoolChangeArrayType::UPD:
            return poolChangeUPD;
        }
    }
#pragma GCC diagnostic pop

    const pool_change_t& getPoolChangeLHS() const noexcept {
        return poolChangeLHS;
    }

    const pool_change_t& getPoolChangeRHS() const noexcept {
        return poolChangeRHS;
    }

    const pool_change_t& getPoolChangeUPD() const noexcept {
        return poolChangeUPD;
    }

    void report(std::ostream& ostr, const mesh::tetrahedron_id_t tet_id) const {
        ostr << "Type: Reaction, ID: " << reaction_id << '\n';

        std::stringstream o_lhs, o_rhs;
        for (auto s: poolChangeRHS.range()) {
            const auto spec_model_idx = pCompdef.getSpecModelIdx(s);
            const auto name = pCompdef.statedef().getSpecID(spec_model_idx);
            report_molecule(o_lhs, name, -poolChangeLHS[s], tet_id);
            report_molecule(o_rhs, name, poolChangeRHS[s], tet_id);
        }

        ostr << o_lhs.str() << " -> " << o_rhs.str();
        ostr << " (kcst: " << kcst << ")\n";
    }

  private:
    const ReacT& reac_;
    const Compdef& pCompdef;
    const container::kproc_id kproc_id;
    const ContReacID reaction_id;
    osh::Real kcst;
    const osh::I64 order;
    osh::I64 extent{0};

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

/**
 * \brief State definition of a complex reaction.
 *
 * In addition to species, that are handled by ReacdefBase, complex reactions involve
 * complex events that consist in the update (UPD) or deletion (DEL) of an existing
 * complex state, or in the creation (CRE) of a new complex state.
 *
 */
class ComplexReacdef: public ReacdefBase<steps::model::ComplexReac> {
  public:
    ComplexReacdef(const Compdef& compdef,
                   container::kproc_id kproc,
                   container::complex_reaction_id reaction,
                   const steps::model::ComplexReac& reac);

    const std::map<model::complex_id, std::set<model::complex_substate_id>>& complexDEPMAP()
        const noexcept {
        return pComplex_DEPMAP;
    }
    const std::map<model::complex_id, std::set<model::complex_substate_id>>& complexUPDMAP()
        const noexcept {
        return pComplex_UPDMAP;
    }
    const std::vector<std::shared_ptr<ComplexUpdateEventdef>>& updEvents() const noexcept {
        return pComplexUPDEvs;
    }
    const std::vector<std::shared_ptr<ComplexDeleteEventdef>>& delEvents() const noexcept {
        return pComplexDELEvs;
    }
    const std::vector<std::shared_ptr<ComplexCreateEventdef>>& creEvents() const noexcept {
        return pComplexCREEvs;
    }
    std::vector<std::shared_ptr<ComplexLHSEventdef>> lhsEvents() const {
        std::vector<std::shared_ptr<ComplexLHSEventdef>> ret;
        ret.reserve(pComplexUPDEvs.size() + pComplexDELEvs.size());
        ret.insert(ret.end(), pComplexUPDEvs.begin(), pComplexUPDEvs.end());
        ret.insert(ret.end(), pComplexDELEvs.begin(), pComplexDELEvs.end());
        return ret;
    }

  private:
    std::vector<std::shared_ptr<ComplexUpdateEventdef>> pComplexUPDEvs;
    std::vector<std::shared_ptr<ComplexDeleteEventdef>> pComplexDELEvs;
    std::vector<std::shared_ptr<ComplexCreateEventdef>> pComplexCREEvs;

    // cmplxIdx -> {sub unit states ind}
    std::map<model::complex_id, std::set<model::complex_substate_id>> pComplex_DEPMAP;
    std::map<model::complex_id, std::set<model::complex_substate_id>> pComplex_UPDMAP;
};

}  // namespace steps::dist
