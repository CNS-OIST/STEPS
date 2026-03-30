#include "reactions.hpp"

#include <Omega_h_map.hpp>
#include <Omega_h_mark.hpp>

#include "geom/dist/distmesh.hpp"
#include "kproc_state.hpp"
#include "math/constants.hpp"
#include "model/complexreac.hpp"
#include "mpi/dist/tetopsplit/definition/fwd.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist::kproc {

void ComplexReactions::addComplexReactions(const Compdef& compartment,
                                           MolState& mol_state,
                                           const mesh::tetrahedron_id_t& tet) {
    for (const auto& reacdef: compartment.reacdefs<ComplexReacdef>()) {
        std::vector<MolStateComplexElementID> comp_deps;
        const auto& deps = reacdef->complexDEPMAP();
        for (const auto& [cId, substates]: deps) {
            for (const auto& sus: substates) {
                comp_deps.emplace_back(tet, cId, sus);
            }
        }
        complex_reactions_deps.push_back(comp_deps);

        std::vector<MolStateComplexElementID> comp_upds;
        const auto& upds = reacdef->complexUPDMAP();
        for (const auto& [cId, substates]: upds) {
            for (const auto& sus: substates) {
                comp_upds.emplace_back(tet, cId, sus);
            }
        }
        complex_reactions_upds.push_back(comp_upds);

        candidates.emplace_back();
        auto& lhscands = candidates.back();
        std::map<model::complex_id, uint> added;
        for (auto ev: reacdef->lhsEvents()) {
            auto it = added.find(ev->complexIdx());
            if (it == added.end()) {
                it = added.insert({ev->complexIdx(), lhscands.size()}).first;
                lhscands.emplace_back(ev->complexIdx(), tet);
            }
            lhscands[it->second].addEvent(ev, mol_state.moleculesOnElements());
        }
    }
}

osh::Real ComplexReactions::computeRate(const MolState& mol_state, size_t index) const {
    auto specRate = ReactionsBase<ComplexReacdef>::computeRate(mol_state, index);
    if (specRate == 0.0) {
        return 0;
    }
    // Get the rates for the complex reaction part
    double cmult = 1.0;
    for (auto& cand: candidates[index]) {
        cmult *= cand.rateMult(mol_state.moleculesOnElements());
    }
    return cmult * specRate;
}

const std::vector<MolStateElementID>& ComplexReactions::updateMolStateAndOccupancy(
    MolState& mol_state,
    rng::RNG& rng,
    size_t index,
    osh::Real event_time) const {
    // Species part of the complex reaction
    auto& upd = ReactionsBase<ComplexReacdef>::updateMolStateAndOccupancy(
        mol_state, rng, index, event_time, 1);

    const auto& entity = ownerPoints_[index];
    auto& entmols = mol_state.moleculesOnElements();
    // Updates
    for (auto& cands: candidates[index]) {
        for (auto& event: cands.selectEvents(entmols, rng)) {
            if (event.first->type() == UPDEvent) {
                const auto ev = std::dynamic_pointer_cast<const ComplexUpdateEventdef>(event.first);
                const auto& state = entmols.complexStates(cands.complexIdx()).at(event.second);
                entmols.updateComplexUpdateOccupancy(cands.complexIdx(),
                                                     event.second,
                                                     ev->getUpdate(state, rng),
                                                     event_time);
            } else {
                // Deletions
                entmols.removeComplexUpdateOccupancy(cands.complexIdx(), event.second, event_time);
            }
        }
    }
    // Creations
    for (auto& ce: reacdefs_[index].get().creEvents()) {
        entmols.createComplexUpdateOccupancy(entity, ce->complexIdx(), ce->init(), event_time);
    }

    // TODO Eventually try to take complex changes into account for the return value, so RSSA can be
    // used
    return upd;
}

}  // namespace steps::dist::kproc
