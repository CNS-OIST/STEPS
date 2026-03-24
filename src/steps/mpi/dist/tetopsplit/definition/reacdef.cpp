#include "reacdef.hpp"

#include <set>

#include "model/complexreac.hpp"
#include "model/reac.hpp"
#include "model/spec.hpp"
#include "solver/fwd.hpp"
#include "statedef.hpp"
#include "util/vocabulary.hpp"


namespace steps::dist {

void report_molecule(std::stringstream& s,
                     const model::species_name& name,
                     const osh::I64 stochiometry,
                     const mesh::tetrahedron_id_t tet_id) {
    if (stochiometry == 0) {
        return;  // discard if not present in the formula
    }
    if (!s.str().empty()) {
        s << " + ";  // add + if there were other molecules before
    }
    if (stochiometry != 1) {
        s << stochiometry << " * ";  // add stochiometric number if relevant
    }
    s << name;
    if (tet_id.valid()) {
        s << "[Tet_" << tet_id << ']';  // add name and tet_id
    }
}

ComplexReacdef::ComplexReacdef(const Compdef& compdef,
                               container::kproc_id kproc,
                               container::complex_reaction_id reaction,
                               const steps::model::ComplexReac& reac)
    : ReacdefBase<steps::model::ComplexReac>::ReacdefBase(compdef, kproc, reaction, reac) {
    // Copy complex events
    const auto& sd = compdef.statedef();
    for (auto* ev: reac.getUPDEvents()) {
        pComplexUPDEvs.push_back(std::make_shared<ComplexUpdateEventdef>(*ev, sd));
    }
    for (auto* ev: reac.getDELEvents()) {
        pComplexDELEvs.push_back(std::make_shared<ComplexDeleteEventdef>(*ev, sd));
    }
    for (auto* ev: reac.getCREEvents()) {
        pComplexCREEvs.push_back(std::make_shared<ComplexCreateEventdef>(*ev, sd));
    }

    // set up deps for complexes
    for (const auto& upd: pComplexUPDEvs) {
        pComplex_DEPMAP[upd->complexIdx()].merge(upd->getDepSet());
        pComplex_UPDMAP[upd->complexIdx()].merge(upd->getUpdSet());
    }
    for (const auto& del: pComplexDELEvs) {
        pComplex_DEPMAP[del->complexIdx()].merge(del->getDepSet());
        pComplex_UPDMAP[del->complexIdx()].merge(del->getUpdSet());
    }
    for (const auto& cre: pComplexCREEvs) {
        pComplex_UPDMAP[cre->complexIdx()].merge(cre->getUpdSet());
    }
}

}  // namespace steps::dist
