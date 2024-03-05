#include "reacdef.hpp"

#include <set>

#include "compdef.hpp"
#include "statedef.hpp"


namespace steps::dist {

Reacdef::Reacdef(const Compdef& compdef,
                 container::kproc_id kproc,
                 container::reaction_id reaction,
                 const std::vector<container::species_id>& reactants,
                 const std::vector<container::species_id>& products,
                 osh::Real t_kcst)
    : pCompdef(compdef)
    , kproc_id(kproc)
    , reaction_id(reaction)
    , kcst(t_kcst)
    , order(static_cast<osh::I64>(reactants.size())) {
    auto num_species = compdef.getNSpecs();

    poolChangeLHS.assign(static_cast<size_t>(num_species), 0);
    poolChangeRHS.assign(static_cast<size_t>(num_species), 0);
    poolChangeUPD.assign(static_cast<size_t>(num_species), 0);

    for (const auto& species: reactants) {
        poolChangeLHS[static_cast<size_t>(species.get())] -= 1;
        poolChangeUPD[static_cast<size_t>(species.get())] -= 1;
    }
    for (const auto& species: products) {
        poolChangeRHS[static_cast<size_t>(species.get())] += 1;
        poolChangeUPD[static_cast<size_t>(species.get())] += 1;
    }

    for (container::species_id species(0); species < num_species; species++) {
        if (poolChangeUPD[static_cast<size_t>(species.get())] != 0) {
            updSpecModelIdxs.push_back(compdef.getSpecModelIdx(species));
        }
    }
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
const Reacdef::pool_change_t& Reacdef::getPoolChangeArray(PoolChangeArrayType type) const noexcept {
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

void Reacdef::report(std::ostream& ostr, const mesh::tetrahedron_id_t tet_id) const {
    ostr << "Type: Reaction, ID: " << reaction_id << '\n';

    std::stringstream o_lhs, o_rhs;
    const auto n_local_specs = pCompdef.getNSpecs();
    for (osh::I32 s = 0; s < n_local_specs; ++s) {
        const auto spec_model_idx = pCompdef.getSpecModelIdx(container::species_id(s));
        const auto name = pCompdef.statedef().getSpecID(spec_model_idx);
        report_molecule(o_lhs, name, -poolChangeLHS[s], tet_id);
        report_molecule(o_rhs, name, poolChangeRHS[s], tet_id);
    }

    ostr << o_lhs.str() << " -> " << o_rhs.str();
    ostr << " (kcst: " << kcst << ")\n";
}

}  // namespace steps::dist
