#include <set>

#include "opsplit/compdef.hpp"
#include "opsplit/reacdef.hpp"
#include "opsplit/statedef.hpp"

namespace zee {

Reacdef::Reacdef(const Compdef& compdef,
                 container::kproc_id kproc,
                 container::reaction_id reaction,
                 const std::vector<container::specie_id>& reactants,
                 const std::vector<container::specie_id>& products,
                 PetscScalar t_kcst)
    : pCompdef(compdef)
    , kproc_id(kproc)
    , reaction_id(reaction)
    , kcst(t_kcst)
    , order(static_cast<PetscInt>(reactants.size())) {
    auto num_species = compdef.getNSpecs();

    poolChangeLHS.assign(static_cast<size_t>(num_species), 0);
    poolChangeRHS.assign(static_cast<size_t>(num_species), 0);
    poolChangeUPD.assign(static_cast<size_t>(num_species), 0);

    for (const auto& specie: reactants) {
        poolChangeLHS[static_cast<size_t>(specie.get())] -= 1;
        poolChangeUPD[static_cast<size_t>(specie.get())] -= 1;
    }
    for (const auto& specie: products) {
        poolChangeRHS[static_cast<size_t>(specie.get())] += 1;
        poolChangeUPD[static_cast<size_t>(specie.get())] += 1;
    }

    for (container::specie_id specie(0); specie < num_species; specie++) {
        if (poolChangeUPD[static_cast<size_t>(specie.get())] != 0) {
            updSpecModelIdxs.push_back(compdef.getSpecModelIdx(specie));
        }
    }
}

const std::vector<PetscInt>& Reacdef::getPoolChangeArray(PoolChangeArrayType type) const {
    switch (type) {
    case PoolChangeArrayType::LHS:
        return poolChangeLHS;
    case PoolChangeArrayType::RHS:
        return poolChangeRHS;
    case PoolChangeArrayType::UPD:
        return poolChangeUPD;
    }
    throw std::logic_error("Unknown pool change type");
}

void Reacdef::report(std::ostream& ostr) const {
    ostr << "Reaction Report" << '\n';
    ostr << "KProc Container Idx: " << kproc_id << " Reac Container Idx: " << reaction_id
         << " Order: " << order << " KCST: " << kcst << '\n';

    const auto n_local_specs = pCompdef.getNSpecs();
    ostr << "Pool Change\n";
    ostr << "Spec:\t[ ";
    for (container::specie_id s(0); s < n_local_specs; s++) {
        const auto spec_model_idx = pCompdef.getSpecModelIdx(s);
        ostr << pCompdef.statedef().getSpecID(spec_model_idx) << '\t';
    }
    ostr << "]\n";

    ostr << "LHS\t[ ";
    for (auto&& value: poolChangeLHS) {
        ostr << value << '\t';
    }
    ostr << "]\n";

    ostr << "RHS:\t[ ";
    for (auto&& value: poolChangeRHS) {
        ostr << value << '\t';
    }
    ostr << "]\n";

    ostr << "UPD:\t[ ";
    for (auto&& value: poolChangeUPD) {
        ostr << value << '\t';
    }
    ostr << "]\n";
    ostr << '\n';
}

}  // namespace zee
