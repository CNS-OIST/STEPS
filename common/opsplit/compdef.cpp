#include <cassert>
#include <set>

#include "../../OpSplitDMPlex/common.hpp"
#include "opsplit/compdef.hpp"
#include "opsplit/diffdef.hpp"
#include "opsplit/reacdef.hpp"

namespace zee {

Compdef::Compdef(const Statedef& statedef,
                 model::compartment_id t_model_compartment,
                 container::compartment_id t_container_compartment)
    : pStatedef(statedef)
    , model_compartment(std::move(t_model_compartment))
    , container_compartment(t_container_compartment) {}

container::specie_id Compdef::addSpec(model::specie_id specie) {
    auto specieIt = specM2C.find(specie);
    if (specieIt != specM2C.end()) {
        return specieIt->second;
    }
    const container::specie_id spec_container_idx(
        static_cast<container::specie_id::value_type>(specC2M.size()));
    specM2C[specie] = spec_container_idx;
    specC2M.push_back(specie);
    return spec_container_idx;
}

container::specie_id Compdef::getSpecContainerIdx(model::specie_id specie) const {
    auto result = specM2C.find(specie);
    if (result != specM2C.end()) {
        return result->second;
    }
    return container::specie_id(UNKNOWN_IDX);
}

model::specie_id Compdef::getSpecModelIdx(container::specie_id specie) const {
    assert(specie < static_cast<container::specie_id::value_type>(specC2M.size()));
    return specC2M[static_cast<size_t>(specie.get())];
}

container::reaction_id Compdef::addReac(const std::vector<container::specie_id>& reactants,
                                        const std::vector<container::specie_id>& products,
                                        PetscScalar kcst) {
    container::kproc_id kproc_id(nKProcs);
    container::reaction_id reac_container_idx(static_cast<PetscInt>(reacdefPtrs.size()));
    reacdefPtrs.emplace_back(
        std::make_unique<Reacdef>(*this, kproc_id, reac_container_idx, reactants, products, kcst));
    nKProcs++;
    return reac_container_idx;
}

Reacdef& Compdef::getReac(container::reaction_id reaction) const {
    assert(reaction < static_cast<PetscInt>(reacdefPtrs.size()));
    return *reacdefPtrs[static_cast<size_t>(reaction.get())];
}


container::diffusion_id Compdef::addDiff(container::specie_id specie, PetscScalar dcst) {
    assert(specie < static_cast<container::specie_id::value_type>(specC2M.size()));
    const container::kproc_id kproc_id(nKProcs);
    const container::diffusion_id diffusion_id(static_cast<PetscInt>(diffdefPtrs.size()));
    diffdefPtrs.emplace_back(
        std::make_unique<Diffdef>(*this, kproc_id, diffusion_id, specie, dcst));
    nKProcs++;
    species_diffused_.insert(specie);
    return diffusion_id;
}

container::specie_id Compdef::getDiffSpecContainerIdx(container::diffusion_id diffusion) {
    return diffdefPtrs[static_cast<size_t>(diffusion.get())]->getSpecContainerIdx();
}

model::specie_id Compdef::getDiffSpecModelIdx(container::diffusion_id diffusion) {
    const auto spec_id = diffdefPtrs[static_cast<size_t>(diffusion.get())]->getSpecContainerIdx();
    return specC2M[static_cast<size_t>(spec_id.get())];
}

Diffdef& Compdef::getDiff(container::diffusion_id diffusion) {
    return *diffdefPtrs[static_cast<size_t>(diffusion.get())];
}

Diffdef& Compdef::getDiffByKProcContainerIdx(container::kproc_id kproc) {
    return *diffdefPtrs[static_cast<size_t>(kproc.get() - getNReacs())];
}

bool Compdef::KProcDepSpec(container::kproc_id kproc,
                           container::specie_id specie,
                           KProcType filter_mask) const {
    const auto type = getKProcType(kproc);
    if (!any(type & filter_mask)) {
        return false;
    }
    switch (type) {
    case KProcType::Reac: {
        return reacdefPtrs[static_cast<size_t>(kproc.get())]->depSpec(specie);
    }
    case KProcType::Diff: {
        return diffdefPtrs[static_cast<size_t>(kproc.get() - getNReacs())]->depSpec(specie);
    }
    case KProcType::SDiff:
    case KProcType::VDepSReac:
    case KProcType::SReac:
        break;
    }
    return false;
}

void Compdef::report(std::ostream& ostr) const {
    ostr << "Compartment ID: " << model_compartment << " Model Idx: " << container_compartment
         << std::endl;
    ostr << "Number of Species: " << specM2C.size() << std::endl;
    ostr << "SpecM2C: [Model Idx, Container Idx]\n";
    for (const auto& spec: specM2C) {
        ostr << "[" << spec.first << ", " << spec.second << "]" << std::endl;
    }
    ostr << "SpecC2M: [Container Idx, Model Idx]\n";
    auto spec_container_idx = 0;
    for (const auto& spec: specC2M) {
        ostr << "[" << spec_container_idx << ", " << spec << "]" << std::endl;
        spec_container_idx++;
    }

    ostr << std::endl;
    ostr << "Number of Kinetic Processes: " << nKProcs << " (Reactions: " << reacdefPtrs.size()
         << ","
         << " Diffusions: " << diffdefPtrs.size() << ")\n";

    ostr << std::endl;
    for (const auto& reacdef: reacdefPtrs) {
        reacdef->report(ostr);
    }

    ostr << std::endl;
    for (const auto& diffdef: diffdefPtrs) {
        diffdef->report(ostr);
    }
    ostr << std::endl;
}

}  // namespace zee
