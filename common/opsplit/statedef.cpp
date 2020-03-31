#include <algorithm>
#include <sstream>

#include <hadoken/format/format.hpp>

#include "opsplit/compdef.hpp"
#include "opsplit/diffdef.hpp"
#include "opsplit/patchdef.hpp"
#include "opsplit/reacdef.hpp"
#include "opsplit/sreacdef.hpp"
#include "opsplit/statedef.hpp"


#include "../../OpSplitDMPlex/common.hpp"

namespace zee {

using hadoken::scat;

model::specie_id Statedef::addSpec(const model::specie_name& name) {
    auto result = specModelIdxs.find(name);
    if (result != specModelIdxs.end()) {
        return result->second;
    }
    auto model_idx = static_cast<PetscInt>(specModelIdxs.size());
    specModelIdxs.emplace(name, model_idx);
    specIDs.push_back(name);
    return specModelIdxs[name];
}

model::specie_id Statedef::getSpecModelIdx(const model::specie_name& name) const {
    auto result = specModelIdxs.find(name);
    if (result == specModelIdxs.end()) {
        return model::specie_id(UNKNOWN_IDX);
    }
    return result->second;
}

container::compartment_id Statedef::addComp(const model::compartment_id& compartment) {
    container::compartment_id model_id(
        static_cast<container::compartment_id::value_type>(compdefPtrs.size()));
    compdefPtrs.emplace_back(std::make_unique<Compdef>(*this, compartment, model_id));
    compModelIdxs[compartment] = model_id;
    return model_id;
}

container::patch_id Statedef::addPatch(
    const model::patch_id& patchId,
    const model::compartment_id& inner_compartment_id,
    const boost::optional<model::compartment_id>& outer_compartment_id) {
    container::patch_id container_id(
        static_cast<container::patch_id::value_type>(patchdefPtrs.size()));
    // we need to specify the inner compartment
    patchdefPtrs.emplace_back(
        std::make_unique<Patchdef>(*this, patchId, inner_compartment_id, outer_compartment_id));
    patchModelIdxs[patchId] = container_id;
    return container_id;
}

container::surface_reaction_id Statedef::addSurfReac(
    const model::patch_id& patchId,
    const std::vector<model::specie_name>& reactants_i,
    const std::vector<model::specie_name>& reactants_s,
    const std::vector<model::specie_name>& reactants_o,
    const std::vector<model::specie_name>& products_i,
    const std::vector<model::specie_name>& products_s,
    const std::vector<model::specie_name>& products_o,
    const PetscScalar kcst) {
    if (!patchModelIdxs.count(patchId)) {
        throw std::invalid_argument("Wrong patch id: " + patchId);
    }
    Patchdef& patch = *patchdefPtrs[static_cast<size_t>(patchModelIdxs[patchId].get())];

    const Compdef& inner_compartment =
        *compdefPtrs[static_cast<size_t>(compModelIdxs[patch.getInnerCompId()].get())];
    boost::optional<const Compdef&> outer_compartment =
        patch.getOuterCompId()
            ? *compdefPtrs[static_cast<size_t>(compModelIdxs[*patch.getOuterCompId()].get())]
            : boost::optional<const Compdef&>();

    auto& specMdlIdxs = this->specModelIdxs;
    auto convert_id_type_patch = [&patch, &specMdlIdxs](const std::vector<model::specie_name>& v) {
        std::vector<container::specie_id> res;
        res.reserve(v.size());
        std::transform(v.begin(),
                       v.end(),
                       std::back_inserter(res),
                       [&patch, &specMdlIdxs](auto& l) {
                           return patch.getSpecPatchIdx(specMdlIdxs[l]);
                       });
        return res;
    };
    auto convert_id_type_comp = [&specMdlIdxs](const Compdef& comp,
                                               const std::vector<model::specie_name>& v) {
        std::vector<container::specie_id> res;
        res.reserve(v.size());
        std::transform(v.begin(), v.end(), std::back_inserter(res), [&specMdlIdxs, &comp](auto& l) {
            auto id = comp.getSpecContainerIdx(specMdlIdxs[l]);
            if (id == container::specie_id(UNKNOWN_IDX)) {
                throw std::logic_error(
                    hadoken::scat("Compartment  ", comp.getID(), ": unknown specie ", l));
            }
            return comp.getSpecContainerIdx(specMdlIdxs[l]);
        });
        return res;
    };
    std::vector<container::specie_id> empty;
    auto id = patch.addSurfaceReac(
        convert_id_type_comp(inner_compartment, reactants_i),
        convert_id_type_patch(reactants_s),
        reactants_o.size() ? convert_id_type_comp(*outer_compartment, reactants_o) : empty,
        convert_id_type_comp(inner_compartment, products_i),
        convert_id_type_patch(products_s),
        products_o.size() ? convert_id_type_comp(*outer_compartment, products_o) : empty,
        kcst);
    return id;
}

container::compartment_id Statedef::getCompModelIdx(const model::compartment_id& compartment) const
    noexcept {
    auto id = compModelIdxs.find(compartment);
    if (id == compModelIdxs.end()) {
        return container::compartment_id(UNKNOWN_IDX);
    }
    return id->second;
}

Compdef& Statedef::getCompdef(container::compartment_id compartment) const {
    return *compdefPtrs[static_cast<size_t>(compartment.get())];
}

Compdef& Statedef::getCompdef(const model::compartment_id& compartment) const {
    return *compdefPtrs[static_cast<size_t>(getCompModelIdx(compartment).get())];
}

Patchdef& Statedef::getPatchdef(const container::patch_id& patchId) const {
    return *patchdefPtrs[static_cast<size_t>(patchId.get())];
}

Patchdef& Statedef::getPatchdef(const model::patch_id& patchId) const {
    return *patchdefPtrs[static_cast<size_t>(patchModelIdxs.at(patchId).get())];
}


model::specie_id Statedef::addCompSpec(const model::compartment_id& compartment,
                                       const model::specie_name& specie) {
    auto spec_id = addSpec(specie);
    auto& compdef = compdefPtrs[static_cast<size_t>(compModelIdxs[compartment].get())];
    compdef->addSpec(spec_id);
    return spec_id;
}

std::vector<model::specie_id> Statedef::addCompSpecs(
    const model::compartment_id& compartment,
    const std::vector<model::specie_name>& species) {
    std::vector<model::specie_id> ids;
    ids.reserve(species.size());
    for (const auto& spec_id: species) {
        ids.push_back(addCompSpec(compartment, spec_id));
    }
    return ids;
}

std::vector<model::specie_id> Statedef::addPatchSpecs(
    const model::patch_id& patchId,
    const std::vector<model::specie_name>& species) {
    std::vector<model::specie_id> ids;
    ids.reserve(species.size());
    for (const auto& spec_id: species) {
        ids.push_back(addPatchSpec(patchId, spec_id));
    }
    return ids;
}

model::specie_id Statedef::addPatchSpec(const model::patch_id& patch_id,
                                        const model::specie_name& specie) {
    auto spec_id = addSpec(specie);
    auto& patchdef = patchdefPtrs[static_cast<size_t>(patchModelIdxs[patch_id].get())];
    patchdef->addSpec(spec_id);
    return spec_id;
}


container::specie_id Statedef::getCompSpecContainerIdx(const model::compartment_id& compartment,
                                                       const model::specie_name& specie) const {
    auto comp_model_idx = getCompModelIdx(compartment);
    if (comp_model_idx == static_cast<model::specie_id::value_type>(UNKNOWN_IDX)) {
        throw std::invalid_argument(scat("Unknown compartment: ", compartment));
    }
    const auto spec_model_idx = getSpecModelIdx(specie);
    if (spec_model_idx == static_cast<model::specie_id::value_type>(UNKNOWN_IDX)) {
        throw std::invalid_argument(scat("Unknown specie: ", specie));
    }

    return compdefPtrs[static_cast<size_t>(comp_model_idx.get())]->getSpecContainerIdx(
        spec_model_idx);
}

container::diffusion_id Statedef::addCompDiff(const model::compartment_id& compartment,
                                              const model::specie_name& specie_name,
                                              PetscScalar dcst) {
    const auto scidx = getCompSpecContainerIdx(compartment, specie_name);
    return compdefPtrs[static_cast<size_t>(compModelIdxs[compartment].get())]->addDiff(scidx, dcst);
}

container::reaction_id Statedef::addCompReac(const model::compartment_id& compartment,
                                             const std::vector<model::specie_name>& reactants,
                                             const std::vector<model::specie_name>& products,
                                             PetscScalar kcst) {
    std::vector<container::specie_id> reactants_ids;
    std::vector<container::specie_id> products_ids;
    reactants_ids.reserve(reactants.size());
    products_ids.reserve(products.size());

    std::transform(reactants.begin(),
                   reactants.end(),
                   std::back_inserter(reactants_ids),
                   [this, &compartment](const auto& specie) {
                       return this->getCompSpecContainerIdx(compartment, specie);
                   });

    std::transform(products.begin(),
                   products.end(),
                   std::back_inserter(products_ids),
                   [this, &compartment](const auto& specie) {
                       return this->getCompSpecContainerIdx(compartment, specie);
                   });

    auto& compartmentPtr = compdefPtrs[static_cast<size_t>(compModelIdxs[compartment].get())];
    return compartmentPtr->addReac(reactants_ids, products_ids, kcst);
}

std::string Statedef::createReport() const {
    std::stringstream report_stream;
    report_stream << "Biochemical Model Report\n";
    report_stream << "Number of Species: " << specModelIdxs.size() << std::endl;
    report_stream << "[ID, ModelIdx]\n";
    for (auto&& spec: specModelIdxs) {
        report_stream << "[" << spec.first << ", " << spec.second << "]" << '\n';
    }
    report_stream << '\n';

    report_stream << "Number of Compartments: " << compModelIdxs.size() << '\n';
    report_stream << "[ID, Model Idx]\n";
    for (auto&& comp: compModelIdxs) {
        report_stream << "[" << comp.first << ", " << comp.second << "]" << '\n';
    }
    report_stream << '\n';

    report_stream << "Number of Patches: " << patchModelIdxs.size() << '\n';
    report_stream << "[ID, Model Idx]\n";
    for (auto&& patch: patchModelIdxs) {
        report_stream << "[" << patch.first << ", " << patch.second << "]" << '\n';
    }
    report_stream << '\n';

    report_stream << "Detail Compartment Report\n";
    for (auto&& compdef: compdefPtrs) {
        compdef->report(report_stream);
    }
    report_stream << '\n';

    return report_stream.str();
}

}  // namespace zee
