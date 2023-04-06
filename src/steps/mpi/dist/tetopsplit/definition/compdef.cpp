#include "compdef.hpp"

#include <cassert>

#include "diffdef.hpp"
#include "reacdef.hpp"

namespace steps {
namespace dist {

Compdef::Compdef(const Statedef& statedef,
                 model::compartment_id t_model_compartment,
                 container::compartment_id t_container_compartment)
    : pStatedef(statedef)
    , model_compartment(std::move(t_model_compartment))
    , container_compartment(t_container_compartment) {}

container::species_id Compdef::addSpec(model::species_id species) {
  auto speciesIt = specM2C.find(species);
  if (speciesIt != specM2C.end()) {
    return speciesIt->second;
  }
  const container::species_id spec_container_idx(
      static_cast<container::species_id::value_type>(specC2M.size()));
  specM2C[species] = spec_container_idx;
  specC2M.push_back(species);
  return spec_container_idx;
}

container::species_id
Compdef::getSpecContainerIdx(model::species_id species) const {
  auto result = specM2C.find(species);
  if (result != specM2C.end()) {
    return result->second;
  }
  return std::nullopt;
}

model::species_id
Compdef::getSpecModelIdx(container::species_id species) const {
  assert(species <
         static_cast<container::species_id::value_type>(specC2M.size()));
  return specC2M[static_cast<size_t>(species.get())];
}

container::reaction_id
Compdef::addReac(const std::vector<container::species_id> &reactants,
                 const std::vector<container::species_id> &products,
                 osh::Real kcst) {
  container::kproc_id kproc_id(nKProcs);
  container::reaction_id reac_container_idx(
      static_cast<osh::I64>(reacdefPtrs.size()));
  reacdefPtrs.emplace_back(std::make_unique<Reacdef>(
      *this, kproc_id, reac_container_idx, reactants, products, kcst));
  nKProcs++;
  return reac_container_idx;
}

Reacdef& Compdef::getReac(container::reaction_id reaction) const {
  assert(reaction < static_cast<osh::I64>(reacdefPtrs.size()));
  return *reacdefPtrs[static_cast<size_t>(reaction.get())];
}

container::diffusion_id Compdef::addDiff(container::species_id species,
                                         osh::Real dcst) {
  assert(species <
         static_cast<container::species_id::value_type>(specC2M.size()));
  const container::kproc_id kproc_id(nKProcs);
  const container::diffusion_id diffusion_id(
      static_cast<osh::I64>(diffdefPtrs.size()));
  diffdefPtrs.emplace_back(
      std::make_unique<Diffdef>(*this, kproc_id, diffusion_id, species, dcst));
  nKProcs++;
  species_diffused_.insert(species);
  return diffusion_id;
}

container::species_id
Compdef::getDiffSpecContainerIdx(container::diffusion_id diffusion) {
  return diffdefPtrs[static_cast<size_t>(diffusion.get())]
      ->getSpecContainerIdx();
}

model::species_id
Compdef::getDiffSpecModelIdx(container::diffusion_id diffusion) {
  const auto spec_id =
      diffdefPtrs[static_cast<size_t>(diffusion.get())]->getSpecContainerIdx();
  return specC2M[static_cast<size_t>(spec_id.get())];
}

Diffdef& Compdef::getDiff(container::diffusion_id diffusion) {
    return *diffdefPtrs[static_cast<size_t>(diffusion.get())];
}

Diffdef& Compdef::getDiffByKProcContainerIdx(container::kproc_id kproc) {
    return *diffdefPtrs[static_cast<size_t>(kproc.get() - getNReacs())];
}

bool Compdef::KProcDepSpec(container::kproc_id kproc,
                           container::species_id species) const {
  const auto type = getKProcType(kproc);
  switch (type) {
  case kproc::KProcType::Reac: {
    return reacdefPtrs[static_cast<size_t>(kproc.get())]->depSpec(species);
  }
  case kproc::KProcType::Diff: {
    return diffdefPtrs[static_cast<size_t>(kproc.get() - getNReacs())]->depSpec(
        species);
  }
  case kproc::KProcType::VDepSReac:
  case kproc::KProcType::GHKSReac:
  case kproc::KProcType::SReac:
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
        reacdef->report(ostr, std::nullopt);
    }

    ostr << std::endl;
    for (const auto& diffdef: diffdefPtrs) {
        diffdef->report(ostr);
    }
    ostr << std::endl;
}

}  // namespace dist
}  // namespace steps
