#include "statedef.hpp"

#include <algorithm>
#include <sstream>

#include "compdef.hpp"
#include "diffdef.hpp"
#include "patchdef.hpp"
#include "reacdef.hpp"
#include "sreacdef.hpp"

#include "geom/dist/distcomp.hpp"
#include "geom/dist/distmemb.hpp"
#include "geom/dist/distmesh.hpp"
#include "geom/dist/distpatch.hpp"
#include "math/constants.hpp"
#include "model/chan.hpp"
#include "model/chanstate.hpp"
#include "model/diff.hpp"
#include "model/ghkcurr.hpp"
#include "model/model.hpp"
#include "model/ohmiccurr.hpp"
#include "model/reac.hpp"
#include "model/spec.hpp"
#include "model/sreac.hpp"
#include "model/vdepsreac.hpp"
#include "util/collections.hpp"

namespace steps::dist {

// Some of the STEPS model / geom objects use std::set<T*> internally to avoid
// duplication. This means the order in which elements are added could be
// different in different cores. To avoid this, we order them by name.
template <typename T>
bool stepsObjComp(T* a, T* b) {
    return a->getID() < b->getID();
}

/**
 * \brief Helper function that augment std::map::at with an error message
 * if \key cannot be found.
 *
 * TMap can be a std::map or a std::unordered_map
 *
 */
template <typename AssociativeContainer>
const typename AssociativeContainer::mapped_type& map_at(
    const AssociativeContainer& container,
    const typename AssociativeContainer::key_type& key,
    const std::string& name) {
    try {
        return container.at(key);
    } catch (const std::out_of_range&) {
        std::stringstream ss;
        ss << "No " << name << " with id " << key << ". Possible values are:\n[";
        for (const auto& i: container) {
            ss << i.first << ", ";
        }
        ss << "].\n";
        throw std::invalid_argument(ss.str());
    }
}

template <typename AssociativeContainer>
typename AssociativeContainer::mapped_type& map_at(
    AssociativeContainer& container,
    const typename AssociativeContainer::key_type& key,
    const std::string& name) {
    try {
        return container.at(key);
    } catch (const std::out_of_range&) {
        std::stringstream ss;
        ss << "No " << name << " with id " << key << ". Possible values are:\n[";
        for (const auto& i: container) {
            ss << i.first << ", ";
        }
        ss << "].\n";
        throw std::invalid_argument(ss.str());
    }
}


Statedef::Statedef(const steps::model::Model& model, const steps::dist::DistMesh& mesh) {
    // Species
    auto allSpecs = model.getAllSpecs();
    std::sort(allSpecs.begin(), allSpecs.end(), stepsObjComp<steps::model::Spec>);
    for (auto* spec: allSpecs) {
        addSpec({spec->getID().c_str()});
    }

    // Compartments
    auto allComps = mesh.getAllComps();
    std::sort(allComps.begin(), allComps.end(), stepsObjComp<steps::dist::DistComp>);
    for (auto* comp: allComps) {
        auto comp_id = comp->getID();
        addComp({comp_id.c_str()});

        addCompartmentConductivity(comp_id.c_str(), comp->getConductivity());

        auto allCompSpecs = comp->getAllSpecs(&model);
        std::sort(allCompSpecs.begin(), allCompSpecs.end(), stepsObjComp<steps::model::Spec>);
        for (auto& spec: allCompSpecs) {
            addCompSpec(comp_id.c_str(), spec->getID().c_str());
        }

        // Reactions
        auto allReacs = comp->getAllReacs(&model);
        std::sort(allReacs.begin(), allReacs.end(), stepsObjComp<steps::model::Reac>);
        for (auto* reac: allReacs) {
            const auto lhs = reac->getLHS();
            std::vector<model::species_name> lhs_str;
            lhs_str.reserve(lhs.size());
            for (const auto* spec: lhs) {
                lhs_str.emplace_back(spec->getID());
            }
            const auto rhs = reac->getRHS();
            std::vector<model::species_name> rhs_str;
            rhs_str.reserve(rhs.size());
            for (const auto* spec: rhs) {
                rhs_str.emplace_back(spec->getID());
            }
            const auto& kcst = reac->getKcst();
            addCompReac(comp_id.c_str(), lhs_str, rhs_str, kcst);
        }

        // Diffusions
        for (const auto* diff: comp->getAllDiffs(&model)) {
            const auto& lig_id = diff->getLig()->getID();
            const auto dcst = diff->getDcst();
            addCompDiff({comp_id.c_str()}, {lig_id.c_str()}, dcst);
        }
    }

    auto allPatches = mesh.getAllPatches();
    std::sort(allPatches.begin(), allPatches.end(), stepsObjComp<steps::dist::DistPatch>);
    for (auto* patch: allPatches) {
        model::patch_id patch_id{patch->getID()};
        const model::compartment_id inner_compartment{patch->getIComp()->getID()};
        const auto* ocomp = patch->getOComp();
        std::optional<model::compartment_id> outer_compartment;
        if (ocomp != nullptr) {
            outer_compartment.emplace(ocomp->getID());
        }
        addPatch(patch_id, inner_compartment, outer_compartment);

        auto allPatchSpecs = patch->getAllSpecs(&model);
        std::sort(allPatchSpecs.begin(), allPatchSpecs.end(), stepsObjComp<steps::model::Spec>);
        for (auto* spec: allPatchSpecs) {
            addPatchSpec(patch_id, model::species_name(spec->getID()));
        }

        // Surface reactions
        auto allSReacs = patch->getAllSReacs(&model);
        std::sort(allSReacs.begin(), allSReacs.end(), stepsObjComp<steps::model::SReac>);
        for (auto* sreac: allSReacs) {
            std::vector<model::species_name> lhs_i, lhs_s, lhs_o, rhs_i, rhs_s, rhs_o;
            for (auto* spec: sreac->getILHS()) {
                lhs_i.emplace_back(spec->getID());
            }
            for (auto* spec: sreac->getSLHS()) {
                lhs_s.emplace_back(spec->getID());
            }
            for (auto* spec: sreac->getOLHS()) {
                lhs_o.emplace_back(spec->getID());
            }
            for (auto* spec: sreac->getIRHS()) {
                rhs_i.emplace_back(spec->getID());
            }
            for (auto* spec: sreac->getSRHS()) {
                rhs_s.emplace_back(spec->getID());
            }
            for (auto* spec: sreac->getORHS()) {
                rhs_o.emplace_back(spec->getID());
            }
            surfReacIdxs[model::surface_reaction_id(sreac->getID())] =
                addSurfReac(patch_id, lhs_i, lhs_s, lhs_o, rhs_i, rhs_s, rhs_o, sreac->getKcst());
        }

        for (auto& ssysName: patch->getSurfsys()) {
            auto* ssys = model.getSurfsys(ssysName);
            // Surface diffusions
            if (not ssys->_getAllDiffs().empty()) {
                throw std::logic_error("Model contains surface diffusion rules.");
            }
            // Voltage dependent surface reactions
            for (auto& vdepsreac: ssys->_getAllVDepSReacs()) {
                std::vector<model::species_name> lhs_i, lhs_s, lhs_o, rhs_i, rhs_s, rhs_o;
                for (auto& spec: vdepsreac.second->getILHS()) {
                    lhs_i.emplace_back(spec->getID());
                }
                for (auto& spec: vdepsreac.second->getSLHS()) {
                    lhs_s.emplace_back(spec->getID());
                }
                for (auto& spec: vdepsreac.second->getOLHS()) {
                    lhs_o.emplace_back(spec->getID());
                }
                for (auto& spec: vdepsreac.second->getIRHS()) {
                    rhs_i.emplace_back(spec->getID());
                }
                for (auto& spec: vdepsreac.second->getSRHS()) {
                    rhs_s.emplace_back(spec->getID());
                }
                for (auto& spec: vdepsreac.second->getORHS()) {
                    rhs_o.emplace_back(spec->getID());
                }
                const auto& kTable = vdepsreac.second->_getK();
                double vmin = vdepsreac.second->_getVMin();
                double vmax = vdepsreac.second->_getVMax();
                double dv = vdepsreac.second->_getDV();
                std::function<osh::Real(osh::Real)> kcst = [vmin, vmax, dv, kTable](double v) {
                    if (v > vmax) {
                        std::ostringstream msg;
                        msg << "Voltage is higher than maximum: " << v << " > " << vmax;
                        throw std::out_of_range(msg.str());
                    }
                    if (v < vmin) {
                        std::ostringstream msg;
                        msg << "Voltage is lower than minimum: " << v << " < " << vmin;
                        throw std::out_of_range(msg.str());
                    }
                    double v2 = ((v - vmin) / dv);
                    double lv = floor(v2);
                    auto lvidx = static_cast<uint>(lv);
                    uint uvidx = static_cast<uint>(ceil(v2));
                    double r = v2 - lv;
                    return (((1.0 - r) * kTable[lvidx]) + (r * kTable[uvidx]));
                };
                addVDepSurfReac(patch_id, lhs_i, lhs_s, lhs_o, rhs_i, rhs_s, rhs_o, kcst);
            }
        }
    }

    for (auto& memb: mesh.membranes()) {
        for (auto& patchId: memb.second->patches()) {
            // Membranes
            addMembrane(memb.first, patchId, memb.second->getCapacitance());
            // Channels
            auto allChans = model.getAllChans();
            std::sort(allChans.begin(), allChans.end(), stepsObjComp<steps::model::Chan>);
            for (auto* chan: allChans) {
                std::vector<model::species_name> states;
                auto allChanStates = chan->getAllChanStates();
                std::sort(allChanStates.begin(),
                          allChanStates.end(),
                          stepsObjComp<steps::model::ChanState>);
                for (auto* state: allChanStates) {
                    states.emplace_back(state->getID());
                }
                addChannel(memb.first, model::channel_id(chan->getID()), states);
            }
            // Currents
            auto* patch = mesh.getPatch(patchId);
            for (auto& ssysName: patch->getSurfsys()) {
                auto* ssys = model.getSurfsys(ssysName);
                // Ohmic currents
                for (auto& curr: ssys->_getAllOhmicCurrs()) {
                    auto* chanstate = curr.second->getChanState();
                    addOhmicCurrent(model::ohmic_current_id(curr.second->getID()),
                                    memb.first,
                                    model::channel_id(chanstate->getChan()->getID()),
                                    model::species_name(chanstate->getID()),
                                    curr.second->getG(),
                                    curr.second->getERev());
                }
                // GHK currents
                for (auto& curr: ssys->_getAllGHKcurrs()) {
                    if (not curr.second->_infosupplied()) {
                        std::ostringstream msg;
                        msg << "GHK current " << curr.first << ": Undefined permeability.";
                        throw std::invalid_argument(msg.str());
                    }
                    if (not curr.second->_realflux()) {
                        throw std::invalid_argument(
                            "GHK currents in distributed STEPS do not support "
                            "the computeflux=False argument.");
                    }
                    if (curr.second->_vshift() != 0.0) {
                        throw std::invalid_argument(
                            "GHK currents in distributed STEPS do not support "
                            "the vshift argument.");
                    }
                    auto* chanstate = curr.second->getChanState();
                    std::optional<double> voconc;
                    if (curr.second->_voconc() >= 0) {
                        voconc = curr.second->_voconc();
                    }
                    addGHKCurrentSurfReac(model::ghk_current_id(curr.second->getID()),
                                          memb.first,
                                          model::channel_id(chanstate->getChan()->getID()),
                                          model::species_name(chanstate->getID()),
                                          model::species_name(curr.second->getIon()->getID()),
                                          curr.second->_P(),
                                          curr.second->_valence(),
                                          voconc);
                }
            }
        }
    }
}

model::species_id Statedef::addSpec(const model::species_name& name) {
    auto result = specModelIdxs.find(name);
    if (result != specModelIdxs.end()) {
        return result->second;
    }
    auto model_idx = static_cast<int>(specModelIdxs.size());
    specModelIdxs.emplace(name, model_idx);
    specIDs.push_back(name);
    return specModelIdxs[name];
}

model::species_id Statedef::getSpecModelIdx(const model::species_name& name) const {
    auto result = specModelIdxs.find(name);
    if (result == specModelIdxs.end()) {
        return {};
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
    const std::optional<model::compartment_id>& outer_compartment_id) {
    container::patch_id container_id(
        static_cast<container::patch_id::value_type>(patchdefPtrs.size()));
    // we need to specify the inner compartment
    patchdefPtrs.emplace_back(std::make_unique<Patchdef>(
        *this, patchId, container_id, inner_compartment_id, outer_compartment_id));
    patchModelIdxs[patchId] = container_id;
    return container_id;
}

container::surface_reaction_id Statedef::addSurfReac(
    const model::patch_id& patchId,
    const std::vector<model::species_name>& reactants_i,
    const std::vector<model::species_name>& reactants_s,
    const std::vector<model::species_name>& reactants_o,
    const std::vector<model::species_name>& products_i,
    const std::vector<model::species_name>& products_s,
    const std::vector<model::species_name>& products_o,
    osh::Real kcst) {
    SReacInfo info{};
    info.kCst = kcst;
    return addSurfReacImpl<SReacInfo>(
        patchId, reactants_i, reactants_s, reactants_o, products_i, products_s, products_o, info);
}

container::surface_reaction_id Statedef::addVDepSurfReac(
    const model::patch_id& patchId,
    const std::vector<model::species_name>& reactants_i,
    const std::vector<model::species_name>& reactants_s,
    const std::vector<model::species_name>& reactants_o,
    const std::vector<model::species_name>& products_i,
    const std::vector<model::species_name>& products_s,
    const std::vector<model::species_name>& products_o,
    const std::function<osh::Real(osh::Real)>& kcst) {
    VDepInfo info;
    info.kCstFun = kcst;
    return addSurfReacImpl(
        patchId, reactants_i, reactants_s, reactants_o, products_i, products_s, products_o, info);
}

void Statedef::addGHKCurrentSurfReac(const model::ghk_current_id& curr_id,
                                     const model::membrane_id& membrane,
                                     const model::channel_id& ion_channel,
                                     const model::species_name& ion_channel_name,
                                     const model::species_name& ion_name,
                                     osh::Real permeability,
                                     osh::I64 valence,
                                     std::optional<osh::Real> outer_conc,
                                     std::optional<osh::Real> inner_conc) {
    auto it = membranePtrs.find(membrane);
    if (it != membranePtrs.end()) {
        Membrane& memb = *it->second;
        auto it_channel = memb.channels().find(ion_channel);
        if (it_channel != memb.channels().end()) {
            auto curr_it = ghkCurrPtrs.find(curr_id);
            if (curr_it == ghkCurrPtrs.end()) {
                curr_it =
                    ghkCurrPtrs
                        .emplace(curr_id,
                                 std::make_unique<GHKCurrent>(ion_channel_name, ion_name, valence))
                        .first;
            }
            it_channel->second.addGHKCurrent(*curr_it->second);

            const auto& channel_states = it_channel->second.channel_states;
            const auto& patch_id = memb.getPatch();
            const auto& patch = getPatchdef(patch_id);
            const auto& ion_channel_state_id = patch.getSpecPatchIdx(
                getSpecModelIdx(ion_channel_name));
            if (std::find(channel_states.begin(), channel_states.end(), ion_channel_state_id) ==
                channel_states.end()) {
                std::ostringstream msg;
                msg << "GHK current : Unknown channel state " << ion_channel_name;
                throw std::invalid_argument(msg.str());
            }
            GHKInfo info;
            info.curr_id = curr_id;
            info.in2out = true;
            info.permeability = permeability;
            info.valence = valence;
            double m3_per_liter = 1.0e-3;
            if (inner_conc) {
                info.inner_conc = *inner_conc / m3_per_liter * math::AVOGADRO;
                /*[m^{-3}]*/
            }
            if (outer_conc) {
                info.outer_conc = *outer_conc / m3_per_liter * math::AVOGADRO;
            }
            if (!inner_conc && !outer_conc) {
                // As the rate of the GHK 'reaction' depends on the concentrations of
                // both reactant and product, we add an ion that serves as a catalyst
                addSurfReacImpl(memb.getPatch(),
                                {ion_name},
                                {ion_channel_name},
                                {ion_name},
                                {},
                                {ion_channel_name},
                                {ion_name, ion_name},
                                info);
                info.in2out = false;
                addSurfReacImpl(memb.getPatch(),
                                {ion_name},
                                {ion_channel_name},
                                {ion_name},
                                {ion_name, ion_name},
                                {ion_channel_name},
                                {},
                                info);
            } else if (inner_conc && !outer_conc) {
                // As the rate of the GHK 'reaction' depends on the concentrations of
                // both reactant and product, we add an ion that serves as a catalyst
                addSurfReacImpl(memb.getPatch(),
                                {},
                                {ion_channel_name},
                                {ion_name},
                                {},
                                {ion_channel_name},
                                {ion_name, ion_name},
                                info);
                info.in2out = false;
                addSurfReacImpl(memb.getPatch(),
                                {},
                                {ion_channel_name},
                                {ion_name},
                                {},
                                {ion_channel_name},
                                {},
                                info);
            } else if (!inner_conc && outer_conc) {
                addSurfReacImpl(memb.getPatch(),
                                {ion_name},
                                {ion_channel_name},
                                {},
                                {},
                                {ion_channel_name},
                                {},
                                info);
                info.in2out = false;
                // As the rate of the GHK 'reaction' depends on the concentrations of
                // both reactant and product, we add an ion that serves as a catalyst
                addSurfReacImpl(memb.getPatch(),
                                {ion_name},
                                {ion_channel_name},
                                {},
                                {ion_name, ion_name},
                                {ion_channel_name},
                                {},
                                info);
            } else {
                addSurfReacImpl(
                    memb.getPatch(), {}, {ion_channel_name}, {}, {}, {ion_channel_name}, {}, info);
                info.in2out = false;
                addSurfReacImpl(
                    memb.getPatch(), {}, {ion_channel_name}, {}, {}, {ion_channel_name}, {}, info);
            }
        } else {
            std::ostringstream msg;
            msg << "Unknown channel " << ion_channel << " in membrane " << membrane;
            throw std::invalid_argument(msg.str());
        }
    } else {
        std::ostringstream msg;
        msg << "Unregistered membrane id " << membrane;
        throw std::invalid_argument(msg.str());
    }
}

template <typename PropensityType>
container::surface_reaction_id Statedef::addSurfReacImpl(
    const model::patch_id& patchId,
    const std::vector<model::species_name>& reactants_i,
    const std::vector<model::species_name>& reactants_s,
    const std::vector<model::species_name>& reactants_o,
    const std::vector<model::species_name>& products_i,
    const std::vector<model::species_name>& products_s,
    const std::vector<model::species_name>& products_o,
    const PropensityType kcst) {
    if (!patchModelIdxs.count(patchId)) {
        throw std::invalid_argument("Wrong patch id: " + patchId);
    }
    Patchdef& patch = *patchdefPtrs[static_cast<size_t>(patchModelIdxs[patchId].get())];

    const Compdef& inner_compartment =
        *compdefPtrs[static_cast<size_t>(compModelIdxs[patch.getInnerCompId()].get())];
    const Compdef* outer_compartment{};
    if (patch.getOuterCompId()) {
        const auto model_index = static_cast<size_t>(compModelIdxs[*patch.getOuterCompId()].get());
        outer_compartment = compdefPtrs[model_index].get();
    };
    auto& specMdlIdxs = this->specModelIdxs;
    auto convert_id_type_patch = [&patch, &specMdlIdxs](const std::vector<model::species_name>& v) {
        std::vector<container::species_id> res;
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
                                               const std::vector<model::species_name>& v) {
        std::vector<container::species_id> res;
        res.reserve(v.size());
        std::transform(v.begin(), v.end(), std::back_inserter(res), [&specMdlIdxs, &comp](auto& l) {
            auto id = comp.getSpecContainerIdx(specMdlIdxs[l]);
            if (id.unknown()) {
                std::ostringstream msg;
                msg << "Compartment " << comp.getID() << ": unknown species " << l;
                throw std::logic_error(msg.str());
            }
            return comp.getSpecContainerIdx(specMdlIdxs[l]);
        });
        return res;
    };
    std::vector<container::species_id> empty;
    auto id = patch.addSurfaceReacImpl<PropensityType>(
        convert_id_type_comp(inner_compartment, reactants_i),
        convert_id_type_patch(reactants_s),
        !reactants_o.empty() ? convert_id_type_comp(*outer_compartment, reactants_o) : empty,
        convert_id_type_comp(inner_compartment, products_i),
        convert_id_type_patch(products_s),
        !products_o.empty() ? convert_id_type_comp(*outer_compartment, products_o) : empty,
        kcst);
    return id;
}

container::compartment_id Statedef::getCompModelIdx(
    const model::compartment_id& compartment) const noexcept {
    const auto id = compModelIdxs.find(compartment);
    if (id == compModelIdxs.end()) {
        return {};
    }
    return id->second;
}

Compdef& Statedef::getCompdef(container::compartment_id compartment) const noexcept {
    return *compdefPtrs[static_cast<size_t>(compartment.get())];
}

Compdef& Statedef::getCompdef(const model::compartment_id& compartment) const noexcept {
    return *compdefPtrs[static_cast<size_t>(getCompModelIdx(compartment).get())];
}

Patchdef& Statedef::getPatchdef(const container::patch_id& patchId) const noexcept {
    return *patchdefPtrs[static_cast<size_t>(patchId.get())];
}

Patchdef& Statedef::getPatchdef(const model::patch_id& patchId) const noexcept {
    return *patchdefPtrs[static_cast<size_t>(patchModelIdxs.at(patchId).get())];
}

model::species_id Statedef::addCompSpec(const model::compartment_id& compartment,
                                        const model::species_name& species) {
    auto spec_id = addSpec(species);
    auto& compdef = compdefPtrs[static_cast<size_t>(compModelIdxs[compartment].get())];
    compdef->addSpec(spec_id);
    return spec_id;
}

std::vector<model::species_id> Statedef::addCompSpecs(
    const model::compartment_id& compartment,
    const std::vector<model::species_name>& species) {
    std::vector<model::species_id> ids;
    ids.reserve(species.size());
    for (const auto& spec_id: species) {
        ids.push_back(addCompSpec(compartment, spec_id));
    }
    return ids;
}

std::vector<model::species_id> Statedef::addPatchSpecs(
    const model::patch_id& patchId,
    const std::vector<model::species_name>& species) {
    std::vector<model::species_id> ids;
    ids.reserve(species.size());
    for (const auto& spec_id: species) {
        ids.push_back(addPatchSpec(patchId, spec_id));
    }
    return ids;
}

model::species_id Statedef::addPatchSpec(const model::patch_id& patch_id,
                                         const model::species_name& species) {
    auto spec_id = addSpec(species);
    auto& patchdef = patchdefPtrs[static_cast<size_t>(patchModelIdxs[patch_id].get())];
    patchdef->addSpec(spec_id);
    return spec_id;
}

container::species_id Statedef::getCompSpecContainerIdx(const model::compartment_id& compartment,
                                                        const model::species_name& specie) const {
    const auto comp_model_idx = getCompModelIdx(compartment);
    if (comp_model_idx.unknown()) {
        std::ostringstream msg;
        msg << "Unknown compartment: " << compartment;
        throw std::invalid_argument(msg.str());
    }
    const auto spec_model_idx = getSpecModelIdx(specie);
    if (spec_model_idx.unknown()) {
        std::ostringstream msg;
        msg << "Unknown species: " << specie;
        throw std::invalid_argument(msg.str());
    }

    return compdefPtrs[static_cast<size_t>(comp_model_idx.get())]->getSpecContainerIdx(
        spec_model_idx);
}

container::diffusion_id Statedef::addCompDiff(const model::compartment_id& compartment,
                                              const model::species_name& species_name,
                                              osh::Real dcst) {
    const auto scidx = getCompSpecContainerIdx(compartment, species_name);
    return compdefPtrs[static_cast<size_t>(compModelIdxs[compartment].get())]->addDiff(scidx, dcst);
}

container::reaction_id Statedef::addCompReac(const model::compartment_id& compartment,
                                             const std::vector<model::species_name>& reactants,
                                             const std::vector<model::species_name>& products,
                                             osh::Real kcst) {
    std::vector<container::species_id> reactants_ids;
    std::vector<container::species_id> products_ids;
    reactants_ids.reserve(reactants.size());
    products_ids.reserve(products.size());

    std::transform(reactants.begin(),
                   reactants.end(),
                   std::back_inserter(reactants_ids),
                   [this, &compartment](const auto& species) {
                       return this->getCompSpecContainerIdx(compartment, species);
                   });

    std::transform(products.begin(),
                   products.end(),
                   std::back_inserter(products_ids),
                   [this, &compartment](const auto& species) {
                       return this->getCompSpecContainerIdx(compartment, species);
                   });

    auto& compartmentPtr = compdefPtrs[static_cast<size_t>(compModelIdxs[compartment].get())];
    return compartmentPtr->addReac(reactants_ids, products_ids, kcst);
}

const container::surface_reaction_id& Statedef::getSReacIdx(
    const model::surface_reaction_id& reac) const {
    return map_at(surfReacIdxs, reac, "surface reaction");
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

void Statedef::addOhmicCurrent(const model::ohmic_current_id& curr_id,
                               const model::membrane_id& membrane,
                               const model::channel_id& channel,
                               const std::optional<model::species_name>& species_name,
                               double conductance,
                               double reversal_potential) {
    auto& memb = map_at(membranePtrs, membrane, "membrane");
    auto& channel_def = map_at(memb->channels(), channel, "channel");
    const auto& patch_id = memb->getPatch();
    const auto& patch = getPatchdef(patch_id);
    auto curr_it = ohmicCurrPtrs.find(curr_id);
    if (curr_it == ohmicCurrPtrs.end()) {
        std::optional<container::species_id> spec;
        if (species_name) {
            spec = patch.getSpecPatchIdx(getSpecModelIdx(*species_name));
        }
        curr_it =
            ohmicCurrPtrs
                .emplace(curr_id,
                         std::make_unique<OhmicCurrent>(conductance, reversal_potential, spec))
                .first;
    }
    channel_def.addOhmicCurrent(*curr_it->second);
}

void Statedef::setStimulus(const model::membrane_id& membrane, osh::Real current) {
    auto& memb = map_at(membranePtrs, membrane, "membrane");
    memb->setStimulus([current](auto) { return current; });
}

void Statedef::setResistance(const model::membrane_id& membrane, osh::Real resistance) {
    if (resistance <= 0) {
        throw std::invalid_argument("resistance must be positive");
    }
    auto& memb = map_at(membranePtrs, membrane, "membrane");
    memb->setConductance(1.0 / resistance);
}

osh::Real Statedef::getResistance(const model::membrane_id& membrane) {
    auto& memb = map_at(membranePtrs, membrane, "membrane");
    return 1.0 / memb->conductance();
}

void Statedef::setReversalPotential(const model::membrane_id& membrane,
                                    osh::Real reversal_potential) {
    auto& memb = map_at(membranePtrs, membrane, "membrane");
    memb->setReversalPotential(reversal_potential);
}

osh::Real Statedef::getReversalPotential(const model::membrane_id& membrane) {
    auto& memb = map_at(membranePtrs, membrane, "membrane");
    return memb->reversal_potential();
}

void Statedef::addCompartmentConductivity(const model::compartment_id& compartment_id,
                                          const osh::Real conductivity) {
    compartment_conductivity_[compartment_id] = conductivity;
}

osh::Real Statedef::getCompartmentConductivity(const model::compartment_id& compartment) const {
    auto cond = map_at(compartment_conductivity_, compartment, "compartment");
    if (cond <= 0) {
        std::ostringstream msg;
        msg << "Invalid or undeclared conductivity for compartment with id " << compartment << ".";
        throw std::invalid_argument(msg.str());
    }
    return cond;
}

void Statedef::addChannel(const model::membrane_id& membrane,
                          const model::channel_id& channel,
                          const std::vector<model::species_name>& channel_states) {
    auto& memb = map_at(membranePtrs, membrane, "membrane");
    const auto& patch_id = memb->getPatch();
    const auto& patch = getPatchdef(patch_id);
    std::vector<container::species_id> specs;
    specs.reserve(channel_states.size());
    for (const auto& channel_state: channel_states) {
        specs.push_back(patch.getSpecPatchIdx(addPatchSpec(patch_id, channel_state)));
    }
    memb->addChannel(channel, Channel(specs));
}

}  // namespace steps::dist
