#include "sreacdef.hpp"

#include "model/chanstate.hpp"
#include "model/complexevents.hpp"
#include "model/complexsreac.hpp"
#include "model/ghkcurr.hpp"
#include "model/model.hpp"
#include "model/surfsys.hpp"
#include "patchdef.hpp"
#include "statedef.hpp"
#include "util/vocabulary.hpp"
#include <type_traits>

namespace steps::dist {

//-------------------------------------------------------

using SL = steps::model::Location;
using SC = SpecieClassifier;


template <typename VDepSReacT, typename>
VDepInfo getPropInfo(const VDepSReacT& reac) {
    VDepInfo info;
    const auto& kTable = reac._getK();
    double vmin = reac._getVMin();
    double vmax = reac._getVMax();
    double dv = reac._getDV();
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
    info.kCstFun = kcst;
    // Surface reaction charge has opposite convention from GHK
    info.charge = -reac.getCharge();
    return info;
}

template <typename CurrType, typename>
GHKInfo getPropInfo(const CurrType& curr, bool in2out) {
    GHKInfo info;
    info.curr_id = model::ghk_current_id(curr.getID());
    info.in2out = in2out;
    info.permeability = curr._P();
    info.charge = curr._valence();
    double m3_per_liter = 1.0e-3;
    if (curr._voconc() >= 0) {
        info.outer_conc = curr._voconc() / m3_per_liter * math::AVOGADRO;
    }
    return info;
}

template <typename SReacT, typename>
SReacInfo getPropInfo(const SReacT& reac) {
    SReacInfo info;
    info.kCst = reac.getKcst();
    // Surface reaction charge has opposite convention from GHK
    info.charge = -reac.getCharge();
    return info;
}

template <typename CurrType>
SurfaceReactionComponents getGHKIonsReactionComponents(const Patchdef& patchdef,
                                                       const CurrType& curr,
                                                       bool in2out) {
    SurfaceReactionComponents components;

    auto ionIdIn = patchdef.getInnerComp().getSpecContainerIdx(curr.getIon());
    // As the rate of the GHK 'reaction' depends on the concentrations of
    // both reactant and product, we add an ion that serves as a catalyst
    components.emplace_back(ionIdIn, SC::Reactant, SL::PATCH_IN);
    if (not in2out) {
        components.emplace_back(ionIdIn, SC::Product, SL::PATCH_IN);
        components.emplace_back(ionIdIn, SC::Product, SL::PATCH_IN);
    }

    if (curr._voconc() < 0) {
        assert(patchdef.getOuterCompId().has_value());
        const Compdef& compdef = patchdef.statedef().getCompdef(patchdef.getOuterCompId().value());
        auto ionIdOut = compdef.getSpecContainerIdx(curr.getIon());
        components.emplace_back(ionIdOut, SC::Reactant, SL::PATCH_OUT);
        if (in2out) {
            components.emplace_back(ionIdOut, SC::Product, SL::PATCH_OUT);
            components.emplace_back(ionIdOut, SC::Product, SL::PATCH_OUT);
        }
    }

    return components;
}

SurfaceReactionComponents getReactionComponents(const Patchdef& patchdef,
                                                const steps::model::GHKcurr& curr,
                                                bool in2out) {
    SurfaceReactionComponents components = getGHKIonsReactionComponents(patchdef, curr, in2out);

    auto chanId = patchdef.getSpecContainerIdx(curr.getChanState());

    components.emplace_back(chanId, SC::Reactant, SL::PATCH_SURF);
    components.emplace_back(chanId, SC::Product, SL::PATCH_SURF);

    return components;
}

SurfaceReactionComponents getReactionComponents(const Patchdef& patchdef,
                                                const steps::model::ComplexGHKcurr& curr,
                                                bool in2out) {
    return getGHKIonsReactionComponents(patchdef, curr, in2out);
}

template <typename SReacT>
SurfaceReactionComponents getReactionComponents(const Patchdef& patchdef, const SReacT& reac) {
    SurfaceReactionComponents components;
    for (auto* spec: reac.getILHS()) {
        auto csid = patchdef.getInnerComp().getSpecContainerIdx(*spec);
        components.emplace_back(csid, SC::Reactant, SL::PATCH_IN);
    }
    for (auto* spec: reac.getIRHS()) {
        container::species_id csid = patchdef.getInnerComp().getSpecContainerIdx(*spec);
        components.emplace_back(csid, SC::Product, SL::PATCH_IN);
    }
    for (auto* spec: reac.getSLHS()) {
        container::species_id csid = patchdef.getSpecContainerIdx(*spec);
        components.emplace_back(csid, SC::Reactant, SL::PATCH_SURF);
    }
    for (auto* spec: reac.getSRHS()) {
        container::species_id csid = patchdef.getSpecContainerIdx(*spec);
        components.emplace_back(csid, SC::Product, SL::PATCH_SURF);
    }
    if (patchdef.getOuterCompId().has_value()) {
        const Compdef& compdef = patchdef.statedef().getCompdef(patchdef.getOuterCompId().value());
        for (auto* spec: reac.getOLHS()) {
            container::species_id csid = compdef.getSpecContainerIdx(*spec);
            components.emplace_back(csid, SC::Reactant, SL::PATCH_OUT);
        }
        for (auto* spec: reac.getORHS()) {
            container::species_id csid = compdef.getSpecContainerIdx(*spec);
            components.emplace_back(csid, SC::Product, SL::PATCH_OUT);
        }
    }
    return components;
}

template <typename PropensityInfo>
SReacdefBase<PropensityInfo>::SReacdefBase(const Patchdef& patchdef,
                                           container::kproc_id kproc,
                                           ContReacID surface_reaction_id,
                                           const SurfaceReactionComponents& reaction_definition,
                                           PropensityInfo info)
    : patchdef_(patchdef)
    , kproc_id_(kproc)
    , surface_reaction_(surface_reaction_id)
    , info_(info) {
    for (auto loc: steps::model::AllPatchLocations) {
        lhs_[loc];
        upd_[loc];
    }
    for (const auto& t: reaction_definition) {
        const auto& [spec_id, spec_class, loc] = t;
        if (spec_class == SC::Reactant) {
            order_++;
            lhs_[loc][spec_id] -= 1;
            upd_[loc][spec_id] -= 1;
        } else {
            upd_[loc][spec_id] += 1;
        }
    }
    // Remove zero elements
    for (auto& [_, d]: this->upd_) {
        auto it = d.begin();
        while (it != d.end()) {
            if (it->second == 0) {
                it = d.erase(it);
            } else {
                it++;
            }
        }
    }

    const bool all_inner = lhs_[SL::PATCH_OUT].empty();
    const bool all_outer = lhs_[SL::PATCH_IN].empty();
    if (!(all_inner || all_outer) &&
        typeid(PropensityInfo).hash_code() != typeid(GHKInfo).hash_code()) {
        throw std::logic_error(
            "A surface reaction involves volume reactants "
            "all in the inner or all in the outer compartment");
    }
    is_inner_ = all_inner;
    is_surf_surf_ = lhs_[SL::PATCH_OUT].empty() && lhs_[SL::PATCH_IN].empty();
}


template <typename SReacT>
ModelSReacdef<SReacT>::ModelSReacdef(const Patchdef& patchdef,
                                     container::kproc_id kproc,
                                     ContReacID surface_reaction_id,
                                     const SReacT& reac)
    : SReacdefBase<SReacT>::template SReacdefBase<SReacT>(patchdef,
                                                          kproc,
                                                          surface_reaction_id,
                                                          getReactionComponents(patchdef, reac),
                                                          getPropInfo(reac))
    , reac_(reac) {}

template <typename SReacT>
void ModelSReacdef<SReacT>::reset() {
    SReacdefBase<SReacT>::reset();
    this->info_ = getPropInfo(reac_);
}

template <typename GHKcurrT>
GHKSReacdefBase<GHKcurrT>::GHKSReacdefBase(const Patchdef& patchdef,
                                           container::kproc_id kproc,
                                           ContReacID surface_reaction_id,
                                           const GHKcurrT& curr,
                                           bool in2out)
    : SReacdefBase<GHKcurrT>::template SReacdefBase<GHKcurrT>(patchdef,
                                                              kproc,
                                                              surface_reaction_id,
                                                              getReactionComponents(patchdef,
                                                                                    curr,
                                                                                    in2out),
                                                              getPropInfo(curr, in2out))
    , curr_(curr)
    , in2out_(in2out) {}

template <typename GHKcurrT>
void GHKSReacdefBase<GHKcurrT>::reset() {
    SReacdefBase<GHKcurrT>::reset();
    this->info_ = getPropInfo(curr_, in2out_);
}

GHKSReacdef::GHKSReacdef(const Patchdef& patchdef,
                         container::kproc_id kproc,
                         ContReacID surface_reaction_id,
                         const steps::model::GHKcurr& curr,
                         bool in2out)
    : GHKSReacdefBase<steps::model::GHKcurr>::template GHKSReacdefBase<steps::model::GHKcurr>(
          patchdef,
          kproc,
          surface_reaction_id,
          curr,
          in2out) {}

ComplexGHKSReacdef::ComplexGHKSReacdef(const Patchdef& patchdef,
                                       container::kproc_id kproc,
                                       ContReacID surface_reaction_id,
                                       const steps::model::ComplexGHKcurr& curr,
                                       bool in2out)
    : GHKSReacdefBase<steps::model::ComplexGHKcurr>::template GHKSReacdefBase<
          steps::model::ComplexGHKcurr>(patchdef, kproc, surface_reaction_id, curr, in2out)
    , pComplexId(patchdef.statedef().getComplexModelIdx(
          model::complex_name(curr.getChanState().complexId))) {
    auto updEv = curr.getUpdEvent();
    pComplexSurfUPDEv = std::make_shared<ComplexUpdateEventdef>(updEv, patchdef.statedef());
    pComplex_DEPSET = pComplexSurfUPDEv->getDepSet();
}

template <typename CSReacT>
ModelComplexSReacdef<CSReacT>::ModelComplexSReacdef(const Patchdef& patchdef,
                                                    container::kproc_id kproc,
                                                    ContReacID surface_reaction_id,
                                                    const CSReacT& reac)
    : ModelSReacdef<CSReacT>::ModelSReacdef(patchdef, kproc, surface_reaction_id, reac) {
    // Copy complex events
    const auto& sd = patchdef.statedef();
    for (auto loc: steps::model::AllPatchLocations) {
        pComplex_DEPMAP[loc];
        pComplex_UPDMAP[loc];
        for (auto* ev: reac.getUPDEvents(loc)) {
            pComplexUPDEvs[loc].push_back(std::make_shared<ComplexUpdateEventdef>(*ev, sd));
            this->order_++;
            this->is_surf_surf_ &= loc == SL::PATCH_SURF;
        }
        for (auto* ev: reac.getDELEvents(loc)) {
            pComplexDELEvs[loc].push_back(std::make_shared<ComplexDeleteEventdef>(*ev, sd));
            this->order_++;
            this->is_surf_surf_ &= loc == SL::PATCH_SURF;
        }
        for (auto* ev: reac.getCREEvents(loc)) {
            pComplexCREEvs[loc].push_back(std::make_shared<ComplexCreateEventdef>(*ev, sd));
        }
    }
    const bool all_inner = this->lhs_[SL::PATCH_OUT].empty() &&
                           pComplexUPDEvs[SL::PATCH_OUT].empty() &&
                           pComplexDELEvs[SL::PATCH_OUT].empty();
    const bool all_outer = this->lhs_[SL::PATCH_IN].empty() &&
                           pComplexUPDEvs[SL::PATCH_IN].empty() &&
                           pComplexDELEvs[SL::PATCH_IN].empty();
    if (!(all_inner || all_outer)) {
        throw std::logic_error(
            "A surface reaction involves volume reactants "
            "all in the inner or all in the outer compartment");
    }
    this->is_inner_ = all_inner;

    // set up deps for complexes
    for (auto loc: steps::model::AllPatchLocations) {
        for (const auto& upd: pComplexUPDEvs[loc]) {
            pComplex_DEPMAP[loc][upd->complexIdx()].merge(upd->getDepSet());
            if (loc == upd->destLoc()) {
                // If the complex stays in the same location
                pComplex_UPDMAP[loc][upd->complexIdx()].merge(upd->getUpdSet());
            } else {
                // If the complex moves to a different location, it is equivalent to a delete and a
                // create
                const auto& filts = upd->filters();
                auto& updset1 = pComplex_UPDMAP[loc][upd->complexIdx()];
                auto& updset2 = pComplex_UPDMAP[upd->destLoc()][upd->complexIdx()];
                for (auto& filt: filts) {
                    for (auto sus: filt.range()) {
                        if (filt[sus].max > 0) {
                            updset1.insert(sus);
                        }
                        for (auto updt: upd->updates()) {
                            if (filt[sus].max + updt.update[sus] > 0) {
                                updset2.insert(sus);
                            }
                        }
                    }
                }
            }
        }
        for (const auto& del: pComplexDELEvs[loc]) {
            pComplex_DEPMAP[loc][del->complexIdx()].merge(del->getDepSet());
            pComplex_UPDMAP[loc][del->complexIdx()].merge(del->getUpdSet());
        }
        for (const auto& cre: pComplexCREEvs[loc]) {
            pComplex_UPDMAP[loc][cre->complexIdx()].merge(cre->getUpdSet());
        }
    }
}


template class ModelSReacdef<steps::model::SReac>;
template class ModelSReacdef<steps::model::ComplexSReac>;
template class ModelSReacdef<steps::model::VDepSReac>;
template class ModelSReacdef<steps::model::VDepComplexSReac>;

template class ModelComplexSReacdef<steps::model::ComplexSReac>;
template class ModelComplexSReacdef<steps::model::VDepComplexSReac>;

template class GHKSReacdefBase<steps::model::GHKcurr>;
template class GHKSReacdefBase<steps::model::ComplexGHKcurr>;

}  // namespace steps::dist
