/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */

#include "comp.hpp"

#include "geom.hpp"
#include "patch.hpp"

#include "model/complexreac.hpp"
#include "model/complexsreac.hpp"
#include "model/diff.hpp"
#include "model/ghkcurr.hpp"
#include "model/model.hpp"
#include "model/raftsreac.hpp"
#include "model/raftsys.hpp"
#include "model/reac.hpp"
#include "model/spec.hpp"
#include "model/sreac.hpp"
#include "model/surfsys.hpp"
#include "model/vdepsreac.hpp"
#include "model/vessreac.hpp"
#include "model/vessurfsys.hpp"
#include "model/volsys.hpp"

#include "util/error.hpp"

namespace steps::wm {

Comp::Comp(std::string id, Geom& container, double vol)
    : pVol(vol)
    , pID(std::move(id))
    , pContainer(container) {
    ArgErrLogIf(pVol < 0.0, "Compartment volume can't be negative.");
    pContainer._handleCompAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

Comp::~Comp() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setID(std::string const& id) {
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pContainer._handleCompIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setVol(double vol) {
    ArgErrLogIf(vol < 0.0, "Compartment volume can't be negative.");
    pVol = vol;
}

////////////////////////////////////////////////////////////////////////////////

void Comp::addVolsys(std::string const& id) {
    // string identifier is only added to set if it is not already included
    pVolsys.insert(id);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::delVolsys(std::string const& id) {
    // string identifier is only removed from set if it is included
    pVolsys.erase(id);
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<model::Spec*> Comp::getAllSpecs(const model::Model& model) const {
    util::flat_set<model::Spec*> specs;
    // Add those from reacs, complex reacs and diffs, in volsys
    for (const auto& volsys_id: pVolsys) {
        model::Volsys& volsys = model.getVolsys(volsys_id);
        const auto& volsys_specs = volsys.getAllSpecs();
        specs.insert(volsys_specs.begin(), volsys_specs.end());
    }

    for (auto const& ipatch: getIPatches()) {
        for (auto const& ssys_id: ipatch->getSurfsys()) {
            model::Surfsys& surfsys = model.getSurfsys(ssys_id);

            // Since this is known as inner patch, patch knows this comp as OUTER
            for (auto const& sreac: surfsys.getAllSReacs()) {
                const auto& sreac_specs_olhs = sreac->getOLHS();
                specs.insert(sreac_specs_olhs.begin(), sreac_specs_olhs.end());
                const auto& sreac_specs_orhs = sreac->getORHS();
                specs.insert(sreac_specs_orhs.begin(), sreac_specs_orhs.end());
            }
            for (auto const& vdepsreac: surfsys.getAllVDepSReacs()) {
                const auto& vdepsreac_specs_olhs = vdepsreac->getOLHS();
                specs.insert(vdepsreac_specs_olhs.begin(), vdepsreac_specs_olhs.end());
                const auto& vdepsreac_specs_orhs = vdepsreac->getORHS();
                specs.insert(vdepsreac_specs_orhs.begin(), vdepsreac_specs_orhs.end());
            }
            for (auto const& cplxsreac: surfsys.getAllComplexSReacs()) {
                const auto& cplxsreac_specs_olhs = cplxsreac->getOLHS();
                specs.insert(cplxsreac_specs_olhs.begin(), cplxsreac_specs_olhs.end());
                const auto& cplxsreac_specs_orhs = cplxsreac->getORHS();
                specs.insert(cplxsreac_specs_orhs.begin(), cplxsreac_specs_orhs.end());
            }
            for (auto const& ghk: surfsys.getAllGHKcurrs()) {
                // This is OUTER so only add if not virtual outer conc
                if (ghk->_voconc() < 0.0) {
                    specs.insert(&ghk->getIon());
                }
            }
        }
    }

    for (auto const& opatch: getOPatches()) {
        for (auto const& ssys_id: opatch->getSurfsys()) {
            model::Surfsys& surfsys = model.getSurfsys(ssys_id);

            // Since this is known as outer patch, patch knows this comp as INNER
            for (auto const& sreac: surfsys.getAllSReacs()) {
                const auto& sreac_specs_ilhs = sreac->getILHS();
                specs.insert(sreac_specs_ilhs.begin(), sreac_specs_ilhs.end());
                const auto& sreac_specs_irhs = sreac->getIRHS();
                specs.insert(sreac_specs_irhs.begin(), sreac_specs_irhs.end());
            }
            for (auto const& vdepsreac: surfsys.getAllVDepSReacs()) {
                const auto& vdepsreac_specs_ilhs = vdepsreac->getILHS();
                specs.insert(vdepsreac_specs_ilhs.begin(), vdepsreac_specs_ilhs.end());
                const auto& vdepsreac_specs_irhs = vdepsreac->getIRHS();
                specs.insert(vdepsreac_specs_irhs.begin(), vdepsreac_specs_irhs.end());
            }
            for (auto const& cplxsreac: surfsys.getAllComplexSReacs()) {
                const auto& cplxsreac_specs_ilhs = cplxsreac->getILHS();
                specs.insert(cplxsreac_specs_ilhs.begin(), cplxsreac_specs_ilhs.end());
                const auto& cplxsreac_specs_irhs = cplxsreac->getIRHS();
                specs.insert(cplxsreac_specs_irhs.begin(), cplxsreac_specs_irhs.end());
            }
            for (auto const& ghk: surfsys.getAllGHKcurrs()) {
                if (ghk->_realflux()) {
                    specs.insert(&ghk->getIon());
                }
            }
        }
    }

    // Comps have all volumetric species from vesicle surface reactions and raft surface reactions
    // added since vesicles and rafts can appear in any compartment or bordering patch
    for (auto const& vessurfsys: model.getAllVesSurfsyss()) {
        for (auto const& vessreac: vessurfsys->getAllVesSReacs()) {
            const auto& vessreac_specs_olhs = vessreac->getOLHS();
            specs.insert(vessreac_specs_olhs.begin(), vessreac_specs_olhs.end());
            const auto& vessreac_specs_orhs = vessreac->getORHS();
            specs.insert(vessreac_specs_orhs.begin(), vessreac_specs_orhs.end());
        }
    }
    for (auto const& raftsys: model.getAllRaftsyss()) {
        for (auto const& raftsreac: raftsys->getAllRaftSReacs()) {
            const auto& raftsreac_specs_olhs = raftsreac->getOLHS();
            specs.insert(raftsreac_specs_olhs.begin(), raftsreac_specs_olhs.end());
            const auto& raftsreac_specs_orhs = raftsreac->getORHS();
            specs.insert(raftsreac_specs_orhs.begin(), raftsreac_specs_orhs.end());
            const auto& raftsreac_specs_ilhs = raftsreac->getILHS();
            specs.insert(raftsreac_specs_ilhs.begin(), raftsreac_specs_ilhs.end());
            const auto& raftsreac_specs_irhs = raftsreac->getIRHS();
            specs.insert(raftsreac_specs_irhs.begin(), raftsreac_specs_irhs.end());
        }
    }

    return specs;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<model::Reac*> Comp::getAllReacs(const model::Model& model) const {
    util::flat_set<model::Reac*> pReacs;
    for (const auto& volsys_id: pVolsys) {
        model::Volsys& volsys = model.getVolsys(volsys_id);
        const auto& reacs = volsys.getAllReacs();
        pReacs.insert(reacs.begin(), reacs.end());
    }
    return pReacs;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<model::Diff*> Comp::getAllDiffs(const model::Model& model) const {
    util::flat_set<model::Diff*> pDiffs;
    for (const auto& volsys_id: pVolsys) {
        model::Volsys& volsys = model.getVolsys(volsys_id);
        const auto& diffs = volsys.getAllDiffs();
        pDiffs.insert(diffs.begin(), diffs.end());
    }

    return pDiffs;
}

////////////////////////////////////////////////////////////////////////////////

void Comp::_addIPatch(Patch& patch) {
    AssertLog(patch.getOComp() == this);
    // patch pointer is only added to set if it is not already included
    pIPatches.insert(&patch);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::_delIPatch(Patch& patch) {
    AssertLog(patch.getOComp() == this);
    pIPatches.erase(&patch);
}
////////////////////////////////////////////////////////////////////////////////

void Comp::_addOPatch(Patch& patch) {
    AssertLog(&patch.getIComp() == this);
    // patch pointer is only added to set if it is not already included
    pOPatches.insert(&patch);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::_delOPatch(Patch& patch) {
    AssertLog(&patch.getIComp() == this);
    pOPatches.erase(&patch);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::_handleSelfDelete() {
    pContainer._handleCompDel(*this);
}

////////////////////////////////////////////////////////////////////////////////

wm::Patch* Comp::_getIPatch(uint lidx) const {
    AssertLog(lidx < pIPatches.size());
    auto pit = pIPatches.begin();
    std::advance(pit, lidx);
    return *pit;
}

////////////////////////////////////////////////////////////////////////////////

wm::Patch* Comp::_getOPatch(uint lidx) const {
    AssertLog(lidx < pOPatches.size());
    auto pit = pOPatches.begin();
    std::advance(pit, lidx);
    return *pit;
}

}  // namespace steps::wm
