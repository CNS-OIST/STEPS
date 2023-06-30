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

#include <cassert>
#include <sstream>

#include <easylogging++.h>

#include "geom.hpp"
#include "patch.hpp"

#include "model/ghkcurr.hpp"
#include "model/model.hpp"
#include "model/sreac.hpp"
#include "model/vdepsreac.hpp"
#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::wm {

Comp::Comp(std::string id, Geom* container, double vol)
    : pVol(vol)
    , pID(std::move(id))
    , pContainer(container) {
    ArgErrLogIf(pContainer == nullptr, "No container provided to Comp initializer function.");
    ArgErrLogIf(pVol < 0.0, "Compartment volume can't be negative.");
    pContainer->_handleCompAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Comp::~Comp() {
    if (pContainer == nullptr) {
        return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setID(std::string const& id) {
    AssertLog(pContainer != nullptr);
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pContainer->_handleCompIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setVol(double vol) {
    AssertLog(pContainer != nullptr);
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

std::vector<model::Spec*> Comp::getAllSpecs(const model::Model* model) const {
    std::set<model::Spec*> specs;
    std::set<std::string>::iterator it;

    // Add those from reacs and diffs, in volsys
    for (it = pVolsys.begin(); it != pVolsys.end(); it++) {
        model::Volsys* volsys = model->getVolsys(*it);
        std::vector<model::Spec*> volsys_specs = volsys->getAllSpecs();
        specs.insert(volsys_specs.begin(), volsys_specs.end());
    }

    for (auto const& ipatch: getIPatches()) {
        for (auto const& ssys_id: ipatch->getSurfsys()) {
            model::Surfsys* surfsys = model->getSurfsys(ssys_id);

            for (auto const& sreac: surfsys->getAllSReacs()) {
                // Since this is known as inner patch, patch knows this comp as OUTER
                const auto& sreac_specs_olhs = sreac->getOLHS();
                specs.insert(sreac_specs_olhs.begin(), sreac_specs_olhs.end());
                const auto& sreac_specs_orhs = sreac->getORHS();
                specs.insert(sreac_specs_orhs.begin(), sreac_specs_orhs.end());
            }
            for (auto const& vdepsreac: surfsys->getAllVDepSReacs()) {
                // Since this is known as inner patch, patch knows this comp as OUTER
                const auto& vdepsreac_specs_olhs = vdepsreac->getOLHS();
                specs.insert(vdepsreac_specs_olhs.begin(), vdepsreac_specs_olhs.end());
                const auto& vdepsreac_specs_orhs = vdepsreac->getORHS();
                specs.insert(vdepsreac_specs_orhs.begin(), vdepsreac_specs_orhs.end());
            }
            for (auto const& ghk: surfsys->getAllGHKcurrs()) {
                // This is OUTER so only add if not virtual outer conc
                if (ghk->_voconc() < 0.0) {
                    specs.insert(ghk->getIon());
                }
            }
        }
    }

    for (auto const& opatch: getOPatches()) {
        for (auto const& ssys_id: opatch->getSurfsys()) {
            model::Surfsys* surfsys = model->getSurfsys(ssys_id);

            for (auto const& sreac: surfsys->getAllSReacs()) {
                // Since this is known as outer patch, patch knows this comp as INNER
                const auto& sreac_specs_ilhs = sreac->getILHS();
                specs.insert(sreac_specs_ilhs.begin(), sreac_specs_ilhs.end());
                const auto& sreac_specs_irhs = sreac->getIRHS();
                specs.insert(sreac_specs_irhs.begin(), sreac_specs_irhs.end());
            }
            for (auto const& vdepsreac: surfsys->getAllVDepSReacs()) {
                // Since this is known as outer patch, patch knows this comp as INNER
                const auto& vdepsreac_specs_ilhs = vdepsreac->getILHS();
                specs.insert(vdepsreac_specs_ilhs.begin(), vdepsreac_specs_ilhs.end());
                const auto& vdepsreac_specs_irhs = vdepsreac->getIRHS();
                specs.insert(vdepsreac_specs_irhs.begin(), vdepsreac_specs_irhs.end());
            }
            for (auto const& ghk: surfsys->getAllGHKcurrs()) {
                if (ghk->_realflux()) {
                    specs.insert(ghk->getIon());
                }
            }
        }
    }

    return {specs.begin(), specs.end()};
}

////////////////////////////////////////////////////////////////////////////////

std::vector<model::Reac*> Comp::getAllReacs(const model::Model* model) const {
    std::set<model::Reac*> pReacs;
    std::set<std::string>::iterator it;
    for (it = pVolsys.begin(); it != pVolsys.end(); it++) {
        model::Volsys* volsys = model->getVolsys(*it);
        std::vector<model::Reac*> reacs = volsys->getAllReacs();
        pReacs.insert(reacs.begin(), reacs.end());
    }

    return {pReacs.begin(), pReacs.end()};
}

////////////////////////////////////////////////////////////////////////////////

std::vector<model::Diff*> Comp::getAllDiffs(const model::Model* model) const {
    std::set<model::Diff*> pDiffs;
    std::set<std::string>::iterator it;
    for (it = pVolsys.begin(); it != pVolsys.end(); it++) {
        model::Volsys* volsys = model->getVolsys(*it);
        std::vector<model::Diff*> diffs = volsys->getAllDiffs();
        pDiffs.insert(diffs.begin(), diffs.end());
    }

    return {pDiffs.begin(), pDiffs.end()};
}

////////////////////////////////////////////////////////////////////////////////

void Comp::_addIPatch(Patch* patch) {
    AssertLog(patch->getOComp() == this);
    // patch pointer is only added to set if it is not already included
    pIPatches.insert(patch);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::_delIPatch(Patch* patch) {
    AssertLog(patch->getOComp() == this);
    pIPatches.erase(patch);
}
////////////////////////////////////////////////////////////////////////////////

void Comp::_addOPatch(Patch* patch) {
    AssertLog(patch->getIComp() == this);
    // patch pointer is only added to set if it is not already included
    pOPatches.insert(patch);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::_delOPatch(Patch* patch) {
    AssertLog(patch->getIComp() == this);
    pOPatches.erase(patch);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::_handleSelfDelete() {
    pContainer->_handleCompDel(this);
    pVol = 0.0;
    pVolsys.clear();
    pIPatches.clear();
    pOPatches.clear();
    pContainer = nullptr;
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
