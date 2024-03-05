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

#include "patch.hpp"

#include "model/chanstate.hpp"
#include "model/complexsreac.hpp"
#include "model/diff.hpp"
#include "model/endocytosis.hpp"
#include "model/ghkcurr.hpp"
#include "model/model.hpp"
#include "model/ohmiccurr.hpp"
#include "model/raftgen.hpp"
#include "model/raftsreac.hpp"
#include "model/raftsys.hpp"
#include "model/sreac.hpp"
#include "model/surfsys.hpp"
#include "model/vdepsreac.hpp"
#include "model/vessreac.hpp"
#include "model/vessurfsys.hpp"
#include "util/error.hpp"

namespace steps::wm {

Patch::Patch(std::string id, Geom& container, Comp& icomp, Comp* ocomp, double area)
    : pArea(area)
    , pID(std::move(id))
    , pContainer(container) {
    _setIComp(icomp);
    _setOComp(ocomp);

    ArgErrLogIf(pArea < 0.0, "Patch area can't be negative.\n");

    pContainer._handlePatchAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

Patch::~Patch() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setID(std::string const& id) {
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pContainer._handlePatchIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setArea(double area) {
    ArgErrLogIf(area < 0.0, "Patch area can't be negative.\n");

    pArea = area;
}

////////////////////////////////////////////////////////////////////////////////

void Patch::addSurfsys(std::string const& id) {
    // string identifier is only added to set if it is not already included
    pSurfsys.insert(id);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::delSurfsys(std::string const& id) {
    // string identifier is only removed from set if it is included
    pSurfsys.erase(id);
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<model::Spec*> Patch::getAllSpecs(const model::Model& model) const {
    util::flat_set<model::Spec*> specs;

    for (const auto& id: pSurfsys) {
        model::Surfsys& surfsys = model.getSurfsys(id);
        for (auto const& sreac: surfsys.getAllSReacs()) {
            const auto& sreac_specs_slhs = sreac->getSLHS();
            specs.insert(sreac_specs_slhs.begin(), sreac_specs_slhs.end());
            const auto& sreac_specs_srhs = sreac->getSRHS();
            specs.insert(sreac_specs_srhs.begin(), sreac_specs_srhs.end());
        }
        for (auto const& sdiff: surfsys.getAllDiffs()) {
            specs.insert(&sdiff->getLig());
        }
        for (auto const& vdepsreac: surfsys.getAllVDepSReacs()) {
            const auto& vdepsreac_specs_slhs = vdepsreac->getSLHS();
            specs.insert(vdepsreac_specs_slhs.begin(), vdepsreac_specs_slhs.end());
            const auto& vdepsreac_specs_srhs = vdepsreac->getSRHS();
            specs.insert(vdepsreac_specs_srhs.begin(), vdepsreac_specs_srhs.end());
        }
        for (auto const& cplxsreac: surfsys.getAllComplexSReacs()) {
            const auto& cplxsreac_specs_slhs = cplxsreac->getSLHS();
            specs.insert(cplxsreac_specs_slhs.begin(), cplxsreac_specs_slhs.end());
            const auto& cplxsreac_specs_srhs = cplxsreac->getSRHS();
            specs.insert(cplxsreac_specs_srhs.begin(), cplxsreac_specs_srhs.end());
        }
        for (auto const& oc: surfsys.getAllOhmicCurrs()) {
            specs.insert(&oc->getChanState());
        }
        for (auto const& ghk: surfsys.getAllGHKcurrs()) {
            specs.insert(&ghk->getChanState());
        }
        for (auto const& raftgen: surfsys.getAllRaftGens()) {
            const auto& raftgen_specs_sig = raftgen->getSpecSignature();
            specs.insert(raftgen_specs_sig.begin(), raftgen_specs_sig.end());
        }
        for (auto const& endo: surfsys.getAllEndocytosis()) {
            const auto& endo_specs = endo->getSpecDeps();
            specs.insert(endo_specs.begin(), endo_specs.end());
        }
    }

    for (auto const& vessurfsys: model.getAllVesSurfsyss()) {
        for (auto const& vessreac: vessurfsys->getAllVesSReacs()) {
            const auto& vessreac_specs_slhs = vessreac->getSLHS();
            specs.insert(vessreac_specs_slhs.begin(), vessreac_specs_slhs.end());
            const auto& vessreac_specs_srhs = vessreac->getSRHS();
            specs.insert(vessreac_specs_srhs.begin(), vessreac_specs_srhs.end());
        }
    }
    for (auto const& raftsys: model.getAllRaftsyss()) {
        for (auto const& raftsreac: raftsys->getAllRaftSReacs()) {
            const auto& raftsreac_specs_slhs = raftsreac->getSLHS();
            specs.insert(raftsreac_specs_slhs.begin(), raftsreac_specs_slhs.end());
            const auto& raftsreac_specs_srhs = raftsreac->getSRHS();
            specs.insert(raftsreac_specs_srhs.begin(), raftsreac_specs_srhs.end());
        }
    }

    return specs;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<model::SReac*> Patch::getAllSReacs(const model::Model& model) const {
    util::flat_set<model::SReac*> pSReacs;
    for (const auto& id: pSurfsys) {
        model::Surfsys& surfsys = model.getSurfsys(id);
        const auto& sreacs = surfsys.getAllSReacs();
        pSReacs.insert(sreacs.begin(), sreacs.end());
    }

    return pSReacs;
}

////////////////////////////////////////////////////////////////////////////////

void Patch::_setIComp(Comp& icomp) {
    ArgErrLogIf(&icomp.getContainer() != &pContainer,
                "Compartment does not belong to same container as patch.\n");

    auto const& ipatches = icomp.getIPatches();

    ArgErrLogIf(ipatches.find(this) != ipatches.end(),
                "Patch is already on inside of compartment.\n");

    // remove the patch if it was already on the outside of some
    // other compartment
    if (pIComp != nullptr) {
        pIComp->_delOPatch(*this);
    }

    pIComp = &icomp;
    pIComp->_addOPatch(*this);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::_setOComp(Comp* ocomp) {
    if (ocomp == nullptr) {
        return;
    }

    ArgErrLogIf(&ocomp->getContainer() != &pContainer,
                "Compartment does not belong to same container as patch.\n");

    auto const& opatches = ocomp->getOPatches();

    ArgErrLogIf(opatches.find(this) != opatches.end(),
                "Patch is already on outside of compartment.\n");

    // remove the patch if it was already on the inside of some
    // other compartment
    if (pOComp != nullptr) {
        pOComp->_delIPatch(*this);
    }

    pOComp = ocomp;
    pOComp->_addIPatch(*this);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::_handleSelfDelete() {
    pContainer._handlePatchDel(*this);
    pArea = 0.0;
    pIComp->_delOPatch(*this);
    if (pOComp != nullptr) {
        pOComp->_delIPatch(*this);
    }
}

}  // namespace steps::wm
