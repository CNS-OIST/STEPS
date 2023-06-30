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

#include "solver/vesbinddef.hpp"

#include "solver/compdef.hpp"
#include "solver/linkspecdef.hpp"
#include "solver/statedef.hpp"
#include "solver/types.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

VesBinddef::VesBinddef(Statedef* sd, vesbind_global_id idx, model::VesBind* vb)
    : pStatedef(sd)
    , pIdx(idx)
    , pKcst()
    , pMaxDistance()
    , pImmobility()
    , pSetupdone(false)
    , pGotVDep1(true)
    , pGotVDep2(true)
    , pGotLDep1(true)
    , pGotLDep2(true) {
    AssertLog(pStatedef != nullptr);
    AssertLog(vb != nullptr);

    pName = vb->getID();
    pKcst = vb->getKcst();

    pReactants1 = vb->getReactants1();
    pReactants2 = vb->getReactants2();
    pProducts1 = vb->getProducts1();
    pProducts2 = vb->getProducts2();

    pMaxDistance = vb->getLengthMax();
    pMinDistance = vb->getLengthMin();

    pVdep1 = vb->getVDeps1();
    pVdep2 = vb->getVDeps2();

    pLdep1 = vb->getLDeps1();
    pLdep2 = vb->getLDeps2();

    pImmobility = vb->getImmobilization();

    if (pVdep1.empty()) {
        pGotVDep1 = false;
    }
    if (pVdep2.empty()) {
        pGotVDep2 = false;
    }
    if (pLdep1.empty()) {
        pGotLDep1 = false;
    }
    if (pLdep2.empty()) {
        pGotLDep2 = false;
    }

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) {
        return;  // Would be weird, but okay.
    }

    pSpec_VDEP1.container().resize(nspecs);
    pSpec_VDEP2.container().resize(nspecs);

    uint nlspecs = pStatedef->countLinkSpecs();

    pSpec_LDEP1.container().resize(nlspecs);
    pSpec_LDEP2.container().resize(nlspecs);
}

////////////////////////////////////////////////////////////////////////////////

void VesBinddef::checkpoint(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void VesBinddef::restore(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void VesBinddef::setup() {
    AssertLog(pSetupdone == false);

    // The vesicles '1' and '2' should match
    AssertLog(pReactants1.first == pProducts1.first);
    AssertLog(pReactants2.first == pProducts2.first);

    pVesicle_1_idx = pStatedef->getVesicleIdx(pReactants1.first);
    pSpec_1_gidx = pStatedef->getSpecIdx(pReactants1.second);

    pVesicle_2_idx = pStatedef->getVesicleIdx(pReactants2.first);
    pSpec_2_gidx = pStatedef->getSpecIdx(pReactants2.second);

    pProduct1_gidx = pStatedef->getLinkSpecIdx(pProducts1.second);
    pProduct2_gidx = pStatedef->getLinkSpecIdx(pProducts2.second);

    pProduct1def = pStatedef->linkspecdef(pProduct1_gidx);
    pProduct2def = pStatedef->linkspecdef(pProduct2_gidx);

    if (pGotVDep1) {
        for (auto const& vdep1: pVdep1) {
            spec_global_id sidx = pStatedef->getSpecIdx(vdep1);
            pSpec_VDEP1[sidx] += 1;
        }
    }

    if (pGotVDep2) {
        for (auto const& vdep2: pVdep2) {
            spec_global_id sidx = pStatedef->getSpecIdx(vdep2);
            pSpec_VDEP2[sidx] += 1;
        }
    }

    if (pGotLDep1) {
        for (auto const& ldep1: pLdep1) {
            linkspec_global_id lsidx = pStatedef->getLinkSpecIdx(ldep1);
            pSpec_LDEP1[lsidx] += 1;
        }
    }

    if (pGotLDep2) {
        for (auto const& ldep2: pLdep2) {
            linkspec_global_id lsidx = pStatedef->getLinkSpecIdx(ldep2);
            pSpec_LDEP2[lsidx] += 1;
        }
    }
    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

const std::string& VesBinddef::name() const {
    return pName;
}

////////////////////////////////////////////////////////////////////////////////

double VesBinddef::kcst() const {
    return pKcst;
}

////////////////////////////////////////////////////////////////////////////////

int VesBinddef::immobility() const {
    return pImmobility;
}

////////////////////////////////////////////////////////////////////////////////

vesicle_global_id VesBinddef::getVes1idx() const {
    AssertLog(pSetupdone == true);
    return pVesicle_1_idx;
}

////////////////////////////////////////////////////////////////////////////////

vesicle_global_id VesBinddef::getVes2idx() const {
    AssertLog(pSetupdone == true);
    return pVesicle_2_idx;
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id VesBinddef::getSpec1gidx() const {
    AssertLog(pSetupdone == true);
    return pSpec_1_gidx;
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id VesBinddef::getSpec2gidx() const {
    AssertLog(pSetupdone == true);
    return pSpec_2_gidx;
}

////////////////////////////////////////////////////////////////////////////////

linkspec_global_id VesBinddef::getLinkSpec1gidx() const {
    AssertLog(pSetupdone == true);
    return pProduct1_gidx;
}

////////////////////////////////////////////////////////////////////////////////

linkspec_global_id VesBinddef::getLinkSpec2gidx() const {
    AssertLog(pSetupdone == true);
    return pProduct2_gidx;
}

////////////////////////////////////////////////////////////////////////////////

LinkSpecdef* VesBinddef::getLinkSpec1def() const {
    AssertLog(pSetupdone == true);
    return pProduct1def;
}

////////////////////////////////////////////////////////////////////////////////

LinkSpecdef* VesBinddef::getLinkSpec2def() const {
    AssertLog(pSetupdone == true);
    return pProduct2def;
}

////////////////////////////////////////////////////////////////////////////////

uint VesBinddef::vdep1(spec_global_id gidx) const {
    return pSpec_VDEP1.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint VesBinddef::vdep2(spec_global_id gidx) const {
    return pSpec_VDEP2.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint VesBinddef::ldep1(linkspec_global_id gidx) const {
    return pSpec_LDEP1.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint VesBinddef::ldep2(linkspec_global_id gidx) const {
    return pSpec_LDEP2.at(gidx);
}

}  // namespace steps::solver
