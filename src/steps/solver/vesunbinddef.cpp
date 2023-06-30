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

#include "solver/vesunbinddef.hpp"

// STEPS headers.
#include "solver/compdef.hpp"
#include "solver/statedef.hpp"
#include "solver/types.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

VesUnbinddef::VesUnbinddef(Statedef* sd, vesunbind_global_id idx, model::VesUnbind* vb)
    : pStatedef(sd)
    , pIdx(idx)
    , pKcst()
    , pImmobility()
    , pSetupdone(false) {
    AssertLog(pStatedef != nullptr);
    AssertLog(vb != nullptr);

    pName = vb->getID();
    pKcst = vb->getKcst();

    pImmobility = vb->getImmobilization();

    pProducts1 = vb->getProducts1();
    pProducts2 = vb->getProducts2();

    pLinks1 = vb->getLinks1();
    pLinks2 = vb->getLinks2();
}

////////////////////////////////////////////////////////////////////////////////

VesUnbinddef::~VesUnbinddef() = default;

////////////////////////////////////////////////////////////////////////////////

void VesUnbinddef::checkpoint(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void VesUnbinddef::restore(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void VesUnbinddef::setup() {
    AssertLog(pSetupdone == false);

    // Make sure the vesicles match up
    AssertLog(pLinks1.first == pProducts1.first);
    AssertLog(pLinks2.first == pProducts2.first);

    pVesicle_1_idx = pStatedef->getVesicleIdx(pProducts1.first);
    pVesicle_2_idx = pStatedef->getVesicleIdx(pProducts2.first);

    pSpec_1_gidx = pStatedef->getSpecIdx(pProducts1.second);
    pSpec_2_gidx = pStatedef->getSpecIdx(pProducts2.second);

    pLinkSpec_1_gidx = pStatedef->getLinkSpecIdx(pLinks1.second);
    pLinkSpec_2_gidx = pStatedef->getLinkSpecIdx(pLinks2.second);

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

vesicle_global_id VesUnbinddef::getVes1idx() const {
    AssertLog(pSetupdone == true);
    return pVesicle_1_idx;
}

////////////////////////////////////////////////////////////////////////////////

vesicle_global_id VesUnbinddef::getVes2idx() const {
    AssertLog(pSetupdone == true);
    return pVesicle_2_idx;
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id VesUnbinddef::getSpec1gidx() const {
    AssertLog(pSetupdone == true);
    return pSpec_1_gidx;
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id VesUnbinddef::getSpec2gidx() const {
    AssertLog(pSetupdone == true);
    return pSpec_2_gidx;
}

////////////////////////////////////////////////////////////////////////////////

linkspec_global_id VesUnbinddef::getLinkSpec1gidx() const {
    AssertLog(pSetupdone == true);
    return pLinkSpec_1_gidx;
}

////////////////////////////////////////////////////////////////////////////////

linkspec_global_id VesUnbinddef::getLinkSpec2gidx() const {
    AssertLog(pSetupdone == true);
    return pLinkSpec_2_gidx;
}

}  // namespace steps::solver
