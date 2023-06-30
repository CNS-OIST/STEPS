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

#include "solver/vessdiffdef.hpp"

// STEPS headers.
#include "model/spec.hpp"
#include "solver/fwd.hpp"
#include "solver/specdef.hpp"
#include "solver/statedef.hpp"
#include "solver/types.hpp"
#include "util/checkpointing.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

VesSDiffdef::VesSDiffdef(Statedef* sd, vessdiff_global_id idx, model::VesSDiff* d)
    : pStatedef(sd)
    , pIdx(idx)
    , pDcst() {
    AssertLog(pStatedef != nullptr);
    AssertLog(d != nullptr);

    pName = d->getID();
    pDcst = d->getDcst();
    pLig = d->getLig()->getID();
}

////////////////////////////////////////////////////////////////////////////////

VesSDiffdef::~VesSDiffdef() = default;

////////////////////////////////////////////////////////////////////////////////

void VesSDiffdef::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pDcst);
}

////////////////////////////////////////////////////////////////////////////////

void VesSDiffdef::restore(std::fstream& cp_file) {
    util::restore(cp_file, pDcst);
}

////////////////////////////////////////////////////////////////////////////////

void VesSDiffdef::setup() {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void VesSDiffdef::setDcst(double d) {
    AssertLog(d >= 0.0);
    pDcst = d;
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id VesSDiffdef::lig() const {
    AssertLog(pStatedef != nullptr);
    return pStatedef->getSpecIdx(pLig);
}

////////////////////////////////////////////////////////////////////////////////

int VesSDiffdef::dep(spec_global_id gidx) const {
    if (gidx == lig()) {
        return 1;
    } else {
        return 0;
    }
}

////////////////////////////////////////////////////////////////////////////////

bool VesSDiffdef::reqspec(spec_global_id gidx) const {
    if (gidx == lig()) {
        return true;
    }
    return false;
}

}  // namespace steps::solver
