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

// STEPS headers.
#include "diffdef.hpp"
#include "model/diff.hpp"
#include "model/spec.hpp"
#include "specdef.hpp"
#include "types.hpp"
#include "util/checkpointing.hpp"
#include "util/error.hpp"
// logging
#include <easylogging++.h>

namespace steps::solver {

////////////////////////////////////////////////////////////////////////////////

template <typename GlobalId>
MetaDiffdef<GlobalId>::MetaDiffdef(Statedef* sd, GlobalId idx, model::Diff* d)
    : pStatedef(sd)
    , pIdx(idx)
    , pDcst()
    , pSetupdone(false) {
    AssertLog(pStatedef != nullptr);
    AssertLog(d != nullptr);

    pName = d->getID();
    pDcst = d->getDcst();
    pLig = d->getLig()->getID();
    ligGIdx = pStatedef->getSpecIdx(pLig);

    const uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) {
        return;
    }
    pSpec_DEP.resize(nspecs, DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

template <typename GlobalId>
void MetaDiffdef<GlobalId>::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pDcst);
}

////////////////////////////////////////////////////////////////////////////////

template <typename GlobalId>

void MetaDiffdef<GlobalId>::restore(std::fstream& cp_file) {
    util::restore(cp_file, pDcst);
}

////////////////////////////////////////////////////////////////////////////////


template <typename GlobalId>
void MetaDiffdef<GlobalId>::setup() {
    AssertLog(pSetupdone == false);

    pSpec_DEP[lig().get()] = DEP_STOICH;

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

template <typename GlobalId>

void MetaDiffdef<GlobalId>::setDcst(double d) {
    AssertLog(d >= 0.0);
    pDcst = d;
}

////////////////////////////////////////////////////////////////////////////////

template <typename GlobalId>
spec_global_id MetaDiffdef<GlobalId>::lig() const {
    AssertLog(pStatedef != nullptr);
    return ligGIdx;
}

////////////////////////////////////////////////////////////////////////////////

template <typename GlobalId>
int MetaDiffdef<GlobalId>::dep(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_DEP.at(gidx.get());
}

////////////////////////////////////////////////////////////////////////////////

template <typename GlobalId>
bool MetaDiffdef<GlobalId>::reqspec(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_DEP.at(gidx.get()) != DEP_NONE;
}

// explicit template instantiation definitions
template class MetaDiffdef<diff_global_id>;
template class MetaDiffdef<surfdiff_global_id>;

}  // namespace steps::solver
