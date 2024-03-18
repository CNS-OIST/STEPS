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

#include "complexsreacdef.hpp"

#include "geom/patch.hpp"
#include "patchdef.hpp"
#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

ComplexSReacdef::ComplexSReacdef(Statedef& sd, complexsreac_global_id idx, model::ComplexSReac& sr)
    : pStatedef(sd)
    , pIdx(idx)
    , pName(sr.getID())
    , pOrder(sr.getOrder())
    , pKcst(sr.getKcst())
    , pIlhs(sr.getILHS())
    , pOlhs(sr.getOLHS())
    , pSlhs(sr.getSLHS())
    , pIrhs(sr.getIRHS())
    , pOrhs(sr.getORHS())
    , pSrhs(sr.getSRHS())
    , pSurface_surface(sr.getSurfSurf()) {
    if (pOrder == 0) {
        std::ostringstream os;
        os << "\nModel contains zero-order surface reaction, which are not "
              "permitted. ";
        os << " Zero-order volume reaction may be used instead\n.";
        ArgErrLog(os.str());
    }

    if (sr.getInner()) {
        pOrient = SReacdef::INSIDE;
    } else {
        pOrient = SReacdef::OUTSIDE;
    }

    for (auto loc: model::AllPatchLocations) {
        for (auto* ev: sr.getUPDEvents(loc)) {
            pComplexUPDEvs[loc].push_back(std::make_shared<ComplexUpdateEventdef>(*ev, sd));
        }
        for (auto* ev: sr.getDELEvents(loc)) {
            pComplexDELEvs[loc].push_back(std::make_shared<ComplexDeleteEventdef>(*ev, sd));
        }
        for (auto* ev: sr.getCREEvents(loc)) {
            pComplexCREEvs[loc].push_back(std::make_shared<ComplexCreateEventdef>(*ev, sd));
        }
    }

    uint nspecs = pStatedef.countSpecs();

    pSpec_S_DEP.container().resize(nspecs, DEP_NONE);
    pSpec_S_LHS.container().resize(nspecs);
    if (pOrient == SReacdef::INSIDE) {
        pSpec_I_DEP.container().resize(nspecs, DEP_NONE);
        pSpec_I_LHS.container().resize(nspecs);
    } else {
        pSpec_O_DEP.container().resize(nspecs, DEP_NONE);
        pSpec_O_LHS.container().resize(nspecs);
    }
    pSpec_I_RHS.container().resize(nspecs);
    pSpec_S_RHS.container().resize(nspecs);
    pSpec_O_RHS.container().resize(nspecs);
    pSpec_I_UPD.container().resize(nspecs);
    pSpec_S_UPD.container().resize(nspecs);
    pSpec_O_UPD.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void ComplexSReacdef::checkpoint(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void ComplexSReacdef::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void ComplexSReacdef::setup() {
    AssertLog(pSetupdone == false);

    if (outside()) {
        AssertLog(pIlhs.empty());
    } else if (inside()) {
        AssertLog(pOlhs.empty());
    } else {
        AssertLog(false);
    }

    for (auto const& ol: pOlhs) {
        spec_global_id sidx = pStatedef.getSpecIdx(*ol);
        pSpec_O_LHS[sidx] += 1;
    }

    for (auto const& il: pIlhs) {
        spec_global_id sidx = pStatedef.getSpecIdx(*il);
        pSpec_I_LHS[sidx] += 1;
    }

    for (auto const& sl: pSlhs) {
        spec_global_id sidx = pStatedef.getSpecIdx(*sl);
        pSpec_S_LHS[sidx] += 1;
    }

    for (auto const& ir: pIrhs) {
        spec_global_id sidx = pStatedef.getSpecIdx(*ir);
        pSpec_I_RHS[sidx] += 1;
    }

    for (auto const& sr: pSrhs) {
        spec_global_id sidx = pStatedef.getSpecIdx(*sr);
        pSpec_S_RHS[sidx] += 1;
    }

    for (auto const& orh: pOrhs) {
        spec_global_id sidx = pStatedef.getSpecIdx(*orh);
        pSpec_O_RHS[sidx] += 1;
    }

    // Now set up the update vector
    uint ngspecs = pStatedef.countSpecs();

    // Deal with surface.
    for (auto s: spec_global_id::range(ngspecs)) {
        const auto lhs = static_cast<int>(pSpec_S_LHS[s]);
        const auto rhs = static_cast<int>(pSpec_S_RHS[s]);
        const int aux = pSpec_S_UPD[s] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_S_DEP[s] |= DEP_STOICH;
        }
        if (aux != 0) {
            pSpec_S_UPD_Coll.push_back(s);
        }
    }

    // Deal with inside.
    for (auto s: spec_global_id::range(ngspecs)) {
        const int lhs = (inside() ? static_cast<int>(pSpec_I_LHS[s]) : 0);
        const auto rhs = static_cast<int>(pSpec_I_RHS[s]);
        const int aux = pSpec_I_UPD[s] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_I_DEP[s] |= DEP_STOICH;
        }
        if (aux != 0) {
            pSpec_I_UPD_Coll.push_back(s);
        }
    }

    // Deal with outside.
    for (auto s: spec_global_id::range(ngspecs)) {
        const int lhs = (outside() ? static_cast<int>(pSpec_O_LHS[s]) : 0);
        const auto rhs = static_cast<int>(pSpec_O_RHS[s]);
        const int aux = pSpec_O_UPD[s] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_O_DEP[s] |= DEP_STOICH;
        }
        if (aux != 0) {
            pSpec_O_UPD_Coll.push_back(s);
        }
    }

    // set up deps for complexes
    for (auto loc: model::AllPatchLocations) {
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

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexSReacdef::reqInside() const {
    AssertLog(pSetupdone);

    // This can be checked by seeing if DEP_I or RHS_I is non-zero
    // for any species.
    for (auto s: spec_global_id::range(pStatedef.countSpecs())) {
        if (reqspec_I(s)) {
            return true;
        }
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexSReacdef::reqOutside() const {
    AssertLog(pSetupdone);

    // This can be checked by seeing if DEP_O or RHS_O is non-zero
    // for any species.
    for (auto s: spec_global_id::range(pStatedef.countSpecs())) {
        if (reqspec_O(s)) {
            return true;
        }
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

uint ComplexSReacdef::lhs_I(spec_global_id gidx) const {
    if (outside()) {
        return 0;
    }
    return pSpec_I_LHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint ComplexSReacdef::lhs_S(spec_global_id gidx) const {
    return pSpec_S_LHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint ComplexSReacdef::lhs_O(spec_global_id gidx) const {
    if (inside()) {
        return 0;
    }
    return pSpec_O_LHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int ComplexSReacdef::dep_I(spec_global_id gidx) const {
    AssertLog(pSetupdone);
    if (outside()) {
        return DEP_NONE;
    }
    return pSpec_I_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int ComplexSReacdef::dep_S(spec_global_id gidx) const {
    AssertLog(pSetupdone);
    return pSpec_S_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int ComplexSReacdef::dep_O(spec_global_id gidx) const {
    AssertLog(pSetupdone);
    if (inside()) {
        return DEP_NONE;
    }
    return pSpec_O_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexSReacdef::complexdep(model::ComplexLocation loc,
                                 complex_global_id gidx,
                                 complex_substate_id sus) const {
    AssertLog(pSetupdone);
    AssertLog(gidx < pStatedef.countComplexes());
    const auto it = pComplex_DEPMAP.find(loc);
    if (it != pComplex_DEPMAP.end()) {
        const auto it2 = it->second.find(gidx);
        return it2 != it->second.end() and it2->second.find(sus) != it2->second.end();
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

uint ComplexSReacdef::rhs_I(spec_global_id gidx) const {
    return pSpec_I_RHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint ComplexSReacdef::rhs_S(spec_global_id gidx) const {
    return pSpec_S_RHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint ComplexSReacdef::rhs_O(spec_global_id gidx) const {
    return pSpec_O_RHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int ComplexSReacdef::upd_I(spec_global_id gidx) const {
    AssertLog(pSetupdone);
    return pSpec_I_UPD.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int ComplexSReacdef::upd_S(spec_global_id gidx) const {
    AssertLog(pSetupdone);
    return pSpec_S_UPD.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int ComplexSReacdef::upd_O(spec_global_id gidx) const {
    AssertLog(pSetupdone);
    return pSpec_O_UPD.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexSReacdef::reqspec_I(spec_global_id gidx) const {
    AssertLog(pSetupdone);
    if (inside()) {
        if (pSpec_I_DEP.at(gidx) != DEP_NONE) {
            return true;
        }
    }
    if (pSpec_I_RHS.at(gidx) != 0) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexSReacdef::reqspec_S(spec_global_id gidx) const {
    AssertLog(pSetupdone);
    if (pSpec_S_DEP.at(gidx) != DEP_NONE) {
        return true;
    }
    if (pSpec_S_RHS.at(gidx) != 0) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexSReacdef::reqspec_O(spec_global_id gidx) const {
    AssertLog(pSetupdone);
    if (outside()) {
        if (pSpec_O_DEP.at(gidx) != DEP_NONE) {
            return true;
        }
    }
    if (pSpec_O_RHS.at(gidx) != 0) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

static std::map<complex_global_id, std::set<complex_substate_id>> empty_map;
const std::map<complex_global_id, std::set<complex_substate_id>>& ComplexSReacdef::complexUPDMAP(
    model::ComplexLocation loc) const {
    const auto it = pComplex_UPDMAP.find(loc);
    if (it != pComplex_UPDMAP.end()) {
        return it->second;
    } else {
        return empty_map;
    }
}

////////////////////////////////////////////////////////////////////////////////

static std::vector<std::shared_ptr<ComplexUpdateEventdef>> empty_upd;
const std::vector<std::shared_ptr<ComplexUpdateEventdef>>& ComplexSReacdef::updEvents(
    model::ComplexLocation loc) const {
    const auto it = pComplexUPDEvs.find(loc);
    if (it != pComplexUPDEvs.end()) {
        return it->second;
    } else {
        return empty_upd;
    }
}

////////////////////////////////////////////////////////////////////////////////

static std::vector<std::shared_ptr<ComplexDeleteEventdef>> empty_del;
const std::vector<std::shared_ptr<ComplexDeleteEventdef>>& ComplexSReacdef::delEvents(
    model::ComplexLocation loc) const {
    const auto it = pComplexDELEvs.find(loc);
    if (it != pComplexDELEvs.end()) {
        return it->second;
    } else {
        return empty_del;
    }
}

////////////////////////////////////////////////////////////////////////////////

static std::vector<std::shared_ptr<ComplexCreateEventdef>> empty_cre;
const std::vector<std::shared_ptr<ComplexCreateEventdef>>& ComplexSReacdef::creEvents(
    model::ComplexLocation loc) const {
    const auto it = pComplexCREEvs.find(loc);
    if (it != pComplexCREEvs.end()) {
        return it->second;
    } else {
        return empty_cre;
    }
}

////////////////////////////////////////////////////////////////////////////////


}  // namespace steps::solver
