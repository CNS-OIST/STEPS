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

#include "model/complexsreac.hpp"

#include "model/model.hpp"
#include "model/spec.hpp"
#include "model/surfsys.hpp"
#include "util/error.hpp"

namespace steps::model {

ComplexSReac::ComplexSReac(std::string const& id,
                           Surfsys& surfsys,
                           std::vector<Spec*> const& ilhs,
                           std::vector<Spec*> const& slhs,
                           std::vector<Spec*> const& olhs,
                           std::vector<Spec*> const& irhs,
                           std::vector<Spec*> const& srhs,
                           std::vector<Spec*> const& orhs,
                           std::vector<ComplexEvent*> const& icompEvs,
                           std::vector<ComplexEvent*> const& scompEvs,
                           std::vector<ComplexEvent*> const& ocompEvs,
                           double kcst)
    : pID(id)
    , pModel(surfsys.getModel())
    , pSurfsys(surfsys)
    , pILHS(ilhs)
    , pSLHS(slhs)
    , pOLHS(olhs)
    , pIRHS(irhs)
    , pSRHS(srhs)
    , pORHS(orhs) {
    setKcst(kcst);

    for (auto const& l: pOLHS) {
        AssertLog(&l->getModel() == &pModel);
    }
    for (auto const& l: pILHS) {
        AssertLog(&l->getModel() == &pModel);
    }
    for (auto const& l: pSLHS) {
        AssertLog(&l->getModel() == &pModel);
    }
    for (auto const& l: pORHS) {
        AssertLog(&l->getModel() == &pModel);
    }
    for (auto const& l: pIRHS) {
        AssertLog(&l->getModel() == &pModel);
    }
    for (auto const& l: pSRHS) {
        AssertLog(&l->getModel() == &pModel);
    }

    pLocOrder[ComplexLocation::PATCH_IN] = pILHS.size();
    pLocOrder[ComplexLocation::PATCH_SURF] = pSLHS.size();
    pLocOrder[ComplexLocation::PATCH_OUT] = pOLHS.size();

    for (auto* ev: icompEvs) {
        _addEvent(ev, ComplexLocation::PATCH_IN);
    }
    for (auto* ev: scompEvs) {
        _addEvent(ev, ComplexLocation::PATCH_SURF);
    }
    for (auto* ev: ocompEvs) {
        _addEvent(ev, ComplexLocation::PATCH_OUT);
    }

    pOuter = pLocOrder[ComplexLocation::PATCH_OUT] > 0;

    ArgErrLogIf(pOuter and pLocOrder[ComplexLocation::PATCH_IN] > 0,
                "Surface reaction cannot contain reactants on both sides of the patch.");

    for (auto ord: pLocOrder) {
        pOrder += ord.second;
    }

    pSurfSurf = pOrder == pLocOrder[ComplexLocation::PATCH_SURF];

    pSurfsys._handleComplexSReacAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

void ComplexSReac::_addEvent(ComplexEvent* ev, ComplexLocation loc) {
    if (dynamic_cast<ComplexUpdateEvent*>(ev) != nullptr) {
        pCompUPD[loc].push_back(dynamic_cast<ComplexUpdateEvent*>(ev));
        pLocOrder[loc]++;
    } else if (dynamic_cast<ComplexDeleteEvent*>(ev) != nullptr) {
        pCompDEL[loc].push_back(dynamic_cast<ComplexDeleteEvent*>(ev));
        pLocOrder[loc]++;
    } else if (dynamic_cast<ComplexCreateEvent*>(ev) != nullptr) {
        pCompCRE[loc].push_back(dynamic_cast<ComplexCreateEvent*>(ev));
    }
}

////////////////////////////////////////////////////////////////////////////////

static std::vector<ComplexUpdateEvent*> empty_upd;

const std::vector<ComplexUpdateEvent*>& ComplexSReac::getUPDEvents(
    ComplexLocation loc) const noexcept {
    const auto it = pCompUPD.find(loc);
    if (it != pCompUPD.end()) {
        return it->second;
    }
    return empty_upd;
}

////////////////////////////////////////////////////////////////////////////////

static std::vector<ComplexDeleteEvent*> empty_del;

const std::vector<ComplexDeleteEvent*>& ComplexSReac::getDELEvents(
    ComplexLocation loc) const noexcept {
    const auto it = pCompDEL.find(loc);
    if (it != pCompDEL.end()) {
        return it->second;
    }
    return empty_del;
}

////////////////////////////////////////////////////////////////////////////////

static std::vector<ComplexCreateEvent*> empty_cre;

const std::vector<ComplexCreateEvent*>& ComplexSReac::getCREEvents(
    ComplexLocation loc) const noexcept {
    const auto it = pCompCRE.find(loc);
    if (it != pCompCRE.end()) {
        return it->second;
    }
    return empty_cre;
}

////////////////////////////////////////////////////////////////////////////////

void ComplexSReac::setKcst(double kcst) {
    ArgErrLogIf(kcst < 0.0, "Surface reaction constant can't be negative");
    pKcst = kcst;
}
////////////////////////////////////////////////////////////////////////////////

util::flat_set<Spec*> ComplexSReac::getAllSpecs() const {
    util::flat_set<Spec*> specSet;
    specSet.insert(pOLHS.begin(), pOLHS.end());
    specSet.insert(pILHS.begin(), pILHS.end());
    specSet.insert(pSLHS.begin(), pSLHS.end());
    specSet.insert(pORHS.begin(), pORHS.end());
    specSet.insert(pIRHS.begin(), pIRHS.end());
    specSet.insert(pSRHS.begin(), pSRHS.end());
    return specSet;
}

////////////////////////////////////////////////////////////////////////////////


}  // namespace steps::model
