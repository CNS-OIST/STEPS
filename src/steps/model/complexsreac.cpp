/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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

ComplexSReacBase::ComplexSReacBase(std::string const& id,
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
                                   int charge)
    : pID(id)
    , pModel(surfsys.getModel())
    , pSurfsys(surfsys)
    , pILHS(ilhs)
    , pSLHS(slhs)
    , pOLHS(olhs)
    , pIRHS(irhs)
    , pSRHS(srhs)
    , pORHS(orhs)
    , pCharge(charge) {
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

    pLocOrder[Location::PATCH_IN] = pILHS.size();
    pLocOrder[Location::PATCH_SURF] = pSLHS.size();
    pLocOrder[Location::PATCH_OUT] = pOLHS.size();

    for (auto* ev: icompEvs) {
        _addEvent(ev, Location::PATCH_IN);
    }
    for (auto* ev: scompEvs) {
        _addEvent(ev, Location::PATCH_SURF);
    }
    for (auto* ev: ocompEvs) {
        _addEvent(ev, Location::PATCH_OUT);
    }

    pOuter = pLocOrder[Location::PATCH_OUT] > 0;

    ArgErrLogIf(pOuter and pLocOrder[Location::PATCH_IN] > 0,
                "Surface reaction cannot contain reactants on both sides of the patch.");

    for (auto ord: pLocOrder) {
        pOrder += ord.second;
    }

    pSurfSurf = pOrder == pLocOrder[Location::PATCH_SURF];
}

////////////////////////////////////////////////////////////////////////////////

void ComplexSReacBase::_addEvent(ComplexEvent* ev, Location loc) {
    if (auto* upd = dynamic_cast<ComplexUpdateEvent*>(ev)) {
        pCompUPD[loc].push_back(upd);
        pLocOrder[loc]++;
    } else if (auto* del = dynamic_cast<ComplexDeleteEvent*>(ev)) {
        pCompDEL[loc].push_back(del);
        pLocOrder[loc]++;
    } else if (auto* cre = dynamic_cast<ComplexCreateEvent*>(ev)) {
        pCompCRE[loc].push_back(cre);
    }
}

////////////////////////////////////////////////////////////////////////////////

static std::vector<ComplexUpdateEvent*> empty_upd;

const std::vector<ComplexUpdateEvent*>& ComplexSReacBase::getUPDEvents(
    Location loc) const noexcept {
    const auto it = pCompUPD.find(loc);
    if (it != pCompUPD.end()) {
        return it->second;
    }
    return empty_upd;
}

////////////////////////////////////////////////////////////////////////////////

static std::vector<ComplexDeleteEvent*> empty_del;

const std::vector<ComplexDeleteEvent*>& ComplexSReacBase::getDELEvents(
    Location loc) const noexcept {
    const auto it = pCompDEL.find(loc);
    if (it != pCompDEL.end()) {
        return it->second;
    }
    return empty_del;
}

////////////////////////////////////////////////////////////////////////////////

static std::vector<ComplexCreateEvent*> empty_cre;

const std::vector<ComplexCreateEvent*>& ComplexSReacBase::getCREEvents(
    Location loc) const noexcept {
    const auto it = pCompCRE.find(loc);
    if (it != pCompCRE.end()) {
        return it->second;
    }
    return empty_cre;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<Spec*> ComplexSReacBase::getAllSpecs() const {
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

void ComplexSReacBase::setCharge(int charge) {
    pCharge = charge;
}

////////////////////////////////////////////////////////////////////////////////

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
                           double kcst,
                           int charge)
    : ComplexSReacBase::ComplexSReacBase(id,
                                         surfsys,
                                         ilhs,
                                         slhs,
                                         olhs,
                                         irhs,
                                         srhs,
                                         orhs,
                                         icompEvs,
                                         scompEvs,
                                         ocompEvs,
                                         charge) {
    setKcst(kcst);

    pSurfsys._handleComplexSReacAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

void ComplexSReac::setKcst(double kcst) {
    ArgErrLogIf(kcst < 0.0, "Surface reaction constant can't be negative");
    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

VDepComplexSReac::VDepComplexSReac(std::string const& id,
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
                                   std::vector<double> ktab,
                                   double vmin,
                                   double vmax,
                                   double dv,
                                   uint tablesize,
                                   int charge)
    : ComplexSReacBase::ComplexSReacBase(id,
                                         surfsys,
                                         ilhs,
                                         slhs,
                                         olhs,
                                         irhs,
                                         srhs,
                                         orhs,
                                         icompEvs,
                                         scompEvs,
                                         ocompEvs,
                                         charge)
    , pVMin(vmin)
    , pVMax(vmax)
    , pDV(dv)
    , pTablesize(tablesize) {
    ArgErrLogIf(ktab.size() != pTablesize, "Table of reaction parameters is not of expected size");
    AssertLog(pDV > 0.0);
    pK = ktab;

    pSurfsys._handleVDepComplexSReacAdd(*this);
}

}  // namespace steps::model
