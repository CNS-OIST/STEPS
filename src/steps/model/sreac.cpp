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

#include "sreac.hpp"

#include "spec.hpp"
#include "surfsys.hpp"
#include "util/error.hpp"

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

SReac::SReac(std::string const& id,
             Surfsys& surfsys,
             std::vector<Spec*> const& olhs,
             std::vector<Spec*> const& ilhs,
             std::vector<Spec*> const& slhs,
             std::vector<Spec*> const& irhs,
             std::vector<Spec*> const& srhs,
             std::vector<Spec*> const& orhs,
             double kcst)
    : pID(id)
    , pModel(surfsys.getModel())
    , pSurfsys(surfsys)
    , pOuter(false)
    , pOrder(0)
    , pKcst(kcst) {
    ArgErrLogIf(pKcst < 0.0, "Surface reaction constant can't be negative");

    // Can't have species on the lhs in the inner and outer compartment
    ArgErrLogIf(!olhs.empty() && !ilhs.empty(),
                "Volume lhs species must belong to either inner or outer "
                "compartment, not both.");

    if (!olhs.empty()) {
        setOLHS(olhs);
    }
    if (!ilhs.empty()) {
        setILHS(ilhs);
    }
    setSLHS(slhs);
    setIRHS(irhs);
    setSRHS(srhs);
    setORHS(orhs);

    pSurfsys._handleSReacAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

SReac::~SReac() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void SReac::_handleSelfDelete() {
    pSurfsys._handleSReacDel(*this);
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setID(std::string const& id) {
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pSurfsys._handleSReacIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setOLHS(std::vector<Spec*> const& olhs) {
    for (auto const& ol: olhs) {
        AssertLog(&ol->getModel() == &pModel);
    }

    if (!pILHS.empty()) {
        CLOG(WARNING, "general_log") << "WARNING: Removing inner compartment species from lhs "
                                        "stoichiometry for SReac "
                                     << getID();
    }
    pILHS.clear();
    pOLHS = olhs;

    pOuter = true;
    pOrder = pOLHS.size() + pSLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setILHS(std::vector<Spec*> const& ilhs) {
    for (auto const& il: ilhs) {
        AssertLog(&il->getModel() == &pModel);
    }

    if (!pOLHS.empty()) {
        CLOG(WARNING, "general_log") << "\nWARNING: Removing outer compartment species from lhs "
                                        "stoichiometry for SReac "
                                     << getID();
    }

    pOLHS.clear();
    pILHS = ilhs;

    pOuter = false;
    pOrder = pILHS.size() + pSLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setSLHS(std::vector<Spec*> const& slhs) {
    for (auto const& sl: slhs) {
        AssertLog(&sl->getModel() == &pModel);
    }
    pSLHS = slhs;

    if (pOuter) {
        pOrder = pOLHS.size() + pSLHS.size();
    } else {
        pOrder = pILHS.size() + pSLHS.size();
    }
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setIRHS(std::vector<Spec*> const& irhs) {
    for (auto const& ir: irhs) {
        AssertLog(&ir->getModel() == &pModel);
    }
    pIRHS = irhs;
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setSRHS(std::vector<Spec*> const& srhs) {
    for (auto const& sr: srhs) {
        AssertLog(&sr->getModel() == &pModel);
    }
    pSRHS = srhs;
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setORHS(std::vector<Spec*> const& orhs) {
    for (auto const& ors: orhs) {
        AssertLog(&ors->getModel() == &pModel);
    }
    pORHS = orhs;
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setKcst(double kcst) {
    ArgErrLogIf(kcst < 0.0, "Surface reaction constant can't be negative");

    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<Spec*> SReac::getAllSpecs() const {
    AssertLog(pOLHS.empty() || pILHS.empty());

    util::flat_set<Spec*> specs;
    specs.insert(getOLHS().begin(), getOLHS().end());
    specs.insert(getILHS().begin(), getILHS().end());
    specs.insert(getSLHS().begin(), getSLHS().end());
    specs.insert(getIRHS().begin(), getIRHS().end());
    specs.insert(getSRHS().begin(), getSRHS().end());
    specs.insert(getORHS().begin(), getORHS().end());
    return specs;
}

}  // namespace steps::model
