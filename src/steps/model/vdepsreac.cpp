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

#include "vdepsreac.hpp"

#include "spec.hpp"
#include "surfsys.hpp"
#include "util/error.hpp"

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

VDepSReac::VDepSReac(std::string const& id,
                     Surfsys& surfsys,
                     std::vector<Spec*> const& olhs,
                     std::vector<Spec*> const& ilhs,
                     std::vector<Spec*> const& slhs,
                     std::vector<Spec*> const& irhs,
                     std::vector<Spec*> const& srhs,
                     std::vector<Spec*> const& orhs,
                     std::vector<double> ktab,
                     double vmin,
                     double vmax,
                     double dv,
                     uint tablesize)
    : pID(id)
    , pModel(surfsys.getModel())
    , pSurfsys(surfsys)
    , pOuter(false)
    , pOrder(0)
    , pVMin(vmin)
    , pVMax(vmax)
    , pDV(dv)
    , pTablesize(tablesize) {
    // Can't have species on the lhs in the inner and outer compartment
    ArgErrLogIf(!olhs.empty() && !ilhs.empty(),
                "Volume lhs species must belong to either inner or outer "
                "compartment, not both.");

    ArgErrLogIf(ktab.size() != pTablesize, "Table of reaction parameters is not of expected size");

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

    AssertLog(pDV > 0.0);

    // Copy the rate information to local vector
    pK = ktab;

    pSurfsys._handleVDepSReacAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac::~VDepSReac() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::_handleSelfDelete() {
    pSurfsys._handleVDepSReacDel(*this);
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setID(std::string const& id) {
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pSurfsys._handleVDepSReacIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setOLHS(std::vector<Spec*> const& olhs) {
    for (auto const& ol: olhs) {
        AssertLog(&ol->getModel() == &pModel);
    }
    if (!pILHS.empty()) {
        CLOG(WARNING, "general_log") << "\nWARNING: Removing inner compartment species from lhs "
                                        "stoichiometry for VDepSreac "
                                     << getID() << "\n";
    }
    pILHS.clear();
    pOLHS = olhs;
    pOuter = true;
    pOrder = pOLHS.size() + pSLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setILHS(std::vector<Spec*> const& ilhs) {
    for (auto const& il: ilhs) {
        AssertLog(&il->getModel() == &pModel);
    }
    if (!pOLHS.empty()) {
        CLOG(WARNING, "general_log") << "\nWARNING: Removing outer compartment species from lhs "
                                        "stoichiometry for VDepSreac "
                                     << getID() << "\n";
    }
    pILHS = ilhs;
    pOLHS.clear();
    pOuter = false;
    pOrder = pILHS.size() + pSLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setSLHS(std::vector<Spec*> const& slhs) {
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

void VDepSReac::setIRHS(std::vector<Spec*> const& irhs) {
    for (auto const& ir: irhs) {
        AssertLog(&ir->getModel() == &pModel);
    }
    pIRHS = irhs;
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setSRHS(std::vector<Spec*> const& srhs) {
    for (auto const& sr: srhs) {
        AssertLog(&sr->getModel() == &pModel);
    }
    pSRHS = srhs;
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setORHS(std::vector<Spec*> const& orhs) {
    for (auto const& ors: orhs) {
        AssertLog(&ors->getModel() == &pModel);
    }
    pORHS = orhs;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<Spec*> VDepSReac::getAllSpecs() const {
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
