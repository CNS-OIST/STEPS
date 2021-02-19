/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

// STL headers.
#include <cassert>
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/model.hpp"
#include "steps/model/spec.hpp"
#include "steps/model/surfsys.hpp"
#include "steps/model/vdepsreac.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

////////////////////////////////////////////////////////////////////////////////

VDepSReac::VDepSReac(std::string const & id, Surfsys * surfsys,
                     vector<Spec *> const & olhs, vector<Spec *> const & ilhs,
                     vector<Spec *> const & slhs,
                     vector<Spec *> const & irhs, vector<Spec *> const & srhs,
                     vector<Spec *> const & orhs,
                     vector<double> ktab, double vmin, double vmax,
                     double dv, uint tablesize)
: pID(id)
, pModel(nullptr)
, pSurfsys(surfsys)
, pOuter(false)
, pOrder(0)
, pK(nullptr)
, pVMin(vmin)
, pVMax(vmax)
, pDV(dv)
, pTablesize(tablesize)
{
    if (pSurfsys == nullptr)
    {
        ostringstream os;
        os << "No surfsys provided to SReac initializer function";
        ArgErrLog(os.str());
    }

    // Can't have species on the lhs in the inner and outer compartment
    if (!olhs.empty() && !ilhs.empty())
    {
        ostringstream os;
        os << "Volume lhs species must belong to either inner or outer ";
        os << "compartment, not both.";
        ArgErrLog(os.str());
    }

    if (ktab.size() != pTablesize)
    {
        ostringstream os;
        os << "Table of reaction parameters is not of expected size";
        ArgErrLog(os.str());
    }
    pModel = pSurfsys->getModel();
    AssertLog(pModel != nullptr);

    if (!olhs.empty()) setOLHS(olhs);
    if (!ilhs.empty()) setILHS(ilhs);
    setSLHS(slhs);
    setIRHS(irhs);
    setSRHS(srhs);
    setORHS(orhs);

    AssertLog(pDV > 0.0);

    // Copy the rate information to local array
    pK = new double[pTablesize];
    std::memcpy(pK, ktab.data(), pTablesize * sizeof(double));

    pSurfsys->_handleVDepSReacAdd(this);

}

////////////////////////////////////////////////////////////////////////////////

VDepSReac::~VDepSReac()
{
    if (pSurfsys == nullptr) { return;
}
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::_handleSelfDelete()
{
    pSurfsys->_handleVDepSReacDel(this);
    delete[] pK;
    pOrder = 0;
    pORHS.clear();
    pSRHS.clear();
    pIRHS.clear();
    pSLHS.clear();
    pILHS.clear();
    pOLHS.clear();
    pSurfsys = nullptr;
    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setID(string const & id)
{
    AssertLog(pSurfsys != nullptr);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pSurfsys->_handleVDepSReacIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setOLHS(vector<Spec *> const & olhs)
{
    AssertLog(pSurfsys != nullptr);

    if (!pILHS.empty())
    {
        ostringstream os;
        os << "\nWARNING: Removing inner compartment species from lhs stoichiometry for VDepSreac " << getID() <<"\n";
    }
    pILHS.clear();
    pOLHS.clear();
    for (auto const& ol: olhs) {
        AssertLog(ol->getModel() == pModel);
        pOLHS.push_back(ol);
    }
    pOuter = true;
    pOrder = pOLHS.size() + pSLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setILHS(vector<Spec *> const & ilhs)
{
    AssertLog(pSurfsys != nullptr);

    if (!pOLHS.empty())
    {
        ostringstream os;
        os << "\nWARNING: Removing outer compartment species from lhs stoichiometry for VDepSreac " << getID() <<"\n";
    }
    pOLHS.clear();
    pILHS.clear();
    for (auto const& il: ilhs) {
        AssertLog(il->getModel() == pModel);
        pILHS.push_back(il);
    }
    pOuter = false;
    pOrder = pILHS.size() + pSLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setSLHS(vector<Spec *> const & slhs)
{
    AssertLog(pSurfsys != nullptr);
    pSLHS.clear();
    for (auto const& sl: slhs) {
        AssertLog(sl->getModel() == pModel);
        pSLHS.push_back(sl);
    }

    if (pOuter) { pOrder = pOLHS.size() + pSLHS.size(); }
    else { pOrder = pILHS.size() + pSLHS.size(); }
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setIRHS(vector<Spec *> const & irhs)
{
    AssertLog(pSurfsys != nullptr);
    pIRHS.clear();
    for (auto const& ir: irhs) {
        AssertLog(ir->getModel() == pModel);
        pIRHS.push_back(ir);
    }
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setSRHS(vector<Spec *> const & srhs)
{
    AssertLog(pSurfsys != nullptr);
    pSRHS.clear();
    pSRHS.reserve(srhs.size());
    for (auto const& sr: srhs) {
        AssertLog(sr->getModel() == pModel);
        pSRHS.push_back(sr);
    }
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setORHS(vector<Spec *> const & orhs)
{
    AssertLog(pSurfsys != nullptr);
    pORHS.clear();
    pORHS.reserve(orhs.size());
    for (auto const& ors: orhs) {
        AssertLog(ors->getModel() == pModel);
        pORHS.push_back(ors);
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> VDepSReac::getK() const
{
    std::vector<double> k(pK, pK + pTablesize);
    return k;
}

////////////////////////////////////////////////////////////////////////////////

vector<Spec *> VDepSReac::getAllSpecs() const
{
    SpecPVec specs;
    bool first_occ;
    AssertLog(pOLHS.empty() || pILHS.empty());

    for (auto const& ol: getOLHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == ol) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) specs.push_back(ol);
    }

    for (auto const& il: getILHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == il) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) specs.push_back(il);
    }

    for (auto const& sl: getSLHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == sl) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) specs.push_back(sl);
    }

    for (auto const& ir: getIRHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == ir) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) specs.push_back(ir);
    }

    for (auto const& sr:getSRHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == sr) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) specs.push_back(sr);
    }

    for (auto const& ors: getORHS()) {
        first_occ = true;
        for (auto const& s: specs) {
            if (s == ors) {
                first_occ = false;
                break;
            }
        }
        if (first_occ) specs.push_back(ors);
    }

    return specs;
}

////////////////////////////////////////////////////////////////////////////////

// END
