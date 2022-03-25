/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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
/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

#include "sreac.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

#include <easylogging++.h>

#include "model.hpp"
#include "spec.hpp"
#include "surfsys.hpp"

#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

////////////////////////////////////////////////////////////////////////////////

SReac::SReac(string const &id, Surfsys *surfsys, vector<Spec *> const &olhs,
             vector<Spec *> const &ilhs, vector<Spec *> const &slhs,
             vector<Spec *> const &irhs, vector<Spec *> const &srhs,
             vector<Spec *> const &orhs, double kcst)
    : pID(id), pModel(nullptr), pSurfsys(surfsys), pOuter(false), pOLHS(),
      pILHS(), pSLHS(), pIRHS(), pSRHS(), pORHS(), pOrder(0), pKcst(kcst) {

  ArgErrLogIf(pSurfsys == nullptr,
              "No surfsys provided to SReac initializer function");
  ArgErrLogIf(pKcst < 0.0, "Surface reaction constant can't be negative");

  // Can't have species on the lhs in the inner and outer compartment
  ArgErrLogIf(!olhs.empty() && !ilhs.empty(),
              "Volume lhs species must belong to either inner or outer "
              "compartment, not both.");

  pModel = pSurfsys->getModel();
  AssertLog(pModel != nullptr);

  if (!olhs.empty())
    setOLHS(olhs);
  if (!ilhs.empty())
    setILHS(ilhs);
  setSLHS(slhs);
  setIRHS(irhs);
  setSRHS(srhs);
  setORHS(orhs);

  pSurfsys->_handleSReacAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

SReac::~SReac() {
  if (pSurfsys == nullptr) {
    return;
  }
  _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void SReac::_handleSelfDelete() {
  pSurfsys->_handleSReacDel(this);
  pKcst = 0.0;
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

void SReac::setID(string const &id) {
  AssertLog(pSurfsys != nullptr);
  // The following might raise an exception, e.g. if the new ID is not
  // valid or not unique. If this happens, we don't catch but simply let
  // it pass by into the Python layer.
  pSurfsys->_handleSReacIDChange(pID, id);
  // This line will only be executed if the previous call didn't raise
  // an exception.
  pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setOLHS(vector<Spec *> const &olhs) {
  AssertLog(pSurfsys != nullptr);

  if (!pILHS.empty()) {
    ostringstream os;
    os << "\nWARNING: Removing inner compartment species from lhs "
          "stoichiometry for SReac "
       << getID() << ".\n";
  }
  pILHS.clear();
  pOLHS.clear();
  pOLHS.reserve(olhs.size());
  for (auto const &ol : olhs) {
    AssertLog(ol->getModel() == pModel);
    pOLHS.push_back(ol);
  }
  pOuter = true;
  pOrder = pOLHS.size() + pSLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setILHS(vector<Spec *> const &ilhs) {
  AssertLog(pSurfsys != nullptr);

  if (!pOLHS.empty()) {
    ostringstream os;
    os << "\nWARNING: Removing outer compartment species from lhs "
          "stoichiometry for SReac "
       << getID() << ".\n";
  }
  pOLHS.clear();
  pILHS.clear();
  pILHS.reserve(ilhs.size());
  for (auto const &il : ilhs) {
    AssertLog(il->getModel() == pModel);
    pILHS.push_back(il);
  }
  pOuter = false;
  pOrder = pILHS.size() + pSLHS.size();
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setSLHS(vector<Spec *> const &slhs) {
  AssertLog(pSurfsys != nullptr);
  pSLHS.clear();
  pSLHS.reserve(slhs.size());
  for (auto const &sl : slhs) {
    AssertLog(sl->getModel() == pModel);
    pSLHS.push_back(sl);
  }

  if (pOuter) {
    pOrder = pOLHS.size() + pSLHS.size();
  } else {
    pOrder = pILHS.size() + pSLHS.size();
  }
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setIRHS(vector<Spec *> const &irhs) {
  AssertLog(pSurfsys != nullptr);
  pIRHS.clear();
  pIRHS.reserve(irhs.size());
  for (auto const &ir : irhs) {
    AssertLog(ir->getModel() == pModel);
    pIRHS.push_back(ir);
  }
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setSRHS(vector<Spec *> const &srhs) {
  AssertLog(pSurfsys != nullptr);
  pSRHS.clear();
  pSRHS.reserve(srhs.size());
  for (auto const &sr : srhs) {
    AssertLog(sr->getModel() == pModel);
    pSRHS.push_back(sr);
  }
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setORHS(vector<Spec *> const &orhs) {
  AssertLog(pSurfsys != nullptr);
  pORHS.clear();
  pORHS.reserve(orhs.size());
  for (auto const &ors : orhs) {
    AssertLog(ors->getModel() == pModel);
    pORHS.push_back(ors);
  }
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setKcst(double kcst) {
  AssertLog(pSurfsys != nullptr);
  ArgErrLogIf(kcst < 0.0, "Surface reaction constant can't be negative");

  pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

vector<Spec *> SReac::getAllSpecs() const {
  SpecPVec specs;
  bool first_occ;
  AssertLog(pOLHS.empty() || pILHS.empty());

  for (auto const &ol : getOLHS()) {
    first_occ = true;
    for (auto const &s : specs) {
      if (s == ol) {
        first_occ = false;
        break;
      }
    }
    if (first_occ)
      specs.push_back(ol);
  }

  for (const auto &il : getILHS()) {
    first_occ = true;
    for (auto const &s : specs) {
      if (s == il) {
        first_occ = false;
        break;
      }
    }
    if (first_occ)
      specs.push_back(il);
  }

  for (auto const &sl : getSLHS()) {
    first_occ = true;
    for (auto const &s : specs) {
      if (s == sl) {
        first_occ = false;
        break;
      }
    }
    if (first_occ)
      specs.push_back(sl);
  }

  for (auto const &ir : getIRHS()) {
    first_occ = true;
    for (auto const &s : specs) {
      if (s == ir) {
        first_occ = false;
        break;
      }
    }
    if (first_occ)
      specs.push_back(ir);
  }

  for (auto const &sr : getSRHS()) {
    first_occ = true;
    for (auto const &s : specs) {
      if (s == sr) {
        first_occ = false;
        break;
      }
    }
    if (first_occ)
      specs.push_back(sr);
  }

  for (auto const &ors : getORHS()) {
    first_occ = true;
    for (auto const &s : specs) {
      if (s == ors) {
        first_occ = false;
        break;
      }
    }
    if (first_occ)
      specs.push_back(ors);
  }

  return specs;
}

////////////////////////////////////////////////////////////////////////////////
