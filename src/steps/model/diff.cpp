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

/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

#include "diff.hpp"

#include <cassert>
#include <sstream>
#include <string>

#include <easylogging++.h>

#include "model.hpp"
#include "spec.hpp"
#include "volsys.hpp"
#include "surfsys.hpp"

#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

////////////////////////////////////////////////////////////////////////////////

Diff::Diff(string const &id, Volsys *volsys, Spec *lig, double dcst)
    : pID(id), pVolsys(volsys), pLig(lig), pDcst(dcst), pIsvolume(true) {

  ArgErrLogIf(pVolsys == nullptr,
              "No volsys provided to Diff initializer function.");
  ArgErrLogIf(pDcst < 0.0, "Diffusion constant can't be negative");

  pModel = pVolsys->getModel();
  AssertLog(pModel != nullptr);

  pVolsys->_handleDiffAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Diff::Diff(string const &id, Surfsys *surfsys, Spec *lig, double dcst)
    : pID(id), pSurfsys(surfsys), pLig(lig), pDcst(dcst), pIsvolume(false) {

  ArgErrLogIf(pSurfsys == nullptr,
              "No surfsys provided to Diff initializer function.");

  ArgErrLogIf(pDcst < 0.0, "Diffusion constant can't be negative");

  pModel = pSurfsys->getModel();
  AssertLog(pModel != nullptr);

  pSurfsys->_handleDiffAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Diff::~Diff() {
  if (pIsvolume) {
    if (pVolsys == nullptr) {
      return;
    }
  } else {
    if (pSurfsys == nullptr) {
      return;
    }
  }
  _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Diff::_handleSelfDelete() {
  if (pIsvolume) {
    pVolsys->_handleDiffDel(this);
    pVolsys = nullptr;
  } else {
    pSurfsys->_handleDiffDel(this);
    pSurfsys = nullptr;
  }
  pDcst = 0.0;
  pLig = nullptr;
  pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setID(string const &id) {
  if (pIsvolume) {
    AssertLog(pVolsys != nullptr);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pVolsys->_handleDiffIDChange(pID, id);
  } else {
    AssertLog(pSurfsys != nullptr);
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pSurfsys->_handleDiffIDChange(pID, id);
  }
  // This line will only be executed if the previous call didn't raise
  // an exception.
  pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setDcst(double dcst) {
  if (pIsvolume) {
    AssertLog(pVolsys != nullptr);
  } else {
    AssertLog(pSurfsys != nullptr);
  }

  ArgErrLogIf(dcst < 0.0, "Diffusion constant can't be negative");

  pDcst = dcst;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setLig(Spec *lig) {
  AssertLog(lig != nullptr);
  pLig = lig;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec *> Diff::getAllSpecs() const { return {pLig}; }

////////////////////////////////////////////////////////////////////////////////
