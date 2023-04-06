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

#include "patch.hpp"

#include <cassert>
#include <sstream>

#include <easylogging++.h>

#include "model/chanstate.hpp"
#include "model/diff.hpp"
#include "model/ghkcurr.hpp"
#include "model/model.hpp"
#include "model/ohmiccurr.hpp"
#include "model/sreac.hpp"
#include "model/vdepsreac.hpp"
#include "model/vdeptrans.hpp"
#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace wm {

Patch::Patch(std::string id, Geom *container, Comp *icomp, Comp *ocomp,
             double area)
    : pArea(area), pID(std::move(id)), pContainer(container) {
  ArgErrLogIf(pContainer == nullptr,
              "No container provided to Patch initializer function.\n");

  _setIComp(icomp);
  if (ocomp != nullptr) {
    _setOComp(ocomp);
  }

  ArgErrLogIf(pArea < 0.0, "Patch area can't be negative.\n");

  pContainer->_handlePatchAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Patch::~Patch() {
  if (pContainer == nullptr) {
    return;
  }
  _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setID(std::string const &id) {
  AssertLog(pContainer != nullptr);
  if (id == pID) {
    return;
  }
  // The following might raise an exception, e.g. if the new ID is not
  // valid or not unique. If this happens, we don't catch but simply let
  // it pass by into the Python layer.
  pContainer->_handlePatchIDChange(pID, id);
  // This line will only be executed if the previous call didn't raise
  // an exception.
  pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setArea(double area) {
  AssertLog(pContainer != nullptr);

  ArgErrLogIf(area < 0.0, "Patch area can't be negative.\n");

  pArea = area;
}

////////////////////////////////////////////////////////////////////////////////

void Patch::addSurfsys(std::string const &id) {
  // string identifier is only added to set if it is not already included
  pSurfsys.insert(id);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::delSurfsys(std::string const &id) {
  // string identifier is only removed from set if it is included
  pSurfsys.erase(id);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::Spec *>
Patch::getAllSpecs(const steps::model::Model *model) const {

  std::set<steps::model::Spec *> specs;

  for (const auto &id : pSurfsys) {
    steps::model::Surfsys *surfsys = model->getSurfsys(id);
    // std::vector<steps::model::Spec *> specs = surfsys->getAllSpecs();
    for (auto const &sreac : surfsys->getAllSReacs()) {
      auto sreac_specs_slhs = sreac->getSLHS();
      specs.insert(sreac_specs_slhs.begin(), sreac_specs_slhs.end());
      auto sreac_specs_srhs = sreac->getSRHS();
      specs.insert(sreac_specs_srhs.begin(), sreac_specs_srhs.end());
    }
    for (auto const &sdiff : surfsys->getAllDiffs()) {
      specs.insert(sdiff->getLig());
    }
    for (auto const &vdeptrans : surfsys->getAllVDepTrans()) {
      specs.insert(vdeptrans->getSrc());
      specs.insert(vdeptrans->getDst());
    }
    for (auto const &vdepsreac : surfsys->getAllVDepSReacs()) {
      auto vdepsreac_specs_slhs = vdepsreac->getSLHS();
      specs.insert(vdepsreac_specs_slhs.begin(), vdepsreac_specs_slhs.end());
      auto vdepsreac_specs_srhs = vdepsreac->getSRHS();
      specs.insert(vdepsreac_specs_srhs.begin(), vdepsreac_specs_srhs.end());
    }
    for (auto const &oc : surfsys->getAllOhmicCurrs()) {
      specs.insert(oc->getChanState());
    }
    for (auto const &ghk : surfsys->getAllGHKcurrs()) {
      specs.insert(ghk->getChanState());
    }
  }

  return {specs.begin(), specs.end()};
} // namespace wm

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::SReac *>
Patch::getAllSReacs(const steps::model::Model *model) const {
  std::set<steps::model::SReac *> pSReacs;
  for (const auto &id : pSurfsys) {
    steps::model::Surfsys *surfsys = model->getSurfsys(id);
    std::vector<steps::model::SReac *> sreacs = surfsys->getAllSReacs();
    pSReacs.insert(sreacs.begin(), sreacs.end());
  }

  return {pSReacs.begin(), pSReacs.end()};
}

////////////////////////////////////////////////////////////////////////////////

void Patch::_setIComp(Comp *icomp) {
  ArgErrLogIf(icomp->getContainer() != pContainer,
              "Compartment does not belong to same container as patch.\n");

  auto const &ipatches = icomp->getIPatches();

  ArgErrLogIf(ipatches.find(this) != ipatches.end(),
              "Patch is already on inside of compartment.\n");

  // remove the patch if it was already on the outside of some
  // other compartment
  if (pIComp != nullptr) {
    pIComp->_delOPatch(this);
  }

  pIComp = icomp;
  pIComp->_addOPatch(this);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::_setOComp(Comp *ocomp) {
  if (ocomp == nullptr) {
    return;
  }

  ArgErrLogIf(ocomp->getContainer() != pContainer,
              "Compartment does not belong to same container as patch.\n");

  auto const &opatches = ocomp->getOPatches();

  ArgErrLogIf(opatches.find(this) != opatches.end(),
              "Patch is already on outside of compartment.\n");

  // remove the patch if it was already on the inside of some
  // other compartment
  if (pOComp != nullptr) {
    pOComp->_delIPatch(this);
  }

  pOComp = ocomp;
  pOComp->_addIPatch(this);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::_handleSelfDelete() {
  pContainer->_handlePatchDel(this);
  pArea = 0.0;
  pSurfsys.clear();
  pIComp = nullptr;
  pOComp = nullptr;
  pContainer = nullptr;
}

} // namespace wm
} // namespace steps
