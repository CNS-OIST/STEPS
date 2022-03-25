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

#include "surfsys.hpp"

#include <cassert>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include <easylogging++.h>

#include "chan.hpp"
#include "chanstate.hpp"
#include "diff.hpp"
#include "ghkcurr.hpp"
#include "model.hpp"
#include "ohmiccurr.hpp"
#include "spec.hpp"
#include "sreac.hpp"
#include "vdepsreac.hpp"
#include "vdeptrans.hpp"

#include "util/error.hpp"
#include "util/checkid.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

using steps::util::checkID;

////////////////////////////////////////////////////////////////////////////////

Surfsys::Surfsys(string const &id, Model *model) : pID(id), pModel(model) {

  ArgErrLogIf(pModel == nullptr,
              "No model provided to Surfsys initializer function");

  pModel->_handleSurfsysAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Surfsys::~Surfsys() {
  if (pModel == nullptr) {
    return;
  }
  _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::setID(string const &id) {
  AssertLog(pModel != nullptr);
  if (id == pID)
    return;
  // The following might raise an exception, e.g. if the new ID is not
  // valid or not unique. If this happens, we don't catch but simply let
  // it pass by into the Python layer.
  pModel->_handleSurfsysIDChange(pID, id);
  // This line will only be executed if the previous call didn't raise
  // an exception.
  pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSelfDelete() {
  for (auto const &sreac : getAllSReacs()) {
    delete sreac;
  }
  for (auto const &vdtrans : getAllVDepTrans()) {
    delete vdtrans;
  }
  for (auto const &vdsreac : getAllVDepSReacs()) {
    delete vdsreac;
  }
  for (auto const &oc : getAllOhmicCurrs()) {
    delete oc;
  }
  for (auto const &ghk : getAllGHKcurrs()) {
    delete ghk;
  }

  for (auto const &diff : getAllDiffs()) {
    delete diff;
  }

  pModel->_handleSurfsysDel(this);

  pSReacs.clear();
  pVDepTrans.clear();
  pVDepSReacs.clear();
  pOhmicCurrs.clear();
  pGHKcurrs.clear();

  pDiffs.clear();

  pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

SReac *Surfsys::getSReac(string const &id) const {
  auto sreac = pSReacs.find(id);

  ArgErrLogIf(sreac == pSReacs.end(), "Model does not contain surface "
                                      "reaction with name '" +
                                          id + "'");

  AssertLog(sreac->second != nullptr);
  return sreac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delSReac(string const &id) {
  SReac *sreac = getSReac(id);
  delete sreac;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<SReac *> Surfsys::getAllSReacs() const {
  SReacPVec sreacs;
  sreacs.reserve(pSReacs.size());
  for (auto const &sr : pSReacs) {
    sreacs.emplace_back(sr.second);
  }
  return sreacs;
}

////////////////////////////////////////////////////////////////////////////////

Diff *Surfsys::getDiff(string const &id) const {
  auto diff = pDiffs.find(id);

  ArgErrLogIf(diff == pDiffs.end(),
              "Model does not contain diffusion with name '" + id + "'");

  AssertLog(diff->second != nullptr);
  return diff->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delDiff(string const &id) {
  Diff *diff = getDiff(id);
  // delete diff object since it is owned by c++, not python
  delete diff;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Diff *> Surfsys::getAllDiffs() const {
  DiffPVec diffs;
  diffs.reserve(pDiffs.size());
  for (auto const &d : pDiffs) {
    diffs.emplace_back(d.second);
  }
  return diffs;
}

////////////////////////////////////////////////////////////////////////////////

SpecPVec Surfsys::getAllSpecs() const {
  std::set<Spec*> specs_set;
  for (auto const &sreac : getAllSReacs()) {
    auto specs = sreac->getAllSpecs();
    specs_set.insert(specs.begin(), specs.end());
  }

  for (auto const &vdepsreac : getAllVDepSReacs()) {
    auto specs = vdepsreac->getAllSpecs();
    specs_set.insert(specs.begin(), specs.end());
  }

  for (auto const &vdeptrans : getAllVDepTrans()) {
    SpecP dst = vdeptrans->getDst();
    SpecP src = vdeptrans->getSrc();
    specs_set.insert(dst);
    specs_set.insert(src);
  }

  for (auto const &ghk : getAllGHKcurrs()) {
    SpecP ghk_spec = ghk->getIon();
    specs_set.insert(ghk_spec);
  }

  for (auto const &oc : getAllOhmicCurrs()) {
    SpecP oc_spec = oc->getChanState();
    specs_set.insert(oc_spec);
  }

  for (auto const &diff : getAllDiffs()) {
    auto specs = diff->getAllSpecs();
    specs_set.insert(specs.begin(), specs.end());
  }

  return {specs_set.begin(), specs_set.end()};
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkSReacID(string const &id) const {
  checkID(id);

  ArgErrLogIf(pSReacs.find(id) != pSReacs.end(),
              "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacIDChange(string const &o, string const &n) {
  SReacPMapCI sr_old = pSReacs.find(o);
  AssertLog(sr_old != pSReacs.end());

  if (o == n)
    return;
  _checkSReacID(n);

  SReac *sr = sr_old->second;
  AssertLog(sr != nullptr);
  pSReacs.erase(sr->getID());
  pSReacs.emplace(n, sr);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacAdd(SReac *sreac) {
  AssertLog(sreac->getSurfsys() == this);
  _checkSReacID(sreac->getID());
  pSReacs.insert(SReacPMap::value_type(sreac->getID(), sreac));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacDel(SReac *sreac) {
  AssertLog(sreac->getSurfsys() == this);
  pSReacs.erase(sreac->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkDiffID(string const &id) const {
  checkID(id);

  ArgErrLogIf(pDiffs.find(id) != pDiffs.end(),
              "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleDiffIDChange(string const &o, string const &n) {
  DiffPMapCI d_old = pDiffs.find(o);
  AssertLog(d_old != pDiffs.end());

  if (o == n)
    return;
  _checkDiffID(n);

  Diff *d = d_old->second;
  AssertLog(d != nullptr);
  AssertLog(pDiffs.erase(d->getID()) == 1);
  pDiffs.insert(DiffPMap::value_type(n, d));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleDiffAdd(Diff *diff) {
  AssertLog(diff->getSurfsys() == this);
  _checkDiffID(diff->getID());
  pDiffs.insert(DiffPMap::value_type(diff->getID(), diff));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleDiffDel(Diff *diff) {
  AssertLog(diff->getSurfsys() == this);
  pDiffs.erase(diff->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSpecDelete(Spec *spec) {
  {
    std::vector<std::string> sreacs_del;
    for (auto const &sreac : pSReacs) {
      SpecPVec specs = sreac.second->getAllSpecs();
      for (auto const &sr_spec : specs) {
        if (sr_spec == spec) {
          sreacs_del.emplace_back(sreac.second->getID());
          break;
        }
      }
    }
    for (auto const &sr_del : sreacs_del) {
      delSReac(sr_del);
    }
  }
  {
    std::vector<std::string> ghks_del;
    for (auto const &ghk : pGHKcurrs) {
      SpecP ion = ghk.second->getIon();
      if (ion == spec) {
        ghks_del.emplace_back(ghk.second->getID());
      }
      // spec may be a channel state
      SpecP cstate = ghk.second->getChanState();
      if (cstate == spec) {
        ghks_del.emplace_back(ghk.second->getID());
      }
    }
    for (auto const &ghk : ghks_del) {
      delGHKcurr(ghk);
    }
  }
  {
    // spec may also be a derived ChanState object -> need to delete any
    // vdeptrans and ohmic currents that include this channel state
    std::vector<std::string> oc_del;
    for (auto const &oc : pOhmicCurrs) {
      SpecP cstate = oc.second->getChanState();
      if (cstate == spec) {
        oc_del.emplace_back(cstate->getID());
      }
    }
    for (auto const &occurr_del : oc_del) {
      delOhmicCurr(occurr_del);
    }
  }
  {
    std::vector<std::string> vdt_del;
    for (auto const &vdt : pVDepTrans) {
      SpecP dst = vdt.second->getDst();
      SpecP src = vdt.second->getSrc();
      if (dst == spec or src == spec) {
        vdt_del.emplace_back(vdt.second->getID());
      }
    }
    for (auto const &vdept_del : vdt_del) {
      delVDepTrans(vdept_del);
    }
  }
  {
    std::vector<std::string> vdepsreacs_del;
    for (auto const &vdepsreac : pVDepSReacs) {
      SpecPVec specs = vdepsreac.second->getAllSpecs();
      for (auto const &sr_spec : specs) {
        if (sr_spec == spec) {
          vdepsreacs_del.emplace_back(vdepsreac.second->getID());
          break;
        }
      }
    }
    for (auto const &vdsr_del : vdepsreacs_del) {
      delVDepSReac(vdsr_del);
    }
  }
  {
    std::vector<std::string> diffs_del;
    for (auto const &diff : pDiffs) {
      SpecPVec specs = diff.second->getAllSpecs();
      for (auto const &d_spec : specs) {
        if (d_spec == spec) {
          diffs_del.emplace_back(diff.second->getID());
          break;
        }
      }
    }
    for (auto const &d_del : diffs_del) {
      delDiff(d_del);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleChanDelete(Chan *chan) {
  {
    std::vector<std::string> vdtrans_del;
    for (auto const &vdtrans : pVDepTrans) {
      ChanP chans = vdtrans.second->getChan();
      if (chans == chan) {
        vdtrans_del.emplace_back(vdtrans.second->getID());
      }
    }
    for (auto const &vd_del : vdtrans_del) {
      delVDepTrans(vd_del);
    }
  }
  {
    std::vector<std::string> ohmcurr_del;
    for (auto const &ohmcurr : pOhmicCurrs) {
      ChanP chans = ohmcurr.second->getChanState()->getChan();
      if (chans == chan) {
        ohmcurr_del.emplace_back(ohmcurr.second->getID());
      }
    }
    for (auto const &oc_del : ohmcurr_del) {
      delOhmicCurr(oc_del);
    }
  }
  {
    std::vector<std::string> ghkcurr_del;
    for (auto const &ghkcurr : pGHKcurrs) {
      ChanP chans = ghkcurr.second->getChanState()->getChan();
      if (chans == chan) {
        ghkcurr_del.emplace_back(ghkcurr.second->getID());
      }
    }
    for (auto const &ghk_del : ghkcurr_del) {
      delGHKcurr(ghk_del);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

VDepTrans *Surfsys::getVDepTrans(std::string const &id) const {
  auto vdeptrans = pVDepTrans.find(id);

  ArgErrLogIf(
      vdeptrans == pVDepTrans.end(),
      "Model does not contain voltage-dependent transition with name '" + id +
          "'");

  AssertLog(vdeptrans->second != nullptr);
  return vdeptrans->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delVDepTrans(std::string const &id) {
  VDepTrans *vdeptrans = getVDepTrans(id);
  delete vdeptrans;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<VDepTrans *> Surfsys::getAllVDepTrans() const {
  VDepTransPVec vdeptrans;
  vdeptrans.reserve(pVDepTrans.size());
  for (auto const &vd : pVDepTrans) {
    vdeptrans.emplace_back(vd.second);
  }
  return vdeptrans;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepTransIDChange(string const &o, string const &n) {
  VDepTransPMapCI vd_old = pVDepTrans.find(o);
  AssertLog(vd_old != pVDepTrans.end());

  if (o == n)
    return;
  _checkVDepTransID(n);

  VDepTrans *vd = vd_old->second;
  AssertLog(vd != nullptr);
  pVDepTrans.erase(vd->getID());
  pVDepTrans.insert(VDepTransPMap::value_type(n, vd));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkVDepTransID(string const &id) const {
  checkID(id);

  ArgErrLogIf(pVDepTrans.find(id) != pVDepTrans.end(),
              "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepTransAdd(VDepTrans *vdeptrans) {
  AssertLog(vdeptrans->getSurfsys() == this);
  _checkVDepTransID(vdeptrans->getID());
  pVDepTrans.insert(VDepTransPMap::value_type(vdeptrans->getID(), vdeptrans));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepTransDel(VDepTrans *vdeptrans) {
  AssertLog(vdeptrans->getSurfsys() == this);
  pVDepTrans.erase(vdeptrans->getID());
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac *Surfsys::getVDepSReac(std::string const &id) const {
  auto vdepsreac = pVDepSReacs.find(id);

  ArgErrLogIf(
      vdepsreac == pVDepSReacs.end(),
      "Model does not contain voltage-dependent surface reaction with name '" +
          id + "'");

  AssertLog(vdepsreac->second != nullptr);
  return vdepsreac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delVDepSReac(std::string const &id) {
  VDepSReac *vdepsreac = getVDepSReac(id);
  delete vdepsreac;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<VDepSReac *> Surfsys::getAllVDepSReacs() const {
  VDepSReacPVec vdepsreac;
  vdepsreac.reserve(pVDepSReacs.size());
  for (auto const &vd : pVDepSReacs) {
    vdepsreac.emplace_back(vd.second);
  }
  return vdepsreac;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacIDChange(string const &o, string const &n) {
  auto vd_old = pVDepSReacs.find(o);
  AssertLog(vd_old != pVDepSReacs.end());

  if (o == n)
    return;
  _checkVDepSReacID(n);

  VDepSReac *vd = vd_old->second;
  AssertLog(vd != nullptr);
  pVDepSReacs.erase(vd->getID());
  pVDepSReacs.insert(VDepSReacPMap::value_type(n, vd));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkVDepSReacID(string const &id) const {
  checkID(id);

  ArgErrLogIf(pVDepSReacs.find(id) != pVDepSReacs.end(),
              "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacAdd(VDepSReac *vdepsreac) {
  AssertLog(vdepsreac->getSurfsys() == this);
  _checkVDepSReacID(vdepsreac->getID());
  pVDepSReacs.insert(VDepSReacPMap::value_type(vdepsreac->getID(), vdepsreac));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacDel(VDepSReac *vdepsreac) {
  AssertLog(vdepsreac->getSurfsys() == this);
  pVDepSReacs.erase(vdepsreac->getID());
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr *Surfsys::getOhmicCurr(std::string const &id) const {
  auto ohmiccurr = pOhmicCurrs.find(id);

  ArgErrLogIf(ohmiccurr == pOhmicCurrs.end(),
              "Model does not contain ohmic current with name '" + id + "'");

  AssertLog(ohmiccurr->second != nullptr);
  return ohmiccurr->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delOhmicCurr(std::string const &id) {
  OhmicCurr *ohmiccurr = getOhmicCurr(id);
  delete ohmiccurr;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<OhmicCurr *> Surfsys::getAllOhmicCurrs() const {
  OhmicCurrPVec ohmiccurr;
  ohmiccurr.reserve(pOhmicCurrs.size());
  for (auto const &oc : pOhmicCurrs) {
    ohmiccurr.emplace_back(oc.second);
  }
  return ohmiccurr;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrIDChange(string const &o, string const &n) {
  OhmicCurrPMapCI oc_old = pOhmicCurrs.find(o);
  AssertLog(oc_old != pOhmicCurrs.end());

  if (o == n)
    return;
  _checkOhmicCurrID(n);

  OhmicCurr *oc = oc_old->second;
  AssertLog(oc != nullptr);
  pOhmicCurrs.erase(oc->getID());
  pOhmicCurrs.insert(OhmicCurrPMap::value_type(n, oc));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkOhmicCurrID(string const &id) const {
  checkID(id);

  ArgErrLogIf(pOhmicCurrs.find(id) != pOhmicCurrs.end(),
              "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrAdd(OhmicCurr *ohmiccurr) {
  AssertLog(ohmiccurr->getSurfsys() == this);
  _checkOhmicCurrID(ohmiccurr->getID());
  pOhmicCurrs.insert(OhmicCurrPMap::value_type(ohmiccurr->getID(), ohmiccurr));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrDel(OhmicCurr *ohmiccurr) {
  AssertLog(ohmiccurr->getSurfsys() == this);
  pOhmicCurrs.erase(ohmiccurr->getID());
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr *Surfsys::getGHKcurr(std::string const &id) const {
  auto ghkcurr = pGHKcurrs.find(id);

  ArgErrLogIf(ghkcurr == pGHKcurrs.end(),
              "Model does not contain ghk current with name '" + id + "'");

  AssertLog(ghkcurr->second != nullptr);
  return ghkcurr->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delGHKcurr(std::string const &id) {
  GHKcurr *ghkcurr = getGHKcurr(id);
  delete ghkcurr;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<GHKcurr *> Surfsys::getAllGHKcurrs() const {
  GHKcurrPVec ghkcurr;
  ghkcurr.reserve(pGHKcurrs.size());
  for (auto const &ghk : pGHKcurrs) {
    ghkcurr.emplace_back(ghk.second);
  }
  return ghkcurr;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrIDChange(string const &o, string const &n) {
  auto ghk_old = pGHKcurrs.find(o);
  AssertLog(ghk_old != pGHKcurrs.end());

  if (o == n)
    return;
  _checkGHKcurrID(n);

  GHKcurr *ghk = ghk_old->second;
  AssertLog(ghk != nullptr);
  pGHKcurrs.erase(ghk->getID());
  pGHKcurrs.insert(GHKcurrPMap::value_type(n, ghk));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkGHKcurrID(string const &id) const {
  checkID(id);

  ArgErrLogIf(pGHKcurrs.find(id) != pGHKcurrs.end(),
              "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrAdd(GHKcurr *ghkcurr) {
  AssertLog(ghkcurr->getSurfsys() == this);
  _checkGHKcurrID(ghkcurr->getID());
  pGHKcurrs.insert(GHKcurrPMap::value_type(ghkcurr->getID(), ghkcurr));
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrDel(GHKcurr *ghkcurr) {
  AssertLog(ghkcurr->getSurfsys() == this);
  pGHKcurrs.erase(ghkcurr->getID());
}

////////////////////////////////////////////////////////////////////////////////

SReac *Surfsys::_getSReac(uint lidx) const {
  AssertLog(lidx < pSReacs.size());
  auto sr_it = pSReacs.begin();
  std::advance(sr_it, lidx);
  return sr_it->second;
}

////////////////////////////////////////////////////////////////////////////////

VDepTrans *Surfsys::_getVDepTrans(uint lidx) const {
  AssertLog(lidx < pVDepTrans.size());
  auto vd_it = pVDepTrans.begin();
  std::advance(vd_it, lidx);
  return vd_it->second;
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac *Surfsys::_getVDepSReac(uint lidx) const {
  AssertLog(lidx < pVDepSReacs.size());
  auto vd_it = pVDepSReacs.begin();
  std::advance(vd_it, lidx);
  return vd_it->second;
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr *Surfsys::_getOhmicCurr(uint lidx) const {
  AssertLog(lidx < pOhmicCurrs.size());
  auto oc_it = pOhmicCurrs.begin();
  std::advance(oc_it, lidx);
  return oc_it->second;
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr *Surfsys::_getGHKcurr(uint lidx) const {
  AssertLog(lidx < pGHKcurrs.size());
  auto ghk_it = pGHKcurrs.begin();
  std::advance(ghk_it, lidx);
  return ghk_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Diff *Surfsys::_getDiff(uint lidx) const {
  AssertLog(lidx < pDiffs.size());
  auto df_it = pDiffs.begin();
  std::advance(df_it, lidx);
  return df_it->second;
}

////////////////////////////////////////////////////////////////////////////////
