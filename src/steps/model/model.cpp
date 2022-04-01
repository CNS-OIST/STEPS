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

#include "model.hpp"

#include <cassert>
#include <sstream>
#include <string>

#include <easylogging++.h>

#include "chan.hpp"
#include "ghkcurr.hpp"
#include "ohmiccurr.hpp"
#include "spec.hpp"
#include "surfsys.hpp"
#include "vdepsreac.hpp"
#include "vdeptrans.hpp"
#include "volsys.hpp"

#include "util/checkid.hpp"
#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

using steps::util::checkID;

////////////////////////////////////////////////////////////////////////////////

Model::~Model() {
  while (!pSpecs.empty()) {
    delete pSpecs.begin()->second;
  }

  while (!pChans.empty()) {
    delete pChans.begin()->second;
  }

  while (!pVolsys.empty()) {
    delete pVolsys.begin()->second;
  }

  while (!pSurfsys.empty()) {
    delete pSurfsys.begin()->second;
  }
}

////////////////////////////////////////////////////////////////////////////////

Spec *Model::getSpec(string const &id) const {
  auto spec = pSpecs.find(id);

  ArgErrLogIf(spec == pSpecs.end(),
              "Model does not contain species with name '" + id + "'");

  AssertLog(spec->second != nullptr);
  return spec->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delSpec(string const &id) {
  Spec *spec = getSpec(id);
  delete spec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec *> Model::getAllSpecs() const {
  SpecPVec specs;
  specs.reserve(pSpecs.size());
  for (auto const &s : pSpecs) {
    specs.push_back(s.second);
  }
  return specs;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Volsys *> Model::getAllVolsyss() const {
  VolsysPVec volsyss;
  volsyss.reserve(pVolsys.size());
  for (auto const &vs : pVolsys) {
    volsyss.push_back(vs.second);
  }
  return volsyss;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Surfsys *> Model::getAllSurfsyss() const {
  SurfsysPVec surfsyss;
  surfsyss.reserve(pSurfsys.size());
  for (auto const &ss : pSurfsys) {
    surfsyss.push_back(ss.second);
  }
  return surfsyss;
}

////////////////////////////////////////////////////////////////////////////////

Chan *Model::getChan(string const &id) const {
  auto chan = pChans.find(id);

  ArgErrLogIf(chan == pChans.end(),
              "Model does not contain channel with name '" + id + "'");

  AssertLog(chan->second != nullptr);
  return chan->second;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Chan *> Model::getAllChans() const {
  ChanPVec chans;
  chans.reserve(pChans.size());
  for (auto const &c : pChans) {
    chans.push_back(c.second);
  }
  return chans;
}

////////////////////////////////////////////////////////////////////////////////

Volsys *Model::getVolsys(string const &id) const {
  auto volsys = pVolsys.find(id);

  ArgErrLogIf(volsys == pVolsys.end(),
              "Model does not contain volume system with name '" + id + "'");

  AssertLog(volsys->second != nullptr);
  return volsys->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delVolsys(string const &id) {
  Volsys *volsys = getVolsys(id);
  delete volsys;
}

////////////////////////////////////////////////////////////////////////////////

Surfsys *Model::getSurfsys(string const &id) const {
  auto surfsys = pSurfsys.find(id);

  ArgErrLogIf(surfsys == pSurfsys.end(),
              "Model does not contain surface system with name '" + id + "'");

  AssertLog(surfsys->second != nullptr);
  return surfsys->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delSurfsys(string const &id) {
  Surfsys *surfsys = getSurfsys(id);
  delete surfsys;
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkSpecID(string const &id) const {
  checkID(id);

  ArgErrLogIf(pSpecs.find(id) != pSpecs.end(),
              "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecIDChange(string const &o, string const &n) {
  auto s_old = pSpecs.find(o);
  AssertLog(s_old != pSpecs.end());

  if (o == n)
    return;
  _checkSpecID(n);

  Spec *s = s_old->second;
  AssertLog(s != nullptr);
  pSpecs.erase(s->getID()); // or s_old->first
  pSpecs.insert(SpecPMap::value_type(n, s));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecAdd(Spec *spec) {
  AssertLog(spec->getModel() == this);
  _checkSpecID(spec->getID());
  pSpecs.insert(SpecPMap::value_type(spec->getID(), spec));
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecDel(Spec *spec) {
  for (auto const &vsys : pVolsys) {
    vsys.second->_handleSpecDelete(spec);
  }
  for (auto const &ssys : pSurfsys) {
    ssys.second->_handleSpecDelete(spec);
  }
  pSpecs.erase(spec->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleChanDel(Chan *chan) {
  for (auto const &ssys : pSurfsys) {
    ssys.second->_handleChanDelete(chan);
  }
  pChans.erase(chan->getID());
}
////////////////////////////////////////////////////////////////////////////////

void Model::_checkChanID(string const &id) const {
  checkID(id);

  ArgErrLogIf(pChans.find(id) != pChans.end(),
              "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleChanIDChange(string const &o, string const &n) {
  auto c_old = pChans.find(o);
  AssertLog(c_old != pChans.end());

  if (o == n)
    return;
  _checkChanID(n);

  Chan *c = c_old->second;
  AssertLog(c != nullptr);
  pChans.erase(c->getID()); // or c_old->first
  pChans.emplace(n, c);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleChanAdd(Chan *chan) {
  AssertLog(chan->getModel() == this);
  _checkChanID(chan->getID());
  pChans.emplace(chan->getID(), chan);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkVolsysID(string const &id) const {
  checkID(id);

  ArgErrLogIf(pVolsys.find(id) != pVolsys.end(),
              "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysIDChange(string const &o, string const &n) {
  auto v_old = pVolsys.find(o);
  AssertLog(v_old != pVolsys.end());

  if (o == n)
    return;
  _checkVolsysID(n);

  Volsys *v = v_old->second;
  AssertLog(v != nullptr);
  pVolsys.erase(v->getID()); // or v_old->first
  pVolsys.emplace(n, v);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysAdd(Volsys *volsys) {
  AssertLog(volsys->getModel() == this);
  _checkVolsysID(volsys->getID());
  pVolsys.emplace(volsys->getID(), volsys);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysDel(Volsys *volsys) {
  AssertLog(volsys->getModel() == this);
  pVolsys.erase(volsys->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkSurfsysID(string const &id) const {
  checkID(id);

  ArgErrLogIf(pSurfsys.find(id) != pSurfsys.end(),
              "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysIDChange(string const &o, string const &n) {
  auto s_old = pSurfsys.find(o);
  AssertLog(s_old != pSurfsys.end());

  if (o == n)
    return;
  _checkSurfsysID(n);

  Surfsys *s = s_old->second;
  AssertLog(s != nullptr);
  pSurfsys.erase(s->getID());
  pSurfsys.emplace(n, s);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysAdd(Surfsys *surfsys) {
  AssertLog(surfsys->getModel() == this);
  _checkSurfsysID(surfsys->getID());
  pSurfsys.emplace(surfsys->getID(), surfsys);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysDel(Surfsys *surfsys) {
  AssertLog(surfsys->getModel() == this);
  pSurfsys.erase(surfsys->getID());
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countReacs() const {
  uint nreacs = 0;

  for (auto const &vs : pVolsys) {
    nreacs += vs.second->_countReacs();
  }
  return nreacs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countSReacs() const {
  uint nsreacs = 0;

  for (auto const &ss : pSurfsys) {
    nsreacs += ss.second->_countSReacs();
  }
  return nsreacs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVDiffs() const {
  uint ndiffs = 0;

  for (auto const &vs : pVolsys) {
    ndiffs += vs.second->_countDiffs();
  }
  return ndiffs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countSDiffs() const {
  uint ndiffs = 0;

  for (auto const &ss : pSurfsys) {
    ndiffs += ss.second->_countDiffs();
  }
  return ndiffs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVDepTrans() const {
  uint nvdts = 0;

  for (auto const &ss : pSurfsys) {
    nvdts += ss.second->_countVDepTrans();
  }
  return nvdts;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVDepSReacs() const {
  uint nvdsrs = 0;

  for (auto const &ss : pSurfsys) {
    nvdsrs += ss.second->_countVDepSReacs();
  }
  return nvdsrs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countOhmicCurrs() const {
  uint nocs = 0;

  for (auto const &ss : pSurfsys) {
    nocs += ss.second->_countOhmicCurrs();
  }
  return nocs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countGHKcurrs() const {
  uint nghks = 0;

  for (auto const &ss : pSurfsys) {
    nghks += ss.second->_countGHKcurrs();
  }
  return nghks;
}

////////////////////////////////////////////////////////////////////////////////

Spec *Model::_getSpec(uint gidx) const {
  AssertLog(gidx < pSpecs.size());
  auto sp_it = pSpecs.begin();
  std::advance(sp_it, gidx);
  return sp_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Chan *Model::_getChan(uint gidx) const {
  AssertLog(gidx < pChans.size());
  auto ch_it = pChans.begin();
  std::advance(ch_it, gidx);
  return ch_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Reac *Model::_getReac(uint gidx) const {
  // first find which volsys this reac (by global index) belongs to
  uint lidx = gidx;
  for (auto const &vs : pVolsys) {
    uint reacs_tot = vs.second->_countReacs();
    if (reacs_tot > lidx) {
      return vs.second->_getReac(lidx);
    }
    lidx -= reacs_tot;
  }

  // we shouldn't have gotten here
  AssertLog(false);
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

SReac *Model::_getSReac(uint gidx) const {
  // first find which surfsys this sreac (by global index) belongs to
  uint lidx = gidx;
  for (auto const &ss : pSurfsys) {
    uint sreacs_tot = ss.second->_countSReacs();
    if (sreacs_tot > lidx) {
      return ss.second->_getSReac(lidx);
    }
    lidx -= sreacs_tot;
  }

  // we shouldn't have gotten here
  AssertLog(false);
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

Diff *Model::_getVDiff(uint gidx) const {
  // first find which volsys this diff (by global index) belongs to
  uint lidx = gidx;
  for (auto const &vs : pVolsys) {
    uint diffs_tot = vs.second->_countDiffs();
    if (diffs_tot > lidx) {
      return vs.second->_getDiff(lidx);
    }
    lidx -= diffs_tot;
  }
  // we shouldn't have gotten to the end
  AssertLog(false);
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

Diff *Model::_getSDiff(uint gidx) const {
  // first find which surfsys this diff (by global index) belongs to
  uint lidx = gidx;
  for (auto const &ss : pSurfsys) {
    uint diffs_tot = ss.second->_countDiffs();
    if (diffs_tot > lidx) {
      return ss.second->_getDiff(lidx);
    }
    lidx -= diffs_tot;
  }
  // we shouldn't have gotten to the end
  AssertLog(false);
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

VDepTrans *Model::_getVDepTrans(uint gidx) const {
  // first find which surfsys this v-dep-trans (by global index) belongs to
  uint lidx = gidx;
  for (auto const &ss : pSurfsys) {
    uint vdts_tot = ss.second->_countVDepTrans();
    if (vdts_tot > lidx) {
      return ss.second->_getVDepTrans(lidx);
    }
    lidx -= vdts_tot;
  }

  // we shouldn't have gotten here
  AssertLog(false);
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac *Model::_getVDepSReac(uint gidx) const {
  // first find which surfsys this v-dep-sreac (by global index) belongs to
  uint lidx = gidx;
  for (auto const &ss : pSurfsys) {
    uint vdsrs_tot = ss.second->_countVDepSReacs();
    if (vdsrs_tot > lidx) {
      return ss.second->_getVDepSReac(lidx);
    }
    lidx -= vdsrs_tot;
  }

  // we shouldn't have gotten here
  AssertLog(false);
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr *Model::_getOhmicCurr(uint gidx) const {
  // first find which surfsys this ohmic current (by global index) belongs to
  uint lidx = gidx;
  for (auto const &ss : pSurfsys) {
    uint ocs_tot = ss.second->_countOhmicCurrs();
    if (ocs_tot > lidx) {
      return ss.second->_getOhmicCurr(lidx);
    }
    lidx -= ocs_tot;
  }

  // we shouldn't have gotten here
  AssertLog(false);
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr *Model::_getGHKcurr(uint gidx) const {
  // first find which surfsys this ghk current (by global index) belongs to
  uint lidx = gidx;
  for (auto const &ss : pSurfsys) {
    uint ghks_tot = ss.second->_countGHKcurrs();
    if (ghks_tot > lidx) {
      return ss.second->_getGHKcurr(lidx);
    }
    lidx -= ghks_tot;
  }

  // we shouldn't have gotten here
  AssertLog(false);
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

// END
