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

#include "model.hpp"

#include <sstream>
#include <string>

#include <easylogging++.h>

#include "chan.hpp"
#include "ghkcurr.hpp"
#include "model/endocytosis.hpp"
#include "model/exocytosis.hpp"
#include "model/linkspec.hpp"
#include "model/raft.hpp"
#include "model/raftdis.hpp"
#include "model/raftendocytosis.hpp"
#include "model/raftgen.hpp"
#include "model/raftsys.hpp"
#include "model/vesicle.hpp"
#include "model/vessurfsys.hpp"
#include "ohmiccurr.hpp"
#include "spec.hpp"
#include "surfsys.hpp"
#include "vdepsreac.hpp"
#include "volsys.hpp"

#include "util/checkid.hpp"
#include "util/error.hpp"

namespace steps::model {

using util::checkID;

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

    while (!pLinkSpecs.empty()) {
        delete pLinkSpecs.begin()->second;
    }

    while (!pVesicles.empty()) {
        delete pVesicles.begin()->second;
    }

    while (!pRafts.empty()) {
        delete pRafts.begin()->second;
    }

    while (!pVesSurfsys.empty()) {
        delete pVesSurfsys.begin()->second;
    }

    while (!pRaftsys.empty()) {
        delete pRaftsys.begin()->second;
    }
}

////////////////////////////////////////////////////////////////////////////////

Spec* Model::getSpec(std::string const& id) const {
    auto spec = pSpecs.find(id);

    ArgErrLogIf(spec == pSpecs.end(), "Model does not contain species with name '" + id + "'");

    AssertLog(spec->second != nullptr);
    return spec->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delSpec(std::string const& id) const {
    Spec* spec = getSpec(id);
    delete spec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec*> Model::getAllSpecs() const {
    SpecPVec specs;
    specs.reserve(pSpecs.size());
    for (auto const& s: pSpecs) {
        specs.push_back(s.second);
    }
    return specs;
}

////////////////////////////////////////////////////////////////////////////////

LinkSpec* Model::getLinkSpec(std::string const& id) const {
    auto lspec = pLinkSpecs.find(id);
    if (lspec == pLinkSpecs.end()) {
        std::ostringstream os;
        os << "Model does not contain link species with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(lspec->second != nullptr);
    return lspec->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delLinkSpec(std::string const& id) const {
    LinkSpec* lspec = getLinkSpec(id);
    delete lspec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<LinkSpec*> Model::getAllLinkSpecs() const {
    LinkSpecPVec specs;
    specs.reserve(pLinkSpecs.size());
    for (auto const& ls: pLinkSpecs) {
        specs.push_back(ls.second);
    }
    return specs;
}

////////////////////////////////////////////////////////////////////////////////

Vesicle* Model::getVesicle(std::string const& id) const {
    auto it = pVesicles.find(id);
    if (it == pVesicles.end()) {
        std::ostringstream os;
        os << "Model does not contain vesicle with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(it->second != nullptr);
    return it->second;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Vesicle*> Model::getAllVesicles() const {
    VesiclePVec res;
    res.reserve(pVesicles.size());
    for (auto const& ls: pVesicles) {
        res.push_back(ls.second);
    }
    return res;
}

////////////////////////////////////////////////////////////////////////////////

Raft* Model::getRaft(std::string const& id) const {
    auto it = pRafts.find(id);
    if (it == pRafts.end()) {
        std::ostringstream os;
        os << "Model does not contain raft with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(it->second != nullptr);
    return it->second;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Raft*> Model::getAllRafts() const {
    RaftPVec res;
    res.reserve(pRafts.size());
    for (auto const& ls: pRafts) {
        res.push_back(ls.second);
    }
    return res;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Volsys*> Model::getAllVolsyss() const {
    VolsysPVec volsyss;
    volsyss.reserve(pVolsys.size());
    for (auto const& vs: pVolsys) {
        volsyss.push_back(vs.second);
    }
    return volsyss;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Surfsys*> Model::getAllSurfsyss() const {
    SurfsysPVec surfsyss;
    surfsyss.reserve(pSurfsys.size());
    for (auto const& ss: pSurfsys) {
        surfsyss.push_back(ss.second);
    }
    return surfsyss;
}

///////////////////////////////////////////////////////////////////////////////

std::vector<VesSurfsys*> Model::getAllVesSurfsyss() const {
    VesSurfsysPVec vessurfsyss;
    vessurfsyss.reserve(pVesSurfsys.size());
    for (auto const& vs: pVesSurfsys) {
        vessurfsyss.push_back(vs.second);
    }
    return vessurfsyss;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Raftsys*> Model::getAllRaftsyss() const {
    RaftsysPVec raftsyss;
    raftsyss.reserve(pRaftsys.size());
    for (auto const& rs: pRaftsys) {
        raftsyss.push_back(rs.second);
    }
    return raftsyss;
}

////////////////////////////////////////////////////////////////////////////////

Chan* Model::getChan(std::string const& id) const {
    auto chan = pChans.find(id);

    ArgErrLogIf(chan == pChans.end(), "Model does not contain channel with name '" + id + "'");

    AssertLog(chan->second != nullptr);
    return chan->second;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Chan*> Model::getAllChans() const {
    ChanPVec chans;
    chans.reserve(pChans.size());
    for (auto const& c: pChans) {
        chans.push_back(c.second);
    }
    return chans;
}

////////////////////////////////////////////////////////////////////////////////

Volsys* Model::getVolsys(std::string const& id) const {
    auto volsys = pVolsys.find(id);

    ArgErrLogIf(volsys == pVolsys.end(),
                "Model does not contain volume system with name '" + id + "'");

    AssertLog(volsys->second != nullptr);
    return volsys->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delVolsys(std::string const& id) const {
    Volsys* volsys = getVolsys(id);
    delete volsys;
}

////////////////////////////////////////////////////////////////////////////////

Surfsys* Model::getSurfsys(std::string const& id) const {
    auto surfsys = pSurfsys.find(id);

    ArgErrLogIf(surfsys == pSurfsys.end(),
                "Model does not contain surface system with name '" + id + "'");

    AssertLog(surfsys->second != nullptr);
    return surfsys->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delSurfsys(std::string const& id) const {
    Surfsys* surfsys = getSurfsys(id);
    delete surfsys;
}

////////////////////////////////////////////////////////////////////////////////

VesSurfsys* Model::getVesSurfsys(std::string const& id) const {
    const auto& vessurfsys = pVesSurfsys.find(id);
    if (vessurfsys == pVesSurfsys.end()) {
        std::ostringstream os;
        os << "Model does not contain vesicle surface system with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(vessurfsys->second != nullptr);
    return vessurfsys->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delVesSurfsys(std::string const& id) const {
    VesSurfsys* vessurfsys = getVesSurfsys(id);
    delete vessurfsys;
}

////////////////////////////////////////////////////////////////////////////////

Raftsys* Model::getRaftsys(std::string const& id) const {
    auto raftsys = pRaftsys.find(id);
    if (raftsys == pRaftsys.end()) {
        std::ostringstream os;
        os << "Model does not contain raft system with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(raftsys->second != nullptr);
    return raftsys->second;
}

////////////////////////////////////////////////////////////////////////////////

void Model::delRaftsys(std::string const& id) const {
    Raftsys* raftsys = getRaftsys(id);
    delete raftsys;
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkSpecID(std::string const& id) const {
    checkID(id);
    if (pSpecs.find(id) != pSpecs.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkLinkSpecID(std::string const& id) const {
    checkID(id);
    if (pLinkSpecs.find(id) != pLinkSpecs.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkVesicleID(std::string const& id) const {
    checkID(id);
    if (pVesicles.find(id) != pVesicles.end()) {
        std::ostringstream os;
        os << '\'' << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkRaftID(std::string const& id) const {
    checkID(id);
    if (pRafts.find(id) != pRafts.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecIDChange(std::string const& o, std::string const& n) {
    auto s_old = pSpecs.find(o);
    AssertLog(s_old != pSpecs.end());

    if (o == n) {
        return;
    }
    _checkSpecID(n);

    Spec* s = s_old->second;
    AssertLog(s != nullptr);
    pSpecs.erase(s->getID());  // or s_old->first
    pSpecs.emplace(n, s);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleLinkSpecIDChange(std::string const& o, std::string const& n) {
    const auto& ls_old = pLinkSpecs.find(o);
    AssertLog(ls_old != pLinkSpecs.end());

    if (o == n) {
        return;
    }
    _checkLinkSpecID(n);

    LinkSpec* ls = ls_old->second;
    AssertLog(ls != nullptr);
    pLinkSpecs.erase(ls->getID());  // or s_old->first
    pLinkSpecs.emplace(n, ls);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVesicleIDChange(std::string const& o, std::string const& n) {
    auto v_old = pVesicles.find(o);
    AssertLog(v_old != pVesicles.end());

    if (o == n) {
        return;
    }
    _checkVesicleID(n);

    Vesicle* v = v_old->second;
    AssertLog(v != nullptr);
    pVesicles.erase(v->getID());
    pVesicles.emplace(n, v);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleRaftIDChange(std::string const& o, std::string const& n) {
    auto r_old = pRafts.find(o);
    AssertLog(r_old != pRafts.end());

    if (o == n) {
        return;
    }
    _checkRaftID(n);

    Raft* r = r_old->second;
    AssertLog(r != nullptr);
    pRafts.erase(r->getID());
    pRafts.emplace(n, r);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecAdd(Spec* spec) {
    AssertLog(spec->getModel() == this);
    _checkSpecID(spec->getID());
    pSpecs.emplace(spec->getID(), spec);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSpecDel(Spec* spec) {
    AssertLog(spec->getModel() == this);

    for (auto const& vsys: pVolsys) {
        vsys.second->_handleSpecDelete(spec);
    }
    for (auto const& ssys: pSurfsys) {
        ssys.second->_handleSpecDelete(spec);
    }
    pSpecs.erase(spec->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleLinkSpecAdd(LinkSpec* lspec) {
    AssertLog(lspec->getModel() == this);
    _checkLinkSpecID(lspec->getID());
    pLinkSpecs.emplace(lspec->getID(), lspec);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleLinkSpecDel(LinkSpec* lspec) {
    AssertLog(lspec->getModel() == this);

    for (auto const& vsys: pVolsys) {
        vsys.second->_handleLinkSpecDelete(lspec);
    }

    for (auto const& vesssys: pVesSurfsys) {
        vesssys.second->_handleLinkSpecDelete(lspec);
    }

    pLinkSpecs.erase(lspec->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVesicleAdd(Vesicle* vesicle) {
    AssertLog(vesicle->getModel() == this);
    _checkVesicleID(vesicle->getID());
    pVesicles.emplace(vesicle->getID(), vesicle);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVesicleDel(Vesicle* vesicle) {
    AssertLog(vesicle->getModel() == this);
    pVesicles.erase(vesicle->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleRaftAdd(Raft* raft) {
    AssertLog(raft->getModel() == this);
    _checkRaftID(raft->getID());
    pRafts.emplace(raft->getID(), raft);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleRaftDel(Raft* raft) {
    AssertLog(raft->getModel() == this);
    pRafts.erase(raft->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleChanDel(Chan* chan) {
    for (auto const& ssys: pSurfsys) {
        ssys.second->_handleChanDelete(chan);
    }
    pChans.erase(chan->getID());
}
////////////////////////////////////////////////////////////////////////////////

void Model::_checkChanID(std::string const& id) const {
    checkID(id);

    ArgErrLogIf(pChans.find(id) != pChans.end(), "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleChanIDChange(std::string const& o, std::string const& n) {
    auto c_old = pChans.find(o);
    AssertLog(c_old != pChans.end());

    if (o == n) {
        return;
    }
    _checkChanID(n);

    Chan* c = c_old->second;
    AssertLog(c != nullptr);
    pChans.erase(c->getID());  // or c_old->first
    pChans.emplace(n, c);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleChanAdd(Chan* chan) {
    AssertLog(chan->getModel() == this);
    _checkChanID(chan->getID());
    pChans.emplace(chan->getID(), chan);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkVolsysID(std::string const& id) const {
    checkID(id);

    ArgErrLogIf(pVolsys.find(id) != pVolsys.end(), "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysIDChange(std::string const& o, std::string const& n) {
    auto v_old = pVolsys.find(o);
    AssertLog(v_old != pVolsys.end());

    if (o == n) {
        return;
    }
    _checkVolsysID(n);

    Volsys* v = v_old->second;
    AssertLog(v != nullptr);
    pVolsys.erase(v->getID());  // or v_old->first
    pVolsys.emplace(n, v);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysAdd(Volsys* volsys) {
    AssertLog(volsys->getModel() == this);
    _checkVolsysID(volsys->getID());
    pVolsys.emplace(volsys->getID(), volsys);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVolsysDel(Volsys* volsys) {
    AssertLog(volsys->getModel() == this);
    pVolsys.erase(volsys->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkSurfsysID(std::string const& id) const {
    checkID(id);

    ArgErrLogIf(pSurfsys.find(id) != pSurfsys.end(), "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysIDChange(std::string const& o, std::string const& n) {
    auto s_old = pSurfsys.find(o);
    AssertLog(s_old != pSurfsys.end());

    if (o == n) {
        return;
    }
    _checkSurfsysID(n);

    Surfsys* s = s_old->second;
    AssertLog(s != nullptr);
    pSurfsys.erase(s->getID());
    pSurfsys.emplace(n, s);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysAdd(Surfsys* surfsys) {
    AssertLog(surfsys->getModel() == this);
    _checkSurfsysID(surfsys->getID());
    pSurfsys.emplace(surfsys->getID(), surfsys);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleSurfsysDel(Surfsys* surfsys) {
    AssertLog(surfsys->getModel() == this);
    pSurfsys.erase(surfsys->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVesSurfsysAdd(VesSurfsys* vessurfsys) {
    AssertLog(vessurfsys->getModel() == this);
    _checkVesSurfsysID(vessurfsys->getID());
    pVesSurfsys.emplace(vessurfsys->getID(), vessurfsys);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVesSurfsysIDChange(std::string const& o, std::string const& n) {
    auto vs_old = pVesSurfsys.find(o);
    AssertLog(vs_old != pVesSurfsys.end());

    if (o == n) {
        return;
    }
    _checkVesSurfsysID(n);

    VesSurfsys* vs = vs_old->second;
    AssertLog(vs != nullptr);
    pVesSurfsys.erase(vs->getID());
    pVesSurfsys.emplace(n, vs);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleVesSurfsysDel(VesSurfsys* vessurfsys) {
    AssertLog(vessurfsys->getModel() == this);
    pVesSurfsys.erase(vessurfsys->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkVesSurfsysID(std::string const& id) const {
    checkID(id);
    if (pVesSurfsys.find(id) != pVesSurfsys.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_checkRaftsysID(std::string const& id) const {
    checkID(id);
    if (pRaftsys.find(id) != pRaftsys.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleRaftsysIDChange(std::string const& o, std::string const& n) {
    auto rs_old = pRaftsys.find(o);
    AssertLog(rs_old != pRaftsys.end());

    if (o == n) {
        return;
    }
    _checkRaftsysID(n);

    Raftsys* rs = rs_old->second;
    AssertLog(rs != nullptr);
    pRaftsys.erase(rs->getID());
    pRaftsys.emplace(n, rs);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleRaftsysAdd(Raftsys* raftsys) {
    AssertLog(raftsys->getModel() == this);
    _checkRaftsysID(raftsys->getID());
    pRaftsys.emplace(raftsys->getID(), raftsys);
}

////////////////////////////////////////////////////////////////////////////////

void Model::_handleRaftsysDel(Raftsys* raftsys) {
    AssertLog(raftsys->getModel() == this);
    pRaftsys.erase(raftsys->getID());
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countReacs() const {
    uint nreacs = 0;

    for (auto const& vs: pVolsys) {
        nreacs += vs.second->_countReacs();
    }
    return nreacs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countSReacs() const {
    uint nsreacs = 0;

    for (auto const& ss: pSurfsys) {
        nsreacs += ss.second->_countSReacs();
    }
    return nsreacs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countRaftGeneses() const {
    uint nraftgens = 0;

    for (auto const& ss: pSurfsys) {
        nraftgens += ss.second->_countRaftGens();
    }
    return nraftgens;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countRaftDiss() const {
    uint nraftdiss = 0;

    for (auto const& rs: pRaftsys) {
        nraftdiss += rs.second->_countRaftDiss();
    }
    return nraftdiss;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countEndocytosis() const {
    uint nendos = 0;

    for (auto const& ss: pSurfsys) {
        nendos += ss.second->_countEndocytosis();
    }
    return nendos;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countRaftEndocytosis() const {
    uint nendos = 0;

    for (auto const& rs: pRaftsys) {
        nendos += rs.second->_countRaftEndocytosis();
    }
    return nendos;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countExocytosis() const {
    uint nexos = 0;

    for (auto const& vs: pVesSurfsys) {
        nexos += vs.second->_countExocytosis();
    }
    return nexos;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVDiffs() const {
    uint ndiffs = 0;

    for (auto const& vs: pVolsys) {
        ndiffs += vs.second->_countDiffs();
    }
    return ndiffs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countSDiffs() const {
    uint ndiffs = 0;

    for (auto const& ss: pSurfsys) {
        ndiffs += ss.second->_countDiffs();
    }
    return ndiffs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVDepSReacs() const {
    uint nvdsrs = 0;

    for (auto const& ss: pSurfsys) {
        nvdsrs += ss.second->_countVDepSReacs();
    }
    return nvdsrs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countOhmicCurrs() const {
    uint nocs = 0;

    for (auto const& ss: pSurfsys) {
        nocs += ss.second->_countOhmicCurrs();
    }
    return nocs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countGHKcurrs() const {
    uint nghks = 0;

    for (auto const& ss: pSurfsys) {
        nghks += ss.second->_countGHKcurrs();
    }
    return nghks;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVesBinds() const {
    uint nvesbinds = 0;

    for (auto const& vs: pVolsys) {
        nvesbinds += vs.second->_countVesBinds();
    }
    return nvesbinds;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVesUnbinds() const {
    uint nvesunbinds = 0;

    for (auto const& vs: pVolsys) {
        nvesunbinds += vs.second->_countVesUnbinds();
    }
    return nvesunbinds;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVesSDiffs() const {
    uint ndiffs = 0;

    for (auto const& ss: pVesSurfsys) {
        ndiffs += ss.second->_countVesSDiffs();
    }
    return ndiffs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countVesSReacs() const {
    uint nreacs = 0;

    for (auto const& ss: pVesSurfsys) {
        nreacs += ss.second->_countVesSReacs();
    }
    return nreacs;
}

////////////////////////////////////////////////////////////////////////////////

uint Model::_countRaftSReacs() const {
    uint nreacs = 0;

    for (auto const& rss: pRaftsys) {
        nreacs += rss.second->_countRaftSReacs();
    }
    return nreacs;
}

////////////////////////////////////////////////////////////////////////////////

Spec* Model::_getSpec(solver::spec_global_id gidx) const {
    AssertLog(gidx.get() < pSpecs.size());
    auto sp_it = pSpecs.begin();
    std::advance(sp_it, gidx.get());
    return sp_it->second;
}

////////////////////////////////////////////////////////////////////////////////

LinkSpec* Model::_getLinkSpec(solver::linkspec_global_id gidx) const {
    AssertLog(gidx.get() < pLinkSpecs.size());
    auto sp_it = pLinkSpecs.begin();
    std::advance(sp_it, gidx.get());
    return sp_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Vesicle* Model::_getVesicle(solver::vesicle_global_id gidx) const {
    AssertLog(gidx.get() < pVesicles.size());
    auto v_it = pVesicles.begin();
    std::advance(v_it, gidx.get());
    return v_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Raft* Model::_getRaft(solver::raft_global_id gidx) const {
    AssertLog(gidx.get() < pRafts.size());
    auto r_it = pRafts.begin();
    std::advance(r_it, gidx.get());
    return r_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Chan* Model::_getChan(solver::chan_global_id gidx) const {
    AssertLog(gidx.get() < pChans.size());
    auto ch_it = pChans.begin();
    std::advance(ch_it, gidx.get());
    return ch_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Reac* Model::_getReac(solver::reac_global_id gidx) const {
    // first find which volsys this reac (by global index) belongs to
    uint lidx = gidx.get();
    for (auto const& vs: pVolsys) {
        uint reacs_tot = vs.second->_countReacs();
        if (reacs_tot > lidx) {
            return vs.second->_getReac(lidx);
        }
        lidx -= reacs_tot;
    }

    // we shouldn't have gotten to the end
    AssertLog(false);
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

SReac* Model::_getSReac(solver::sreac_global_id gidx) const {
    // first find which surfsys this sreac (by global index) belongs to
    uint lidx = gidx.get();
    for (auto const& ss: pSurfsys) {
        uint sreacs_tot = ss.second->_countSReacs();
        if (sreacs_tot > lidx) {
            return ss.second->_getSReac(lidx);
        }
        lidx -= sreacs_tot;
    }

    // we shouldn't have gotten to the end
    AssertLog(false);
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

RaftGen* Model::_getRaftGen(solver::raftgen_global_id gidx) const {
    // first find which surfsys this raftgen (by global index) belongs to
    uint lidx = gidx.get();

    for (auto const& ss: pSurfsys) {
        uint raftgens_tot = ss.second->_countRaftGens();
        if (raftgens_tot > lidx) {
            return ss.second->_getRaftGen(lidx);
        }
        lidx -= raftgens_tot;
    }

    // we shouldn't have gotten here
    AssertLog(false);
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

RaftDis* Model::_getRaftDis(solver::raftdis_global_id gidx) const {
    // first find which surfsys this raftdis (by global index) belongs to
    uint lidx = gidx.get();

    for (auto const& rs: pRaftsys) {
        uint raftdiss_tot = rs.second->_countRaftDiss();
        if (raftdiss_tot > lidx) {
            return rs.second->_getRaftDis(lidx);
        }
        lidx -= raftdiss_tot;
    }

    // we shouldn't have gotten here
    AssertLog(false);
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

Endocytosis* Model::_getEndocytosis(solver::endocytosis_global_id gidx) const {
    // first find which surfsys this endocyt (by global index) belongs to
    uint lidx = gidx.get();

    for (auto const& ss: pSurfsys) {
        uint endos_tot = ss.second->_countEndocytosis();
        if (endos_tot > lidx) {
            return ss.second->_getEndocytosis(lidx);
        }
        lidx -= endos_tot;
    }

    // we shouldn't have gotten here
    AssertLog(false);
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

RaftEndocytosis* Model::_getRaftEndocytosis(solver::raftendocytosis_global_id gidx) const {
    // first find which raftsys this endocyt (by global index) belongs to
    uint lidx = gidx.get();

    for (auto const& rs: pRaftsys) {
        uint endos_tot = rs.second->_countRaftEndocytosis();
        if (endos_tot > lidx) {
            return rs.second->_getRaftEndocytosis(lidx);
        }
        lidx -= endos_tot;
    }

    // we shouldn't have gotten here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

Exocytosis* Model::_getExocytosis(solver::exocytosis_global_id gidx) const {
    // first find which surfsys this exocyt (by global index) belongs to
    uint lidx = gidx.get();

    for (auto const& vs: pVesSurfsys) {
        uint exos_tot = vs.second->_countExocytosis();
        if (exos_tot > lidx) {
            return vs.second->_getExocytosis(lidx);
        }
        lidx -= exos_tot;
    }

    // we shouldn't have gotten here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

Diff* Model::_getVDiff(solver::diff_global_id gidx) const {
    // first find which volsys this diff (by global index) belongs to
    uint lidx = gidx.get();
    for (auto const& vs: pVolsys) {
        uint diffs_tot = vs.second->_countDiffs();
        if (diffs_tot > lidx) {
            return vs.second->_getDiff(lidx);
        }
        lidx -= diffs_tot;
    }
    // we shouldn't have gotten to the end
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

Diff* Model::_getSDiff(solver::surfdiff_global_id gidx) const {
    // first find which surfsys this diff (by global index) belongs to
    uint lidx = gidx.get();
    for (auto const& ss: pSurfsys) {
        uint diffs_tot = ss.second->_countDiffs();
        if (diffs_tot > lidx) {
            return ss.second->_getDiff(lidx);
        }
        lidx -= diffs_tot;
    }
    // we shouldn't have gotten to the end
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

VesBind* Model::_getVesBind(solver::vesbind_global_id gidx) const {
    // first find which volsys this reac (by global index) belongs to
    uint lidx = gidx.get();

    for (auto const& vs: pVolsys) {
        uint vesbinds_tot = vs.second->_countVesBinds();
        if (vesbinds_tot > lidx) {
            return vs.second->_getVesBind(lidx);
        }
        lidx -= vesbinds_tot;
    }

    // we shouldn't have gotten here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

VesUnbind* Model::_getVesUnbind(solver::vesunbind_global_id gidx) const {
    // first find which volsys this reac (by global index) belongs to
    uint lidx = gidx.get();

    for (auto const& vs: pVolsys) {
        uint vesbinds_tot = vs.second->_countVesUnbinds();
        if (vesbinds_tot > lidx) {
            return vs.second->_getVesUnbind(lidx);
        }
        lidx -= vesbinds_tot;
    }

    // we shouldn't have gotten here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

VesSDiff* Model::_getVesSDiff(solver::vessdiff_global_id gidx) const {
    // first find which surfsys this diff (by global index) belongs to
    uint lidx = gidx.get();

    for (auto const& vs: pVesSurfsys) {
        uint diffs_tot = vs.second->_countVesSDiffs();
        if (diffs_tot > lidx) {
            return vs.second->_getVesSDiff(lidx);
        }
        lidx -= diffs_tot;
    }
    // we shouldn't have gotten to the end
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

VesSReac* Model::_getVesSReac(solver::vessreac_global_id gidx) const {
    // first find which surfsys this reac (by global index) belongs to
    uint lidx = gidx.get();

    for (auto const& vs: pVesSurfsys) {
        uint reacs_tot = vs.second->_countVesSReacs();
        if (reacs_tot > lidx) {
            return vs.second->_getVesSReac(lidx);
        }
        lidx -= reacs_tot;
    }
    // we shouldn't have gotten to the end
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

RaftSReac* Model::_getRaftSReac(solver::raftsreac_global_id gidx) const {
    // first find which raftsys this sreac (by global index) belongs to
    uint lidx = gidx.get();

    for (auto const& rss: pRaftsys) {
        uint reacs_tot = rss.second->_countRaftSReacs();
        if (reacs_tot > lidx) {
            return rss.second->_getRaftSReac(lidx);
        }
        lidx -= reacs_tot;
    }
    // we shouldn't have gotten to the end
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac* Model::_getVDepSReac(solver::vdepsreac_global_id gidx) const {
    // first find which surfsys this v-dep-sreac (by global index) belongs to
    uint lidx = gidx.get();
    for (auto const& ss: pSurfsys) {
        uint vdsrs_tot = ss.second->_countVDepSReacs();
        if (vdsrs_tot > lidx) {
            return ss.second->_getVDepSReac(lidx);
        }
        lidx -= vdsrs_tot;
    }

    // we shouldn't have gotten here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr* Model::_getOhmicCurr(solver::ohmiccurr_global_id gidx) const {
    // first find which surfsys this ohmic current (by global index) belongs to
    uint lidx = gidx.get();
    for (auto const& ss: pSurfsys) {
        uint ocs_tot = ss.second->_countOhmicCurrs();
        if (ocs_tot > lidx) {
            return ss.second->_getOhmicCurr(lidx);
        }
        lidx -= ocs_tot;
    }

    // we shouldn't have gotten here
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr* Model::_getGHKcurr(solver::ghkcurr_global_id gidx) const {
    // first find which surfsys this ghk current (by global index) belongs to
    uint lidx = gidx.get();
    for (auto const& ss: pSurfsys) {
        uint ghks_tot = ss.second->_countGHKcurrs();
        if (ghks_tot > lidx) {
            return ss.second->_getGHKcurr(lidx);
        }
        lidx -= ghks_tot;
    }

    // we shouldn't have gotten here
    AssertLog(false);
}

}  // namespace steps::model
