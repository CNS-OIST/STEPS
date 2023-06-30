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

#include "surfsys.hpp"

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
#include "model/endocytosis.hpp"
#include "model/raftgen.hpp"
#include "ohmiccurr.hpp"
#include "spec.hpp"
#include "sreac.hpp"
#include "vdepsreac.hpp"

#include "util/checkid.hpp"
#include "util/error.hpp"

namespace steps::model {

using util::checkID;

////////////////////////////////////////////////////////////////////////////////

Surfsys::Surfsys(std::string const& id, Model* model)
    : pID(id)
    , pModel(model) {
    ArgErrLogIf(pModel == nullptr, "No model provided to Surfsys initializer function");

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

void Surfsys::setID(std::string const& id) {
    AssertLog(pModel != nullptr);
    if (id == pID) {
        return;
    }
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
    for (auto const& sreac: getAllSReacs()) {
        delete sreac;
    }
    for (auto const& vdsreac: getAllVDepSReacs()) {
        delete vdsreac;
    }
    for (auto const& oc: getAllOhmicCurrs()) {
        delete oc;
    }
    for (auto const& ghk: getAllGHKcurrs()) {
        delete ghk;
    }

    for (auto const& diff: getAllDiffs()) {
        delete diff;
    }

    for (auto const& endo: getAllEndocytosis()) {
        delete endo;
    }

    for (auto const& raftgen: getAllRaftGens()) {
        delete raftgen;
    }

    pModel->_handleSurfsysDel(this);

    pSReacs.clear();
    pVDepSReacs.clear();
    pOhmicCurrs.clear();
    pGHKcurrs.clear();
    pDiffs.clear();

    pRaftGens.clear();
    pEndocytosis.clear();

    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

SReac* Surfsys::getSReac(std::string const& id) const {
    auto sreac = pSReacs.find(id);

    ArgErrLogIf(sreac == pSReacs.end(),
                "Model does not contain surface "
                "reaction with name '" +
                    id + "'");

    AssertLog(sreac->second != nullptr);
    return sreac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delSReac(std::string const& id) {
    SReac* sreac = getSReac(id);
    pSReacs.erase(id);
    delete sreac;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<SReac*> Surfsys::getAllSReacs() const {
    SReacPVec sreacs;
    sreacs.reserve(pSReacs.size());
    for (auto const& sr: pSReacs) {
        sreacs.emplace_back(sr.second);
    }
    return sreacs;
}

////////////////////////////////////////////////////////////////////////////////

RaftGen* Surfsys::getRaftGen(std::string const& id) const {
    auto raftgen = pRaftGens.find(id);
    if (raftgen == pRaftGens.end()) {
        std::ostringstream os;
        os << "Model does not contain raft "
              "genesis with name '"
           << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(raftgen->second != nullptr);
    return raftgen->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delRaftGen(std::string const& id) {
    RaftGen* raftgen = getRaftGen(id);
    pRaftGens.erase(id);
    delete raftgen;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<RaftGen*> Surfsys::getAllRaftGens() const {
    RaftGenPVec raftgens;
    raftgens.reserve(pRaftGens.size());

    for (auto const& rg: pRaftGens) {
        raftgens.push_back(rg.second);
    }
    return raftgens;
}

////////////////////////////////////////////////////////////////////////////////

Endocytosis* Surfsys::getEndocytosis(std::string const& id) const {
    auto endo = pEndocytosis.find(id);
    if (endo == pEndocytosis.end()) {
        std::ostringstream os;
        os << "Model does not contain endocytotic "
              "reaction with name '"
           << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(endo->second != nullptr);
    return endo->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delEndocytosis(std::string const& id) {
    Endocytosis* endo = getEndocytosis(id);
    pEndocytosis.erase(id);
    delete endo;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Endocytosis*> Surfsys::getAllEndocytosis() const {
    EndocytosisPVec endos;
    endos.reserve(pEndocytosis.size());
    for (auto const& endo: pEndocytosis) {
        endos.push_back(endo.second);
    }
    return endos;
}

////////////////////////////////////////////////////////////////////////////////

Diff* Surfsys::getDiff(std::string const& id) const {
    auto diff = pDiffs.find(id);
    if (diff == pDiffs.end()) {
        std::ostringstream os;
        os << "Model does not contain diffusion with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(diff->second != nullptr);
    return diff->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delDiff(std::string const& id) {
    Diff* diff = getDiff(id);
    pDiffs.erase(id);
    // delete diff object since it is owned by c++, not python
    delete diff;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Diff*> Surfsys::getAllDiffs() const {
    DiffPVec diffs;
    diffs.reserve(pDiffs.size());
    for (auto const& d: pDiffs) {
        diffs.emplace_back(d.second);
    }
    return diffs;
}

////////////////////////////////////////////////////////////////////////////////

SpecPVec Surfsys::getAllSpecs() const {
    std::set<Spec*> specs_set;
    for (auto const& sreac: getAllSReacs()) {
        const auto& specs = sreac->getAllSpecs();
        specs_set.insert(specs.begin(), specs.end());
    }

    for (auto const& diff: getAllDiffs()) {
        const auto& specs = diff->getAllSpecs();
        specs_set.insert(specs.begin(), specs.end());
    }


    for (auto const& vdepsreac: getAllVDepSReacs()) {
        const auto& specs = vdepsreac->getAllSpecs();
        specs_set.insert(specs.begin(), specs.end());
    }

    for (auto const& ghk: getAllGHKcurrs()) {
        SpecP ghk_spec = ghk->getIon();
        specs_set.insert(ghk_spec);
    }

    for (auto const& oc: getAllOhmicCurrs()) {
        SpecP oc_spec = oc->getChanState();
        specs_set.insert(oc_spec);
    }

    for (auto const& endo: getAllEndocytosis()) {
        const auto& specs = endo->getAllSpecs();
        specs_set.insert(specs.begin(), specs.end());
    }

    for (auto const& raftgen: getAllRaftGens()) {
        const auto& specs = raftgen->getAllSpecs();
        specs_set.insert(specs.begin(), specs.end());
    }

    return {specs_set.begin(), specs_set.end()};
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkSReacID(std::string const& id) const {
    checkID(id);

    ArgErrLogIf(pSReacs.find(id) != pSReacs.end(), "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacIDChange(std::string const& o, std::string const& n) {
    auto sr_old = pSReacs.find(o);
    AssertLog(sr_old != pSReacs.end());

    if (o == n) {
        return;
    }
    _checkSReacID(n);

    SReac* sr = sr_old->second;
    AssertLog(sr != nullptr);
    pSReacs.erase(sr->getID());
    pSReacs.emplace(n, sr);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacAdd(SReac* sreac) {
    AssertLog(sreac->getSurfsys() == this);
    _checkSReacID(sreac->getID());
    pSReacs.emplace(sreac->getID(), sreac);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSReacDel(SReac* sreac) {
    AssertLog(sreac->getSurfsys() == this);
    pSReacs.erase(sreac->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkRaftGenID(std::string const& id) const {
    checkID(id);
    if (pRaftGens.find(id) != pRaftGens.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleRaftGenIDChange(std::string const& o, std::string const& n) {
    if (o == n) {
        return;
    }
    auto rg_old = pRaftGens.find(o);
    AssertLog(rg_old != pRaftGens.end());

    _checkRaftGenID(n);

    RaftGen* rg = rg_old->second;
    AssertLog(rg != nullptr);
    pRaftGens.erase(rg->getID());
    pRaftGens.emplace(n, rg);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleRaftGenAdd(RaftGen* raftgen) {
    AssertLog(raftgen->getSurfsys() == this);
    _checkRaftGenID(raftgen->getID());
    pRaftGens.emplace(raftgen->getID(), raftgen);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleRaftGenDel(RaftGen* raftgen) {
    AssertLog(raftgen->getSurfsys() == this);
    pRaftGens.erase(raftgen->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkEndocytosisID(std::string const& id) const {
    checkID(id);
    if (pEndocytosis.find(id) != pEndocytosis.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleEndocytosisIDChange(std::string const& o, std::string const& n) {
    auto e_old = pEndocytosis.find(o);
    AssertLog(e_old != pEndocytosis.end());

    if (o == n) {
        return;
    }
    _checkEndocytosisID(n);

    Endocytosis* e = e_old->second;
    AssertLog(e != nullptr);
    pEndocytosis.erase(e->getID());
    pEndocytosis.emplace(n, e);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleEndocytosisAdd(Endocytosis* endocyt) {
    AssertLog(endocyt->getSurfsys() == this);
    _checkEndocytosisID(endocyt->getID());
    pEndocytosis.emplace(endocyt->getID(), endocyt);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleEndocytosisDel(Endocytosis* endocyt) {
    AssertLog(endocyt->getSurfsys() == this);
    pEndocytosis.erase(endocyt->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkDiffID(std::string const& id) const {
    checkID(id);
    if (pDiffs.find(id) != pDiffs.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleDiffIDChange(std::string const& o, std::string const& n) {
    auto d_old = pDiffs.find(o);
    AssertLog(d_old != pDiffs.end());

    if (o == n) {
        return;
    }
    _checkDiffID(n);

    Diff* d = d_old->second;
    AssertLog(d != nullptr);
    AssertLog(pDiffs.erase(d->getID()) == 1);
    pDiffs.emplace(n, d);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleDiffAdd(Diff* diff) {
    AssertLog(diff->getSurfsys() == this);
    _checkDiffID(diff->getID());
    pDiffs.emplace(diff->getID(), diff);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleDiffDel(Diff* diff) {
    AssertLog(diff->getSurfsys() == this);
    pDiffs.erase(diff->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleSpecDelete(Spec* spec) {
    {
        std::vector<std::string> sreacs_del;
        for (auto const& sreac: pSReacs) {
            for (auto const& sr_spec: sreac.second->getAllSpecs()) {
                if (sr_spec == spec) {
                    sreacs_del.push_back(sreac.second->getID());
                    break;
                }
            }
        }
        for (auto const& sr_del: sreacs_del) {
            delSReac(sr_del);
        }
    }
    {
        std::vector<std::string> ghks_del;
        for (auto const& ghk: pGHKcurrs) {
            SpecP ion = ghk.second->getIon();
            if (ion == spec) {
                ghks_del.push_back(ghk.second->getID());
            }
            // spec may be a channel state
            SpecP cstate = ghk.second->getChanState();
            if (cstate == spec) {
                ghks_del.push_back(ghk.second->getID());
            }
        }
        for (auto const& ghk: ghks_del) {
            delGHKcurr(ghk);
        }
    }
    {
        // spec may also be a derived ChanState object -> need to delete any
        // ohmic currents that include this channel state
        std::vector<std::string> oc_del;
        for (auto const& oc: pOhmicCurrs) {
            SpecP cstate = oc.second->getChanState();
            if (cstate == spec) {
                oc_del.push_back(cstate->getID());
            }
        }
        for (auto const& occurr_del: oc_del) {
            delOhmicCurr(occurr_del);
        }
    }
    {
        std::vector<std::string> vdepsreacs_del;
        for (auto const& vdepsreac: pVDepSReacs) {
            for (auto const& sr_spec: vdepsreac.second->getAllSpecs()) {
                if (sr_spec == spec) {
                    vdepsreacs_del.push_back(vdepsreac.second->getID());
                    break;
                }
            }
        }
        for (auto const& vdsr_del: vdepsreacs_del) {
            delVDepSReac(vdsr_del);
        }
    }
    {
        std::vector<std::string> diffs_del;
        for (auto const& diff: pDiffs) {
            for (auto const& d_spec: diff.second->getAllSpecs()) {
                if (d_spec == spec) {
                    diffs_del.push_back(diff.second->getID());
                    break;
                }
            }
        }
        for (auto const& d_del: diffs_del) {
            delDiff(d_del);
        }
    }

    std::vector<std::string> endos_del;
    for (auto const& endo: pEndocytosis) {
        for (auto const& e_spec: endo.second->getAllSpecs()) {
            if (e_spec == spec) {
                endos_del.push_back(endo.second->getID());
                break;
            }
        }
    }
    for (auto const& e_del: endos_del) {
        delEndocytosis(e_del);
    }

    std::vector<std::string> raftgens_del;
    for (auto const& raftgen: pRaftGens) {
        for (auto const& rg_spec: raftgen.second->getAllSpecs()) {
            if (rg_spec == spec) {
                raftgens_del.push_back(raftgen.second->getID());
                break;
            }
        }
    }
    for (auto const& rg_del: raftgens_del) {
        delRaftGen(rg_del);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleChanDelete(Chan* chan) {
    {
        std::vector<std::string> ohmcurr_del;
        for (auto const& ohmcurr: pOhmicCurrs) {
            ChanP chans = ohmcurr.second->getChanState()->getChan();
            if (chans == chan) {
                ohmcurr_del.push_back(ohmcurr.second->getID());
            }
        }
        for (auto const& oc_del: ohmcurr_del) {
            delOhmicCurr(oc_del);
        }
    }
    {
        std::vector<std::string> ghkcurr_del;
        for (auto const& ghkcurr: pGHKcurrs) {
            ChanP chans = ghkcurr.second->getChanState()->getChan();
            if (chans == chan) {
                ghkcurr_del.push_back(ghkcurr.second->getID());
            }
        }
        for (auto const& ghk_del: ghkcurr_del) {
            delGHKcurr(ghk_del);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac* Surfsys::getVDepSReac(std::string const& id) const {
    auto vdepsreac = pVDepSReacs.find(id);

    ArgErrLogIf(vdepsreac == pVDepSReacs.end(),
                "Model does not contain voltage-dependent surface reaction with name '" + id + "'");

    AssertLog(vdepsreac->second != nullptr);
    return vdepsreac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delVDepSReac(std::string const& id) {
    VDepSReac* vdepsreac = getVDepSReac(id);
    pVDepSReacs.erase(id);
    delete vdepsreac;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<VDepSReac*> Surfsys::getAllVDepSReacs() const {
    VDepSReacPVec vdepsreac;
    vdepsreac.reserve(pVDepSReacs.size());
    for (auto const& vd: pVDepSReacs) {
        vdepsreac.emplace_back(vd.second);
    }
    return vdepsreac;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacIDChange(std::string const& o, std::string const& n) {
    auto vd_old = pVDepSReacs.find(o);
    AssertLog(vd_old != pVDepSReacs.end());

    if (o == n) {
        return;
    }
    _checkVDepSReacID(n);

    VDepSReac* vd = vd_old->second;
    AssertLog(vd != nullptr);
    pVDepSReacs.erase(vd->getID());
    pVDepSReacs.emplace(n, vd);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkVDepSReacID(std::string const& id) const {
    checkID(id);

    ArgErrLogIf(pVDepSReacs.find(id) != pVDepSReacs.end(), "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacAdd(VDepSReac* vdepsreac) {
    AssertLog(vdepsreac->getSurfsys() == this);
    _checkVDepSReacID(vdepsreac->getID());
    pVDepSReacs.emplace(vdepsreac->getID(), vdepsreac);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleVDepSReacDel(VDepSReac* vdepsreac) {
    AssertLog(vdepsreac->getSurfsys() == this);
    pVDepSReacs.erase(vdepsreac->getID());
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr* Surfsys::getOhmicCurr(std::string const& id) const {
    auto ohmiccurr = pOhmicCurrs.find(id);

    ArgErrLogIf(ohmiccurr == pOhmicCurrs.end(),
                "Model does not contain ohmic current with name '" + id + "'");

    AssertLog(ohmiccurr->second != nullptr);
    return ohmiccurr->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delOhmicCurr(std::string const& id) {
    OhmicCurr* ohmiccurr = getOhmicCurr(id);
    pOhmicCurrs.erase(id);
    delete ohmiccurr;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<OhmicCurr*> Surfsys::getAllOhmicCurrs() const {
    OhmicCurrPVec ohmiccurr;
    ohmiccurr.reserve(pOhmicCurrs.size());
    for (auto const& oc: pOhmicCurrs) {
        ohmiccurr.emplace_back(oc.second);
    }
    return ohmiccurr;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrIDChange(std::string const& o, std::string const& n) {
    auto oc_old = pOhmicCurrs.find(o);
    AssertLog(oc_old != pOhmicCurrs.end());

    if (o == n) {
        return;
    }
    _checkOhmicCurrID(n);

    OhmicCurr* oc = oc_old->second;
    AssertLog(oc != nullptr);
    pOhmicCurrs.erase(oc->getID());
    pOhmicCurrs.emplace(n, oc);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkOhmicCurrID(std::string const& id) const {
    checkID(id);

    ArgErrLogIf(pOhmicCurrs.find(id) != pOhmicCurrs.end(), "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrAdd(OhmicCurr* ohmiccurr) {
    AssertLog(ohmiccurr->getSurfsys() == this);
    _checkOhmicCurrID(ohmiccurr->getID());
    pOhmicCurrs.emplace(ohmiccurr->getID(), ohmiccurr);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleOhmicCurrDel(OhmicCurr* ohmiccurr) {
    AssertLog(ohmiccurr->getSurfsys() == this);
    pOhmicCurrs.erase(ohmiccurr->getID());
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr* Surfsys::getGHKcurr(std::string const& id) const {
    auto ghkcurr = pGHKcurrs.find(id);

    ArgErrLogIf(ghkcurr == pGHKcurrs.end(),
                "Model does not contain ghk current with name '" + id + "'");

    AssertLog(ghkcurr->second != nullptr);
    return ghkcurr->second;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::delGHKcurr(std::string const& id) {
    GHKcurr* ghkcurr = getGHKcurr(id);
    pGHKcurrs.erase(id);
    delete ghkcurr;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<GHKcurr*> Surfsys::getAllGHKcurrs() const {
    GHKcurrPVec ghkcurr;
    ghkcurr.reserve(pGHKcurrs.size());
    for (auto const& ghk: pGHKcurrs) {
        ghkcurr.emplace_back(ghk.second);
    }
    return ghkcurr;
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrIDChange(std::string const& o, std::string const& n) {
    auto ghk_old = pGHKcurrs.find(o);
    AssertLog(ghk_old != pGHKcurrs.end());

    if (o == n) {
        return;
    }
    _checkGHKcurrID(n);

    GHKcurr* ghk = ghk_old->second;
    AssertLog(ghk != nullptr);
    pGHKcurrs.erase(ghk->getID());
    pGHKcurrs.emplace(n, ghk);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_checkGHKcurrID(std::string const& id) const {
    checkID(id);

    ArgErrLogIf(pGHKcurrs.find(id) != pGHKcurrs.end(), "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrAdd(GHKcurr* ghkcurr) {
    AssertLog(ghkcurr->getSurfsys() == this);
    _checkGHKcurrID(ghkcurr->getID());
    pGHKcurrs.emplace(ghkcurr->getID(), ghkcurr);
}

////////////////////////////////////////////////////////////////////////////////

void Surfsys::_handleGHKcurrDel(GHKcurr* ghkcurr) {
    AssertLog(ghkcurr->getSurfsys() == this);
    pGHKcurrs.erase(ghkcurr->getID());
}

////////////////////////////////////////////////////////////////////////////////

SReac* Surfsys::_getSReac(uint lidx) const {
    AssertLog(lidx < pSReacs.size());
    auto sr_it = pSReacs.begin();
    std::advance(sr_it, lidx);
    return sr_it->second;
}

////////////////////////////////////////////////////////////////////////////////

RaftGen* Surfsys::_getRaftGen(uint lidx) const {
    AssertLog(lidx < pRaftGens.size());
    auto rg_it = pRaftGens.begin();
    std::advance(rg_it, lidx);
    return rg_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Endocytosis* Surfsys::_getEndocytosis(uint lidx) const {
    AssertLog(lidx < pEndocytosis.size());
    auto endo_it = pEndocytosis.begin();
    std::advance(endo_it, lidx);
    return endo_it->second;
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac* Surfsys::_getVDepSReac(uint lidx) const {
    AssertLog(lidx < pVDepSReacs.size());
    auto vd_it = pVDepSReacs.begin();
    std::advance(vd_it, lidx);
    return vd_it->second;
}

////////////////////////////////////////////////////////////////////////////////

OhmicCurr* Surfsys::_getOhmicCurr(uint lidx) const {
    AssertLog(lidx < pOhmicCurrs.size());
    auto oc_it = pOhmicCurrs.begin();
    std::advance(oc_it, lidx);
    return oc_it->second;
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr* Surfsys::_getGHKcurr(uint lidx) const {
    AssertLog(lidx < pGHKcurrs.size());
    auto ghk_it = pGHKcurrs.begin();
    std::advance(ghk_it, lidx);
    return ghk_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Diff* Surfsys::_getDiff(uint lidx) const {
    AssertLog(lidx < pDiffs.size());
    auto df_it = pDiffs.begin();
    std::advance(df_it, lidx);
    return df_it->second;
}

}  // namespace steps::model
