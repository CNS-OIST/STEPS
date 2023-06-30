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


#include "volsys.hpp"

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cassert>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include <easylogging++.h>

#include "diff.hpp"
#include "model.hpp"
#include "model/linkspec.hpp"
#include "model/vesbind.hpp"
#include "model/vesicle.hpp"
#include "model/vesunbind.hpp"
#include "reac.hpp"
#include "spec.hpp"

#include "util/checkid.hpp"
#include "util/error.hpp"

namespace steps::model {

using util::checkID;

////////////////////////////////////////////////////////////////////////////////

Volsys::Volsys(std::string const& id, Model* model)
    : pID(id)
    , pModel(model) {
    ArgErrLogIf(pModel == nullptr, "No model provided to Volsys initializer function");

    pModel->_handleVolsysAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Volsys::~Volsys() {
    if (pModel == nullptr) {
        return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleSelfDelete() {
    for (auto const& reac: getAllReacs()) {
        delete reac;
    }

    for (auto const& diff: getAllDiffs()) {
        delete diff;
    }

    for (auto const& vb: getAllVesBinds()) {
        delete vb;
    }

    for (auto const& vub: getAllVesUnbinds()) {
        delete vub;
    }

    pModel->_handleVolsysDel(this);
    pReacs.clear();
    pDiffs.clear();
    pVesBinds.clear();
    pVesUnbinds.clear();
    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::setID(std::string const& id) {
    AssertLog(pModel != nullptr);
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel->_handleVolsysIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

Reac* Volsys::getReac(std::string const& id) const {
    auto reac = pReacs.find(id);

    ArgErrLogIf(reac == pReacs.end(), "Model does not contain reaction with name '" + id + "'");

    AssertLog(reac->second != nullptr);
    return reac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::delReac(std::string const& id) const {
    auto reac = getReac(id);
    // Delete reac object since it is owned by c++, not python
    delete reac;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Reac*> Volsys::getAllReacs() const {
    ReacPVec reacs;
    reacs.reserve(pReacs.size());
    for (auto const& r: pReacs) {
        reacs.push_back(r.second);
    }
    return reacs;
}

////////////////////////////////////////////////////////////////////////////////

VesBind* Volsys::getVesBind(std::string const& id) const {
    auto vesbind = pVesBinds.find(id);
    if (vesbind == pVesBinds.end()) {
        std::ostringstream os;
        os << "Model does not contain vesicle binding reaction with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(vesbind->second != nullptr);
    return vesbind->second;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::delVesBind(std::string const& id) const {
    VesBind* vesbind = getVesBind(id);
    // Delete VesBind object since it is owned by c++, not python
    delete vesbind;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<VesBind*> Volsys::getAllVesBinds() const {
    VesBindPVec vesbinds;
    vesbinds.reserve(pVesBinds.size());
    for (auto const& vb: pVesBinds) {
        vesbinds.push_back(vb.second);
    }
    return vesbinds;
}

////////////////////////////////////////////////////////////////////////////////

VesUnbind* Volsys::getVesUnbind(std::string const& id) const {
    auto vesunbind = pVesUnbinds.find(id);
    if (vesunbind == pVesUnbinds.end()) {
        std::ostringstream os;
        os << "Model does not contain vesicle unbinding reaction with name '" << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(vesunbind->second != nullptr);
    return vesunbind->second;
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::delVesUnbind(std::string const& id) const {
    VesUnbind* vesunbind = getVesUnbind(id);
    // Delete VesUnbind object since it is owned by c++, not python
    delete vesunbind;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<VesUnbind*> Volsys::getAllVesUnbinds() const {
    VesUnbindPVec vesunbinds;
    vesunbinds.reserve(pVesUnbinds.size());
    for (auto const& vub: pVesUnbinds) {
        vesunbinds.push_back(vub.second);
    }
    return vesunbinds;
}

////////////////////////////////////////////////////////////////////////////////

Diff* Volsys::getDiff(std::string const& id) const {
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

void Volsys::delDiff(std::string const& id) const {
    auto diff = getDiff(id);
    // delete diff object since it is owned by c++, not python
    delete diff;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Diff*> Volsys::getAllDiffs() const {
    DiffPVec diffs;
    diffs.reserve(pDiffs.size());
    for (auto const& d: pDiffs) {
        diffs.push_back(d.second);
    }
    return diffs;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec*> Volsys::getAllSpecs() const {
    std::set<Spec*> specs_set;

    for (auto const& reac: getAllReacs()) {
        const auto& specs = reac->getAllSpecs();
        specs_set.insert(specs.begin(), specs.end());
    }

    for (auto const& diff: getAllDiffs()) {
        const auto& specs = diff->getAllSpecs();
        specs_set.insert(specs.begin(), specs.end());
    }

    return {specs_set.begin(), specs_set.end()};
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_checkReacID(std::string const& id) const {
    checkID(id);

    ArgErrLogIf(pReacs.find(id) != pReacs.end(), "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleReacIDChange(std::string const& o, std::string const& n) {
    auto r_old = pReacs.find(o);
    AssertLog(r_old != pReacs.end());

    if (o == n) {
        return;
    }
    _checkReacID(n);

    Reac* r = r_old->second;
    AssertLog(r != nullptr);
    pReacs.erase(r->getID());
    pReacs.emplace(n, r);
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleReacAdd(Reac* reac) {
    AssertLog(reac->getVolsys() == this);
    _checkReacID(reac->getID());
    pReacs.emplace(reac->getID(), reac);
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleReacDel(Reac* reac) {
    AssertLog(reac->getVolsys() == this);
    pReacs.erase(reac->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_checkDiffID(std::string const& id) const {
    checkID(id);

    ArgErrLogIf(pDiffs.find(id) != pDiffs.end(), "'" + id + "' is already in use");
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleDiffIDChange(std::string const& o, std::string const& n) {
    auto d_old = pDiffs.find(o);
    AssertLog(d_old != pDiffs.end());

    if (o == n) {
        return;
    }
    _checkDiffID(n);

    Diff* d = d_old->second;
    AssertLog(d != nullptr);
    pDiffs.erase(d->getID());
    pDiffs.emplace(n, d);
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleDiffAdd(Diff* diff) {
    AssertLog(diff->getVolsys() == this);
    _checkDiffID(diff->getID());
    pDiffs.emplace(diff->getID(), diff);
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleDiffDel(Diff* diff) {
    AssertLog(diff->getVolsys() == this);
    pDiffs.erase(diff->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_checkVesBindID(std::string const& id) const {
    checkID(id);
    if (pVesBinds.find(id) != pVesBinds.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleVesBindIDChange(std::string const& o, std::string const& n) {
    auto vb_old = pVesBinds.find(o);
    AssertLog(vb_old != pVesBinds.end());

    if (o == n) {
        return;
    }
    _checkVesBindID(n);

    VesBind* vb = vb_old->second;
    AssertLog(vb != nullptr);
    pVesBinds.erase(vb->getID());
    pVesBinds.emplace(n, vb);
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleVesBindAdd(VesBind* vesbind) {
    AssertLog(vesbind->getVolsys() == this);
    _checkVesBindID(vesbind->getID());
    pVesBinds.emplace(vesbind->getID(), vesbind);
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleVesBindDel(VesBind* vesbind) {
    AssertLog(vesbind->getVolsys() == this);
    pVesBinds.erase(vesbind->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_checkVesUnbindID(std::string const& id) const {
    checkID(id);
    if (pVesUnbinds.find(id) != pVesUnbinds.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleVesUnbindIDChange(std::string const& o, std::string const& n) {
    auto vub_old = pVesUnbinds.find(o);
    AssertLog(vub_old != pVesUnbinds.end());

    if (o == n) {
        return;
    }
    _checkVesUnbindID(n);

    VesUnbind* vub = vub_old->second;
    AssertLog(vub != nullptr);
    pVesUnbinds.erase(vub->getID());
    pVesUnbinds.emplace(n, vub);
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleVesUnbindAdd(VesUnbind* vesunbind) {
    AssertLog(vesunbind->getVolsys() == this);
    _checkVesUnbindID(vesunbind->getID());
    pVesUnbinds.emplace(vesunbind->getID(), vesunbind);
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleVesUnbindDel(VesUnbind* vesunbind) {
    AssertLog(vesunbind->getVolsys() == this);
    pVesUnbinds.erase(vesunbind->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleSpecDelete(Spec* spec) {
    {
        std::vector<std::string> reacs_del;
        for (auto const& reac: pReacs) {
            for (auto const& r_spec: reac.second->getAllSpecs()) {
                if (r_spec == spec) {
                    reacs_del.push_back(reac.second->getID());
                    break;
                }
            }
        }
        for (auto const& r_del: reacs_del) {
            delReac(r_del);
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

    {
        std::vector<std::string> vesbinds_del;
        for (auto const& vesbind: pVesBinds) {
            for (auto const& vb_spec: vesbind.second->getAllSpecs()) {
                if (vb_spec == spec) {
                    vesbinds_del.push_back(vesbind.second->getID());
                    break;
                }
            }
        }
        for (auto const& vb_del: vesbinds_del) {
            delVesBind(vb_del);
        }
    }

    {
        std::vector<std::string> vesunbinds_del;
        for (auto const& vesunbind: pVesUnbinds) {
            for (auto const& vub_spec: vesunbind.second->getAllSpecs()) {
                if (vub_spec == spec) {
                    vesunbinds_del.push_back(vesunbind.second->getID());
                    break;
                }
            }
        }
        for (auto const& vub_del: vesunbinds_del) {
            delVesBind(vub_del);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Volsys::_handleLinkSpecDelete(LinkSpec* spec) {
    {
        std::vector<std::string> vesbinds_del;
        for (auto const& vesbind: pVesBinds) {
            if (vesbind.second->getProducts1().second == spec ||
                vesbind.second->getProducts2().second == spec) {
                vesbinds_del.push_back(vesbind.second->getID());
                break;
            }
        }

        for (auto const& vb_del: vesbinds_del) {
            delVesBind(vb_del);
        }
    }

    {
        std::vector<std::string> vesunbinds_del;
        ;
        for (auto const& vesunbind: pVesUnbinds) {
            if (vesunbind.second->getLinks1().second == spec ||
                vesunbind.second->getLinks2().second == spec) {
                vesunbinds_del.push_back(vesunbind.second->getID());
                break;
            }
        }

        for (auto const& vub_del: vesunbinds_del) {
            delVesUnbind(vub_del);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////

Reac* Volsys::_getReac(uint lidx) const {
    AssertLog(lidx < pReacs.size());
    auto rc_it = pReacs.begin();
    std::advance(rc_it, lidx);
    return rc_it->second;
}

////////////////////////////////////////////////////////////////////////////////

Diff* Volsys::_getDiff(uint lidx) const {
    AssertLog(lidx < pDiffs.size());
    auto df_it = pDiffs.begin();
    std::advance(df_it, lidx);
    return df_it->second;
}

////////////////////////////////////////////////////////////////////////////////

VesBind* Volsys::_getVesBind(uint lidx) const {
    AssertLog(lidx < pVesBinds.size());
    auto vbc_it = pVesBinds.begin();
    std::advance(vbc_it, lidx);
    return vbc_it->second;
}

////////////////////////////////////////////////////////////////////////////////

VesUnbind* Volsys::_getVesUnbind(uint lidx) const {
    AssertLog(lidx < pVesUnbinds.size());
    auto vubc_it = pVesUnbinds.begin();
    std::advance(vubc_it, lidx);
    return vubc_it->second;
}

}  // namespace steps::model
