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

#include "raftsys.hpp"

#include "model.hpp"
#include "raftdis.hpp"
#include "raftendocytosis.hpp"
#include "raftsreac.hpp"
#include "spec.hpp"

#include "util/checkid.hpp"
#include "util/error.hpp"

namespace steps::model {

using util::checkID;

////////////////////////////////////////////////////////////////////////////////

Raftsys::Raftsys(std::string const& id, Model& model)
    : pID(id)
    , pModel(model) {
    pModel._handleRaftsysAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

Raftsys::~Raftsys() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::setID(std::string const& id) {
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel._handleRaftsysIDChange(pID, id);

    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_handleSelfDelete() {
    for (auto raftsreac: getAllRaftSReacs()) {
        delete raftsreac;
    }

    for (auto raftendo: getAllRaftEndocytosiss()) {
        delete raftendo;
    }

    for (auto raftdis: getAllRaftDiss()) {
        delete raftdis;
    }

    pModel._handleRaftsysDel(*this);
}

////////////////////////////////////////////////////////////////////////////////

RaftSReac& Raftsys::getRaftSReac(std::string const& id) const {
    auto raftsreac = pRaftSReacs.find(id);
    if (raftsreac == pRaftSReacs.end()) {
        std::ostringstream os;
        os << "Model does not contain raft surface "
              "reaction with name '"
           << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(raftsreac->second != nullptr);
    return *raftsreac->second;
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::delRaftSReac(std::string const& id) const {
    RaftSReac& raftsreac = getRaftSReac(id);
    delete &raftsreac;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<RaftSReac*> Raftsys::getAllRaftSReacs() const {
    std::vector<RaftSReac*> raftsreacs;
    raftsreacs.reserve(pRaftSReacs.size());
    for (const auto& [_, rsr]: pRaftSReacs) {
        raftsreacs.push_back(rsr);
    }
    return raftsreacs;
}

////////////////////////////////////////////////////////////////////////////////

util::flat_set<Spec*> Raftsys::getAllSpecs() const {
    util::flat_set<Spec*> specs;

    for (auto const& rsreac: getAllRaftSReacs()) {
        const auto& rsr_specs = rsreac->getAllSpecs();
        specs.insert(rsr_specs.begin(), rsr_specs.end());
    }

    for (auto const& raftendo: getAllRaftEndocytosiss()) {
        const auto& raftendo_specs = raftendo->getAllSpecs();
        specs.insert(raftendo_specs.begin(), raftendo_specs.end());
    }

    for (auto const& raftdis: getAllRaftDiss()) {
        const auto& rd_specs = raftdis->getAllSpecs();
        specs.insert(rd_specs.begin(), rd_specs.end());
    }

    return specs;
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_checkRaftSReacID(std::string const& id) const {
    checkID(id);
    if (pRaftSReacs.find(id) != pRaftSReacs.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_handleRaftSReacIDChange(std::string const& o, std::string const& n) {
    auto rsr_old = pRaftSReacs.find(o);
    AssertLog(rsr_old != pRaftSReacs.end());

    if (o == n) {
        return;
    }
    _checkRaftSReacID(n);

    RaftSReac* rsr = rsr_old->second;
    AssertLog(rsr != nullptr);
    pRaftSReacs.erase(rsr_old);
    pRaftSReacs.emplace(n, rsr);
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_handleRaftSReacAdd(RaftSReac& rsreac) {
    AssertLog(&rsreac.getRaftsys() == this);
    _checkRaftSReacID(rsreac.getID());
    pRaftSReacs.emplace(rsreac.getID(), &rsreac);
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_handleRaftSReacDel(RaftSReac& rsreac) {
    AssertLog(&rsreac.getRaftsys() == this);
    pRaftSReacs.erase(rsreac.getID());
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_handleSpecDelete(Spec& spec) {
    std::vector<std::string> rsreacs_del;

    for (auto const& rsreac: pRaftSReacs) {
        for (auto const& rsr_spec: rsreac.second->getAllSpecs()) {
            if (rsr_spec == &spec) {
                rsreacs_del.push_back(rsreac.second->getID());
                break;
            }
        }
    }

    for (auto rsr_del: rsreacs_del) {
        delRaftSReac(rsr_del);
    }

    std::vector<std::string> raftendos_del;

    for (auto const& raftendo: pRaftEndocytosis) {
        for (auto const& e_spec: raftendo.second->getAllSpecs()) {
            if (e_spec == &spec) {
                raftendos_del.push_back(raftendo.second->getID());
                break;
            }
        }
    }

    for (auto e_del: raftendos_del) {
        delRaftEndocytosis(e_del);
    }

    std::vector<std::string> raftdiss_del;

    for (auto const& raftdis: pRaftDiss) {
        for (auto const& rd_spec: raftdis.second->getAllSpecs()) {
            if (rd_spec == &spec) {
                raftdiss_del.push_back(raftdis.second->getID());
                break;
            }
        }
    }

    for (auto rd_del: raftdiss_del) {
        delRaftDis(rd_del);
    }
}

////////////////////////////////////////////////////////////////////////////////

RaftSReac& Raftsys::_getRaftSReac(uint lidx) const {
    AssertLog(lidx < pRaftSReacs.size());
    auto rsr_it = pRaftSReacs.begin();
    std::advance(rsr_it, lidx);
    return *rsr_it->second;
}

////////////////////////////////////////////////////////////////////////////////

RaftEndocytosis& Raftsys::getRaftEndocytosis(std::string const& id) const {
    auto raftendo = pRaftEndocytosis.find(id);
    if (raftendo == pRaftEndocytosis.end()) {
        std::ostringstream os;
        os << "Model does not contain raftendocytotic "
              "reaction with name '"
           << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(raftendo->second != nullptr);
    return *raftendo->second;
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::delRaftEndocytosis(std::string const& id) const {
    RaftEndocytosis& raftendo = getRaftEndocytosis(id);
    delete &raftendo;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<RaftEndocytosis*> Raftsys::getAllRaftEndocytosiss() const {
    std::vector<RaftEndocytosis*> raftendos;

    for (auto const& raftendo: pRaftEndocytosis) {
        raftendos.push_back(raftendo.second);
    }
    return raftendos;
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_checkRaftEndocytosisID(std::string const& id) const {
    checkID(id);
    if (pRaftEndocytosis.find(id) != pRaftEndocytosis.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_handleRaftEndocytosisIDChange(std::string const& o, std::string const& n) {
    auto e_old = pRaftEndocytosis.find(o);
    AssertLog(e_old != pRaftEndocytosis.end());

    if (o == n) {
        return;
    }
    _checkRaftEndocytosisID(n);

    RaftEndocytosis* e = e_old->second;
    AssertLog(e != nullptr);
    pRaftEndocytosis.erase(e_old);
    pRaftEndocytosis.emplace(n, e);
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_handleRaftEndocytosisAdd(RaftEndocytosis& raftendocyt) {
    AssertLog(&raftendocyt.getRaftsys() == this);
    _checkRaftEndocytosisID(raftendocyt.getID());
    pRaftEndocytosis.emplace(raftendocyt.getID(), &raftendocyt);
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_handleRaftEndocytosisDel(RaftEndocytosis& raftendocyt) {
    AssertLog(&raftendocyt.getRaftsys() == this);
    pRaftEndocytosis.erase(raftendocyt.getID());
}

////////////////////////////////////////////////////////////////////////////////

RaftEndocytosis& Raftsys::_getRaftEndocytosis(uint lidx) const {
    AssertLog(lidx < pRaftEndocytosis.size());
    auto raftendo_it = pRaftEndocytosis.begin();
    std::advance(raftendo_it, lidx);
    return *raftendo_it->second;
}

////////////////////////////////////////////////////////////////////////////////

RaftDis& Raftsys::getRaftDis(std::string const& id) const {
    auto raftdis = pRaftDiss.find(id);
    if (raftdis == pRaftDiss.end()) {
        std::ostringstream os;
        os << "Model does not contain raft "
              "dissolution with name '"
           << id << "'";
        ArgErrLog(os.str());
    }
    AssertLog(raftdis->second != nullptr);
    return *raftdis->second;
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::delRaftDis(std::string const& id) const {
    RaftDis& raftdis = getRaftDis(id);
    delete &raftdis;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<RaftDis*> Raftsys::getAllRaftDiss() const {
    std::vector<RaftDis*> raftdiss;
    raftdiss.reserve(pRaftDiss.size());

    for (auto const& rd: pRaftDiss) {
        raftdiss.push_back(rd.second);
    }
    return raftdiss;
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_checkRaftDisID(std::string const& id) const {
    checkID(id);
    if (pRaftDiss.find(id) != pRaftDiss.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_handleRaftDisIDChange(std::string const& o, std::string const& n) {
    auto rd_old = pRaftDiss.find(o);
    AssertLog(rd_old != pRaftDiss.end());

    if (o == n) {
        return;
    }
    _checkRaftDisID(n);

    RaftDis* rd = rd_old->second;
    AssertLog(rd != nullptr);
    pRaftDiss.erase(rd_old);
    pRaftDiss.emplace(n, rd);
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_handleRaftDisAdd(RaftDis& raftdis) {
    AssertLog(&raftdis.getRaftsys() == this);
    _checkRaftDisID(raftdis.getID());
    pRaftDiss.emplace(raftdis.getID(), &raftdis);
}

////////////////////////////////////////////////////////////////////////////////

void Raftsys::_handleRaftDisDel(RaftDis& raftdis) {
    AssertLog(&raftdis.getRaftsys() == this);
    pRaftDiss.erase(raftdis.getID());
}

////////////////////////////////////////////////////////////////////////////////

RaftDis& Raftsys::_getRaftDis(uint lidx) const {
    AssertLog(lidx < pRaftDiss.size());
    auto rd_it = pRaftDiss.begin();
    std::advance(rd_it, lidx);
    return *rd_it->second;
}

}  // namespace steps::model
