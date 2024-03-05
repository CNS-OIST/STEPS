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

#include "geom.hpp"

#include "comp.hpp"
#include "patch.hpp"
#include "util/checkid.hpp"
#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::wm {

using util::checkID;

////////////////////////////////////////////////////////////////////////////////

Geom::~Geom() {
    while (!pComps.empty()) {
        auto comp = pComps.begin();
        delete comp->second;
    }
    while (!pPatches.empty()) {
        auto patch = pPatches.begin();
        delete patch->second;
    }
}

////////////////////////////////////////////////////////////////////////////////

Comp& Geom::getComp(std::string const& id) const {
    auto comp = pComps.find(id);
    if (comp == pComps.end()) {
        std::ostringstream os;
        os << "Container does not contain compartment with name '" << id << "'\n";
        ArgErrLog(os.str());
    }
    AssertLog(comp->second != nullptr);
    return *comp->second;
}

////////////////////////////////////////////////////////////////////////////////

void Geom::delComp(std::string const& id) const {
    Comp& comp = getComp(id);
    delete &comp;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Comp*> Geom::getAllComps() const {
    std::vector<Comp*> comps;
    comps.reserve(pComps.size());

    for (const auto& c: pComps) {
        comps.push_back(c.second);
    }
    return comps;
}

////////////////////////////////////////////////////////////////////////////////

Patch& Geom::getPatch(std::string const& id) const {
    auto patch = pPatches.find(id);

    if (patch == pPatches.end()) {
        std::ostringstream os;
        os << "Container does not contain patch with name '" << id << "'\n";
        ArgErrLog(os.str());
    }
    AssertLog(patch->second != nullptr);
    return *patch->second;
}

////////////////////////////////////////////////////////////////////////////////

void Geom::delPatch(std::string const& id) const {
    Patch& patch = getPatch(id);
    delete &patch;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Patch*> Geom::getAllPatches() const {
    PatchPVec patches;
    patches.reserve(pPatches.size());

    for (const auto& p: pPatches) {
        patches.push_back(p.second);
    }
    return patches;
}

////////////////////////////////////////////////////////////////////////////////

void Geom::_checkCompID(std::string const& id) const {
    checkID(id);
    if (pComps.find(id) != pComps.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use.\n";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Geom::_handleCompIDChange(std::string const& o, std::string const& n) {
    auto c_old = pComps.find(o);
    AssertLog(c_old != pComps.end());

    if (o == n) {
        return;
    }
    _checkCompID(n);

    Comp* c = c_old->second;
    AssertLog(c != nullptr);
    pComps.erase(c_old);
    pComps.emplace(n, c);
}

////////////////////////////////////////////////////////////////////////////////

void Geom::_handleCompAdd(Comp& comp) {
    AssertLog(&comp.getContainer() == this);
    _checkCompID(comp.getID());
    pComps.emplace(comp.getID(), &comp);
}

////////////////////////////////////////////////////////////////////////////////

void Geom::_handleCompDel(Comp& comp) {
    AssertLog(&comp.getContainer() == this);
    pComps.erase(comp.getID());
}

////////////////////////////////////////////////////////////////////////////////

void Geom::_checkPatchID(std::string const& id) const {
    checkID(id);
    if (pPatches.find(id) != pPatches.end()) {
        std::ostringstream os;
        os << "'" << id << "' is already in use.\n";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Geom::_handlePatchIDChange(std::string const& o, std::string const& n) {
    if (o == n) {
        return;
    }

    auto p_old = pPatches.find(o);
    AssertLog(p_old != pPatches.end());

    _checkPatchID(n);

    Patch* p = p_old->second;
    AssertLog(p != nullptr);
    pPatches.erase(p_old);
    pPatches.emplace(n, p);
}

////////////////////////////////////////////////////////////////////////////////

void Geom::_handlePatchAdd(Patch& patch) {
    AssertLog(&patch.getContainer() == this);
    checkID(patch.getID());
    if (!pPatches.emplace(patch.getID(), &patch).second) {
        std::ostringstream os;
        os << "'" << patch.getID() << "' is already in use.\n";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void Geom::_handlePatchDel(Patch& patch) {
    AssertLog(&patch.getContainer() == this);
    pPatches.erase(patch.getID());
}

////////////////////////////////////////////////////////////////////////////////

wm::Comp& Geom::_getComp(solver::comp_global_id gidx) const {
    AssertLog(gidx.get() < pComps.size());
    auto cp_it = pComps.begin();
    std::advance(cp_it, gidx.get());
    return *cp_it->second;
}

////////////////////////////////////////////////////////////////////////////////

wm::Patch& Geom::_getPatch(solver::patch_global_id gidx) const {
    AssertLog(gidx.get() < pPatches.size());
    auto pt_it = pPatches.begin();
    std::advance(pt_it, gidx.get());
    return *pt_it->second;
}

}  // namespace steps::wm
