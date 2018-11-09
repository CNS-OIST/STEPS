/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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

// STL headers.
#include <cassert>
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/geom/comp.hpp"
#include "steps/geom/geom.hpp"
#include "steps/geom/patch.hpp"
#include "steps/util/checkid.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace swm = steps::wm;
using steps::util::checkID;

////////////////////////////////////////////////////////////////////////////////

swm::Geom::~Geom()
{
    while (!pComps.empty())
    {
        CompPMapCI comp = pComps.begin();
        delete(comp->second);
    }
    while (!pPatches.empty())
    {
        PatchPMapCI patch = pPatches.begin();
        delete(patch->second);
    }
}

////////////////////////////////////////////////////////////////////////////////

swm::Comp * swm::Geom::getComp(std::string const & id) const
{
    auto comp = pComps.find(id);
    if (comp == pComps.end())
    {
        std::ostringstream os;
        os << "Container does not contain compartment with name '" << id << "'\n";
        ArgErrLog(os.str());
    }
    AssertLog(comp->second != nullptr);
    return comp->second;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::delComp(std::string const & id)
{
    swm::Comp * comp = getComp(id);
    delete(comp);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<swm::Comp *> swm::Geom::getAllComps() const
{
    CompPVec comps = CompPVec();
    for (const auto& c: pComps) {
        comps.push_back(c.second);
    }
    return comps;
}

////////////////////////////////////////////////////////////////////////////////

swm::Patch * swm::Geom::getPatch(std::string const & id) const
{
    auto patch = pPatches.find(id);
    
    if (patch == pPatches.end())
    {
        std::ostringstream os;
        os << "Container does not contain patch with name '" << id << "'\n";
        ArgErrLog(os.str());
    }
    AssertLog(patch->second != nullptr);
    return patch->second;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::delPatch(std::string const & id)
{
    swm::Patch * patch = getPatch(id);
    delete(patch);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<swm::Patch *> swm::Geom::getAllPatches() const
{
    PatchPVec patches = PatchPVec();
    patches.reserve(pPatches.size());

    for (const auto&  p: pPatches) {
        patches.push_back(p.second);
    }
    return patches;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_checkCompID(std::string const & id) const
{
    checkID(id);
    if (pComps.find(id) != pComps.end())
    {
        std::ostringstream os;
        os << "'" << id << "' is already in use.\n";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handleCompIDChange(std::string const & o, std::string const & n)
{
    CompPMapCI c_old = pComps.find(o);
    AssertLog(c_old != pComps.end());

    if (o == n) {
      return;
    }
    _checkCompID(n);

    swm::Comp * c = c_old->second;
    AssertLog(c != nullptr);
    pComps.erase(c->getID());                        // or s_old->first
    pComps.insert(CompPMap::value_type(n, c));
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handleCompAdd(swm::Comp * comp)
{
    AssertLog(comp->getContainer() == this);
    _checkCompID(comp->getID());
    pComps.insert(CompPMap::value_type(comp->getID(), comp));
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handleCompDel(swm::Comp * comp)
{
    AssertLog(comp->getContainer() == this);
    pComps.erase(comp->getID());
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_checkPatchID(std::string const & id) const
{
    checkID(id);
    if (pPatches.find(id) != pPatches.end())
    {
        std::ostringstream os;
        os << "'" << id << "' is already in use.\n";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handlePatchIDChange(std::string const & o, std::string const & n)
{
    if (o == n) {
      return;
    }

    PatchPMapCI p_old = pPatches.find(o);
    AssertLog(p_old != pPatches.end());

    _checkPatchID(n);

    swm::Patch * p = p_old->second;
    AssertLog(p != nullptr);
    pPatches.erase(p->getID());                        // or s_old->first
    pPatches.insert(PatchPMap::value_type(n, p));
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handlePatchAdd(swm::Patch * patch)
{
    AssertLog(patch->getContainer() == this);
    checkID(patch->getID());
    if (!pPatches.insert(std::make_pair(patch->getID(), patch)).second) {
        std::ostringstream os;
        os << "'" << patch->getID() << "' is already in use.\n";
        ArgErrLog(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handlePatchDel(swm::Patch * patch)
{
    AssertLog(patch->getContainer() == this);
    pPatches.erase(patch->getID());
}

////////////////////////////////////////////////////////////////////////////////

steps::wm::Comp * swm::Geom::_getComp(uint gidx) const
{
    AssertLog(gidx < pComps.size());
    auto cp_it = pComps.begin();
    std::advance(cp_it, gidx);
    return cp_it->second;
}

////////////////////////////////////////////////////////////////////////////////

steps::wm::Patch * swm::Geom::_getPatch(uint gidx) const
{
    AssertLog(gidx < pPatches.size());
    auto pt_it = pPatches.begin();
    std::advance(pt_it, gidx);
    return pt_it->second;
}

////////////////////////////////////////////////////////////////////////////////


// END
