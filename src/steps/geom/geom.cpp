/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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
#include "steps/geom/geom.hpp"
#include "steps/geom/comp.hpp"
#include "steps/geom/patch.hpp"
#include "steps/util/checkid.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
namespace swm = steps::wm;

using steps::util::checkID;

////////////////////////////////////////////////////////////////////////////////

swm::Geom::Geom(void)
: pComps()
, pPatches()
{
}

////////////////////////////////////////////////////////////////////////////////

swm::Geom::~Geom(void)
{
    while (pComps.empty() == false)
    {
        CompPMapCI comp = pComps.begin();
        delete(comp->second);
    }
    while(pPatches.empty() == false)
    {
        PatchPMapCI patch = pPatches.begin();
        delete(patch->second);
    }
}

////////////////////////////////////////////////////////////////////////////////

swm::Comp * swm::Geom::getComp(string const & id) const
{
    CompPMapCI comp = pComps.find(id);
    if (comp == pComps.end())
    {
        ostringstream os;
        os << "Container does not contain compartment with name '" << id << "'\n";
        throw steps::ArgErr(os.str());
    }
    assert(comp->second != 0);
    return comp->second;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::delComp(string const & id)
{
    swm::Comp * comp = getComp(id);
    delete(comp);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<swm::Comp *> swm::Geom::getAllComps(void) const
{
    CompPVec comps = CompPVec();
    CompPMapCI comp_end = pComps.end();
    for (CompPMapCI c = pComps.begin(); c != comp_end; ++c)
    {
        comps.push_back(c->second);
    }
    return comps;
}

////////////////////////////////////////////////////////////////////////////////

swm::Patch * swm::Geom::getPatch(string const & id) const
{
    PatchPMapCI patch = pPatches.find(id);
    if (patch == pPatches.end())
    {
        ostringstream os;
        os << "Container does not contain patch with name '" << id << "'\n";
        throw steps::ArgErr(os.str());
    }
    assert(patch->second != 0);
    return patch->second;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::delPatch(string const & id)
{
    swm::Patch * patch = getPatch(id);
    delete(patch);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<swm::Patch *> swm::Geom::getAllPatches(void) const
{
    PatchPVec patches = PatchPVec();
    PatchPMapCI patch_end = pPatches.end();
    for (PatchPMapCI p = pPatches.begin(); p != patch_end; ++p)
    {
        patches.push_back(p->second);
    }
    return patches;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_checkCompID(string const & id) const
{
    checkID(id);
    if (pComps.find(id) != pComps.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use.\n";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handleCompIDChange(string const & o, string const & n)
{
    CompPMapCI c_old = pComps.find(o);
    assert(c_old != pComps.end());

    if (o == n) return;
    _checkCompID(n);

    swm::Comp * c = c_old->second;
    assert(c != 0);
    pComps.erase(c->getID());                        // or s_old->first
    pComps.insert(CompPMap::value_type(n, c));
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handleCompAdd(swm::Comp * comp)
{
    assert(comp->getContainer() == this);
    _checkCompID(comp->getID());
    pComps.insert(CompPMap::value_type(comp->getID(), comp));
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handleCompDel(swm::Comp * comp)
{
    assert(comp->getContainer() == this);
    pComps.erase(comp->getID());
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_checkPatchID(string const & id) const
{
    checkID(id);
    if (pPatches.find(id) != pPatches.end())
    {
        ostringstream os;
        os << "'" << id << "' is already in use.\n";
        throw steps::ArgErr(os.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handlePatchIDChange(string const & o, string const & n)
{
    PatchPMapCI p_old = pPatches.find(o);
    assert(p_old != pPatches.end());

    if (o == n) return;
    _checkPatchID(n);

    swm::Patch * p = p_old->second;
    assert(p != 0);
    pPatches.erase(p->getID());                        // or s_old->first
    pPatches.insert(PatchPMap::value_type(n, p));
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handlePatchAdd(swm::Patch * patch)
{
    assert(patch->getContainer() == this);
    _checkPatchID(patch->getID());
    pPatches.insert(PatchPMap::value_type(patch->getID(), patch));
}

////////////////////////////////////////////////////////////////////////////////

void swm::Geom::_handlePatchDel(swm::Patch * patch)
{
    assert(patch->getContainer() == this);
    pPatches.erase(patch->getID());
}

////////////////////////////////////////////////////////////////////////////////

steps::wm::Comp * swm::Geom::_getComp(uint gidx) const
{
    assert (gidx < pComps.size());
    std::map<std::string, Comp *>::const_iterator cp_it = pComps.begin();
    for (uint i=0; i< gidx; ++i) ++cp_it;
    return cp_it->second;
}

////////////////////////////////////////////////////////////////////////////////

steps::wm::Patch * swm::Geom::_getPatch(uint gidx) const
{
    assert(gidx < pPatches.size());
    std::map<std::string, Patch *>::const_iterator pt_it = pPatches.begin();
    for (uint i=0; i< gidx; ++i) ++pt_it;
    return pt_it->second;
}

////////////////////////////////////////////////////////////////////////////////


// END
