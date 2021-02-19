/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

#include "steps/model/model.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace wm {

Comp::Comp(std::string id, Geom * container, double vol)
: pVol(vol)
, pID(std::move(id))
, pContainer(container)
{
    if (pContainer == nullptr)
    {
        ArgErrLog("No container provided to Comp initializer function.");
    }

    if (pVol < 0.0)
    {
        ArgErrLog("Compartment volume can't be negative.");
    }
    pContainer->_handleCompAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Comp::~Comp()
{
    if (pContainer == nullptr) {
      return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setID(std::string const & id)
{
    AssertLog(pContainer != nullptr);
    if (id == pID) {
      return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pContainer->_handleCompIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setVol(double vol)
{
    AssertLog(pContainer != nullptr);
    if (vol < 0.0)
    {
        ArgErrLog("Compartment volume can't be negative.");
    }
    pVol = vol;
}

////////////////////////////////////////////////////////////////////////////////

void Comp::addVolsys(std::string const & id)
{
    // string identifier is only added to set if it is not already included
    pVolsys.insert(id);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::delVolsys(std::string const & id)
{
    // string identifier is only removed from set if it is included
    pVolsys.erase(id);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::Spec*> Comp::getAllSpecs(steps::model::Model* model) const
{
    std::set<steps::model::Spec*> pSpecs;
    std::set<std::string>::iterator it;
    for (it = pVolsys.begin(); it != pVolsys.end(); it++) {
        steps::model::Volsys* volsys = model->getVolsys(*it);
        std::vector<steps::model::Spec*> specs = volsys->getAllSpecs();
        pSpecs.insert(specs.begin(), specs.end());
    }

    return {pSpecs.begin(), pSpecs.end()};
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::Reac*> Comp::getAllReacs(steps::model::Model* model) const
{
    std::set<steps::model::Reac*> pReacs;
    std::set<std::string>::iterator it;
    for (it = pVolsys.begin(); it != pVolsys.end(); it++) {
        steps::model::Volsys* volsys = model->getVolsys(*it);
        std::vector<steps::model::Reac*> reacs = volsys->getAllReacs();
        pReacs.insert(reacs.begin(), reacs.end());
    }

    return {pReacs.begin(), pReacs.end()};
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::Diff*> Comp::getAllDiffs(steps::model::Model* model) const
{
    std::set<steps::model::Diff*> pDiffs;
    std::set<std::string>::iterator it;
    for (it = pVolsys.begin(); it != pVolsys.end(); it++) {
        steps::model::Volsys* volsys = model->getVolsys(*it);
        std::vector<steps::model::Diff*> diffs = volsys->getAllDiffs();
        pDiffs.insert(diffs.begin(), diffs.end());
    }

    return {pDiffs.begin(), pDiffs.end()};
}


////////////////////////////////////////////////////////////////////////////////

void Comp::_addIPatch(Patch * patch)
{
    AssertLog(patch->getOComp() == this);
    // patch pointer is only added to set if it is not already included
    pIPatches.insert(patch);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::_delIPatch(Patch * patch)
{
    AssertLog(patch->getOComp() == this);
    pIPatches.erase(patch);
}
////////////////////////////////////////////////////////////////////////////////

void Comp::_addOPatch(Patch * patch)
{
    AssertLog(patch->getIComp() == this);
    // patch pointer is only added to set if it is not already included
    pOPatches.insert(patch);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::_delOPatch(Patch * patch)
{
    AssertLog(patch->getIComp() == this);
    pOPatches.erase(patch);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::_handleSelfDelete()
{
    pContainer->_handleCompDel(this);
    pVol = 0.0;
    pVolsys.clear();
    pIPatches.clear();
    pOPatches.clear();
    pContainer = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

steps::wm::Patch * Comp::_getIPatch(uint lidx) const
{
    AssertLog(lidx < pIPatches.size());
    auto pit = pIPatches.begin();
    std::advance(pit, lidx);
    return *pit;
}

////////////////////////////////////////////////////////////////////////////////

steps::wm::Patch * Comp::_getOPatch(uint lidx) const
{
    AssertLog(lidx < pOPatches.size());
    auto pit = pOPatches.begin();
    std::advance(pit, lidx);
    return *pit;
}

} // namespace wm
} // namespace steps

