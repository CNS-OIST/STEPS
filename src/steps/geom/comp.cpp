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

#include "steps/model/model.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace swm = steps::wm;

////////////////////////////////////////////////////////////////////////////////

swm::Comp::Comp(std::string id, swm::Geom * container, double vol)
: pID(std::move(id))
, pContainer(container)
, pVolsys()
, pVol(vol)
, pIPatches()
, pOPatches()
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

swm::Comp::~Comp()
{
    if (pContainer == nullptr) {
      return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::setID(std::string const & id)
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

void swm::Comp::setVol(double vol)
{
    AssertLog(pContainer != nullptr);
    if (vol < 0.0)
    {
        ArgErrLog("Compartment volume can't be negative.");
    }
    pVol = vol;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::addVolsys(std::string const & id)
{
    // string identifier is only added to set if it is not already included
    pVolsys.insert(id);
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::delVolsys(std::string const & id)
{
    // string identifier is only removed from set if it is included
    pVolsys.erase(id);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::Spec*> swm::Comp::getAllSpecs(steps::model::Model* model)
{
    std::set<steps::model::Spec*> pSpecs;
    std::set<std::string>::iterator it;
    for (it = pVolsys.begin(); it != pVolsys.end(); it++) {
        steps::model::Volsys* volsys = model->getVolsys(*it);
        std::vector<steps::model::Spec*> specs = volsys->getAllSpecs();
        pSpecs.insert(specs.begin(), specs.end());
    }

    std::vector<steps::model::Spec*> spec_vec(pSpecs.begin(), pSpecs.end());
    return spec_vec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::Reac*> swm::Comp::getAllReacs(steps::model::Model* model)
{
    std::set<steps::model::Reac*> pReacs;
    std::set<std::string>::iterator it;
    for (it = pVolsys.begin(); it != pVolsys.end(); it++) {
        steps::model::Volsys* volsys = model->getVolsys(*it);
        std::vector<steps::model::Reac*> reacs = volsys->getAllReacs();
        pReacs.insert(reacs.begin(), reacs.end());
    }

    std::vector<steps::model::Reac*> reac_vec(pReacs.begin(), pReacs.end());
    return reac_vec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::Diff*> swm::Comp::getAllDiffs(steps::model::Model* model)
{
    std::set<steps::model::Diff*> pDiffs;
    std::set<std::string>::iterator it;
    for (it = pVolsys.begin(); it != pVolsys.end(); it++) {
        steps::model::Volsys* volsys = model->getVolsys(*it);
        std::vector<steps::model::Diff*> diffs = volsys->getAllDiffs();
        pDiffs.insert(diffs.begin(), diffs.end());
    }

    std::vector<steps::model::Diff*> diff_vec(pDiffs.begin(), pDiffs.end());
    return diff_vec;
}


////////////////////////////////////////////////////////////////////////////////

void swm::Comp::_addIPatch(swm::Patch * patch)
{
    AssertLog(patch->getOComp() == this);
    // patch pointer is only added to set if it is not already included
    pIPatches.insert(patch);
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::_delIPatch(swm::Patch * patch)
{
    AssertLog(patch->getOComp() == this);
    pIPatches.erase(patch);
}
////////////////////////////////////////////////////////////////////////////////

void swm::Comp::_addOPatch(swm::Patch * patch)
{
    AssertLog(patch->getIComp() == this);
    // patch pointer is only added to set if it is not already included
    pOPatches.insert(patch);
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::_delOPatch(swm::Patch * patch)
{
    AssertLog(patch->getIComp() == this);
    pOPatches.erase(patch);
}

////////////////////////////////////////////////////////////////////////////////

void swm::Comp::_handleSelfDelete()
{
    pContainer->_handleCompDel(this);
    pVol = 0.0;
    pVolsys.clear();
    pIPatches.clear();
    pOPatches.clear();
    pContainer = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

steps::wm::Patch * swm::Comp::_getIPatch(uint lidx) const
{
    AssertLog(lidx < pIPatches.size());
    auto pit = pIPatches.begin();
    std::advance(pit, lidx);
    return *pit;
}

////////////////////////////////////////////////////////////////////////////////

steps::wm::Patch * swm::Comp::_getOPatch(uint lidx) const
{
    AssertLog(lidx < pOPatches.size());
    auto pit = pOPatches.begin();
    std::advance(pit, lidx);
    return *pit;
}

////////////////////////////////////////////////////////////////////////////////

/// END
