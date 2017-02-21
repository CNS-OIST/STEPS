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
#include "steps/geom/patch.hpp"

#include "steps/model/model.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
namespace swm = steps::wm;

////////////////////////////////////////////////////////////////////////////////

swm::Patch::Patch(std::string const & id, swm::Geom * container, swm::Comp* icomp,
        swm::Comp* ocomp, double area)
: pID(id)
, pContainer(container)
, pIComp()
, pOComp()
, pSurfsys()
, pArea(area)
{
    if (pContainer == 0)
    {
        ostringstream os;
        os << "No container provided to Patch initializer function.\n";
        throw steps::ArgErr(os.str());
    }

    _setIComp(icomp);
    if (ocomp != 0) _setOComp(ocomp);

    if (pArea < 0.0)
    {
        ostringstream os;
        os << "Patch area can't be negative.\n";
        throw steps::ArgErr(os.str());
    }
    pContainer->_handlePatchAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

swm::Patch::~Patch(void)
{
    if (pContainer == 0) return;
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void swm::Patch::setID(std::string const & id)
{
    assert(pContainer != 0);
    if (id == pID) return;
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pContainer->_handlePatchIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Patch::setArea(double area)
{
    assert(pContainer != 0);
    if (area < 0.0)
    {
        ostringstream os;
        os << "Patch area can't be negative.\n";
        throw steps::ArgErr(os.str());
    }
    pArea = area;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Patch::addSurfsys(std::string const & id)
{
    // string identifier is only added to set if it is not already included
    pSurfsys.insert(id);
}

////////////////////////////////////////////////////////////////////////////////

void swm::Patch::delSurfsys(std::string const & id)
{
    // string identifier is only removed from set if it is included
    pSurfsys.erase(id);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::Spec*> swm::Patch::getAllSpecs(steps::model::Model* model)
{
    std::set<steps::model::Spec*> pSpecs;
    std::set<std::string>::iterator it;
    for (it = pSurfsys.begin(); it != pSurfsys.end(); it++) {
        steps::model::Surfsys* surfsys = model->getSurfsys(*it);
        std::vector<steps::model::Spec*> specs = surfsys->getAllSpecs();
        pSpecs.insert(specs.begin(), specs.end());
    }

    std::vector<steps::model::Spec*> spec_vec(pSpecs.begin(), pSpecs.end());
    return spec_vec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::SReac*> swm::Patch::getAllSReacs(steps::model::Model* model)
{
    std::set<steps::model::SReac*> pSReacs;
    std::set<std::string>::iterator it;
    for (it = pSurfsys.begin(); it != pSurfsys.end(); it++) {
        steps::model::Surfsys* surfsys = model->getSurfsys(*it);
        std::vector<steps::model::SReac*> sreacs = surfsys->getAllSReacs();
        pSReacs.insert(sreacs.begin(), sreacs.end());
    }

    std::vector<steps::model::SReac*> sreac_vec(pSReacs.begin(), pSReacs.end());
    return sreac_vec;
}

////////////////////////////////////////////////////////////////////////////////

void swm::Patch::_setIComp(swm::Comp* icomp)
{
    if (icomp->getContainer() != pContainer)
    {
        ostringstream os;
        os << "Compartment does not belong to same container as patch.\n";
        throw steps::ArgErr(os.str());
    }
    std::set<swm::Patch *> ipatches  = icomp->getIPatches();
    if (ipatches.find(this) != ipatches.end())
    {
        ostringstream os;
        os << "Patch is already on inside of compartment.\n";
        throw steps::ArgErr(os.str());
    }
    // remove the patch if it was already on the outside of some
    // other compartment
    if (pIComp != 0)
    {
        pIComp->_delOPatch(this);
    }

    pIComp = icomp;
    pIComp->_addOPatch(this);

}

////////////////////////////////////////////////////////////////////////////////

void swm::Patch::_setOComp(swm::Comp* ocomp)
{
    if (ocomp == 0) return;

    if (ocomp->getContainer() != pContainer)
    {
        ostringstream os;
           os << "Compartment does not belong to same container as patch.\n";
           throw steps::ArgErr(os.str());
    }
    std::set<swm::Patch *> opatches  = ocomp->getOPatches();
    if (opatches.find(this) != opatches.end())
    {
           ostringstream os;
          os << "Patch is already on outside of compartment.\n";
           throw steps::ArgErr(os.str());
    }
    // remove the patch if it was already on the inside of some
    // other compartment
    if (pOComp != 0)
    {
        pOComp->_delIPatch(this);
    }

    pOComp = ocomp;
    pOComp->_addIPatch(this);
}

////////////////////////////////////////////////////////////////////////////////

void swm::Patch::_handleSelfDelete(void)
{
    pContainer->_handlePatchDel(this);
    pArea = 0.0;
    pSurfsys.clear();
    pIComp = 0;
    pOComp = 0;
    pContainer = 0;
}

////////////////////////////////////////////////////////////////////////////////

/// END

