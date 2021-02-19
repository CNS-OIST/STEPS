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
#include "steps/geom/patch.hpp"

#include "steps/model/model.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace wm {

Patch::Patch(std::string id, Geom * container, Comp* icomp,
        Comp* ocomp, double area)
: pID(std::move(id))
, pContainer(container)
, pArea(area)
{
    if (pContainer == nullptr)
    {
        std::ostringstream os;
        os << "No container provided to Patch initializer function.\n";
        ArgErrLog(os.str());
    }

    _setIComp(icomp);
    if (ocomp != nullptr) {
        _setOComp(ocomp);
    }

    if (pArea < 0.0)
    {
        std::ostringstream os;
        os << "Patch area can't be negative.\n";
        ArgErrLog(os.str());
    }
    pContainer->_handlePatchAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Patch::~Patch()
{
    if (pContainer == nullptr) {
        return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setID(std::string const & id)
{
    AssertLog(pContainer != nullptr);
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pContainer->_handlePatchIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setArea(double area)
{
    AssertLog(pContainer != nullptr);
    if (area < 0.0)
    {
        std::ostringstream os;
        os << "Patch area can't be negative.\n";
        ArgErrLog(os.str());
    }
    pArea = area;
}

////////////////////////////////////////////////////////////////////////////////

void Patch::addSurfsys(std::string const & id)
{
    // string identifier is only added to set if it is not already included
    pSurfsys.insert(id);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::delSurfsys(std::string const & id)
{
    // string identifier is only removed from set if it is included
    pSurfsys.erase(id);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::Spec*> Patch::getAllSpecs(steps::model::Model* model) const
{
    std::set<steps::model::Spec*> pSpecs;
    for (const auto& id : pSurfsys) {
        steps::model::Surfsys* surfsys = model->getSurfsys(id);
        std::vector<steps::model::Spec*> specs = surfsys->getAllSpecs();
        pSpecs.insert(specs.begin(), specs.end());
    }

    return {pSpecs.begin(), pSpecs.end()};
}

////////////////////////////////////////////////////////////////////////////////

std::vector<steps::model::SReac*> Patch::getAllSReacs(steps::model::Model* model) const
{
    std::set<steps::model::SReac*> pSReacs;
    for (const auto& id : pSurfsys) {
        steps::model::Surfsys* surfsys = model->getSurfsys(id);
        std::vector<steps::model::SReac*> sreacs = surfsys->getAllSReacs();
        pSReacs.insert(sreacs.begin(), sreacs.end());
    }

    return {pSReacs.begin(), pSReacs.end()};
}

////////////////////////////////////////////////////////////////////////////////

void Patch::_setIComp(Comp* icomp)
{
    if (icomp->getContainer() != pContainer)
    {
        std::ostringstream os;
        os << "Compartment does not belong to same container as patch.\n";
        ArgErrLog(os.str());
    }
    auto const& ipatches  = icomp->getIPatches();
    if (ipatches.find(this) != ipatches.end())
    {
        std::ostringstream os;
        os << "Patch is already on inside of compartment.\n";
        ArgErrLog(os.str());
    }
    // remove the patch if it was already on the outside of some
    // other compartment
    if (pIComp != nullptr)
    {
        pIComp->_delOPatch(this);
    }

    pIComp = icomp;
    pIComp->_addOPatch(this);

}

////////////////////////////////////////////////////////////////////////////////

void Patch::_setOComp(Comp* ocomp)
{
    if (ocomp == nullptr) {
        return;
    }

    if (ocomp->getContainer() != pContainer)
    {
        std::ostringstream os;
           os << "Compartment does not belong to same container as patch.\n";
           ArgErrLog(os.str());
    }
    auto const& opatches  = ocomp->getOPatches();
    if (opatches.find(this) != opatches.end())
    {
          std::ostringstream os;
          os << "Patch is already on outside of compartment.\n";
           ArgErrLog(os.str());
    }
    // remove the patch if it was already on the inside of some
    // other compartment
    if (pOComp != nullptr)
    {
        pOComp->_delIPatch(this);
    }

    pOComp = ocomp;
    pOComp->_addIPatch(this);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::_handleSelfDelete()
{
    pContainer->_handlePatchDel(this);
    pArea = 0.0;
    pSurfsys.clear();
    pIComp = nullptr;
    pOComp = nullptr;
    pContainer = nullptr;
}

} // namespace wm
} // namespace steps

