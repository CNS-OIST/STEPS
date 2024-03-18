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

#include "diff.hpp"

#include "model.hpp"
#include "spec.hpp"
#include "surfsys.hpp"
#include "volsys.hpp"

#include "util/error.hpp"

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

Diff::Diff(std::string const& id, Volsys& volsys, Spec& lig, double dcst)
    : pID(id)
    , pModel(volsys.getModel())
    , pVolsys(&volsys)
    , pLig(&lig)
    , pDcst(dcst)
    , pIsvolume(true) {
    ArgErrLogIf(pDcst < 0.0, "Diffusion constant can't be negative");
    pVolsys->_handleDiffAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

Diff::Diff(std::string const& id, Surfsys& surfsys, Spec& lig, double dcst)
    : pID(id)
    , pModel(surfsys.getModel())
    , pSurfsys(&surfsys)
    , pLig(&lig)
    , pDcst(dcst)
    , pIsvolume(false) {
    ArgErrLogIf(pDcst < 0.0, "Diffusion constant can't be negative");

    pSurfsys->_handleDiffAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

Diff::~Diff() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Diff::_handleSelfDelete() {
    if (pIsvolume) {
        pVolsys->_handleDiffDel(*this);
    } else {
        pSurfsys->_handleDiffDel(*this);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setID(std::string const& id) {
    if (pIsvolume) {
        AssertLog(pVolsys != nullptr);
        // The following might raise an exception, e.g. if the new ID is not
        // valid or not unique. If this happens, we don't catch but simply let
        // it pass by into the Python layer.
        pVolsys->_handleDiffIDChange(pID, id);
    } else {
        AssertLog(pSurfsys != nullptr);
        // The following might raise an exception, e.g. if the new ID is not
        // valid or not unique. If this happens, we don't catch but simply let
        // it pass by into the Python layer.
        pSurfsys->_handleDiffIDChange(pID, id);
    }
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setDcst(double dcst) {
    if (pIsvolume) {
        AssertLog(pVolsys != nullptr);
    } else {
        AssertLog(pSurfsys != nullptr);
    }

    ArgErrLogIf(dcst < 0.0, "Diffusion constant can't be negative");

    pDcst = dcst;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setLig(Spec& lig) {
    pLig = &lig;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec*> Diff::getAllSpecs() const {
    return {pLig};
}

}  // namespace steps::model
