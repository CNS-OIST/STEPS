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

// STL headers.
#include <sstream>
#include <string>

// STEPS headers.
#include "model/model.hpp"
#include "model/spec.hpp"
#include "model/vessdiff.hpp"
#include "model/vessurfsys.hpp"
#include "model/volsys.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

VesSDiff::VesSDiff(std::string const& id, VesSurfsys* vessurfsys, Spec* lig, double dcst)
    : pID(id)
    , pModel(nullptr)
    , pVesSurfsys(vessurfsys)
    , pLig(lig)
    , pDcst(dcst) {
    if (pVesSurfsys == nullptr) {
        std::ostringstream os;
        os << "No surfsys provided to Diff initializer function.";
        ArgErrLog(os.str());
    }

    if (pDcst < 0.0) {
        std::ostringstream os;
        os << "Diffusion constant can't be negative";
        ArgErrLog(os.str());
    }

    pModel = pVesSurfsys->getModel();
    AssertLog(pModel != nullptr);

    pVesSurfsys->_handleVesSDiffAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

VesSDiff::~VesSDiff() {
    if (pVesSurfsys == nullptr) {
        return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void VesSDiff::_handleSelfDelete() {
    pVesSurfsys->_handleVesSDiffDel(this);
    pVesSurfsys = nullptr;
    pDcst = 0.0;
    pLig = nullptr;
    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void VesSDiff::setID(std::string const& id) {
    AssertLog(pVesSurfsys != nullptr);

    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pVesSurfsys->_handleVesSDiffIDChange(pID, id);

    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void VesSDiff::setDcst(double dcst) {
    AssertLog(pVesSurfsys != nullptr);

    if (dcst < 0.0) {
        std::ostringstream os;
        os << "Diffusion constant can't be negative";
        ArgErrLog(os.str());
    }
    pDcst = dcst;
}

////////////////////////////////////////////////////////////////////////////////

void VesSDiff::setLig(Spec* lig) {
    AssertLog(lig != nullptr);
    pLig = lig;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec*> VesSDiff::getAllSpecs() const {
    return {pLig};
}

}  // namespace steps::model
