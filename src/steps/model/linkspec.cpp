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
#include "model/linkspec.hpp"
#include "model/model.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

LinkSpec::LinkSpec(std::string const& id, Model* model, double dcst)
    : pID(id)
    , pModel(model) {
    if (pModel == nullptr) {
        std::ostringstream os;
        os << "No model provided to LinkSpec initializer function";
        ArgErrLog(os.str());
    }

    if (dcst < 0.0) {
        std::ostringstream os;
        os << "Diffusion coefficient must not be negative!";
        ArgErrLog(os.str());
    }

    pDcst = dcst;

    pModel->_handleLinkSpecAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

LinkSpec::~LinkSpec() {
    if (pModel == nullptr) {
        return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpec::_handleSelfDelete() {
    pModel->_handleLinkSpecDel(this);
    pModel = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpec::setID(std::string const& id) {
    AssertLog(pModel != nullptr);
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel->_handleLinkSpecIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

}  // namespace steps::model
