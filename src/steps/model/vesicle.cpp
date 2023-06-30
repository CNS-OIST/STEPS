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
#include "model/vesicle.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

Vesicle::Vesicle(std::string const& id, Model* model, double diameter, double dcst)
    : pID(id)
    , pModel(model)
    , pDiameter(diameter)
    , pDcst(dcst) {
    if (pModel == nullptr) {
        std::ostringstream os;
        os << "No model provided to Vesicle initializer function";
        ArgErrLog(os.str());
    }
    if (pDcst < 0.0) {
        std::ostringstream os;
        os << "Vesicle diffusion constant can't be negative";
        ArgErrLog(os.str());
    }
    pModel->_handleVesicleAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Vesicle::~Vesicle() {
    if (pModel == nullptr) {
        return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::_handleSelfDelete() {
    pModel->_handleVesicleDel(this);
    pModel = nullptr;
    pDcst = 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::addVesSurfsys(std::string const& id) {
    // string identifier is only added to set if it is not already included
    pVesSurfsys.insert(id);
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::setDcst(double dcst) {
    if (dcst < 0.0) {
        std::ostringstream os;
        os << "Vesicle diffusion constant can't be negative";
        ArgErrLog(os.str());
    }
    pDcst = dcst;
}

////////////////////////////////////////////////////////////////////////////////

void Vesicle::setID(std::string const& id) {
    AssertLog(pModel != nullptr);
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel->_handleVesicleIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

}  // namespace steps::model
