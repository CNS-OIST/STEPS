/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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

#include "model/linkspec.hpp"

#include "model.hpp"
#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

LinkSpec::LinkSpec(std::string const& id, Model& model, double dcst, double max_angle)
    : pID(id)
    , pModel(model) {
    if (dcst < 0.0) {
        std::ostringstream os;
        os << "Diffusion coefficient must not be negative!";
        ArgErrLog(os.str());
    }

    if (max_angle > math::PI or max_angle < -math::PI) {
        std::ostringstream os;
        os << "Maximum angle must not be above pi radians!";
        ArgErrLog(os.str());
    }

    pDcst = dcst;
    pMaxAngle = max_angle;

    pModel._handleLinkSpecAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

LinkSpec::~LinkSpec() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpec::_handleSelfDelete() {
    pModel._handleLinkSpecDel(*this);
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpec::setID(std::string const& id) {
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel._handleLinkSpecIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

double LinkSpec::getMinCos2() const noexcept {
    const double cos_angle = cos(pMaxAngle);
    if (cos_angle < 0.0) {
        return -cos_angle * cos_angle;
    }
    return cos_angle * cos_angle;
}


}  // namespace steps::model
