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

#include "model/vessdiff.hpp"

#include "model.hpp"
#include "spec.hpp"
#include "vessurfsys.hpp"
#include "volsys.hpp"

#include "util/error.hpp"

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

VesSDiff::VesSDiff(std::string const& id, VesSurfsys& vessurfsys, Spec& lig, double dcst)
    : pID(id)
    , pModel(vessurfsys.getModel())
    , pVesSurfsys(vessurfsys)
    , pLig(&lig)
    , pDcst(dcst) {
    if (pDcst < 0.0) {
        std::ostringstream os;
        os << "Diffusion constant can't be negative";
        ArgErrLog(os.str());
    }
    pVesSurfsys._handleVesSDiffAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

VesSDiff::~VesSDiff() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void VesSDiff::_handleSelfDelete() {
    pVesSurfsys._handleVesSDiffDel(*this);
}

////////////////////////////////////////////////////////////////////////////////

void VesSDiff::setID(std::string const& id) {
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pVesSurfsys._handleVesSDiffIDChange(pID, id);

    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void VesSDiff::setDcst(double dcst) {
    if (dcst < 0.0) {
        std::ostringstream os;
        os << "Diffusion constant can't be negative: " << dcst;
        ArgErrLog(os.str());
    }
    pDcst = dcst;
}

////////////////////////////////////////////////////////////////////////////////

void VesSDiff::setLig(Spec& lig) {
    pLig = &lig;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Spec*> VesSDiff::getAllSpecs() const {
    return {pLig};
}

}  // namespace steps::model
