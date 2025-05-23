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

#include "spec.hpp"

#include "model.hpp"

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////

Spec::Spec(std::string const& id, Model& model, int valence)
    : pID(id)
    , pModel(model)
    , pValence(valence) {
    pModel._handleSpecAdd(*this);
}

////////////////////////////////////////////////////////////////////////////////

Spec::~Spec() {
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Spec::_handleSelfDelete() {
    pModel._handleSpecDel(*this);
}

////////////////////////////////////////////////////////////////////////////////

void Spec::setID(std::string const& id) {
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel._handleSpecIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Spec::setValence(int valence) {
    pValence = valence;
}

}  // namespace steps::model
