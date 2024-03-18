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

#include "complexdef.hpp"

#include "model/complex.hpp"
#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

Complexdef::Complexdef(complex_global_id idx, steps::model::Complex& s)
    : pIdx(idx)
    , pNbSubStates(s.getNbSubStates())
    , pName(s.getID()) {}

////////////////////////////////////////////////////////////////////////////////

Complexdef::~Complexdef() = default;

////////////////////////////////////////////////////////////////////////////////

void Complexdef::checkpoint(std::fstream& /*cp_file*/) const {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void Complexdef::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void Complexdef::setup() {
    pSetupdone = true;
}

}  // namespace steps::solver
