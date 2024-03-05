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

#pragma once

#include "model.hpp"

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////
/// Multi-state complex reactant.
/// Component that represents a complex that can be referred to from
/// volume and surface systems.
class Complex {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Complex(std::string const& id, Model& model, const uint& nbSub, const uint& nbSt);
    virtual ~Complex();

    ////////////////////////////////////////////////////////////////////////
    // COMPLEX PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    const std::string& getID() const noexcept {
        return pID;
    }

    Model& getModel() const noexcept {
        return pModel;
    }

    unsigned int getNbSubStates() const noexcept {
        return pnbSubStates;
    }

    unsigned int getNbSubUnits() const noexcept {
        return pnbSubunits;
    }

  private:
    const std::string pID;
    Model& pModel;
    const unsigned int pnbSubunits;
    const unsigned int pnbSubStates;
};

}  // namespace steps::model
