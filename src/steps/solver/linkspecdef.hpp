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

// STL headers.
#include <fstream>
#include <string>

// STEPS headers.
#include "model/linkspec.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"

namespace steps::solver {

class Statedef;

/// Defined LinkSpecies
class LinkSpecdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the linkspecies.
    /// \param l Reference to the associated LinkSpec object.
    LinkSpecdef(Statedef* sd, linkspec_global_id idx, model::LinkSpec* l);

    /// Destructor
    ~LinkSpecdef(void);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: LINKSPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this link species.
    inline linkspec_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the species.
    const std::string& name() const noexcept {
        return pName;
    }

    inline double dcst() const noexcept {
        return pDcst;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
    ///
    /// This method is included for consistency with other def objects,
    /// but currently does very litte.
    void setup();

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef* pStatedef;
    linkspec_global_id pIdx;
    std::string pName;
    double pDcst;
    bool pSetupdone;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::solver
