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
#include <string>

// STEPS headers.
#include "geom/endocyticzone.hpp"
#include "solver/fwd.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

/// Defined Raft
class EndocyticZonedef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param z Pointer to the associated EndocyticZone object.
    EndocyticZonedef(Statedef* sd, tetmesh::EndocyticZone* z);

    /// Destructor
    ~EndocyticZonedef();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    /// Return the name of the raft.
    inline std::string const& name() const noexcept {
        return pName;
    }

    inline const std::vector<triangle_global_id>& tris() const noexcept {
        return pTris;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef* pStatedef;
    std::string pName;

    ////////////////////////////////////////////////////////////////////////
    // DATA
    ////////////////////////////////////////////////////////////////////////

    // Triangles
    std::vector<triangle_global_id> pTris;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::solver
