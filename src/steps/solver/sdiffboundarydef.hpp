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
#include "geom/sdiffboundary.hpp"
#include "statedef.hpp"
#include "util/common.hpp"

namespace steps::solver {

////////////////////////////////////////////////////////////////////////////////
/// Defined surface diffusion boundary object.
class SDiffBoundarydef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the object.
    /// \param sdb Pointer to the associated Surface Diffusion boundary object.
    SDiffBoundarydef(Statedef* sd, sdiffboundary_global_id idx, tetmesh::SDiffBoundary* sdb);

    /// Destructor
    ~SDiffBoundarydef();

    void setup();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION BOUNDARY
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this diffusion boundary.
    inline sdiffboundary_global_id gidx() const {
        return pIdx;
    }

    /// Return the name of this diffusion boundary.
    const std::string& name() const noexcept {
        return pName;
    }

    inline const std::vector<index_t>& bars() const {
        return pBars;
    }

    inline patch_global_id patcha() const {
        return pPatchA;
    }
    inline patch_global_id patchb() const {
        return pPatchB;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef* pStatedef;

    bool pSetupdone;

    // The global index of this diffusion boundary
    sdiffboundary_global_id pIdx;

    // The string identifier of this diffusion rule
    std::string pName;

    // List of all the bars

    std::vector<index_t> pBars;

    patch_global_id pPatchA;
    patch_global_id pPatchB;

    // The pointer to the well-mixed comps is stored, but should not be used
    // only here so it's available during setup.
    wm::Patch* pPatchA_temp;
    wm::Patch* pPatchB_temp;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::solver
