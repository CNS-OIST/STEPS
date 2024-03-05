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

#include <iosfwd>
#include <string>

#include "fwd.hpp"
#include "geom/fwd.hpp"
#include "util/vocabulary.hpp"

namespace steps::solver {

/// Defined diffusion boundary object.
class DiffBoundarydef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the object.
    /// \param db Reference to the associated Diff boundary object.
    DiffBoundarydef(Statedef& sd, diffboundary_global_id idx, tetmesh::DiffBoundary& db);

    DiffBoundarydef(const DiffBoundarydef&) = delete;
    DiffBoundarydef& operator=(const DiffBoundarydef&) = delete;

    void setup(Statedef& sd);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION BOUNDARY
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this diffusion boundary.
    inline diffboundary_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of this diffusion boundary.
    const std::string& name() const noexcept {
        return pName;
    }

    inline const std::vector<triangle_global_id>& tris() const noexcept {
        return pTris;
    }

    inline comp_global_id compa() const noexcept {
        return pCompA;
    }

    inline comp_global_id compb() const noexcept {
        return pCompB;
    }

  private:
    /// The global index of this diffusion boundary
    const diffboundary_global_id pIdx;

    /// The string identifier of this diffusion rule
    const std::string pName;

    /// List of all the triangles
    const std::vector<triangle_global_id>& pTris;

    bool pSetupdone{false};

    // Diffboundarydef will have a setup() to copy the Compdef pointers
    // Lets store pointers
    comp_global_id pCompA;
    comp_global_id pCompB;

    // The pointer to the well-mixed comps is stored, but should not be used
    // only here so it's available during setup.
    wm::Comp* pCompA_temp{nullptr};
    wm::Comp* pCompB_temp{nullptr};
};

}  // namespace steps::solver
