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

#ifndef STEPS_SOLVER_SDIFFBOUNDARYDEF_HPP
#define STEPS_SOLVER_SDIFFBOUNDARYDEF_HPP 1


// STL headers.
#include <string>
#include <fstream>

// STEPS headers.
#include "util/common.h"
#include "statedef.hpp"
#include "geom/sdiffboundary.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace solver {

// Forwards declarations
class   SDiffBoundarydef;

// Auxiliary declarations.
typedef SDiffBoundarydef *                      SDiffBoundaryDefP;
typedef std::vector<SDiffBoundaryDefP>          SDiffBoundaryDefPVec;
typedef SDiffBoundaryDefPVec::iterator          SDiffBoundaryDefPVecI;
typedef SDiffBoundaryDefPVec::const_iterator    SDiffBoundaryDefPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Defined surface diffusion boundary object.
class SDiffBoundarydef
{

public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the object.
    /// \param sdb Pointer to the associated Surface Diffusion boundary object.
    SDiffBoundarydef(Statedef * sd, uint idx, steps::tetmesh::SDiffBoundary * sdb);

    /// Destructor
    ~SDiffBoundarydef();

    void setup();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION BOUNDARY
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this diffusion boundary.
    inline uint gidx() const
    { return pIdx; }

    /// Return the name of this diffusion boundary.
    const std::string& name() const noexcept {
        return pName;
    }

    inline const std::vector<index_t>& bars() const
    { return pBars; }

    inline uint patcha() const
    { return pPatchA; }
    inline uint patchb() const
    { return pPatchB; }


    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    Statedef                           * pStatedef;

    bool                                 pSetupdone;

    // The global index of this diffusion boundary
    uint                                 pIdx;

    // The string identifier of this diffusion rule
    std::string                          pName;

    // List of all the bars

    std::vector<index_t>    pBars;

    uint                                 pPatchA;
    uint                                 pPatchB;

    // The pointer to the well-mixed comps is stored, but should not be used
    // only here so it's available during setup.
    steps::wm::Patch *                   pPatchA_temp;
    steps::wm::Patch *                   pPatchB_temp;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_SOLVER_SDIFFBOUNDARYDEF_HPP

// END
