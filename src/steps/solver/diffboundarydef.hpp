/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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



#ifndef STEPS_SOLVER_DIFFBOUNDARYDEF_HPP
#define STEPS_SOLVER_DIFFBOUNDARYDEF_HPP 1


// STL headers.
#include <string>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/statedef.hpp"
#include "steps/geom/diffboundary.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace solver {

// Forwards declarations
class DiffBoundarydef;

// Auxiliary declarations.
typedef DiffBoundarydef *                      DiffBoundaryDefP;
typedef std::vector<DiffBoundaryDefP>          DiffBoundaryDefPVec;
typedef DiffBoundaryDefPVec::iterator          DiffBoundaryDefPVecI;
typedef DiffBoundaryDefPVec::const_iterator    DiffBoundaryDefPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Defined diffusion boundary object.
class DiffBoundarydef
{

public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the object.
    /// \param d Pointer to the associated Diff boundary object.
    DiffBoundarydef(Statedef * sd, uint idx, steps::tetmesh::DiffBoundary * db);

    /// Destructor
    ~DiffBoundarydef();

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
    inline uint gidx() const noexcept
    { return pIdx; }

    /// Return the name of this diffusion boundary.
    const std::string& name() const noexcept
    { return pName; }

    inline const std::vector<triangle_id_t>& tris() const noexcept
    { return pTris; }

    inline uint compa() const noexcept
    { return pCompA; }

    inline uint compb() const noexcept
    { return pCompB; }


    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    Statedef                          * pStatedef;

    bool                                 pSetupdone{false};

    // The global index of this diffusion boundary
    uint                                pIdx;

    // The string identifier of this diffusion rule
    std::string                         pName;

    // List of all the triangles

    std::vector<triangle_id_t>                     pTris;

    // Diffboundarydef will have a setup() to copy the Compdef pointers
    //Lets store pointers
    uint                                 pCompA{0};
    uint                                 pCompB{0};

    // The pointer to the well-mixed comps is stored, but should not be used
    // only here so it's available during setup.
    steps::wm::Comp *                     pCompA_temp{nullptr};
    steps::wm::Comp *                     pCompB_temp{nullptr};

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_SOLVER_DIFFBOUNDARYDEF_HPP

// END
