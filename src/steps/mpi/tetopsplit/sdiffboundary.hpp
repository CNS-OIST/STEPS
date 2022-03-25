////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2014 Okinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006 University of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////


#ifndef STEPS_MPI_TETOPSPLIT_SDIFFBOUNDARY_HPP
#define STEPS_MPI_TETOPSPLIT_SDIFFBOUNDARY_HPP 1


// STL headers.
#include <cassert>
#include <vector>

// STEPS headers.
#include "patch.hpp"
#include "util/common.h"
#include "solver/types.hpp"
#include "solver/sdiffboundarydef.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace mpi {
namespace tetopsplit {

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;

////////////////////////////////////////////////////////////////////////////////


// Forward declarations.
class SDiffBoundary;

// Auxiliary declarations.
typedef SDiffBoundary *                         SDiffBoundaryP;
typedef std::vector<SDiffBoundaryP>             SDiffBoundaryPVec;
typedef SDiffBoundaryPVec::iterator             SDiffBoundaryPVecI;
typedef SDiffBoundaryPVec::const_iterator       SDiffBoundaryPVecCI;


////////////////////////////////////////////////////////////////////////////////

class SDiffBoundary
{
public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    explicit SDiffBoundary(steps::solver::SDiffBoundarydef * sdbdef);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::SDiffBoundarydef * def() const noexcept
    { return pSDiffBoundarydef; }

    // We need access to the compartments so as to check if species are defined
    smtos::Patch * patchA();

    smtos::Patch * patchB();

    void setPatches(smtos::Patch * patcha, smtos::Patch * patchb);


    void setTriDirection(triangle_id_t tri, uint direction);

    const std::vector<triangle_id_t>& getTris() const noexcept
    { return pTris; }

    const std::vector<uint>& getTriDirection() const noexcept
    { return pTriDirection; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::SDiffBoundarydef    *  pSDiffBoundarydef;

    // Bool to check if patches have been specified
    bool                                  pSetPatches;

    // Patch arbitrarily labelled 'A'
    smtos::Patch                         * pPatchA;
    // Compartment arbitrarily labelled 'B'
    smtos::Patch                         * pPatchB;

    // A big vector of all the tris connected to this diffusion boundary
    std::vector<triangle_id_t>             pTris;

    // Directions have to be stored here - a tri could be connected to 2
    // different diffusion boundaries for example
    // This will be the same length as the pTets vector
    std::vector<uint>                     pTriDirection;

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

#endif
// STEPS_MPI_TETOPSPLIT_SDIFFBOUNDARY_HPP

// END
