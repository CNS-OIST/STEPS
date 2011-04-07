////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2011�Okinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006�University of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPS�is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPS�is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.�If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_TETEXACT_DIFFBOUNDARY_HPP
#define STEPS_TETEXACT_DIFFBOUNDARY_HPP 1


// STL headers.
#include <cassert>
#include <vector>
#include <fstream>s

// STEPS headers.
#include "../common.h"
#include "comp.hpp"
//#include "tri.hpp"
#include "../solver/types.hpp"
#include "../solver/diffboundarydef.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetexact)

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);

////////////////////////////////////////////////////////////////////////////////


// Forward declarations.
class DiffBoundary;

// Auxiliary declarations.
typedef DiffBoundary *                         DiffBoundaryP;
typedef std::vector<DiffBoundaryP>             DiffBoundaryPVec;
typedef DiffBoundaryPVec::iterator             DiffBoundaryPVecI;
typedef DiffBoundaryPVec::const_iterator       DiffBoundaryPVecCI;


////////////////////////////////////////////////////////////////////////////////

class DiffBoundary
{
public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    DiffBoundary(steps::solver::DiffBoundarydef * dbdef);
    ~DiffBoundary(void);

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

    inline steps::solver::DiffBoundarydef * def(void) const
    { return pDiffBoundarydef; }

    // We need access to the compartments so as to check if species are defined
    stex::Comp * compA(void);

    stex::Comp * compB(void);

    void setComps(stex::Comp * compa, stex::Comp * compb);

    // Other data we need is the TETRAHEDRONS (not triangles) affected by
    // this diff boundary- if a solver method like setDiffBoundarySpec is called
    // then this object provides a list of tetrahedrons and importantly
    // THE DIRECTION OF DIFFUSION (the direction that is to the next compartment)
    // the solver can then simply loop over the tets and activate diffusion
    // in that direction (or we could do it here)


    // This information is the tetrahedron connected to this diffusion boundary
    // and the direction of the diffusion boundary for that tetrahedron (0 to 3)
    void setTetDirection(uint tet, uint direction);


    std::vector<uint> getTets(void) const
    { return pTets; }

    std::vector<uint> getTetDirection(void) const
    { return pTetDirection; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::DiffBoundarydef    * pDiffBoundarydef;

    // Bool to check if compartments have been specified
    bool 								pSetComps;

    // Compartment arbitrarily labelled 'A'
    stex::Comp 						  * pCompA;
    // Compartment arbitrarily labelled 'B'
    stex::Comp 						  * pCompB;

    // A big vector of all the tetrahedrons connected to this diffusion boundary
    // If the diff boundary allows passage of an ion these tets will tell
    // their diffs for that ion to allow diffusion
    // The direction information I'm thinking will come from WHERE??
    std::vector<uint> 					pTets;

    // Directions have to be stored here - a tet could be connected to 2
    // different diffusion boundaries for example
    // This will be the same length as the pTets vector
    std::vector<uint> 					pTetDirection;

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetexact)
END_NAMESPACE(steps)

#endif
// STEPS_TETEXACT_DIFFBOUNDARY_HPP

// END
