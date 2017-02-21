/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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

// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/diffboundarydef.hpp"
#include "steps/mpi/tetopsplit/diffboundary.hpp"
#include "steps/mpi/tetopsplit/comp.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

smtos::DiffBoundary::DiffBoundary(steps::solver::DiffBoundarydef * dbdef)
: pDiffBoundarydef(dbdef)
, pCompA(0)
, pCompB(0)
, pTets()
, pTetDirection()
, pSetComps(false)
{
    assert(dbdef != 0);
}

////////////////////////////////////////////////////////////////////////////////

smtos::DiffBoundary::~DiffBoundary(void)
{

}

////////////////////////////////////////////////////////////////////////////////

void smtos::DiffBoundary::checkpoint(std::fstream & cp_file)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void smtos::DiffBoundary::restore(std::fstream & cp_file)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void smtos::DiffBoundary::setComps(smtos::Comp * compa, smtos::Comp * compb)
{
    assert (pSetComps == false);
    assert(compa != 0);
    assert(compb != 0);
    assert(compa != compb);

    pCompA = compa;
    pCompB = compb;
    pSetComps = true;
}

////////////////////////////////////////////////////////////////////////////////

smtos::Comp * smtos::DiffBoundary::compA(void)
{
    assert(pSetComps == true);
    return pCompA;
}

////////////////////////////////////////////////////////////////////////////////

smtos::Comp * smtos::DiffBoundary::compB(void)
{
    assert(pSetComps == true);
    return pCompB;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::DiffBoundary::setTetDirection(uint tet, uint direction)
{
    assert(direction < 4);

    pTets.push_back(tet);
    pTetDirection.push_back(direction);
}

////////////////////////////////////////////////////////////////////////////////

// END

