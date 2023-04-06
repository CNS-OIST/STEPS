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

// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include "util/common.h"
#include "util/error.hpp"
#include "mpi/tetopsplit/comp.hpp"
#include "mpi/tetopsplit/diffboundary.hpp"
#include "solver/diffboundarydef.hpp"

// logging
#include <easylogging++.h>
////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

smtos::DiffBoundary::DiffBoundary(steps::solver::DiffBoundarydef * dbdef)
: pDiffBoundarydef(dbdef)
{
    AssertLog(dbdef != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::DiffBoundary::checkpoint(std::fstream & /*cp_file*/)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void smtos::DiffBoundary::restore(std::fstream & /*cp_file*/)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void smtos::DiffBoundary::setComps(smtos::Comp * compa, smtos::Comp * compb)
{
    AssertLog(pSetComps == false);
    AssertLog(compa != nullptr);
    AssertLog(compb != nullptr);
    AssertLog(compa != compb);

    pCompA = compa;
    pCompB = compb;
    pSetComps = true;
}

////////////////////////////////////////////////////////////////////////////////

smtos::Comp * smtos::DiffBoundary::compA()
{
    AssertLog(pSetComps);
    return pCompA;
}

////////////////////////////////////////////////////////////////////////////////

smtos::Comp * smtos::DiffBoundary::compB()
{
    AssertLog(pSetComps);
    return pCompB;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::DiffBoundary::setTetDirection(tetrahedron_id_t tet, uint direction)
{
    AssertLog(direction < 4);

    pTets.push_back(tet);
    pTetDirection.push_back(direction);
}

////////////////////////////////////////////////////////////////////////////////

// END

