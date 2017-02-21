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
#include "steps/solver/compdef.hpp"
#include "steps/tetexact/comp.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/tetexact/reac.hpp"
#include "steps/tetexact/tet.hpp"
#include "steps/tetexact/wmvol.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

stex::Comp::Comp(steps::solver::Compdef * compdef)
: pCompdef(compdef)
, pVol(0.0)
, pTets()
{
    assert(pCompdef != 0);
}

////////////////////////////////////////////////////////////////////////////////

stex::Comp::~Comp(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stex::Comp::checkpoint(std::fstream & cp_file)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void stex::Comp::restore(std::fstream & cp_file)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void stex::Comp::addTet(stex::WmVol * tet)
{
    assert (tet->compdef() == def());
    pTets.push_back(tet);
    pVol += tet->vol();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Comp::modCount(uint slidx, double count)
{
    assert (slidx < def()->countSpecs());
    double newcount = (def()->pools()[slidx] + count);
    assert (newcount >= 0.0);
    def()->setCount(slidx, newcount);
}

////////////////////////////////////////////////////////////////////////////////

stex::WmVol * stex::Comp::pickTetByVol(double rand01) const
{
    if (countTets() == 0) return 0;
    if (countTets() == 1) return pTets[0];

    double accum = 0.0;
    double selector = rand01 * vol();
    WmVolPVecCI t_end = endTet();
    for (WmVolPVecCI t = bgnTet(); t != t_end; ++t)
    {
        accum += (*t)->vol();
        if (selector < accum) return (*t);
    }
    assert(false);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
