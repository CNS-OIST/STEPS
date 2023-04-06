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
#include "comp.hpp"
#include "kproc.hpp"
#include "model/reac.hpp"
#include "tet.hpp"
#include "wmvol.hpp"
#include "solver/compdef.hpp"

// logging
#include "util/error.hpp"
#include <easylogging++.h>
////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

smtos::Comp::Comp(steps::solver::Compdef * compdef)
: pCompdef(compdef)
{
    AssertLog(pCompdef != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Comp::checkpoint(std::fstream & /*cp_file*/)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Comp::restore(std::fstream & /*cp_file*/)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Comp::addTet(smtos::WmVol * tet)
{
    AssertLog(tet->compdef() == def());
    pTets.push_back(tet);
    pVol += tet->vol();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Comp::modCount(uint slidx, double count)
{
    AssertLog(slidx < def()->countSpecs());
    double newcount = (def()->pools()[slidx] + count);
    AssertLog(newcount >= 0.0);
    def()->setCount(slidx, newcount);
}

////////////////////////////////////////////////////////////////////////////////

smtos::WmVol * smtos::Comp::pickTetByVol(double rand01) const
{
    if (countTets() == 0) { return nullptr;
}
    if (countTets() == 1) return pTets[0];

    double accum = 0.0;
    double selector = rand01 * vol();
    for (auto const& t: pTets) {
        accum += t->vol();
        if (selector < accum) {
            return t;
        }
    }
    AssertLog(false);
}

///////////////////////////////////////////////////////////////////////////

// END
