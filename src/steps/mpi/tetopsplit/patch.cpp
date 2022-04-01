/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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



/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include "kproc.hpp"
#include "patch.hpp"
#include "model/reac.hpp"
#include "solver/compdef.hpp"

// logging
#include <easylogging++.h>
#include "util/error.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

smtos::Patch::Patch(ssolver::Patchdef * patchdef)
: pPatchdef(patchdef)
{
    AssertLog(pPatchdef != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Patch::checkpoint(std::fstream & /*cp_file*/)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Patch::restore(std::fstream & /*cp_file*/)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Patch::addTri(smtos::Tri * tri)
{
    AssertLog(tri->patchdef() == def());
    pTris.push_back(tri);
    pArea += tri->area();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Patch::modCount(uint slidx, double count)
{
    AssertLog(slidx < def()->countSpecs());
    double newcount = def()->pools()[slidx] + count;
    AssertLog(newcount >= 0.0);
    def()->setCount(slidx, newcount);
}

////////////////////////////////////////////////////////////////////////////////

smtos::Tri * smtos::Patch::pickTriByArea(double rand01) const
{
    if (countTris() == 0) { return nullptr;
}
    if (countTris() == 1) return pTris[0];

    double accum = 0.0;
    double selector = rand01 * area();
    for (auto const& t: pTris) {
        accum += t->area();
        if (selector <= accum) {
            return t;
        }
    }

    return *(endTri() - 1);
}

////////////////////////////////////////////////////////////////////////////////
/*
void smtos::Patch::setArea(double a)
{
    AssertLog(a > 0.0);
    pArea = a;
}
*/
////////////////////////////////////////////////////////////////////////////////

// END
