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
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>

// STEPS headers.
#include "steps/common.h"
#include "steps/tetode/tetode.hpp"
#include "steps/tetode/tet.hpp"
#include "steps/tetode/tri.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace stode = steps::tetode;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

stode::Tet::Tet
(
    uint idx, solver::Compdef * cdef, double vol,
    double a0, double a1, double a2, double a3,
    double d0, double d1, double d2, double d3,
    int tet0, int tet1, int tet2, int tet3
)
: pCompdef(cdef)
, pIdx(idx)
, pVol(vol)
, pTets()
//, pTris()
, pNextTet()
, pNextTri()
, pAreas()
, pDist()
{
    assert (a0 > 0.0 && a1 > 0.0 && a2 > 0.0 && a3 > 0.0);
    assert (d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0);


    // At this point we don't have neighbouring tet pointers,
    // but we can store their indices
    for (uint i=0; i <= 3; ++i)
    {
        pNextTet[i] = 0;
        pNextTri[i] = 0;
    }
    pTets[0] = tet0;
    pTets[1] = tet1;
    pTets[2] = tet2;
    pTets[3] = tet3;

    pAreas[0] = a0;
    pAreas[1] = a1;
    pAreas[2] = a2;
    pAreas[3] = a3;

    pDist[0] = d0;
    pDist[1] = d1;
    pDist[2] = d2;
    pDist[3] = d3;


}

////////////////////////////////////////////////////////////////////////////////

 stode::Tet::~Tet(void)
{

}

////////////////////////////////////////////////////////////////////////////////

void stode::Tet::setNextTet(uint i, stode::Tet * t)
{

    if (t->compdef() != compdef())
    {
        pNextTet[i] = 0;
    }
    else
    {
        pNextTet[i] = t;
        if (pNextTri[i] != 0) std::cout << "WARNING: writing over nextTri index " << i;
        pNextTri[i] = 0;
    }

}

////////////////////////////////////////////////////////////////////////////////
/*
void stode::Tet::setNextTri(stex::Tri *t)
{
    uint index = pNextTris.size();
    pNextTris.push_back(t);
}
*/
////////////////////////////////////////////////////////////////////////////////

void stode::Tet::setNextTri(uint i, stode::Tri * t)
{


    // This is too common now to include this message- for any internal patch this happens
    //if (pNextTet[i] != 0) std::cout << "WARNING: writing over nextTet index " << i;

    pNextTet[i] = 0;
    pNextTri[i]= t;
}


////////////////////////////////////////////////////////////////////////////////

void stode::Tet::checkpoint(std::fstream & cp_file)
{
}

////////////////////////////////////////////////////////////////////////////////

void stode::Tet::restore(std::fstream & cp_file)
{
}

////////////////////////////////////////////////////////////////////////////////
// END
