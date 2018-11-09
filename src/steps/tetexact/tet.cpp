/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/reacdef.hpp"
#include "steps/tetexact/diff.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/tetexact/reac.hpp"
#include "steps/tetexact/tet.hpp"
#include "steps/tetexact/tetexact.hpp"
#include "steps/tetexact/tri.hpp"
#include "steps/tetexact/wmvol.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

stex::Tet::Tet
(
    uint idx, solver::Compdef * cdef, double vol,
    double a0, double a1, double a2, double a3,
    double d0, double d1, double d2, double d3,
    int tet0, int tet1, int tet2, int tet3
)
: WmVol(idx, cdef, vol)
, pTets()
//, pTris()
, pNextTet()
, pAreas()
, pDist()
{
    AssertLog(a0 > 0.0 && a1 > 0.0 && a2 > 0.0 && a3 > 0.0);
    AssertLog(d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0);

    pNextTris.resize(4);

    // At this point we don't have neighbouring tet pointers,
    // but we can store their indices
    for (uint i=0; i <= 3; ++i)
    {
        pNextTet[i] = nullptr;
        pNextTris[i] = 0;
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

    std::fill_n(pDiffBndDirection, 4, false);
    kprocs().resize(compdef()->countDiffs() + compdef()->countReacs());

}

////////////////////////////////////////////////////////////////////////////////

stex::Tet::~Tet()
= default;

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)pDiffBndDirection, sizeof(bool) * 4);
    WmVol::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::restore(std::fstream & cp_file)
{
    cp_file.read((char*)pDiffBndDirection, sizeof(bool) * 4);
    WmVol::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setNextTet(uint i, stex::Tet * t)
{

    // Now adding all tets, even those from other compartments, due to the diffusion boundaries
    pNextTet[i] = t;

    //if (pNextTris[i] != 0) CLOG(INFO, "general_log") << "WARNING: writing over nextTri index " << i;
    pNextTris[i] = 0;

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setDiffBndDirection(uint i)
{
    AssertLog(i < 4);

    pDiffBndDirection[i] = true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setNextTri(stex::Tri *t)
{
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setNextTri(uint i, stex::Tri * t)
{
    AssertLog(pNextTris.size() == 4);
    AssertLog(i <= 3);

    pNextTet[i] = nullptr;
    pNextTris[i]= t;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setupKProcs(stex::Tetexact * tex)
{
    uint j = 0;

    // Create reaction kproc's.
    uint nreacs = compdef()->countReacs();
    for (uint i = 0; i < nreacs; ++i)
    {
        ssolver::Reacdef * rdef = compdef()->reacdef(i);
        stex::Reac * r = new stex::Reac(rdef, this);
        kprocs()[j++] = r;
        tex->addKProc(r);
    }

    // Create diffusion kproc's.
    // NOTE: The order is important here- diffs should come after reacs,
    // because diffs will not be stored in WmVols and the Comp will call the
    // parent method often.
    uint ndiffs = compdef()->countDiffs();
    for (uint i = 0; i < ndiffs; ++i)
    {
        ssolver::Diffdef * ddef = compdef()->diffdef(i);
        auto * d = new stex::Diff(ddef, this);
        kprocs()[j++] = d;
        tex->addKProc(d);
    }
}

////////////////////////////////////////////////////////////////////////////////

stex::Diff * stex::Tet::diff(uint lidx) const
{
    AssertLog(lidx < compdef()->countDiffs());
    return dynamic_cast<stex::Diff*>(pKProcs[compdef()->countReacs() + lidx]);
}

////////////////////////////////////////////////////////////////////////////////

 int stex::Tet::getTetDirection(uint tidx)
{
    for (uint i = 0; i < 4; i++) {
        if (pTets[i] == tidx) {
            return i;
        }
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////////

// END
