////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
// $Id:tet.cpp 64 2007-08-20 06:25:41Z stefan $
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/diffdef.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/tetexact/solver_core/diff.hpp>
#include <steps/tetexact/solver_core/reac.hpp>
#include <steps/tetexact/solver_core/sched.hpp>
#include <steps/tetexact/solver_core/tet.hpp>
#include <steps/tetexact/solver_core/tri.hpp>

NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

Tet::Tet
(
    ssim::CompDef * cdef, double vol, 
    double a0, double a1, double a2, double a3, 
    double d0, double d1, double d2, double d3
)
{
    // Copy all this stuff.
    pCompDef = cdef;
    pNextTet[0] = 0;
    pNextTet[1] = 0;
    pNextTet[2] = 0;
    pNextTet[3] = 0;
    pNextTri[0] = 0;
    pNextTri[1] = 0;
    pNextTri[2] = 0;
    pNextTri[3] = 0;
    // Tetrahedral volumes.
    assert(vol >= 0.0);
    pVol = vol;
    // Areas of boundary triangles.
    assert(a0 >= 0.0);
    pAreas[0] = a0;
    assert(a1 >= 0.0);
    pAreas[1] = a1;
    assert(a2 >= 0.0);
    pAreas[2] = a2;
    assert(a3 >= 0.0);
    pAreas[3] = a3;
    // Distances to neighbouring tetrahedrons.
    assert(d0 >= 0.0);
    pDist[0] = d0;
    assert(d1 >= 0.0);
    pDist[1] = d1;
    assert(d2 >= 0.0);
    pDist[2] = d2;
    assert(d3 >= 0.0);
    pDist[3] = d3;
    
    // Based on compartment definition, build other structures.
    uint nspecs = compdef()->countSpecs();
    pPoolCount = new uint[nspecs];
    pPoolFlags = new uint[nspecs];
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);
    pKProcs.resize(cdef->countDiffs() + cdef->countReacs());
}

////////////////////////////////////////////////////////////////////////////////

Tet::~Tet(void)
{
    // Delete species pool information.
    delete[] pPoolCount;
    delete[] pPoolFlags;
    
    // Delete diffusion rules.
    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) delete *i;
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setNextTet(uint i, Tet * t)
{
    if (t->compdef() != compdef()) 
    {
        pNextTet[i] = 0;
    }
    else
    {
        pNextTet[i] = t;
        pNextTri[i] = 0;
    }
}


////////////////////////////////////////////////////////////////////////////////

void Tet::setNextTri(uint i, Tri * t)
{
    pNextTet[i] = 0;
    pNextTri[i] = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setupKProcs(Sched * sched)
{
    uint j = 0;
    // Create diffusion kproc's.
    uint ndiffs = compdef()->countDiffs();
    for (uint i = 0; i < ndiffs; ++i)
    {
        ssim::DiffDef * ddef = compdef()->diff(i);
        Diff * d = new Diff(ddef, this);
        pKProcs[j++] = d;
        sched->addKProc(d);
    }
    
    // Create reaction kproc's.
    uint nreacs = compdef()->countReacs();
    for (uint i = 0; i < nreacs; ++i)
    {
        ssim::ReacDef * rdef = compdef()->reac(i);
        Reac * r = new Reac(rdef, this);
        pKProcs[j++] = r;
        sched->addKProc(r);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setupDeps(void)
{
    std::for_each(pKProcs.begin(), pKProcs.end(), 
        std::mem_fun(&KProc::setupDeps));
}

////////////////////////////////////////////////////////////////////////////////

void Tet::reset(void)
{
    uint nspecs = compdef()->countSpecs();
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);
    std::for_each(pKProcs.begin(), pKProcs.end(), std::mem_fun(&KProc::reset));
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setClamped(uint lidx, bool clamp)
{
    if (clamp == true) pPoolFlags[lidx] |= CLAMPED;
    else pPoolFlags[lidx] &= ~CLAMPED;
}

////////////////////////////////////////////////////////////////////////////////

Diff * Tet::diff(uint lidx) const
{
    assert(lidx < compdef()->countDiffs());
    return dynamic_cast<Diff*>(pKProcs[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

Reac * Tet::reac(uint lidx) const
{
    assert(lidx < compdef()->countReacs());
    return dynamic_cast<Reac*>(pKProcs[compdef()->countDiffs() + lidx]);
}
    
////////////////////////////////////////////////////////////////////////////////

// END
