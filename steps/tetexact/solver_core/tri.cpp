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
// $Id$
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/shared/sreacdef.hpp>
#include <steps/tetexact/solver_core/kproc.hpp>
#include <steps/tetexact/solver_core/sched.hpp>
#include <steps/tetexact/solver_core/sreac.hpp>
#include <steps/tetexact/solver_core/tet.hpp>
#include <steps/tetexact/solver_core/tri.hpp>

NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

Tri::Tri(ssim::PatchDef * pdef, double area)
: pPatchDef(pdef)
, pInnerTet(0)
, pOuterTet(0)
, pArea(area)
, pPoolCount(0)
, pPoolFlags(0)
, pKProcs()
{
    assert(pPatchDef != 0);
    assert(pArea >= 0.0);
    uint nspecs = patchdef()->countSpecs();
    pPoolCount = new uint[nspecs];
    pPoolFlags = new uint[nspecs];
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);
    pKProcs.resize(pPatchDef->countSReacs());
}

////////////////////////////////////////////////////////////////////////////////

Tri::~Tri(void)
{
    delete[] pPoolCount;
    delete[] pPoolFlags;
    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) delete *i;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setInnerTet(Tet * t)
{
    pInnerTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setOuterTet(Tet * t)
{
    pOuterTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setupKProcs(Sched * sched)
{
    uint j = 0;
    // Create surface reaction kproc's.
    uint nsreacs = patchdef()->countSReacs();
    for (uint i = 0; i < nsreacs; ++i)
    {
        ssim::SReacDef * srdef = patchdef()->sreac(i);
        SReac * sr = new SReac(srdef, this);
        assert(sr != 0);
        pKProcs[j++] = sr;
        sched->addKProc(sr);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setupDeps(void)
{
    std::for_each(pKProcs.begin(), pKProcs.end(), 
        std::mem_fun(&KProc::setupDeps));
}

////////////////////////////////////////////////////////////////////////////////

void Tri::reset(void)
{
    uint nspecs = patchdef()->countSpecs();
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);
    std::for_each(pKProcs.begin(), pKProcs.end(), 
        std::mem_fun(&KProc::reset));
}

////////////////////////////////////////////////////////////////////////////////

SReac * Tri::sreac(uint lidx) const
{
    assert(lidx < patchdef()->countSReacs());
    return dynamic_cast<SReac*>(pKProcs[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

// END
