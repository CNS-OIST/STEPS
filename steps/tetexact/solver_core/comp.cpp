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
#include <cassert>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
#include <steps/tetexact/solver_core/comp.hpp>
#include <steps/tetexact/solver_core/tet.hpp>

NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

Comp::Comp(ssim::CompDef * compdef)
: pCompDef(compdef)
, pVol(compdef->vol())
, pTets()
{
    assert(pCompDef != 0);
}

////////////////////////////////////////////////////////////////////////////////

Comp::~Comp(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void Comp::addTet(Tet * tet)
{
    assert(tet->compdef() == def());
    pTets.push_back(tet);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::computeVol(void)
{
    TetPVecCI tet_end = endTet();
    pVol = 0.0;
    for (TetPVecCI tet = bgnTet(); tet != tet_end; ++tet)
    {
        pVol += (*tet)->vol();
    }
}

////////////////////////////////////////////////////////////////////////////////

Tet * Comp::pickTetByVol(double rand01) const
{
    if (countTets() == 0) return 0;
    if (countTets() == 1) return pTets[0];
    
    double accum = 0.0;
    double selector = rand01 * vol();
    TetPVecCI t_end = endTet();
    for (TetPVecCI t = bgnTet(); t != t_end; ++t)
    {
        accum += (*t)->vol();
        if (selector < accum) return *t;
    }
    
    assert(0);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
