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
#include <steps/tetexact/solver_core/patch.hpp>
#include <steps/tetexact/solver_core/tri.hpp>

NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

Patch::Patch(ssim::PatchDef * patchdef)
: pPatchDef(patchdef)
, pTris()
{
    assert(pPatchDef != 0);
}

////////////////////////////////////////////////////////////////////////////////

Patch::~Patch(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void Patch::addTri(Tri * tri)
{
    assert(tri->patchdef() == def());
    pTris.push_back(tri);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::computeArea(void)
{
    TriPVecCI tri_end = endTri();
    pArea = 0.0;
    for (TriPVecCI tri = bgnTri(); tri != tri_end; ++tri)
    {
        pArea += (*tri)->area();
    }
}

////////////////////////////////////////////////////////////////////////////////

// END
