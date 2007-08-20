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
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/specdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/sim/swiginf/func_ssa.hpp>
#include <steps/sim/wmdirect/state.hpp>

////////////////////////////////////////////////////////////////////////////////

double siStep(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

uint siGetNSteps(State * s)
{
    assert(s != 0);
    return s->fNSteps;
}

////////////////////////////////////////////////////////////////////////////////

double siGetA0(State * s)
{
    assert(s != 0);
    return s->computeZeroProp();
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompReacC(State * s, uint cidx, uint ridx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidReac(ridx) == true);
    
    uint l_ridx = s->def()->comp(cidx)->reacG2L(ridx);
    if (l_ridx == 0xFFFF) return 0.0;
    return s->fReacCcsts[cidx][l_ridx];
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompReacH(State * s, uint cidx, uint ridx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidReac(ridx) == true);
    
    uint l_ridx = s->def()->comp(cidx)->reacG2L(ridx);
    if (l_ridx == 0xFFFF) return 0.0;
    return s->fReacHs[cidx][l_ridx];
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompReacA(State * s, uint cidx, uint ridx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidReac(ridx) == true);
    
    uint l_ridx = s->def()->comp(cidx)->reacG2L(ridx);
    if (l_ridx == 0xFFFF) return 0.0;
    return s->fReacProps[cidx][l_ridx];
}

////////////////////////////////////////////////////////////////////////////////

uint siGetCompReacExtent(State * s, uint cidx, uint ridx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidReac(ridx) == true);
    
    uint l_ridx = s->def()->comp(cidx)->reacG2L(ridx);
    if (l_ridx == 0xFFFF) return 0;
    return s->fReacExtents[cidx][l_ridx];
}

////////////////////////////////////////////////////////////////////////////////

void siResetCompReacExtent(State * s, uint cidx, uint ridx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidReac(ridx) == true);
    
    uint l_ridx = s->def()->comp(cidx)->reacG2L(ridx);
    if (l_ridx == 0xFFFF) return;
    s->fReacExtents[cidx][l_ridx] = 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
