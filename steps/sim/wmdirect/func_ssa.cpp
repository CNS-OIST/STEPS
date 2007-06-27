////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

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
    if (l_ridx == 0xFFFF) return 0.0;
    return s->fReacExtents[cidx][l_ridx];
}

////////////////////////////////////////////////////////////////////////////////

void siResetCompReacExtent(State * s, uint cidx, uint ridx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidReac(ridx) == true);
    
    uint l_ridx = s->def()->comp(cidx)->reacG2L(ridx);
    if (l_ridx == 0xFFFF) return 0.0;
    s->fReacExtents[cidx][l_ridx] = 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
