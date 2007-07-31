////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

%module tetexact_core

%import "steps/common.h"

%{
#include <steps/rng/rng.hpp>
#include <steps/sim/swiginf/func_core.hpp>
#include <steps/sim/swiginf/func_ssa.hpp>
#include <steps/sim/swiginf/func_tetmesh.hpp>
%}

%include swiginf/func_core.i
%include swiginf/func_ssa.i
%include swiginf/func_tetmesh.i

////////////////////////////////////////////////////////////////////////////////

// END
