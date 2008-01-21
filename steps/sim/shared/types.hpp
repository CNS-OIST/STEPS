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

#ifndef STEPS_SIM_SHARED_TYPES_HPP
#define STEPS_SIM_SHARED_TYPES_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <vector>

START_NAMESPACE(steps)
START_NAMESPACE(sim)

////////////////////////////////////////////////////////////////////////////////

typedef uint                            gidxT;

typedef std::vector<gidxT>              gidxTVec;
typedef gidxTVec::iterator              gidxTVecI;
typedef gidxTVec::const_iterator        gidxTVecCI;

typedef std::vector<gidxT*>             gidxTPVec;
typedef gidxTPVec::iterator             gidxTPVecI;
typedef gidxTPVec::const_iterator       gidxTPVecCI;

static const gidxT GIDX_UNDEFINED = 0xFFFF;

////////////////////////////////////////////////////////////////////////////////

typedef uint                            lidxT;

typedef std::vector<lidxT>              lidxTVec;
typedef lidxTVec::iterator              lidxTVecI;
typedef lidxTVec::const_iterator        lidxTVecCI;

typedef std::vector<lidxT*>             lidxTPVec;
typedef lidxTPVec::iterator             lidxTPVecI;
typedef lidxTPVec::const_iterator       lidxTPVecCI;

static const lidxT LIDX_UNDEFINED = 0xFFFF;

////////////////////////////////////////////////////////////////////////////////

typedef int                             depT;
typedef std::vector<depT>               depTVec;
typedef depTVec::iterator               depTVecI;
typedef depTVec::const_iterator         depTVecCI;

static const depT DEP_NONE              = 0;
static const depT DEP_STOICH            = 1;
static const depT DEP_RATE              = 2;

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(sim)
END_NAMESPACE(steps)

#endif 
// STEPS_SIM_SHARED_TYPES_HPP

// END
