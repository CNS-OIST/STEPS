////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

#ifndef STEPS_SIM_SHARED_TYPES_HPP
#define STEPS_SIM_SHARED_TYPES_HPP 1


// STL headers.
#include <vector>

START_NAMESPACE(steps)
START_NAMESPACE(solver)

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

END_NAMESPACE(solver)
END_NAMESPACE(steps)

#endif
// STEPS_SIM_SHARED_TYPES_HPP

// END

