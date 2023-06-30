/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#pragma once

// STL headers.
#include <limits>
#include <vector>

/*
namespace steps::solver {

typedef unsigned int gidxT;

static const uint GIDX_UNDEFINED = std::numeric_limits<uint>::max();

struct gidxT {
    uint value;
    gidxT(uint value = GIDX_UNDEFINED)
        : value(value) {}
};

typedef std::vector<gidxT> gidxTVec;
typedef gidxTVec::iterator gidxTVecI;
typedef gidxTVec::const_iterator gidxTVecCI;

typedef std::vector<gidxT*> gidxTPVec;
typedef gidxTPVec::iterator gidxTPVecI;
typedef gidxTPVec::const_iterator gidxTPVecCI;

typedef unsigned int lidxT;

////////////////////////////////////////////////////////////////////////////////

static const uint LIDX_UNDEFINED = std::numeric_limits<uint>::max();

struct lidxT {
    uint value;
    lidxT(uint value = LIDX_UNDEFINED)
        : value(value) {}
};

typedef std::vector<lidxT> lidxTVec;
typedef lidxTVec::iterator lidxTVecI;
typedef lidxTVec::const_iterator lidxTVecCI;

typedef std::vector<lidxT*> lidxTPVec;
typedef lidxTPVec::iterator lidxTPVecI;
typedef lidxTPVec::const_iterator lidxTPVecCI;

////////////////////////////////////////////////////////////////////////////////

typedef int depT;
typedef std::vector<depT> depTVec;
typedef depTVec::iterator depTVecI;
typedef depTVec::const_iterator depTVecCI;

static const depT DEP_NONE = 0;
static const depT DEP_STOICH = 1;
static const depT DEP_RATE = 2;

class Compdef;

}  // namespace steps::solver
*/
