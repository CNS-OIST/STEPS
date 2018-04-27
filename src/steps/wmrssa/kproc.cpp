/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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


// Standard library & STL headers.
// #include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/wmrssa/kproc.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace swmrssa = steps::wmrssa;

////////////////////////////////////////////////////////////////////////////////

swmrssa::KProc::KProc(void)
: rExtent(0)
// , pFlags(0)
, pSchedIDX(0)
{

}

////////////////////////////////////////////////////////////////////////////////

swmrssa::KProc::~KProc(void)
{

}

////////////////////////////////////////////////////////////////////////////////

uint swmrssa::KProc::getExtent(void) const
{
    return rExtent;
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::KProc::resetExtent(void)
{
    rExtent = 0;
}

////////////////////////////////////////////////////////////////////////////////

steps::solver::Reacdef * swmrssa::KProc::defr(void) const
{
    // Should only be called on derived object
    assert (false);
}

////////////////////////////////////////////////////////////////////////////////

steps::solver::SReacdef * swmrssa::KProc::defsr(void) const
{
    // Should olny be called on derived object
    assert (false);
}

////////////////////////////////////////////////////////////////////////////////

// END

