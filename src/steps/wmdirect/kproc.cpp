/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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
#include "steps/wmdirect/kproc.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace swmd = steps::wmdirect;

////////////////////////////////////////////////////////////////////////////////

swmd::KProc::KProc(void)
: rExtent(0)
// , pFlags(0)
, pSchedIDX(0)
{

}

////////////////////////////////////////////////////////////////////////////////

swmd::KProc::~KProc(void)
{

}

////////////////////////////////////////////////////////////////////////////////

uint swmd::KProc::getExtent(void) const
{
    return rExtent;
}

////////////////////////////////////////////////////////////////////////////////

void swmd::KProc::resetExtent(void)
{
    rExtent = 0;
}

////////////////////////////////////////////////////////////////////////////////

steps::solver::Reacdef * swmd::KProc::defr(void) const
{
    // Should only be called on derived object
    assert (false);
}

////////////////////////////////////////////////////////////////////////////////

steps::solver::SReacdef * swmd::KProc::defsr(void) const
{
    // Should olny be called on derived object
    assert (false);
}

////////////////////////////////////////////////////////////////////////////////

// END

