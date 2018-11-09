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
#include "steps/error.hpp"
#include "steps/wmdirect/kproc.hpp"
// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace swmd = steps::wmdirect;

////////////////////////////////////////////////////////////////////////////////

swmd::KProc::KProc()
: rExtent(0)
// , pFlags(0)
, pSchedIDX(0)
{

}

////////////////////////////////////////////////////////////////////////////////

swmd::KProc::~KProc()
= default;

////////////////////////////////////////////////////////////////////////////////

uint swmd::KProc::getExtent() const
{
    return rExtent;
}

////////////////////////////////////////////////////////////////////////////////

void swmd::KProc::resetExtent()
{
    rExtent = 0;
}

////////////////////////////////////////////////////////////////////////////////

steps::solver::Reacdef * swmd::KProc::defr() const
{
    // Should only be called on derived object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

steps::solver::SReacdef * swmd::KProc::defsr() const
{
    // Should olny be called on derived object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

// END

