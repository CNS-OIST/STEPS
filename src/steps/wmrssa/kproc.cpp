/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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
#include "kproc.hpp"

// logging
#include <easylogging++.h>
#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace swmrssa = steps::wmrssa;

////////////////////////////////////////////////////////////////////////////////

steps::solver::Reacdef * swmrssa::KProc::defr() const
{
    // Should only be called on derived object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

steps::solver::SReacdef * swmrssa::KProc::defsr() const
{
    // Should only be called on derived object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

// END

