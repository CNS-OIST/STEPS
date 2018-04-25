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


/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */


// Standard library & STL headers.
#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;

////////////////////////////////////////////////////////////////////////////////

smtos::KProc::KProc(void)
: rExtent(0)
, pFlags(0)
, pSchedIDX(0)
, crData()
{
}

////////////////////////////////////////////////////////////////////////////////

smtos::KProc::~KProc(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void smtos::KProc::setActive(bool active)
{
    if (active == true) pFlags &= ~INACTIVATED;
    else pFlags |= INACTIVATED;
}

////////////////////////////////////////////////////////////////////////////////

uint smtos::KProc::getExtent(void) const
{
    return rExtent;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::KProc::resetExtent(void)
{
    rExtent = 0;
}
////////////////////////////////////////////////////////////////////////////////

void smtos::KProc::resetCcst(void) const
{
    // This should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::KProc::c(void) const
{
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::KProc::h(void)
{
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

int smtos::KProc::apply(steps::rng::RNG * rng)
{
    // Should never get called on base object
	AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

int smtos::KProc::apply(steps::rng::RNG * rng, uint nmolcs)
{
    // Should never get called on base object
	AssertLog(false);
}


////////////////////////////////////////////////////////////////////////////////

void smtos::KProc::apply(steps::rng::RNG * rng, double dt, double simtime, double period)
{
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::KProc::resetOccupancies(void)
{
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<smtos::KProc*> const & smtos::KProc::getLocalUpdVec(int direction)
{
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & smtos::KProc::getRemoteUpdVec(int direction)
{
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////
/*
std::vector<smtos::KProc*> const & smtos::KProc::getSharedUpd(void)
{
    // Should never get called on base object
	AssertLog(false);
}
*/

////////////////////////////////////////////////////////////////////////////////
// END
