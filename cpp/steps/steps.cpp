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
// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>

// STEPS headers.
#include <steps/common.h>
#include <steps/steps.hpp>
#include <steps/console.hpp>

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::console, scons);

////////////////////////////////////////////////////////////////////////////////

static bool initdone = false;

////////////////////////////////////////////////////////////////////////////////

void steps::init(void)
{
    // Pre-initialization code.
    if (initdone == true) return;
    steps::console::init();
    scons::dbg() << "Initializing STEPS" << scons::endm;

    // Initialize other stuff.
    // ...

    // Post-initialization code.
    scons::info() << "STEPS - STochastic Engine for Pathway Simulation";
    scons::info() << scons::endm;
    scons::info() << "Copyright (C) 2005-2008 Stefan Wils. ";
    scons::info() << "All rights reserved." << scons::endm;
    initdone = false;
}

////////////////////////////////////////////////////////////////////////////////

void steps::finish(void)
{
    // Pre-finishing code.
    if (initdone == false) return;
    steps::dbg() << "Cleaning up STEPS library" << scons::endm;

    // Clean up other stuff.
    // ...

    // Post-initialization code.
    scons::info() << "STEPS library successfully cleaned up" << scons::endm;
    steps::console::finish();
    initdone = true;
}

////////////////////////////////////////////////////////////////////////////////

// END
*/
