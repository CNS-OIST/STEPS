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

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cassert>
//#include <cunistd>
#include <string>

// STEPS headers.
#include <steps/common.h>
#include <steps/compat/hostname.hpp>

////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);

////////////////////////////////////////////////////////////////////////////////

string steps::getHostname(void)
{
    // Fetch hostname.
    char n[255];
#ifdef NDEBUG
    gethostname(n, 255);
#else
    int res = gethostname(n, 255);
    assert(res == 0);
#endif
    // Copy hostname.
    return string(n);
}

////////////////////////////////////////////////////////////////////////////////

void steps::getHostname(string & name)
{
    // Fetch hostname.
    char n[255];
#ifdef NDEBUG
    gethostname(n, 255);
#else
    int res = gethostname(n, 255);
    assert(res == 0);
#endif
    // Copy hostname.
    name = n;
}

////////////////////////////////////////////////////////////////////////////////

// END
