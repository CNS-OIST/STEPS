////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
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

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cassert>
#include <ctime>
#include <string>
#include <sys/time.h>

// STEPS headers.
#include <steps/common.h>
#include <steps/compat/time.hpp>

////////////////////////////////////////////////////////////////////////////////

USING_NAMESPACE(std);

////////////////////////////////////////////////////////////////////////////////

void steps::breakdownLocalTime
(
    time_t t, 
    uint & year, 
    uint & month, 
    uint & day,
    uint & hour,
    uint & min,
    uint & sec
)
{
    // Fetch the time.
    struct tm t2;
    localtime_r(& t, & t2);
    // Copy.
    year = t2.tm_year + 1900;
    month = t2.tm_mon + 1;
    day = t2.tm_mday;
    hour = t2.tm_hour;
    min = t2.tm_min;
    sec = t2.tm_sec;
}

////////////////////////////////////////////////////////////////////////////////

time_t steps::getLocalTime(void)
{
    return time(0);
}

////////////////////////////////////////////////////////////////////////////////

double steps::getPreciseTime(void)
{
    struct timeval tval;
    gettimeofday(&tval, 0);
    return (tval.tv_sec + tval.tv_usec / 1000000.0);
}

////////////////////////////////////////////////////////////////////////////////

// END
