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
