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

#ifndef STEPS_COMPAT_TIME_HPP
#define STEPS_COMPAT_TIME_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cassert>
#include <string>

// STEPS headers.
#include <steps/common.h>

START_NAMESPACE(steps)

////////////////////////////////////////////////////////////////////////////////

/// Converts a given time (t) into its components.
///
STEPS_EXTERN
void breakdownLocalTime
(
    time_t t,
    uint & year,
    uint & month,
    uint & day,
    uint & hour,
    uint & min,
    uint & sec
);

/// Fetch the local time in seconds since 0 hours, 0 minutes, 0
/// seconds, January 1, 1970, Coordinated Universal Time, without
/// including leap seconds. (I.e. using time().) This function is
/// handy for use in date stamps.
///
STEPS_EXTERN
time_t getLocalTime(void);

/// Fetch the current time in seconds and microseconds (fractional part)
/// since 0 hours, 0 minutes, 0 seconds, January 1, 1970, Coordinated
/// Universal Time, without including leap seconds. The actual resolution
/// depends on the system (typically milli- or microseconds) and is
/// returned in GMT. This function is handy for use in timers.
///
STEPS_EXTERN
double getPreciseTime(void);

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(steps)

#endif
// STEPS_COMPAT_TIME_HPP

// END
