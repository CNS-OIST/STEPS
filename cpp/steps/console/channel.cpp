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
#include <sstream>
#include <string>
#include <cstring>
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/compat/time.hpp>
#include <steps/console/channel.hpp>

USING_NAMESPACE(std);
USING_NAMESPACE(steps::console);

////////////////////////////////////////////////////////////////////////////////

Channel::Channel(void)
: ostringstream()
, pStream(&cerr)
{
}

////////////////////////////////////////////////////////////////////////////////

Channel::Channel(ostream * os)
: ostringstream()
, pStream(os)
{
}

////////////////////////////////////////////////////////////////////////////////

Channel::~Channel(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void Channel::setStream(ostream * os)
{
    pStream = os;
}

////////////////////////////////////////////////////////////////////////////////

void Channel::commit(void)
{
    // Discard.
    if (pStream == 0)
    {
        str("");
        return;
    }

    // Convert time.
    time_t timestamp = steps::getLocalTime();
    uint year, month, day, hour, min, sec;
    steps::breakdownLocalTime(timestamp, year, month, day, hour, min, sec);

    // Output date.
    pStream->fill('0');
    pStream->width(4); * pStream << year; pStream->width(0);
    * pStream << "/";
    pStream->width(2); * pStream << month; pStream->width(0);
    * pStream << "/";
    pStream->width(2); * pStream << day; pStream->width(0);
    * pStream << " ";

    // Output time.
    pStream->width(2); * pStream << hour; pStream->width(0);
    * pStream << ":";
    pStream->width(2); * pStream << min; pStream->width(0);
    * pStream << ":";
    pStream->width(2); * pStream << sec; pStream->width(0);
    * pStream << " ";

    // Output message and reset the stringstream.
    * pStream << str() << endl;
    str("");
}

////////////////////////////////////////////////////////////////////////////////

Channel & operator<< (Channel & c, EndMsg const & e)
{
    c.commit();
    return c;
}

////////////////////////////////////////////////////////////////////////////////

// END
