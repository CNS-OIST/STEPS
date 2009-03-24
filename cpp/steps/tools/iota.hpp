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

#ifndef STEPS_TOOLS_IOTA_HPP
#define STEPS_TOOLS_IOTA_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>

START_NAMESPACE(steps)
START_NAMESPACE(tools)

////////////////////////////////////////////////////////////////////////////////

/// Provides a very handy algorithm, taken from the non-standard 
/// additions of the SGI/STL.
///
template <typename ForwardIterator, typename ValueType >
void iota(ForwardIterator begin, ForwardIterator end, ValueType value)
{
    while ( begin != end ) 
    {
        *begin = value;
        ++begin;
        ++value;
    }
}

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tools)
END_NAMESPACE(steps)

#endif
// STEPS_TOOLSIOTA_HPP

// END
