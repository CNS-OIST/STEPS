////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_MATH_IOTA_HPP
#define STEPS_MATH_IOTA_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>

START_NAMESPACE(steps)
START_NAMESPACE(math)

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

END_NAMESPACE(math)
END_NAMESPACE(steps)

#endif
// STEPS_MATH_IOTA_HPP

// END
