////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
////////////////////////////////////////////////////////////////////////////////

// $Id$

#ifndef STEPS_MATH_ACCUMULATE_HPP
#define STEPS_MATH_ACCUMULATE_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>

START_NAMESPACE(steps)
START_NAMESPACE(math)

////////////////////////////////////////////////////////////////////////////////

template <typename InputIterator, typename ValueType >
ValueType accumulate(InputIterator begin, InputIterator end, ValueType init)
{
    while (begin != end) 
    {
        init += *begin;
        ++begin;
    }
    return init;
}

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(math)
END_NAMESPACE(steps)

#endif
// STEPS_MATH_ACCUMULATE_HPP

// END
