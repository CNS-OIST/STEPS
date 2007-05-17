////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/normal3.hpp>
#include <steps/math/vector3.hpp>
#include <steps/math/tools.hpp>

////////////////////////////////////////////////////////////////////////////////

// STEPS library.
NAMESPACE_ALIAS(steps::math, smath);
USING(smath, Normal3);

////////////////////////////////////////////////////////////////////////////////

std::ostream & smath::operator<< (std::ostream & os, Normal3 const & n)
{
    os << "(" << n.getX() << ", "
              << n.getY() << ", "
              << n.getZ() << ")";
    return os;
}

////////////////////////////////////////////////////////////////////////////////

// END
