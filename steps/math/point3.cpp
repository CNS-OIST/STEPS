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
#include <steps/math/point3.hpp>
#include <steps/math/vector3.hpp>
#include <steps/math/tools.hpp>

////////////////////////////////////////////////////////////////////////////////

// STEPS library.
NAMESPACE_ALIAS(steps::math, smath);
USING(smath, Point3);

////////////////////////////////////////////////////////////////////////////////

std::ostream & smath::operator<< (std::ostream & os, Point3 const & pnt)
{
    os << "(" << pnt.getX() << ", "
              << pnt.getY() << ", "
              << pnt.getZ() << ")";
    return os;
}

////////////////////////////////////////////////////////////////////////////////

// END
