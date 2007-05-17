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

// STEPS headers.
#include <steps/common.h>
#include <steps/math/aabbox3.hpp>
#include <steps/math/point3.hpp>

////////////////////////////////////////////////////////////////////////////////

// STEPS library.
NAMESPACE_ALIAS(steps::math, smath);
USING(smath, AABBox3);
USING(smath, Point3);

////////////////////////////////////////////////////////////////////////////////

void AABBox3::reset(void)
{
    pXMin = 1000000; pXMax = -1000000;
    pYMin = 1000000; pYMax = -1000000;
    pZMin = 1000000; pZMax = -1000000;
}

////////////////////////////////////////////////////////////////////////////////

void AABBox3::grow(AABBox3 const & box)
{
    if (box.pXMin < pXMin) pXMin = box.pXMin;
    if (box.pXMax > pXMax) pXMax = box.pXMax;
    if (box.pYMin < pYMin) pYMin = box.pYMin;
    if (box.pYMax > pYMax) pYMax = box.pYMax;
    if (box.pZMin < pZMin) pZMin = box.pZMin;
    if (box.pZMax > pZMax) pZMax = box.pZMax;
}

////////////////////////////////////////////////////////////////////////////////

void AABBox3::grow(Point3 const & point)
{
    if (point.getX() < pXMin) pXMin = point.getX();
    if (point.getX() > pXMax) pXMax = point.getX();
    if (point.getY() < pYMin) pYMin = point.getY();
    if (point.getY() > pYMax) pYMax = point.getY();
    if (point.getZ() < pZMin) pZMin = point.getZ();
    if (point.getZ() > pZMax) pZMax = point.getZ();
}

////////////////////////////////////////////////////////////////////////////////

bool AABBox3::isInside(AABBox3 const & box) const
{
    if (box.pXMin < pXMin) return false;
    if (box.pXMax > pXMax) return false;
    if (box.pYMin < pYMin) return false;
    if (box.pYMax > pYMax) return false;
    if (box.pZMin < pZMin) return false;
    if (box.pZMax > pZMax) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

bool AABBox3::isInside(Point3 const & point) const
{
    if (point.getX() < pXMin) return false;
    if (point.getX() > pXMax) return false;
    if (point.getY() < pYMin) return false;
    if (point.getY() > pYMax) return false;
    if (point.getZ() < pZMin) return false;
    if (point.getZ() > pZMax) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

std::ostream & smath::operator<< (std::ostream & os, AABBox3 const & b)
{
    os << "(" << b.pXMin << ", " << b.pXMax << "; ";
    os <<        b.pYMin << ", " << b.pYMax << "; ";
    os <<        b.pZMin << ", " << b.pZMax << ")";
    return os;
}

////////////////////////////////////////////////////////////////////////////////

// END
