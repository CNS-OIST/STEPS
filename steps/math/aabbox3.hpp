////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_MATH_AABBOX3_HPP
#define STEPS_MATH_AABBOX3_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cassert>
#include <ostream>

// Boost headers.
#include <steps/boost/scoped_ptr.hpp>
#include <steps/boost/shared_ptr.hpp>
#include <steps/boost/weak_ptr.hpp>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/point3.hpp>
#include <steps/math/vector3.hpp>

START_NAMESPACE(steps)
START_NAMESPACE(math)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.

class AABBox3;

// Auxiliary declarations.

typedef boost::scoped_ptr<AABBox3>              AABBox3ScPtr;
typedef boost::shared_ptr<AABBox3>              AABBox3ShPtr;
typedef boost::weak_ptr<AABBox3>                AABBox3WkPtr;

////////////////////////////////////////////////////////////////////////////////

/// Implements a 3-dimensional, axis-aligned bounding box. Typically used
/// for speeding up geometric comparisons.
///
class AABBox3
{

private:

    ////////////////////////////////////////////////////////////////////////
    // CONSTRAINTS: X
    ////////////////////////////////////////////////////////////////////////
    
    double                                      pXMin;
    
    double                                      pXMax;
    
    ////////////////////////////////////////////////////////////////////////
    // CONSTRAINTS: Y
    ////////////////////////////////////////////////////////////////////////
    
    double                                      pYMin;
    
    double                                      pYMax;
    
    ////////////////////////////////////////////////////////////////////////
    // CONSTRAINTS: Z
    ////////////////////////////////////////////////////////////////////////
    
    double                                      pZMin;
    
    double                                      pZMax;

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Default constructor.
    ///
    AABBox3(void)
    : pXMin(10000000)
    , pXMax(-1000000)
    , pYMin(10000000)
    , pYMax(-1000000)
    , pZMin(10000000)
    , pZMax(-1000000)
    {
    }
    
    /// Constructor.
    ///
    AABBox3(double xmin, double xmax,
            double ymin, double ymax,
            double zmin, double zmax)
    : pXMin(xmin)
    , pXMax(xmax)
    , pYMin(ymin)
    , pYMax(ymax)
    , pZMin(zmin)
    , pZMax(zmax)
    {
        assert(pXMin <= pXMax);
        assert(pYMin <= pYMax);
        assert(pZMin <= pZMax);
    }
    
    /// Destructor.
    ///
    ~AABBox3(void)
    {
    }

    /// Resets the bounding box to include the entire space.
    ///
    void reset(void);
    
    ////////////////////////////////////////////////////////////////////////
    // BOUNDING BOX GROWING
    ////////////////////////////////////////////////////////////////////////
    
    void grow(AABBox3 const & box);
    
    void grow(Point3 const & point);

    ////////////////////////////////////////////////////////////////////////
    // ACCESS TO CONSTRAINTS
    ////////////////////////////////////////////////////////////////////////

    double getXMin(void) const
    {
        return pXMin;
    }
    
    double getXMax(void) const
    {
        return pXMax;
    }
    
    double getYMin(void) const
    {
        return pYMin;
    }
    
    double getYMax(void) const
    {
        return pYMax;
    }
    
    double getZMin(void) const
    {
        return pZMin;
    }
    
    double getZMax(void) const
    {
        return pZMax;
    }
    
    void set(double xmin, double xmax,
             double ymin, double ymax,
             double zmin, double zmax)
    {
        pXMin = xmin;
        pXMax = xmax;
        pYMin = ymin;
        pYMax = ymax;
        pZMin = zmin;
        pZMax = zmax;
    }
    
    void setX(double xmin, double xmax)
    {
        pXMin = xmin;
        pXMax = xmax;
    }
    
    void setXMin(double xmin)
    {
        pXMin = xmin;
    }
    
    void setXMax(double xmax)
    {
        pXMax = xmax;
    }
    
    void setY(double ymin, double ymax)
    {
        pYMin = ymin;
        pYMax = ymax;
    }
    
    void setYMin(double ymin)
    {
        pYMin = ymin;
    }
    
    void setYMax(double ymax)
    {
        pYMax = ymax;
    }
    
    void setZ(double zmin, double zmax)
    {
        pZMin = zmin;
        pZMax = zmax;
    }
    
    void setZMin(double zmin)
    {
        pZMin = zmin;
    }
    
    void setZMax(double zmax)
    {
        pZMax = zmax;
    }
    
    ////////////////////////////////////////////////////////////////////////
    // PREDICATES
    ////////////////////////////////////////////////////////////////////////
    
    /// Checks whether another bounding box is inside this box.
    ///
    bool isInside(AABBox3 const & box) const;
    
    /// Checks whether a point is inside this box.
    ///
    bool isInside(Point3 const & point) const;
    
    ////////////////////////////////////////////////////////////////////////
    // STREAM OUTPUT
    ////////////////////////////////////////////////////////////////////////
    
    friend std::ostream & operator<< (std::ostream & os, AABBox3 const & b);
    
};

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN std::ostream & operator<< (std::ostream & os, AABBox3 const & b);

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(math)
END_NAMESPACE(steps)

#endif
// STEPS_MATH_AABBOX3_HPP

// END
