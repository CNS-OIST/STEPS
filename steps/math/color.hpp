////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_MATH_COLOR_HPP
#define STEPS_MATH_COLOR_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cmath>

// Boost headers.
#include <steps/boost/scoped_ptr.hpp>
#include <steps/boost/shared_ptr.hpp>
#include <steps/boost/weak_ptr.hpp>

// STEPS headers.
#include <steps/math/tools.hpp>

START_NAMESPACE(steps)
START_NAMESPACE(math)

////////////////////////////////////////////////////////////////////////////////

class RGBColor
{

private:

    double                                      pR;
    
    double                                      pG;
    
    double                                      pB;

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Default constructor.
    ///
    RGBColor(void) { }
    
    /// Constructor.
    ///
    RGBColor(double r, double g, double b)
    : pR(r), pG(g), pB(b) { }

    /// Destructor.
    ///
    ~RGBColor(void) { }
    
    ////////////////////////////////////////////////////////////////////////
    // COMPONENT ACCESS
    ////////////////////////////////////////////////////////////////////////
    
    double red(void) const
    {
        return pR;
    }
    
    double green(void) const
    {
        return pG;
    }
    
    double blue(void) const
    {
        return pB;
    }

    void red(double r)
    {
        pR = r;
    }
    
    void green(double g)
    {
        pG = g;
    }
    
    void blue(double b)
    {
        pB = b;
    }

    void set(double r, double g, double b)
    {
        pR = r;
        pG = g;
        pB = b;
    }
    
    ////////////////////////////////////////////////////////////////////////
    // ASSIGNMENT
    ////////////////////////////////////////////////////////////////////////

    RGBColor & operator= (float c[3])
    {
        pR = c[0];
        pG = c[1];
        pB = c[2];
        return (*this);
    }

    RGBColor & operator= (double c[3])
    {
        pR = c[0];
        pG = c[1];
        pB = c[2];
        return (*this);
    }
    
    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

inline void clampRed(RGBColor & color, double min, double max)
{
    if (color.red() < min) color.red(min);
    if (color.red() > max) color.red(max);
}

inline void clampGreen(RGBColor & color, double min, double max)
{
    if (color.green() < min) color.green(min);
    if (color.green() > max) color.green(max);
}

inline void clampBlue(RGBColor & color, double min, double max)
{
    if (color.blue() < min) color.blue(min);
    if (color.blue() > max) color.blue(max);
}

inline void clamp(RGBColor & color, double min, double max)
{
    clampRed(color, min, max);
    clampGreen(color, min, max);
    clampBlue(color, min, max);
}

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(math)
END_NAMESPACE(steps)

#endif
// STEPS_MATH_COLOR_HPP

// END
