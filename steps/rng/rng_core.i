////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
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

%module rng_core

%import "steps/common.h"

%{
// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

#include <steps/rng/rng.hpp>
%}

START_NAMESPACE(steps)
START_NAMESPACE(rng)

////////////////////////////////////////////////////////////////////////////////

class RNG
{
    
public:

    RNG(uint bufsize);
    virtual ~RNG(void);

    void initialize(ulong const & seed);

    uint get(void);
    
    double getUnfII(void);
    double getUnfIE(void);
    double getUnfEE(void);
    double getUnfIE53(void);
    
    float getStdExp(void);
    virtual double getExp(double lambda);
    long getPsn(float lambda);
    float getStdNrm(void);

protected:
    // Mark this class as abstract.
    virtual void concreteInitialize(ulong seed) = 0;
    virtual void concreteFillBuffer(void) = 0;

};

////////////////////////////////////////////////////////////////////////////////

RNG * create_mt19937(uint bufsize);

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(rng)
END_NAMESPACE(steps)

// END
