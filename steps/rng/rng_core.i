////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

%module rng_core

%import "steps/common.h"

%{
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
