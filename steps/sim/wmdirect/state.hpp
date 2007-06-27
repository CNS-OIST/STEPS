////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_WMDIRECT_STATE_HPP
#define STEPS_SIM_WMDIRECT_STATE_HPP 1

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/statedef.hpp>

////////////////////////////////////////////////////////////////////////////////

class State
{

public:

    enum PoolFlags
    {
        CLAMPED_POOLFLAG       = 1
    };
    static const uint PoolFlagDefault = 0;
    
    enum ReacFlags
    {
        ACTIVE_REACFLAG         = 1,
        KCONST_REACFLAG         = 2
    };
    static const uint ReacFlagDefault = ACTIVE_REACFLAG | KCONST_REACFLAG;
    
    ////////////////////////////////////////////////////////////////////////
    
    State(void);
    ~State(void);
    
    StateDef * def(void) const
    { return pStateDef; }

    ////////////////////////////////////////////////////////////////////////

    /// Should this be a method, or be implemented in wmdirect.cpp? I chose
    /// the first, for unclear reasons ;-)
    void setup(void);

    /// Should this be a method, or be implemented in wmdirect.cpp? I chose
    /// the first, for unclear reasons ;-)
    void reset(void);

    ////////////////////////////////////////////////////////////////////////

    /// Computes a mesoscopic reaction constant (quantity 'c_mu' in 
    /// Gillespie's papers) from the macroscopic reaction constant K. 
    /// The scaling depends on the order of the reaction.
    ///
    /// This method does not automatically recompute the reaction 
    /// propensity of the corresponding reaction.
    void computeCcst(uint cidx, uint l_ridx);
    
    /// If the reaction is not active, this method just sets a propensity
    /// of zero and returns immediately. 
    ///
    void computeReacProp(uint cidx, uint l_ridx);

    double computeZeroProp(void) const;

    ////////////////////////////////////////////////////////////////////////

    double                      fTime;
    
    uint                        fNSteps;
    
    uint **                     fPoolCount;
    uint **                     fPoolFlags;
    /// Forward reaction constants.
    double **                   fReacKcsts;
    /// Scaled forward reaction constants.
    double **                   fReacCcsts;
    double **                   fReacHs;
    double **                   fReacProps;
    uint **                     fReacFlags;
    uint **                     fReacExtents;
    double *                    fCompVols;
    
    steps::rng::RNG *           fRNG;

private:

    StateDef *                  pStateDef;

};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_WMDIRECT_STATE_HPP

// END
