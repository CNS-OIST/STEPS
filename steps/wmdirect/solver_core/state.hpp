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
// $Id:state.hpp 64 2007-08-20 06:25:41Z stefan $
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_WMDIRECT_SOLVER_CORE_STATE_HPP
#define STEPS_WMDIRECT_SOLVER_CORE_STATE_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

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

    void setupState(void);
    
    /// Should this be a method, or be implemented in wmdirect.cpp? I chose
    /// the first, for unclear reasons ;-)
    ///
    void setupTetmesh(void);

    /// Should this be a method, or be implemented in wmdirect.cpp? I chose
    /// the first, for unclear reasons ;-)
    ///
    void reset(void);

    ////////////////////////////////////////////////////////////////////////

    /// Computes a mesoscopic reaction constant (quantity 'c_mu' in 
    /// Gillespie's papers) from the macroscopic reaction constant K. 
    /// The scaling depends on the order of the reaction.
    ///
    /// This method does not automatically recompute the reaction 
    /// propensity of the corresponding reaction.
    ///
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
// STEPS_WMDIRECT_SOLVER_CORE_STATE_HPP

// END
