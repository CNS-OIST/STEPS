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

#ifndef STEPS_TETEXACT_SOLVER_CORE_STATE_HPP
#define STEPS_TETEXACT_SOLVER_CORE_STATE_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/tetexact/solver_core/sched.hpp>
#include <steps/tetexact/solver_core/tet.hpp>

////////////////////////////////////////////////////////////////////////////////

class State
{
    
public:
    
    ////////////////////////////////////////////////////////////////////////
    
    /// Default constructor.
    ///
    State(void);
    
    /// Destructor.
    ///
    ~State(void);
    
    /// Return the state definition.
    ///
    inline StateDef * def(void) const
    { return pStateDef; }

    inline Sched * sched(void) const
    { return pSched; }
    
    ////////////////////////////////////////////////////////////////////////
    
    void setupState(void);
    
    void setupTetmesh(void);
    
    void reset(void);
    
    ////////////////////////////////////////////////////////////////////////
    
    void step(void);
    
    void run(double maxt);
    
    ////////////////////////////////////////////////////////////////////////
    
    inline void setRNG(steps::rng::RNG * rng)
    { pRNG = rng; }
    
    inline steps::rng::RNG * rng(void) const
    { return pRNG; }
    
    ////////////////////////////////////////////////////////////////////////
    
    inline double time(void) const
    { return pTime; }
    
    inline uint nsteps(void) const
    { return pNSteps; }
    
    ////////////////////////////////////////////////////////////////////////
    
    uint addTet
    (
        CompDef * cdef, double vol, 
        double a1, double a2, double a3, double a4,
        double d1, double d2, double d3, double d4
    );
    
    inline Tet * tet(uint tidx) const
    { return pTets[tidx]; }
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    
    void executeStep(KProc * kp, double dt);
    
    ////////////////////////////////////////////////////////////////////////
    
    StateDef *                  pStateDef;
    
    steps::rng::RNG *           pRNG;
    
    ////////////////////////////////////////////////////////////////////////
    // TIME
    ////////////////////////////////////////////////////////////////////////
    
    inline void incTime(double dt)
    { pTime += dt; }
    
    inline void setTime(double t)
    { pTime = t; }
    
    inline void resetTime(void)
    { pTime = 0.0; }
    
    double                      pTime;
    
    ////////////////////////////////////////////////////////////////////////
    // DISCRETE STEP COUNTER
    ////////////////////////////////////////////////////////////////////////
    
    inline void incNSteps(uint i = 1)
    { pNSteps += i; }
    
    inline void resetNSteps(void)
    { pNSteps = 0; }
    
    uint                        pNSteps;
    
    ////////////////////////////////////////////////////////////////////////
    
    Sched *                     pSched;
    
    ////////////////////////////////////////////////////////////////////////
    // THE MESH ELEMENTS
    ////////////////////////////////////////////////////////////////////////
    
    std::vector<Tet *>          pTets;
    
    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_TETEXACT_SOLVER_CORE_STATE_HPP

// END
