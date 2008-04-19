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

#ifndef STEPS_TETEXACT_SOLVER_CORE_SREAC_HPP
#define STEPS_TETEXACT_SOLVER_CORE_SREAC_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/sim/shared/sreacdef.hpp>
#include <steps/tetexact/solver_core/kproc.hpp>
#include <steps/tetexact/solver_core/sched.hpp>

// Forward declarations.
class State;
class Tri;

////////////////////////////////////////////////////////////////////////////////

class SReac
: public KProc
{

public:
    
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    
    SReac(steps::sim::SReacDef * srdef, Tri * tri);
    virtual ~SReac(void);

    ////////////////////////////////////////////////////////////////////////
    
    inline steps::sim::SReacDef * def(void) const
    { return pSReacDef; }
    
    inline Tri * tri(void) const
    { return pTri; }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////
    
    virtual void setupDeps(void);
    virtual bool depSpecTet(uint gidx, Tet * tet);
    virtual bool depSpecTri(uint gidx, Tri * triangle);
    virtual void reset(void);
    virtual double rate(void) const;
    virtual SchedIDXVec const & apply(State * s);

    ////////////////////////////////////////////////////////////////////////
    
    double c(void) const
    { return pCcst; }
    
    ////////////////////////////////////////////////////////////////////////
    
private:

    ////////////////////////////////////////////////////////////////////////
    
    steps::sim::SReacDef *      pSReacDef;
    Tri *                       pTri;
    SchedIDXVec                 pUpdVec;
    
    /// Properly scaled reaction constant.
    double                      pCcst;
    
};

////////////////////////////////////////////////////////////////////////////////

#endif 
// STEPS_TETEXACT_SOLVER_CORE_SREAC_HPP

////////////////////////////////////////////////////////////////////////////////

// END
