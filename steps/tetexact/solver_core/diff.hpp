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
// $Id:diff.hpp 64 2007-08-20 06:25:41Z stefan $
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_TETEXACT_SOLVER_CORE_DIFF_HPP
#define STEPS_TETEXACT_SOLVER_CORE_DIFF_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>
#include <steps/tetexact/solver_core/kproc.hpp>
#include <steps/tetexact/solver_core/sched.hpp>

// Forward declarations.
class DiffDef;
class State;
class Tet;

////////////////////////////////////////////////////////////////////////////////

class Diff
: public KProc
{

public:
    
    ////////////////////////////////////////////////////////////////////////

    Diff(DiffDef * ddef, Tet * tet);
    virtual ~Diff(void);

    ////////////////////////////////////////////////////////////////////////
    
    inline DiffDef * def(void) const
    { return pDiffDef; }
        
    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////
    
    virtual void setupDeps(void);
    virtual bool depSpecTet(uint gidx, Tet * tet);
    virtual void reset(void);
    virtual double rate(void) const;
    virtual SchedIDXVec const & apply(State * s);

    ////////////////////////////////////////////////////////////////////////
    
private:

    DiffDef *                   pDiffDef;
    Tet *                       pTet;
    SchedIDXVec                 pUpdVec[4];
    /// Properly scaled diffusivity constant.
    double                      pScaledDcst;
    /// Used in selecting which directory the molecule should go.
    double                      pCDFSelector[3];

};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_TETEXACT_SOLVER_CORE_DIFF_HPP

// END
