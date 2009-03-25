////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
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
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_WMDIRECT_SREAC_HPP
#define STEPS_WMDIRECT_SREAC_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.

// STEPS headers.
#include <steps/common.h>
#include <steps/wmdirect/kproc.hpp>
#include <steps/solver/sreacdef.hpp>
#include <steps/solver/types.hpp>


////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(wmdirect)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
class Comp;
class Patch;

////////////////////////////////////////////////////////////////////////////////

class SReac
: public steps::wmdirect::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SReac(steps::solver::SReacdef * srdef, Patch * patch);
    ~SReac(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    static const int INACTIVATED = 1;

    bool active(void) const;

    bool inactive(void) const
    { return (! active()); }


    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    void setupDeps(void);
    bool depSpecComp(uint gidx, Comp * comp);
    bool depSpecPatch(uint gidx, Patch * patch);
    void reset(void);
    double rate(void) const;
    std::vector<uint> const & apply(void);

    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::SReacdef * defsr(void) const
    { return pSReacdef; }

    void resetCcst(void);

    double c(void) const
    { return pCcst; }

    double h(void) const
    { return (rate()/pCcst); }

	uint updVecSize(void) const
    { return pUpdVec.size(); }
	
    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::SReacdef           * pSReacdef;
    Patch                             * pPatch;
    std::vector<uint>                   pUpdVec;

    /// Properly scaled reaction constant.
    double                              pCcst;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(wmdirect)
END_NAMESPACE(steps)

#endif
// STEPS_WMDIRECT_SREAC_HPP
