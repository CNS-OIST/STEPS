/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */


#ifndef STEPS_WMDIRECT_SREAC_HPP
#define STEPS_WMDIRECT_SREAC_HPP 1


// STL headers.

// STEPS headers.
#include "steps/common.h"
#include "steps/wmdirect/kproc.hpp"
#include "steps/solver/sreacdef.hpp"
#include "steps/solver/types.hpp"


////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace wmdirect {

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
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

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

}
}

#endif
// STEPS_WMDIRECT_SREAC_HPP
