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


#ifndef STEPS_WMRSSA_SREAC_HPP
#define STEPS_WMRSSA_SREAC_HPP 1


// STL headers.

// STEPS headers.
#include "steps/common.h"
#include "steps/wmrssa/kproc.hpp"
#include "steps/solver/sreacdef.hpp"
#include "steps/solver/types.hpp"


////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace wmrssa {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
class Comp;
class Patch;

////////////////////////////////////////////////////////////////////////////////

class SReac
: public steps::wmrssa::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SReac(steps::solver::SReacdef * srdef, Patch * patch);
    ~SReac();

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

    bool active() const;

    bool inactive() const
    { return (! active()); }


    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    void setupDeps();
    bool depSpecComp(uint gidx, Comp * comp);
    bool depSpecPatch(uint gidx, Patch * patch);
    void reset();
    double rate(steps::wmrssa::PropensityRSSA prssa = steps::wmrssa::CURRENT);
    std::vector<uint> const & apply();

    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::SReacdef * defsr() const
    { return pSReacdef; }

    void resetCcst();

    double c() const
    { return pCcst; }

    double propensityLB() const
    { return pPropensityLB; }

    double h()
    { return (rate()/pCcst); }

    uint updVecSize() const
    { return pUpdVec.size(); }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::SReacdef           * pSReacdef;
    Patch                             * pPatch;
    std::vector<uint>                   pUpdVec;
    std::vector<uint>                   emptyVec;

    /// Properly scaled reaction constant.
    double                              pCcst;
    double                              pPropensityLB;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_WMRSSA_SREAC_HPP
