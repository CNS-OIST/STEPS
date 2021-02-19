/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

#ifndef STEPS_WMRSSA_REAC_HPP
#define STEPS_WMRSSA_REAC_HPP 1


// STL headers.

// STEPS headers.
#include "steps/common.h"
#include "steps/wmrssa/kproc.hpp"
#include "steps/solver/reacdef.hpp"
#include "steps/solver/types.hpp"


////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace wmrssa {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
class Patch;
class Comp;

////////////////////////////////////////////////////////////////////////////////

class Reac
: public steps::wmrssa::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Reac(steps::solver::Reacdef * rdef, Comp * comp);
    ~Reac() override;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file) override;

    /// restore data
    void restore(std::fstream & cp_file) override;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    static const int INACTIVATED = 1;

    bool active() const;

    inline bool inactive() const
    { return (! active()); }


    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    void setupDeps() override;
    bool depSpecComp(uint gidx, Comp * comp) override;
    bool depSpecPatch(uint gidx, Patch * patch) override;
    void reset() override;
    double rate(steps::wmrssa::PropensityRSSA prssa = steps::wmrssa::CURRENT) override;
    std::vector<uint> const & apply() override;

    inline uint updVecSize() const noexcept override
    { return pUpdVec.size(); }

    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Reacdef * defr() const noexcept override
    { return pReacdef; }

    void resetCcst() override;

    inline double c() const noexcept override
    { return pCcst; }

    inline double propensityLB() const noexcept override
    { return pPropensityLB; }

    inline double h() noexcept override
    { return rate() / pCcst; }

    ////////////////////////////////////////////////////////////////////////

private:

  ////////////////////////////////////////////////////////////////////////

    steps::solver::Reacdef            * pReacdef;
    Comp                              * pComp;
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
// STEPS_WMRSSA_REAC_HPP

// END
