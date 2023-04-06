/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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
#include "util/common.h"
#include "kproc.hpp"
#include "solver/sreacdef.hpp"
#include "solver/types.hpp"


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
    ~SReac() override;

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

    inline bool inactive() const noexcept
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

    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::SReacdef * defsr() const noexcept override
    { return pSReacdef; }

    void resetCcst() override;

    inline double c() const noexcept override
    { return pCcst; }

    inline double propensityLB() const noexcept override
    { return pPropensityLB; }

    inline double h() noexcept override
    { return (rate()/pCcst); }

    uint updVecSize() const override
    { return pUpdVec.size(); }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::SReacdef           * pSReacdef;
    Patch                             * pPatch;
    std::vector<uint>                   pUpdVec;

    /// Properly scaled reaction constant.
    double                              pCcst{0.0};
    double                              pPropensityLB{0.0};

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_WMRSSA_SREAC_HPP
