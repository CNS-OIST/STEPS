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

#pragma once

// STEPS headers.
#include "kproc.hpp"
#include "solver/sreacdef.hpp"

namespace steps::wmrssa {

// Forward declarations
class Comp;
class Patch;

class SReac: public KProc {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SReac(solver::SReacdef* srdef, Patch* patch);
    ~SReac() override;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) override;

    /// restore data
    void restore(std::fstream& cp_file) override;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    static const int INACTIVATED = 1;

    bool active() const;

    inline bool inactive() const noexcept {
        return !active();
    }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    void setupDeps() override;
    bool depSpecComp(solver::spec_global_id gidx, Comp* comp) override;
    bool depSpecPatch(solver::spec_global_id gidx, Patch* patch) override;
    void reset() override;
    double rate(wmrssa::PropensityRSSA prssa = wmrssa::CURRENT) override;
    std::vector<solver::kproc_global_id> const& apply() override;

    ////////////////////////////////////////////////////////////////////////

    inline solver::SReacdef* defsr() const noexcept override {
        return pSReacdef;
    }

    void resetCcst() override;

    inline double c() const noexcept override {
        return pCcst;
    }

    inline double propensityLB() const noexcept override {
        return pPropensityLB;
    }

    inline double h() noexcept override {
        return rate() / pCcst;
    }

    uint updVecSize() const override {
        return pUpdVec.size();
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::SReacdef* pSReacdef;
    Patch* pPatch;
    std::vector<solver::kproc_global_id> pUpdVec;

    /// Properly scaled reaction constant.
    double pCcst{0.0};
    double pPropensityLB{0.0};

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::wmrssa
