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

#include "solver/complexreacdef.hpp"
#include "util/common.hpp"
#include "wmdirect/complexevents.hpp"
#include "wmdirect/kproc.hpp"

namespace steps::wmdirect {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
class Patch;
class Comp;
class ComplexReac;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class ComplexReac: public KProc {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    ComplexReac(solver::ComplexReacdef& rdef, Comp& comp);
    ~ComplexReac() override;

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

    bool active() const;

    bool inactive() const {
        return !active();
    }


    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    void setupDeps() override;
    bool depSpecComp(solver::spec_global_id gidx, Comp* comp) override;
    bool depComplexComp(solver::complex_global_id gidx,
                        solver::complex_substate_id sus,
                        Comp* comp) override;
    bool depSpecPatch(solver::spec_global_id gidx, Patch* patch) override;
    void reset() override;
    double rate() const override;
    std::vector<solver::kproc_global_id> const& apply() override;

    uint updVecSize() const noexcept override {
        return pUpdVec.size();
    }

    ////////////////////////////////////////////////////////////////////////

    solver::ComplexReacdef& defcr() const noexcept override {
        return pComplexReacdef;
    }

    void resetCcst() override;

    double c() const noexcept override {
        return pCcst;
    }

    double h() const noexcept override {
        return rate() / pCcst;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::ComplexReacdef& pComplexReacdef;
    Comp& pComp;
    std::vector<solver::kproc_global_id> pUpdVec;
    /// Properly scaled reaction constant.
    double pCcst{0.0};

    // Needs to be mutable because ComplexLHSCandidates::rateMult is not const
    mutable std::map<solver::complex_global_id, ComplexLHSCandidates> candidates;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::wmdirect
