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

#include <vector>

#include "kproc.hpp"
#include "wmvol.hpp"

#include "math/constants.hpp"
#include "solver/reacdef.hpp"
#include "util/common.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::mpi::tetopsplit {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class WmVol;
class Tri;
class TetOpSplitP;

////////////////////////////////////////////////////////////////////////////////

class Reac: public KProc {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Reac(solver::Reacdef* rdef, mpi::tetopsplit::WmVol* tet);
    Reac(const Reac&) = delete;
    ~Reac() override;

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

    double c() const override {
        return pCcst;
    }
    void resetCcst() override;

    inline double kcst() const {
        return pKcst;
    }
    void setKcst(double k);

    double h() noexcept override {
        return (rate() / pCcst);
    }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    void setupDeps() override;
    /*
    bool depSpecTet(solver::spec_global_id gidx, mpi::tetopsplit::WmVol *
    tet) override; bool depSpecTri(solver::spec_global_id gidx,
    mpi::tetopsplit::Tri * tri) override;
        */
    void reset() override;
    double rate(mpi::tetopsplit::TetOpSplitP* solver = nullptr) override;
    inline double getScaledDcst(
        mpi::tetopsplit::TetOpSplitP* /*solver*/ = nullptr) const noexcept override {
        return 0.0;
    }

    // at the moment we assume that reactions are applied globally so no sync is
    // required
    using KProc::apply;
    void apply(const rng::RNGptr& rng, double dt, double simtime, double period) override;

    std::vector<KProc*> const& getLocalUpdVec(int direction = -1) const override;
    std::vector<solver::kproc_global_id> const& getRemoteUpdVec(int direction = -1) const override;

    void resetOccupancies() override;

    /// MPI
    inline bool getInHost() const noexcept override {
        return pTet->getInHost();
    }

    inline int getHost() const noexcept override {
        return pTet->getHost();
    }

    inline mpi::tetopsplit::WmVol* container() const noexcept {
        return pTet;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::Reacdef* pReacdef;
    mpi::tetopsplit::WmVol* pTet;

    std::vector<KProc*> localUpdVec;
    std::vector<solver::kproc_global_id> remoteUpdVec;

    /// Properly scaled reaction constant.
    double pCcst;
    // Also store the K constant for convenience
    double pKcst;

    ////////////////////////////////////////////////////////////////////////
};

////////////////////////////////////////////////////////////////////////////////

}  // namespace steps::mpi::tetopsplit
