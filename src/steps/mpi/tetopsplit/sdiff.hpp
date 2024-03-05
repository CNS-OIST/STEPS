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

// Standard library & STL headers.
#include <array>
#include <fstream>
#include <map>
#include <random>
#include <string>
#include <vector>

// STEPS headers.
#include "kproc.hpp"
#include "math/constants.hpp"
#include "solver/diffdef.hpp"
#include "tetopsplit.hpp"

namespace steps::mpi::tetopsplit {

// Forward declarations.
class Tri;

////////////////////////////////////////////////////////////////////////////////

class SDiff: public KProc {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SDiff(solver::SurfDiffdef* sdef, Tri* tri);
    SDiff(const SDiff&) = delete;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) override;

    /// restore data
    void restore(std::fstream& cp_file) override;

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    inline solver::SurfDiffdef* def() const noexcept {
        return pSDiffdef;
    }

    double dcst(int direction = -1);
    void setDcst(double d);
    void setDirectionDcst(int direction, double dcst);

    void setupDeps() override;

    /*
    bool depSpecTet(solver::spec_global_id gidx, WmVol * tet) override;
    bool depSpecTri(solver::spec_global_id gidx, Tri * tri) override;
        */
    void reset() override;
    double rate(TetOpSplitP* solver = nullptr) override;
    inline double getScaledDcst(
        mpi::tetopsplit::TetOpSplitP* /*solver */ = nullptr) const noexcept override {
        return pScaledDcst;
    }

    using KProc::apply;
    int apply(const rng::RNGptr& rng) override;
    int apply(const rng::RNGptr& rng, uint nmolcs) override;

    std::vector<KProc*> const& getLocalUpdVec(int direction = -1) const override;
    std::vector<solver::kproc_global_id> const& getRemoteUpdVec(int direction = -1) const override;

    bool getInHost() const noexcept override {
        return pTri->getInHost();
    }

    int getHost() const noexcept override {
        return pTri->getHost();
    }
    ////////////////////////////////////////////////////////////////////////

    void setSDiffBndActive(uint i, bool active);

    bool getSDiffBndActive(uint i) const;

    ////////////////////////////////////////////////////////////////////////

    inline solver::spec_local_id getLigLidx() const noexcept {
        return lidxTri;
    }

    ////////////////////////////////////////////////////////////////////////

    inline Tri* getTri() const noexcept {
        return pTri;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::spec_local_id lidxTri;
    solver::SurfDiffdef* pSDiffdef;
    Tri* pTri;

    std::vector<KProc*> localUpdVec[3];
    std::vector<KProc*> localAllUpdVec;

    std::vector<solver::kproc_global_id> remoteUpdVec[3];
    std::vector<solver::kproc_global_id> remoteAllUpdVec;

    // empty vec to return if no update occurs

    std::vector<KProc*> pEmptyvec;
    std::vector<solver::kproc_global_id> idxEmptyvec;

    // Storing the species local index for each neighbouring tri: Needed
    // because neighbours may belong to different patches for
    // diffusion boundary for surfaces
    std::array<solver::spec_local_id, 3> pNeighbPatchLidx{std::nullopt, std::nullopt, std::nullopt};

    /// Properly scaled diffusivity constant.
    double pScaledDcst{};
    // Compartmental dcst. Stored for convenience
    double pDcst{};
    std::map<uint, double> directionalDcsts;

    std::array<double, 3> pNonCDFSelector{0.0, 0.0, 0.0};
    std::vector<uint> pDirections;
    uint pNdirections{0};

    // A flag to see if the species can move between compartments
    std::array<bool, 3> pSDiffBndActive{false, false, false};

    // Flags to store if a direction is a diffusion boundary direction
    std::array<bool, 3> pSDiffBndDirection{false, false, false};

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::mpi::tetopsplit
