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
#include <fstream>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "math/constants.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "solver/diffdef.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations.
class TetRDEF;
class TriRDEF;

////////////////////////////////////////////////////////////////////////////////

class Diff: public KProc {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Diff(solver::Diffdef* ddef, TetRDEF* tet);
    Diff(const Diff&) = delete;
    ~Diff();

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

    inline solver::Diffdef* def() const noexcept {
        return pDiffdef;
    }

    double dcst(int direction = -1);
    void setDcst(double d);
    void setDirectionDcst(int direction, double dcst);

    // This is necessary for the TetVesicle solver- after vesicle diffusion the
    // rates must be recalculated because volumes may have changed by occupancies.
    void recalcDcst();

    void setupDeps() override;

    void reset() override;
#if defined(__GNUC__) or defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif

    double rate(TetVesicleRDEF* solver = nullptr) override;

#if defined(__GNUC__) or defined(__GNUG__)
#pragma GCC diagnostic pop
#endif
    using KProc::apply;
    int apply(const rng::RNGptr& rng) override;
    int apply(const rng::RNGptr& rng, uint nmolcs) override;

    ////////////////////////////////////////////////////////////////////////

    void setDiffBndActive(uint i, bool active);

    bool getDiffBndActive(uint i) const;

    //////////////// ADDED FOR MPI STEPS ////////////////////

    inline double getScaledDcst() const noexcept {
        return pScaledDcst;
    }

    std::vector<KProc*> const& getLocalUpdVec(int direction = -1) const override;
    std::vector<solver::kproc_global_id> const& getRemoteUpdVec(int direction = -1) const override;

    void resetOccupancies() override{};  // not necessary for this special kproc

    inline solver::spec_local_id getLigLidx() const noexcept {
        return lidxTet;
    }

    inline TetRDEF* getTet() noexcept {
        return pTet;
    }

    // MPI STUFF
    inline bool getInHost() const noexcept override {
        return pTet->getInHost();
    }

    inline int getHost() const noexcept override {
        return pTet->getHost();
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::Diffdef* pDiffdef;
    TetRDEF* pTet;

    /////////////////////// SPEC INDICES ///////////////////////

    // Global index of the diffusing species
    solver::spec_global_id ligGIdx;
    // Local index (to tet) of the diffusing species
    solver::spec_local_id lidxTet;

    // Storing the species local index for each neighbouring tet: Needed
    // because neighbours may belong to different compartments
    // and therefore have different spec indices
    std::array<solver::spec_local_id, 4> pNeighbCompLidx{std::nullopt,
                                                         std::nullopt,
                                                         std::nullopt,
                                                         std::nullopt};

    //////////////////// DIFFUSION COEFFICIENT  ////////////////

    // Compartmental dcst. Stored for convenience
    double pDcst;
    /// Properly scaled diffusivity constant.
    double pScaledDcst;

    std::map<uint, double> directionalDcsts;

    /// Used in selecting which directory the molecule should go.
    std::array<double, 4> pNonCDFSelector{0.0, 0.0, 0.0, 0.0};

    //////////////////// DIFFUSION BOUNDARY  ////////////////

    // A flag to see if the species can move between compartments
    std::array<bool, 4> pDiffBndActive{false, false, false, false};

    // Flags to store if a direction is a diffusion boundary direction
    std::array<bool, 4> pDiffBndDirection{false, false, false, false};

    /////////////////////// OPSPLIT/MPI ////////////////////

    std::vector<KProc*> localUpdVec[4];
    std::vector<KProc*> localAllUpdVec;

    std::vector<solver::kproc_global_id> remoteUpdVec[4];
    std::vector<solver::kproc_global_id> remoteAllUpdVec;

    // empty vec to return if no update occurs
    const std::vector<KProc*> pEmptyvec;
    const std::vector<solver::kproc_global_id> idxEmptyvec;

    std::vector<uint> pDirections;
    uint pNdirections{};

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::mpi::tetvesicle
