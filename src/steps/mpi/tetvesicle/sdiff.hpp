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
class TriRDEF;

////////////////////////////////////////////////////////////////////////////////

class SDiff: public KProc {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SDiff(solver::SurfDiffdef* sdef, TriRDEF* tri);
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

    inline solver::SurfDiffdef* sdef() const noexcept {
        return pSDiffdef;
    }

    double dcst(int direction = -1);
    void setDcst(double d);
    void setDirectionDcst(int direction, double dcst);

#if defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif
    void setupDeps() override;
    void reset() override;
    double rate(TetVesicleRDEF* solver = nullptr) override;
    int apply(const rng::RNGptr& rng) override;
    int apply(const rng::RNGptr& rng, uint nmolcs) override;
#if defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

    ////////////////////////////////////////////////////////////////////////

    void setSDiffBndActive(uint i, bool active);

    bool getSDiffBndActive(uint i) const;

    //////////////// ADDED FOR MPI STEPS ////////////////////

    inline double getScaledDcst() const noexcept {
        return pScaledDcst;
    }

    std::vector<KProc*> const& getLocalUpdVec(int direction = -1) const override;
    std::vector<solver::kproc_global_id> const& getRemoteUpdVec(int direction = -1) const override;

    void resetOccupancies() override{};  // not necessary for this special kproc

    inline solver::spec_local_id getLigLidx() const noexcept {
        return lidxTri;
    }

    inline TriRDEF* getTri() const noexcept {
        return pTri;
    }

    bool getInHost() const noexcept override {
        return pTri->getInHost();
    }

    int getHost() const noexcept override {
        return pTri->getHost();
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::SurfDiffdef* pSDiffdef;
    TriRDEF* pTri;

    /////////////////////// SPEC INDICES ///////////////////////

    solver::spec_global_id ligGIdx;
    solver::spec_local_id lidxTri;

    // Storing the species local index for each neighbouring tri: Needed
    // because neighbours may belong to different patches
    // and therefore have different spec indices
    std::array<solver::spec_local_id, 3> pNeighbPatchLidx{std::nullopt, std::nullopt, std::nullopt};

    //////////////////// DIFFUSION COEFFICIENT  ////////////////

    // Compartmental dcst. Stored for convenience
    double pDcst;

    /// Properly scaled diffusivity constant.
    double pScaledDcst;

    std::map<uint, double> directionalDcsts;

    std::array<double, 3> pNonCDFSelector{0.0, 0.0, 0.0};

    //////////////////// DIFFUSION BOUNDARY  ////////////////

    // A flag to see if the species can move between compartments
    std::array<bool, 3> pSDiffBndActive{false, false, false};

    // Flags to store if a direction is a diffusion boundary direction
    std::array<bool, 3> pSDiffBndDirection{false, false, false};

    /////////////////////// OPSPLIT/MPI ////////////////////

    std::vector<KProc*> localUpdVec[3];
    std::vector<KProc*> localAllUpdVec;

    std::vector<solver::kproc_global_id> remoteUpdVec[3];
    std::vector<solver::kproc_global_id> remoteAllUpdVec;

    // empty vec to return if no update occurs
    std::vector<KProc*> pEmptyvec;
    std::vector<solver::kproc_global_id> idxEmptyvec;

    std::vector<uint> pDirections;
    uint pNdirections{};

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::mpi::tetvesicle
