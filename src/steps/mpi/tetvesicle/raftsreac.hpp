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
#include "mpi/tetvesicle/raftproxy.hpp"
#include "mpi/tetvesicle/tri_rdef.hpp"
#include "solver/raftsreacdef.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations.
class TetVesicleRDEF;
class TetRDEF;

////////////////////////////////////////////////////////////////////////////////

class RaftSReac: public KProc {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    RaftSReac(solver::RaftSReacdef* rsrdef, TriRDEF* tri);
    ~RaftSReac();

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

    inline solver::RaftSReacdef* def() const noexcept {
        return pRaftSReacdef;
    }

    // This is now read from the def object to allow control of kcsts during
    // simulation for all objects of this type
    inline double kcst() const noexcept {
        return pRaftSReacdef->kcst();
    }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

#if defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif

    // This must be implmented but does nothing since they are setup on the fly
    // for this object
    void setupDeps() override;

    void reset() override;
    double rate(TetVesicleRDEF* solver = nullptr) override;

    void apply(const rng::RNGptr& rng,
               double dt,
               double simtime,
               double period,
               TetVesicleRDEF* solver = nullptr) override;

#if defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

    //////////////// ADDED FOR MPI STEPS ////////////////////

    std::vector<KProc*> const& getLocalUpdVec(int /*direction = -1*/) const override {
        return localUpdVec;
    }
    std::vector<solver::kproc_global_id> const& getRemoteUpdVec(
        int /*direction = -1*/) const override {
        return remoteUpdVec;
    }

    void resetOccupancies() override;

    inline bool getInHost() const noexcept override {
        return pTri->getInHost();
    }

    inline int getHost() const noexcept override {
        return pTri->getHost();
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::RaftSReacdef* pRaftSReacdef;
    TriRDEF* pTri;

    std::vector<KProc*> localUpdVec;
    std::vector<solver::kproc_global_id> remoteUpdVec;

    ////////////////////////////////////////////////////////////////////////

    // Store the current rate as per raft within this tri
    std::map<solver::raft_individual_id, double> pRate_per_raft;

    // For convenience
    double pTotal_rate;
};

}  // namespace steps::mpi::tetvesicle
