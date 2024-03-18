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
#include "mpi/tetvesicle/vesproxy.hpp"
#include "solver/vessreacdef.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations.
class TetRDEF;
class TriRDEF;
class TetVesicleRDEF;

////////////////////////////////////////////////////////////////////////////////

// NEW VesReac needs to be one per tet, but stores rate of rections per vesicle.
// In this way one object needs to be created at the beginning of the simulation
// and rate updated depending on vesicle occupancy
class VesReac: public KProc {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    VesReac(solver::VesSReacdef* vsrdef, TetRDEF* tet);
    ~VesReac();

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

    inline solver::VesSReacdef* def() const noexcept {
        return pVesSReacdef;
    }

    // This is now read from the def object to allow control of kcsts during
    // simulation for all objects of this type
    inline double kcst() const noexcept {
        return pVesSReacdef->kcst();
    }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    // These used to be dynamic but are now constant. HOWEVER DOES NOT CURRENTLY
    // include KPROCS THAT USED TO BE DYNAMIC SUCH AS EXO, RAFTSREAC, RAFTENDO,
    // RAFTDIS,

#if defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif
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
    ////////////////////////////////////////////////////////////////////////

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
        return pTet->getInHost();
    }

    inline int getHost() const noexcept override {
        return pTet->getHost();
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::VesSReacdef* pVesSReacdef;

    TetRDEF* pTet;

    std::vector<math::point3d> pTri_barycentres;

    bool pRate_zero;

    std::vector<KProc*> localUpdVec;
    std::vector<solver::kproc_global_id> remoteUpdVec;

    // Store the current rate as per vesicle within this tet
    std::map<solver::vesicle_individual_id, double> pRate_per_ves;

    // For convenience- the sum of the rates held in pRate_per_ves
    double pTotal_rate;

    // Better to save this bool than keep calling reqSurface in def with loop over
    // all species
    bool pReqSurface;
};

}  // namespace steps::mpi::tetvesicle
