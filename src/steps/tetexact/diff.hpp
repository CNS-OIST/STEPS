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
#include <string>
#include <vector>

// STEPS headers.
#include "kproc.hpp"
#include "util/common.hpp"
// #include "tetexact.hpp"
#include "math/constants.hpp"
#include "solver/diffdef.hpp"
#include "solver/types.hpp"

namespace steps::tetexact {

// Forward declarations.
class Tet;
class Tri;
class WmVol;

class Diff: public KProc {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Diff(solver::Diffdef* ddef, steps::tetexact::Tet* tet);
    Diff(const Diff&) = delete;
    virtual ~Diff();

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

    void setupDeps() override;
    bool depSpecTet(solver::spec_global_id gidx, steps::tetexact::WmVol* tet) override;
    bool depSpecTri(solver::spec_global_id gidx, steps::tetexact::Tri* tri) override;
    void reset() override;
    double rate(steps::tetexact::Tetexact* solver = nullptr) override;
    std::vector<KProc*> const& apply(const rng::RNGptr& rng, double dt, double simtime) override;

    uint updVecSize() const override;

    ////////////////////////////////////////////////////////////////////////

    void setDiffBndActive(uint i, bool active);

    bool getDiffBndActive(uint i) const;

    ////////////////////////////////////////////////////////////////////////

    // inline solver::Reacdef * defr() const
    //{ return pReacdef; }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::spec_local_id lidxTet;
    solver::Diffdef* pDiffdef;
    steps::tetexact::Tet* pTet;
    std::vector<KProc*> pUpdVec[4];
    std::map<uint, double> directionalDcsts;

    /// Properly scaled diffusivity constant.
    double pScaledDcst{};
    // Compartmental dcst. Stored for convenience
    double pDcst{};

    // A flag to see if the species can move between compartments
    std::array<bool, 4> pDiffBndActive{false, false, false, false};

    // Flags to store if a direction is a diffusion boundary direction
    std::array<bool, 4> pDiffBndDirection{false, false, false, false};

    std::array<double, 3> pCDFSelector{0.0, 0.0, 0.0};

    // Storing the species local index for each neighbouring tet: Needed
    // because neighbours may belong to different compartments
    // and therefore have different spec indices
    std::array<solver::spec_local_id, 4> pNeighbCompLidx{std::nullopt,
                                                         std::nullopt,
                                                         std::nullopt,
                                                         std::nullopt};
};

}  // namespace steps::tetexact
