/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

#ifndef STEPS_MPI_TETOPSPLIT_DIFF_HPP
#define STEPS_MPI_TETOPSPLIT_DIFF_HPP 1


// Standard library & STL headers.
#include <array>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <random>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/constants.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/types.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/mpi/tetopsplit/tetopsplit.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace mpi{
namespace tetopsplit{

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Tet;
class Tri;
class WmVol;

////////////////////////////////////////////////////////////////////////////////

class Diff
: public steps::mpi::tetopsplit::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Diff(steps::solver::Diffdef * ddef, steps::mpi::tetopsplit::Tet * tet);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file) override;

    /// restore data
    void restore(std::fstream & cp_file) override;

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Diffdef * def() const
    { return pDiffdef; }

    double dcst(int direction = -1);
    void setDcst(double d);
    void setDirectionDcst(int direction, double dcst);

    void setupDeps() override;
    bool depSpecTet(uint gidx, steps::mpi::tetopsplit::WmVol * tet) override;
    bool depSpecTri(uint gidx, steps::mpi::tetopsplit::Tri * tri) override;
    void reset() override;
    double rate(steps::mpi::tetopsplit::TetOpSplitP * solver = nullptr) override;
    inline double getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * /*solver */ = nullptr) const noexcept override {
        return pScaledDcst;
    }

    using KProc::apply;
    int apply(const rng::RNGptr &rng) override;
    int apply(const rng::RNGptr &rng, uint nmolcs) override;
    
    std::vector<KProc*> const & getLocalUpdVec(int direction = -1) const override;
    std::vector<uint> const & getRemoteUpdVec(int direction = -1) const override;

    ////////////////////////////////////////////////////////////////////////

    void setDiffBndActive(uint i, bool active);

    bool getDiffBndActive(uint i) const;

    ////////////////////////////////////////////////////////////////////////

    //inline steps::solver::Reacdef * defr() const
    //{ return pReacdef; }

    inline uint getLigLidx() const noexcept {return lidxTet;}

    ////////////////////////////////////////////////////////////////////////

    inline steps::mpi::tetopsplit::Tet* getTet() noexcept {return pTet;}
    ////////////////////////////////////////////////////////////////////////

    // MPI STUFF
    inline bool getInHost() const noexcept override {
        return pTet->getInHost();
    }
    
    inline int getHost() const noexcept override {
        return pTet->getHost();
    }
    
private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::Diffdef            * pDiffdef;
    steps::mpi::tetopsplit::Tet       * pTet;
    std::map<uint, double>              directionalDcsts;
    
    std::vector<KProc*>                 localUpdVec[4];
    std::vector<KProc*>					localAllUpdVec;
    
    std::vector<uint>                   remoteUpdVec[4];
    std::vector<uint>					remoteAllUpdVec;

    /// Properly scaled diffusivity constant.
    double                              pScaledDcst{};
    // Compartmental dcst. Stored for convenience
    double                              pDcst{};

    // A flag to see if the species can move between compartments
    std::array<bool, 4>     pDiffBndActive{false, false, false, false};

    // Flags to store if a direction is a diffusion boundary direction
    std::array<bool, 4>     pDiffBndDirection{false, false, false, false};

    std::array<double, 4>   pNonCDFSelector{0.0, 0.0, 0.0, 0.0};

    // Storing the species local index for each neighbouring tet: Needed
    // because neighbours may belong to different compartments
    // and therefore have different spec indices
    std::array<solver::lidxT, 4>  pNeighbCompLidx{solver::LIDX_UNDEFINED, solver::LIDX_UNDEFINED, solver::LIDX_UNDEFINED, solver::LIDX_UNDEFINED};


    // empty vec to return if no update occurs

    const std::vector<KProc*>					pEmptyvec;
    const std::vector<uint>					idxEmptyvec;


    std::vector<uint> 					pDirections;
    uint								pNdirections{};

    uint                                lidxTet;


    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif // STEPS_MPI_TETOPSPLIT_DIFF_HPP
