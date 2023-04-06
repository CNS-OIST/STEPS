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


#ifndef STEPS_MPI_TETOPSPLIT_SDIFF_HPP
#define STEPS_MPI_TETOPSPLIT_SDIFF_HPP 1


// Standard library & STL headers.
#include <array>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <random>

// STEPS headers.
#include "kproc.hpp"
#include "tetopsplit.hpp"
#include "util/common.h"
#include "math/constants.hpp"
#include "solver/diffdef.hpp"
#include "solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace mpi{
namespace tetopsplit{

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Tri;

////////////////////////////////////////////////////////////////////////////////

class SDiff: public KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SDiff(steps::solver::Diffdef * sdef, Tri * tri);

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

    inline steps::solver::Diffdef * def() const noexcept
    { return pSDiffdef; }

    double dcst(int direction = -1);
    void setDcst(double d);
    void setDirectionDcst(int direction, double dcst);

    void setupDeps() override;

    bool depSpecTet(uint gidx, WmVol * tet) override;
    bool depSpecTri(uint gidx, Tri * tri) override;

    void reset() override;
    double rate(TetOpSplitP * solver = nullptr) override;
    double getScaledDcst(TetOpSplitP * solver = nullptr) const override;

    using KProc::apply;
    int apply(const rng::RNGptr &rng) override;
    int apply(const rng::RNGptr &rng, uint nmolcs) override;

    std::vector<KProc*> const & getLocalUpdVec(int direction = -1) const override;
    std::vector<uint> const & getRemoteUpdVec(int direction = -1) const override;

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

    inline uint getLigLidx() const noexcept {return lidxTri;}

    ////////////////////////////////////////////////////////////////////////

    inline Tri* getTri() const noexcept {return pTri;}

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    uint                                lidxTri;
    steps::solver::Diffdef              * pSDiffdef;
    Tri         * pTri;

    std::vector<KProc*>                 localUpdVec[3];
    std::vector<KProc*>					localAllUpdVec;

    std::vector<uint>                   remoteUpdVec[3];
    std::vector<uint>					remoteAllUpdVec;

    // empty vec to return if no update occurs

    std::vector<KProc*>					pEmptyvec;
    std::vector<uint>					idxEmptyvec;


    // Storing the species local index for each neighbouring tri: Needed
    // because neighbours may belong to different patches for
    // diffusion boundary for surfaces
    std::array<solver::lidxT, 3>  pNeighbPatchLidx{solver::LIDX_UNDEFINED, solver::LIDX_UNDEFINED, solver::LIDX_UNDEFINED};

    /// Properly scaled diffusivity constant.
    double                              pScaledDcst{};
    // Compartmental dcst. Stored for convenience
    double                              pDcst{};
    std::map<uint, double>              directionalDcsts;

    std::array<double, 3>  				pNonCDFSelector{0.0, 0.0, 0.0};
    std::vector<uint> 					pDirections;
    uint				 				pNdirections{0};

    // A flag to see if the species can move between compartments
    std::array<bool, 3>     pSDiffBndActive{false, false, false};

    // Flags to store if a direction is a diffusion boundary direction
    std::array<bool, 3>     pSDiffBndDirection{false, false, false};

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif // STEPS_MPI_TETOPSPLIT_SDIFF_HPP
