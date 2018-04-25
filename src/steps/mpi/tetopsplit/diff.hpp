/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <random>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/constants.hpp"
#include "steps/solver/diffdef.hpp"
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
    ~Diff(void);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Diffdef * def(void) const
    { return pDiffdef; }

    double dcst(int direction = -1);
    void setDcst(double d);
    void setDirectionDcst(int direction, double dcst);

    void setupDeps(void);
    bool depSpecTet(uint gidx, steps::mpi::tetopsplit::WmVol * tet);
    bool depSpecTri(uint gidx, steps::mpi::tetopsplit::Tri * tri);
    void reset(void);
    double rate(steps::mpi::tetopsplit::TetOpSplitP * solver = 0);
    inline double getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * solver = 0) {
        return pScaledDcst;
    }

    int apply(steps::rng::RNG * rng);
    int apply(steps::rng::RNG * rng, uint nmolcs);
    
    std::vector<KProc*> const & getLocalUpdVec(int direction = -1);
    std::vector<uint> const & getRemoteUpdVec(int direction = -1);

    ////////////////////////////////////////////////////////////////////////

    void setDiffBndActive(uint i, bool active);

    bool getDiffBndActive(uint i) const;

    ////////////////////////////////////////////////////////////////////////

    //inline steps::solver::Reacdef * defr(void) const
    //{ return pReacdef; }

    inline uint getLigLidx(void) {return lidxTet;}

    ////////////////////////////////////////////////////////////////////////

    inline steps::mpi::tetopsplit::Tet* getTet(void) {return pTet;}
    ////////////////////////////////////////////////////////////////////////

    // MPI STUFF
    bool getInHost(void) {
        return pTet->getInHost();
    }
    
    int getHost(void) {
        return pTet->getHost();
    }
    
private:

    ////////////////////////////////////////////////////////////////////////

    uint                                ligGIdx;
    uint                                lidxTet;
    steps::solver::Diffdef            * pDiffdef;
    steps::mpi::tetopsplit::Tet       * pTet;
    std::map<uint, double>              directionalDcsts;
    
    std::vector<KProc*>                 localUpdVec[4];
    std::vector<KProc*>					localAllUpdVec;
    
    std::vector<uint>                   remoteUpdVec[4];
    std::vector<uint>					remoteAllUpdVec;

    // empty vec to return if no update occurs
    
    std::vector<KProc*>					pEmptyvec;
    std::vector<uint>					idxEmptyvec;

    // Storing the species local index for each neighbouring tet: Needed
    // because neighbours may belong to different compartments
    // and therefore have different spec indices
    int                                   pNeighbCompLidx[4];

    /// Properly scaled diffusivity constant.
    double                              pScaledDcst;
    // Compartmental dcst. Stored for convenience
    double                              pDcst;

    double                              pNonCDFSelector[4];

    // A flag to see if the species can move between compartments
    bool                                 pDiffBndActive[4];

    // Flags to store if a direction is a diffusion boundary direction
    bool                                 pDiffBndDirection[4];


    std::vector<uint> 					pDirections;
    uint								pNdirections;



    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif // STEPS_MPI_TETOPSPLIT_DIFF_HPP
