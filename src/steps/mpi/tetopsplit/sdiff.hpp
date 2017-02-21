/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_MPI_TETOPSPLIT_SDIFF_HPP
#define STEPS_MPI_TETOPSPLIT_SDIFF_HPP 1


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
class Tri;

////////////////////////////////////////////////////////////////////////////////

class SDiff
: public steps::mpi::tetopsplit::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SDiff(steps::solver::Diffdef * sdef, steps::mpi::tetopsplit::Tri * tri);
    ~SDiff(void);

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

    inline steps::solver::Diffdef * sdef(void) const
    { return pSDiffdef; }

    double dcst(int direction = -1);
    void setDcst(double d);
    void setDirectionDcst(int direction, double dcst);

    void setupDeps(void);

    bool depSpecTet(uint gidx, steps::mpi::tetopsplit::WmVol * tet);
    bool depSpecTri(uint gidx, steps::mpi::tetopsplit::Tri * tri);

    void reset(void);
    double rate(steps::mpi::tetopsplit::TetOpSplitP * solver = 0);
    double getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * solver = 0);
    
    
    int apply(steps::rng::RNG * rng);
    int apply(steps::rng::RNG * rng, uint nmolcs);

    std::vector<KProc*> const & getLocalUpdVec(int direction = -1);
    std::vector<uint> const & getRemoteUpdVec(int direction = -1);

    bool getInHost(void) {
        return pTri->getInHost();
    }
    
    int getHost(void) {
        return pTri->getHost();
    }
    ////////////////////////////////////////////////////////////////////////

    //inline steps::solver::Reacdef * defr(void) const
    //{ return pReacdef; }

    ////////////////////////////////////////////////////////////////////////


    inline uint getLigLidx(void) {return lidxTri;}

    ////////////////////////////////////////////////////////////////////////

    inline steps::mpi::tetopsplit::Tri* getTri(void) {return pTri;}
private:

    ////////////////////////////////////////////////////////////////////////

    uint                                ligGIdx;
    uint                                lidxTri;
    steps::solver::Diffdef              * pSDiffdef;
    steps::mpi::tetopsplit::Tri         * pTri;

    std::vector<KProc*>                 localUpdVec[3];
    std::vector<KProc*>					localAllUpdVec;
    
    std::vector<uint>                   remoteUpdVec[3];
    std::vector<uint>					remoteAllUpdVec;
    
    // empty vec to return if no update occurs
    
    std::vector<KProc*>					pEmptyvec;
    std::vector<uint>					idxEmptyvec;

    /*
    // Storing the species local index for each neighbouring tri: Needed
    // because neighbours may belong to different patches if we ever
    // implement diffusion boundary for surfaces
    // and therefore have different spec indices
    int                                   pNeighbPatchLidx[3];

     */


    /// Properly scaled diffusivity constant.
    double                              pScaledDcst;
    // Compartmental dcst. Stored for convenience
    double                              pDcst;
    std::map<uint, double>              directionalDcsts;

    double 								pNonCDFSelector[3];
    std::vector<uint> 					pDirections;
    uint				 				pNdirections;

    // A flag to see if the species can move between compartments
    bool                                 pDiffBndActive[3];

    // Flags to store if a direction is a diffusion boundary direction
    bool                                 pDiffBndDirection[3];

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif // STEPS_MPI_TETOPSPLIT_SDIFF_HPP
