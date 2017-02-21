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


#ifndef STEPS_MPI_TETOPSPLIT_VDEPSREAC_HPP
#define STEPS_MPI_TETOPSPLIT_VDEPSREAC_HPP 1

////////////////////////////////////////////////////////////////////////////////

// Standard library & STL headers.
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/constants.hpp"
#include "steps/solver/vdepsreacdef.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
//#include "tetopsplit.hpp"
#include "steps/mpi/tetopsplit/tri.hpp"
////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace mpi {
 namespace tetopsplit {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Tri;
class TetOpSplitP;

////////////////////////////////////////////////////////////////////////////////

class VDepSReac
: public steps::mpi::tetopsplit::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    VDepSReac(steps::solver::VDepSReacdef * vdsrdef, steps::mpi::tetopsplit::Tri * tri);
    ~VDepSReac(void);

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

    void setupDeps(void);
    bool depSpecTet(uint gidx, steps::mpi::tetopsplit::WmVol * tet);
    bool depSpecTri(uint gidx, steps::mpi::tetopsplit::Tri * tri);
    void reset(void);

    double rate(steps::mpi::tetopsplit::TetOpSplitP * solver = 0);
    double getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * solver = 0)
    {return 0.0;}
    
    void apply(steps::rng::RNG * rng, double dt, double simtime, double period);

    std::vector<KProc*> const & getLocalUpdVec(int direction = -1);
    std::vector<uint> const & getRemoteUpdVec(int direction = -1);
    
    void resetOccupancies(void);
    
    bool getInHost(void) {
        return pTri->getInHost();
    }
    
    int getHost(void) {
        return pTri->getHost();
    }
    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::VDepSReacdef       * pVDepSReacdef;
    steps::mpi::tetopsplit::Tri       * pTri;

    std::vector<KProc*>                 localUpdVec;
    std::vector<uint>                   remoteUpdVec;

    // The information about the size of the comaprtment or patch, and the
    // dimensions. Important for scaling the constant.
    // As volumes and areas currently don't change this can be stored as
    // a constant.
    double                                 pScaleFactor;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_MPI_TETOPSPLIT_VDEPSREAC_HPP
