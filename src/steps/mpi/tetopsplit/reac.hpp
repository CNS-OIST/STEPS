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


#ifndef STEPS_MPI_TETOPSPLIT_REAC_HPP
#define STEPS_MPI_TETOPSPLIT_REAC_HPP 1


// Standard library & STL headers.
#include <map>
#include <string>
#include <vector>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/constants.hpp"
#include "steps/solver/reacdef.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "wmvol.hpp"


////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace mpi {
 namespace tetopsplit {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class WmVol;
class Tri;
class TetOpSplitP;

////////////////////////////////////////////////////////////////////////////////

class Reac
: public steps::mpi::tetopsplit::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Reac(steps::solver::Reacdef * rdef, steps::mpi::tetopsplit::WmVol * tet);
    ~Reac();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    double c() const
    { return pCcst; }
    void resetCcst();

    inline double kcst() const
    { return pKcst; }
    void setKcst(double k);

    double h()
    { return (rate()/pCcst); }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    void setupDeps();
    bool depSpecTet(uint gidx, steps::mpi::tetopsplit::WmVol * tet);
    bool depSpecTri(uint gidx, steps::mpi::tetopsplit::Tri * tri);

    void reset();
    double rate(steps::mpi::tetopsplit::TetOpSplitP * solver = 0);
    double getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * solver = 0)
    {return 0.0;}
    
	// at the moment we assume that reactions are applied globally so no sync is required
    void apply(steps::rng::RNG * rng, double dt, double simtime, double period);
    
    std::vector<KProc*> const & getLocalUpdVec(int direction = -1);
    std::vector<uint> const & getRemoteUpdVec(int direction = -1);
	
    void resetOccupancies();
	
    /// MPI
    bool getInHost() {
        return pTet->getInHost();
    }

    int getHost() {
        return pTet->getHost();
    }
    
    steps::mpi::tetopsplit::WmVol* container() {
        return pTet;
    }
    ////////////////////////////////////////////////////////////////////////


private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::Reacdef                              * pReacdef;
    steps::mpi::tetopsplit::WmVol                       * pTet;
    
    std::vector<KProc*>                 localUpdVec;
    std::vector<uint>                   remoteUpdVec;
  
    /// Properly scaled reaction constant.
    double                                                pCcst;
    // Also store the K constant for convenience
    double                                                pKcst;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_MPI_TETOPSPLIT_REAC_HPP
