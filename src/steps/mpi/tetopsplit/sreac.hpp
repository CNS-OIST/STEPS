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


#ifndef STEPS_MPI_TETOPSPLIT_SREAC_HPP
#define STEPS_MPI_TETOPSPLIT_SREAC_HPP 1

////////////////////////////////////////////////////////////////////////////////


// Standard library & STL headers.
#include <map>
#include <string>
#include <vector>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/constants.hpp"
#include "steps/solver/sreacdef.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/mpi/tetopsplit/tri.hpp"
//#include "tetopsplit.hpp"


////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace mpi {
 namespace tetopsplit {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Tri;
class TetOpSplitP;

////////////////////////////////////////////////////////////////////////////////

class SReac
: public steps::mpi::tetopsplit::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SReac(steps::solver::SReacdef * srdef, steps::mpi::tetopsplit::Tri * tri);
    ~SReac();

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
    
    // We need the works here: dt and simtime needed for Ohmic Currents
    void apply(steps::rng::RNG * rng, double dt, double simtime, double period);

    std::vector<KProc*> const & getLocalUpdVec(int direction = -1);
    std::vector<uint> const & getRemoteUpdVec(int direction = -1);

    void resetOccupancies();

    ////////////////////////////////////////////////////////////////////////

    //inline steps::solver::Reacdef * defr() const
    //{ return pReacdef; }

    ////////////////////////////////////////////////////////////////////////
    // mpi
    bool getInHost() {
        return pTri->getInHost();
    }
    
    int getHost() {
        return pTri->getHost();
    }
    
private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::SReacdef           * pSReacdef;
    steps::mpi::tetopsplit::Tri              * pTri;

    std::vector<KProc*>                 localUpdVec;
    std::vector<uint>                   remoteUpdVec;
    
    /// Properly scaled reaction constant.
    double                              pCcst;
    // Store the kcst for convenience
    double                              pKcst;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_MPI_TETOPSPLIT_SREAC_HPP
