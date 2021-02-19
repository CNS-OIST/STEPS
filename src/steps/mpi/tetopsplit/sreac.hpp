/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

class SReac: public KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SReac(steps::solver::SReacdef * srdef, steps::mpi::tetopsplit::Tri * tri);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file) override;

    /// restore data
    void restore(std::fstream & cp_file) override;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline double c() const noexcept override
    { return pCcst; }
    void resetCcst() override;

    inline double kcst() const noexcept
    { return pKcst; }
    void setKcst(double k);

    inline double h() override
    { return rate() / pCcst; }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    void setupDeps() override;
    bool depSpecTet(uint gidx, steps::mpi::tetopsplit::WmVol * tet) override;
    bool depSpecTri(uint gidx, steps::mpi::tetopsplit::Tri * tri) override;
    void reset() override;
    double rate(steps::mpi::tetopsplit::TetOpSplitP * solver = nullptr) override;
    inline double getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * /*solver*/ = nullptr) const override
    {return 0.0;}
    
    // We need the works here: dt and simtime needed for Ohmic Currents
    using KProc::apply;
    void apply(const rng::RNGptr &rng, double dt, double simtime, double period) override;

    std::vector<KProc*> const & getLocalUpdVec(int direction = -1) const override;
    std::vector<uint> const & getRemoteUpdVec(int direction = -1) const override;

    void resetOccupancies() override;

    ////////////////////////////////////////////////////////////////////////

    //inline steps::solver::Reacdef * defr() const
    //{ return pReacdef; }

    ////////////////////////////////////////////////////////////////////////
    // mpi
    inline bool getInHost() const noexcept override {
        return pTri->getInHost();
    }
    
    inline int getHost() const noexcept override {
        return pTri->getHost();
    }
    
private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::SReacdef           * pSReacdef;
    steps::mpi::tetopsplit::Tri              * pTri;

    std::vector<KProc*>                 localUpdVec;
    std::vector<uint>                   remoteUpdVec;
    
    /// Properly scaled reaction constant.
    double                              pCcst{0.0};
    // Store the kcst for convenience
    double                              pKcst{0.0};

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_MPI_TETOPSPLIT_SREAC_HPP
