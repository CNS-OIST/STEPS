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


#ifndef STEPS_MPI_TETOPSPLIT_VDEPTRANS_HPP
#define STEPS_MPI_TETOPSPLIT_VDEPTRANS_HPP 1

////////////////////////////////////////////////////////////////////////////////

// Standard library & STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/constants.hpp"
#include "steps/solver/vdeptransdef.hpp"
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

class VDepTrans
: public steps::mpi::tetopsplit::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    VDepTrans(steps::solver::VDepTransdef * vdtdef, steps::mpi::tetopsplit::Tri * tri);

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

    void setupDeps() override;
    bool depSpecTet(uint gidx, steps::mpi::tetopsplit::WmVol * tet) override;
    bool depSpecTri(uint gidx, steps::mpi::tetopsplit::Tri * tri) override;
    void reset() override;

    double rate(steps::mpi::tetopsplit::TetOpSplitP * solver) override;
    double getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * /*solver*/ = nullptr) const noexcept override
    {return 0.0;}

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
    void apply(steps::rng::RNG * rng, double dt,double simtime, double period);
#pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
    void apply(steps::rng::RNG * rng, double dt,double simtime, double period);
#pragma GCC diagnostic pop
#endif

    std::vector<KProc*> const & getLocalUpdVec(int direction = -1) const override;
    std::vector<uint> const & getRemoteUpdVec(int direction = -1) const override;

    void resetOccupancies() override;
    
    inline bool getInHost() const noexcept override
    { return pTri->getInHost(); }
    
    inline int getHost() const noexcept override
    { return pTri->getHost(); }
    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::VDepTransdef       * pVDepTransdef;
    steps::mpi::tetopsplit::Tri       * pTri;
    
    std::vector<KProc*>                 localUpdVec;
    std::vector<uint>                   remoteUpdVec;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_MPI_TETOPSPLIT_VDEPTRANS_HPP

