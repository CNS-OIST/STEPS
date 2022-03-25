/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_MPI_TETOPSPLIT_GHKCURR_HPP
#define STEPS_MPI_TETOPSPLIT_GHKCURR_HPP 1

////////////////////////////////////////////////////////////////////////////////

// Standard library & STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include "kproc.hpp"
#include "tri.hpp"
#include "util/common.h"
#include "math/constants.hpp"
#include "solver/ghkcurrdef.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace mpi {
namespace tetopsplit {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Tri;
class TetOpSplitP;

////////////////////////////////////////////////////////////////////////////////

class GHKcurr
: public KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    GHKcurr(steps::solver::GHKcurrdef * ghkdef, steps::mpi::tetopsplit::Tri * tri);

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
    inline double getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * /*solver*/ = nullptr) const noexcept override
    {return 0.0;}

    using KProc::apply;
    void apply(const rng::RNGptr &rng, double dt, double simtime, double period) override;
    std::vector<KProc*> const & getLocalUpdVec(int direction = -1) const override;
    std::vector<uint> const & getRemoteUpdVec(int direction = -1) const override;

    void resetOccupancies() override;

    inline bool efflux() const noexcept
    { return pEffFlux; }

    void setEffFlux(bool efx) noexcept
    { pEffFlux = efx; }

    bool getInHost() const noexcept override {
        return pTri->getInHost();
    }

    int getHost() const noexcept override {
        return pTri->getHost();
    }
    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::GHKcurrdef         * pGHKcurrdef;
    steps::mpi::tetopsplit::Tri              * pTri;

    std::vector<KProc*>                 localUpdVec;
    std::vector<uint>                   remoteUpdVec;
    // Flag if flux is outward, positive flux (true) or inward, negative flux (false)
    bool                                pEffFlux;

    ////////////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_MPI_TETOPSPLIT_GHKCURR_HPP
