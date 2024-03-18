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

#pragma once

// Standard library & STL headers.
#include <fstream>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "math/constants.hpp"
#include "solver/raftendocytosisdef.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations.
class Raft;

////////////////////////////////////////////////////////////////////////////////

class RaftEndocytosis {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    RaftEndocytosis(solver::RaftEndocytosisdef* endodef, Raft* raft);
    RaftEndocytosis(solver::RaftEndocytosisdef* endodef, Raft* raft, std::fstream& cp_file);
    ~RaftEndocytosis();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);
    // restore done by 2nd constructor

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    double c() const noexcept {
        return pCcst;
    }

    void resetCcst();

    inline double kcst() const noexcept {
        return pKcst;
    }

    void setKcst(double k);

    double h() {
        return rate() / pCcst;
    }

    inline bool active() const noexcept {
        return pActive;
    }
    inline bool inactive() const noexcept {
        return !pActive;
    }
    void setActive(bool active);

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    void reset();

    double rate();

    void apply();

    ////////////////////////////////////////////////////////////////////////

    inline solver::RaftEndocytosisdef* endodef() const noexcept {
        return pRaftEndocytosisdef;
    }

    inline bool inner() const noexcept {
        return endodef()->inner();
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::RaftEndocytosisdef* pRaftEndocytosisdef;
    Raft* pRaft;

    /// Properly scaled reaction constant.
    double pCcst;

    /// Store the kcst for convenience
    double pKcst;

    bool pActive;

    uint pFails;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::mpi::tetvesicle
