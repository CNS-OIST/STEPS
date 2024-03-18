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
#include "mpi/tetvesicle/tri_vesraft.hpp"
#include "solver/endocytosisdef.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations.
class TetVesicleVesRaft;

////////////////////////////////////////////////////////////////////////////////

class Endocytosis {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Endocytosis(solver::Endocytosisdef* endodef, std::vector<TriVesRaft*>& tri);
    ~Endocytosis();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    void resetCcst();

    inline double kcst() const noexcept {
        return pKcst;
    }

    void setKcst(double k);

    inline double c() const noexcept {
        return pCcst;
    }

    inline double h() const noexcept {
        return rate() / pCcst;
    }

    inline unsigned long getExtent() const noexcept {
        return rExtent;
    }
    inline void resetExtent() noexcept {
        rExtent = 0;
        pEvents.clear();
    }

    std::vector<solver::EndocytosisEvent> getEvents();

    inline void addEvent(double time,
                         triangle_global_id tidx,
                         solver::vesicle_individual_id vidx) noexcept {
        rExtent++;
        pEvents.emplace_back(time, tidx, vidx);
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
    double rate() const;

    void apply(TetVesicleVesRaft* solver);

    ////////////////////////////////////////////////////////////////////////

    inline solver::Endocytosisdef* endodef() const noexcept {
        return pEndocytosisdef;
    }

    inline bool inner() const noexcept {
        return endodef()->inner();
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::Endocytosisdef* pEndocytosisdef;
    std::vector<TriVesRaft*> pTris;

    /// Properly scaled reaction constant.
    double pCcst;
    // Store the kcst for convenience
    double pKcst;

    unsigned long rExtent;

    std::vector<solver::EndocytosisEvent> pEvents;

    bool pActive;

    std::vector<math::point3d> pPos;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::mpi::tetvesicle
