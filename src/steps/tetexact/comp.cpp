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

// STEPS headers.
#include "comp.hpp"
#include "model/reac.hpp"

// logging
#include "util/error.hpp"

namespace steps::tetexact {


Comp::Comp(solver::Compdef* compdef)
    : pCompdef(compdef)
    , pVol(0.0) {
    AssertLog(pCompdef != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

Comp::~Comp() = default;

////////////////////////////////////////////////////////////////////////////////

void Comp::checkpoint(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void Comp::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void Comp::addTet(WmVol* tet) {
    AssertLog(tet->compdef() == def());
    pTets.push_back(tet);
    pVol += tet->vol();
}

////////////////////////////////////////////////////////////////////////////////

WmVol* Comp::pickTetByVol(double rand01) const {
    if (countTets() == 0) {
        return nullptr;
    }
    if (countTets() == 1) {
        return pTets[0];
    }

    double accum = 0.0;
    double selector = rand01 * vol();
    auto t_end = endTet();
    for (auto t = bgnTet(); t != t_end; ++t) {
        accum += (*t)->vol();
        if (selector < accum) {
            return *t;
        }
    }
    AssertLog(false);
    return nullptr;
}

}  // namespace steps::tetexact
