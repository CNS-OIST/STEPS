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


// Standard library & STL headers.
#include <map>
#include <vector>

// STEPS headers.
#include "comp.hpp"
// logging
#include "util/error.hpp"
#include <easylogging++.h>
////////////////////////////////////////////////////////////////////////////////

namespace stode = steps::tetode;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

stode::Comp::Comp(steps::solver::Compdef * compdef)
: pCompdef(compdef)
{
    AssertLog(pCompdef != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

stode::Comp::~Comp()
= default;

////////////////////////////////////////////////////////////////////////////////

void stode::Comp::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&pVol), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stode::Comp::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&pVol), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stode::Comp::addTet(stode::Tet * tet)
{
    AssertLog(tet->compdef() == def());
    auto lidx = static_cast<index_t>(pTets.size());
    pTets.push_back(tet);
    pTets_GtoL.emplace(tet->idx(), lidx);
    pVol+=tet->vol();
}

////////////////////////////////////////////////////////////////////////////////

stode::Tet * stode::Comp::getTet(tetrahedron_id_t lidx)
{
    AssertLog(lidx < static_cast<index_t>(pTets.size()));
    return pTets[lidx.get()];
}


////////////////////////////////////////////////////////////////////////////////

steps::tetrahedron_id_t stode::Comp::getTet_GtoL(tetrahedron_id_t gidx)
{
    auto lidx_it = pTets_GtoL.find(gidx);
    AssertLog(lidx_it != pTets_GtoL.end());
    return lidx_it->second;
}

////////////////////////////////////////////////////////////////////////////////

// END
