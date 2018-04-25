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


// STL headers.
#include <string>
#include <cassert>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/types.hpp"
#include "steps/error.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/chandef.hpp"
#include "steps/model/chan.hpp"
#include "steps/model/chanstate.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

ssolver::Chandef::Chandef(Statedef * sd, uint idx, steps::model::Chan * c)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pChanStates(0)
, pNChanStates(0)
, pChanStatesVec()
, pSetupdone(false)
{
    AssertLog(pStatedef != 0);
    AssertLog(c != 0);
    pName = c->getID();

    pChanStatesVec = c->getAllChanStates();
    pNChanStates = pChanStatesVec.size();
    if (pNChanStates == 0) return;

    pChanStates = new uint[pNChanStates];
    std::fill_n(pChanStates, pNChanStates, GIDX_UNDEFINED);

    ////// anything else????
}

////////////////////////////////////////////////////////////////////////////////

ssolver::Chandef::~Chandef(void)
{
    if (pNChanStates > 0) delete[] pChanStates;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Chandef::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&pNChanStates, sizeof(uint));
    cp_file.write((char*)pChanStates, sizeof(uint) * pNChanStates);
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Chandef::restore(std::fstream & cp_file)
{
    if (pNChanStates > 0) delete[] pChanStates;

    cp_file.read((char*)&pNChanStates, sizeof(uint));
    pChanStates = new uint[pNChanStates];
    cp_file.read((char*)pChanStates, sizeof(uint) * pNChanStates);
}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::Chandef::name(void) const
{
    return pName;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Chandef::setup(void)
{

    AssertLog(pSetupdone == false);
    AssertLog(pChanStatesVec.size() == nchanstates());
    for (uint i = 0; i < nchanstates(); ++i)
    {
        uint gidx = pStatedef->getSpecIdx(pChanStatesVec[i]);
        pChanStates[i] = gidx;
    }

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////


// END



