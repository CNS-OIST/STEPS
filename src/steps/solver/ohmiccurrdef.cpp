/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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
#include <iostream>
#include <sstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/types.hpp"
#include "steps/error.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/ohmiccurrdef.hpp"
#include "steps/model/chanstate.hpp"
#include "steps/model/ohmiccurr.hpp"

namespace ssolver = steps::solver;
namespace smod = steps::model;

////////////////////////////////////////////////////////////////////////////////

ssolver::OhmicCurrdef::OhmicCurrdef(Statedef * sd, uint gidx, smod::OhmicCurr * oc)
: pStatedef(sd)
, pIdx(gidx)
, pName()
, pSetupdone(false)
, pChanState()
, pG(0.0)
, pERev(0.0)
, pSpec_DEP(0)
, pSpec_CHANSTATE(GIDX_UNDEFINED)
{
    assert(pStatedef != 0);
    assert(oc != 0);

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) return; // Would be weird, but okay.
    pSpec_DEP = new int[nspecs];
    std::fill_n(pSpec_DEP, nspecs, DEP_NONE);

    pName = oc->getID();
    pChanState = oc->getChanState()->getID();
    pG = oc->getG();
    assert (pG >= 0.0);
    pERev = oc->getERev();
}

////////////////////////////////////////////////////////////////////////////////

ssolver::OhmicCurrdef::~OhmicCurrdef(void)
{
    if (pStatedef->countSpecs() > 0) delete[] pSpec_DEP;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::OhmicCurrdef::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&pG, sizeof(double));
    cp_file.write((char*)&pERev, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::OhmicCurrdef::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&pG, sizeof(double));
    cp_file.read((char*)&pERev, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::OhmicCurrdef::setup(void)
{
    assert(pSetupdone == false);

    uint chidx = pStatedef->getSpecIdx(pChanState);

    pSpec_CHANSTATE = chidx;
    pSpec_DEP[chidx] |= DEP_STOICH;

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::OhmicCurrdef::chanstate(void) const
{
    assert(pSetupdone == true);
    return pSpec_CHANSTATE;
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::OhmicCurrdef::dep(uint gidx) const
{
    assert(pSetupdone == true);
    assert(gidx < pStatedef->countSpecs());
    return pSpec_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::OhmicCurrdef::req(uint gidx) const
{
    assert(pSetupdone == true);
    assert(gidx < pStatedef->countSpecs());
    if (pSpec_DEP[gidx] != DEP_NONE) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

// END
