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


// STL headers.
#include <cassert>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/spec.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/specdef.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/types.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace ssolver = steps::solver;
namespace smod = steps::model;

////////////////////////////////////////////////////////////////////////////////

ssolver::Diffdef::Diffdef(Statedef * sd, uint idx, steps::model::Diff * d)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pDcst()
, pLig()
, pSetupdone(false)
, pSpec_DEP(nullptr)
{
    AssertLog(pStatedef != nullptr);
    AssertLog(d != nullptr);

    pName = d->getID();
    pDcst = d->getDcst();
    pLig = d->getLig()->getID();
    ligGIdx = pStatedef->getSpecIdx(pLig);

    const uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) { return;
    }
    pSpec_DEP = new int[nspecs];
    std::fill_n(pSpec_DEP, nspecs, DEP_NONE);

}

////////////////////////////////////////////////////////////////////////////////

ssolver::Diffdef::~Diffdef()
{
    if (pStatedef->countSpecs() > 0) delete[] pSpec_DEP;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Diffdef::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&pDcst), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Diffdef::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&pDcst), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Diffdef::setup()
{
    AssertLog(pSetupdone == false);

    pSpec_DEP[lig()] = DEP_STOICH;

    pSetupdone = true;

}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

void ssolver::Diffdef::setDcst(double d)
{
    AssertLog(d >= 0.0);
    pDcst = d;
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Diffdef::lig() const
{
    AssertLog(pStatedef != nullptr);
    return ligGIdx;
}

////////////////////////////////////////////////////////////////////////////////
/*
void ssolver::Diffdef::setLig(uint gidx)
{
    AssertLog(gidx < pStatedef->countSpecs());
    ssolver::Specdef * spec = pStatedef->specdef(gidx);
    pLig = spec->name();
}
*/
////////////////////////////////////////////////////////////////////////////////

int ssolver::Diffdef::dep(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::Diffdef::reqspec(uint gidx) const
{
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_DEP[gidx] != DEP_NONE;
}

////////////////////////////////////////////////////////////////////////////////

// END
