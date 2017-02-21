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

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/types.hpp"
#include "steps/error.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/specdef.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/model/spec.hpp"

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
, pSpec_DEP(0)
{
    assert(pStatedef != 0);
    assert(d != 0);

    pName = d->getID();
    pDcst = d->getDcst();
    pLig = d->getLig()->getID();

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) return;
    pSpec_DEP = new int[nspecs];
    std::fill_n(pSpec_DEP, nspecs, DEP_NONE);

}

////////////////////////////////////////////////////////////////////////////////

ssolver::Diffdef::~Diffdef(void)
{
    if (pStatedef->countSpecs() > 0) delete[] pSpec_DEP;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Diffdef::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&pDcst, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Diffdef::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&pDcst, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Diffdef::setup(void)
{
    assert (pSetupdone == false);

    pSpec_DEP[lig()] = DEP_STOICH;

    pSetupdone = true;

}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::Diffdef::name(void) const
{
    return pName;
}

////////////////////////////////////////////////////////////////////////////////

double ssolver::Diffdef::dcst(void) const
{
    return pDcst;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Diffdef::setDcst(double d)
{
    assert (d >= 0.0);
    pDcst = d;
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Diffdef::lig(void) const
{
    assert (pStatedef != 0);
    return pStatedef->getSpecIdx(pLig);
}

////////////////////////////////////////////////////////////////////////////////
/*
void ssolver::Diffdef::setLig(uint gidx)
{
    assert (gidx < pStatedef->countSpecs());
    ssolver::Specdef * spec = pStatedef->specdef(gidx);
    pLig = spec->name();
}
*/
////////////////////////////////////////////////////////////////////////////////

int ssolver::Diffdef::dep(uint gidx) const
{
    assert(pSetupdone == true);
    assert(gidx < pStatedef->countSpecs());
    return pSpec_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::Diffdef::reqspec(uint gidx) const
{
    assert(pSetupdone == true);
    assert(gidx < pStatedef->countSpecs());
    if (pSpec_DEP[gidx] != DEP_NONE) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

// END
