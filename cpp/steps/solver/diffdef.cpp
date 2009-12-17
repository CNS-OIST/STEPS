////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>
#include <cassert>

// STEPS headers.
#include <steps/common.h>
#include <steps/solver/types.hpp>
#include <steps/error.hpp>
#include <steps/solver/statedef.hpp>
#include <steps/solver/specdef.hpp>
#include <steps/solver/diffdef.hpp>
#include <steps/model/spec.hpp>

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::solver, ssolver);
NAMESPACE_ALIAS(steps::model, smod);

////////////////////////////////////////////////////////////////////////////////

ssolver::Diffdef::Diffdef(Statedef * sd, uint idx, steps::model::Diff * d)
: pStatedef(sd)
, pIdx(idx)
, pDiff(d)
, pSetupdone(false)
, pSpec_DEP(0)
{
    assert(pStatedef != 0);
    assert(pDiff != 0);

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) return;
    pSpec_DEP = new int[nspecs];
    std::fill_n(pSpec_DEP, nspecs, DEP_NONE);

}

////////////////////////////////////////////////////////////////////////////////

ssolver::Diffdef::~Diffdef(void)
{
	delete[] pSpec_DEP;
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
	assert (pDiff != 0);
	return pDiff->getID();
}

////////////////////////////////////////////////////////////////////////////////

double ssolver::Diffdef::dcst(void) const
{
	assert (pDiff != 0);
	return pDiff->getDcst();
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Diffdef::setDcst(double d)
{
	assert (d >= 0.0);
	pDiff->setDcst(d);
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Diffdef::lig(void) const
{
	assert (pDiff != 0);
	assert (pStatedef != 0);
	smod::Spec * spec =  pDiff->getLig();
	return pStatedef->getSpecIdx(spec);
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Diffdef::setLig(uint gidx)
{
	assert (pDiff != 0);
	assert (gidx < pStatedef->countSpecs());
	ssolver::Specdef * spec = pStatedef->specdef(gidx);
	pDiff->setLig(spec->spec());
}

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
