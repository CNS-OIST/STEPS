////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
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
#include <steps/solver/compdef.hpp>
#include <steps/solver/reacdef.hpp>
#include <steps/geom/comp.hpp>
#include <steps/model/spec.hpp>

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::solver, ssolver);
NAMESPACE_ALIAS(steps::model, smod);

////////////////////////////////////////////////////////////////////////////////

ssolver::Reacdef::Reacdef(Statedef * sd, uint idx, steps::model::Reac * r)
: pStatedef(sd)
, pIdx(idx)
, pReac(r)
, pSetupdone(false)
, pSpec_DEP(0)
, pSpec_LHS(0)
, pSpec_RHS(0)
, pSpec_UPD(0)
, pSpec_UPD_Coll()
{
    assert(pStatedef != 0);
    assert(pReac != 0);

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) return; // Would be weird, but okay.
    pSpec_DEP = new int[nspecs];
    std::fill_n(pSpec_DEP, nspecs, DEP_NONE);
    pSpec_LHS = new uint[nspecs];
    std::fill_n(pSpec_LHS, nspecs, 0);
    pSpec_RHS = new uint[nspecs];
    std::fill_n(pSpec_RHS, nspecs, 0);
    pSpec_UPD = new int[nspecs];
    std::fill_n(pSpec_UPD, nspecs, 0);

}

////////////////////////////////////////////////////////////////////////////////

ssolver::Reacdef::~Reacdef(void)
{
	delete[] pSpec_DEP;
	delete[] pSpec_LHS;
	delete[] pSpec_RHS;
	delete[] pSpec_UPD;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Reacdef::setup(void)
{
	assert (pSetupdone == false);
	assert (pReac != 0);

	// first copy the information about the reaction stoichiometry from Reac object
	smod::SpecPVec lhs = pReac->getLHS();
	smod::SpecPVecCI l_end = lhs.end();
	for (smod::SpecPVecCI l = lhs.begin(); l != l_end; ++l)
	{
		uint sidx = pStatedef->getSpecIdx(*l);
		pSpec_LHS[sidx] += 1;
	}
	smod::SpecPVec rhs = pReac->getRHS();
	smod::SpecPVecCI r_end = rhs.end();
	for (smod::SpecPVecCI r = rhs.begin(); r != r_end; ++r)
	{
		uint sidx = pStatedef->getSpecIdx(*r);
		pSpec_RHS[sidx] += 1;
	}

	// Now set up the update vector
	uint nspecs = pStatedef->countSpecs();
	for (uint i = 0; i < nspecs; ++i)
	{
	    int lhs = static_cast<int>(pSpec_LHS[i]);
	    int rhs = static_cast<int>(pSpec_RHS[i]);
	    int aux = pSpec_UPD[i] = (rhs - lhs);
	    if (lhs != 0) pSpec_DEP[i] |= DEP_STOICH;
	    if (aux != 0) pSpec_UPD_Coll.push_back(i);
	}

	pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::Reacdef::name(void) const
{
	assert (pReac != 0);
	return pReac->getID();
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Reacdef::order(void) const
{
	assert (pReac != 0);
	return pReac->getOrder();
}

////////////////////////////////////////////////////////////////////////////////

double ssolver::Reacdef::kcst(void) const
{
	assert (pReac != 0);
	return pReac->getKcst();
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Reacdef::lhs(uint gidx) const
{
    assert(gidx < pStatedef->countSpecs());
    return pSpec_LHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Reacdef::dep(uint gidx) const
{
	assert(pSetupdone == true);
    assert(gidx < pStatedef->countSpecs());
    return pSpec_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ssolver::Reacdef::rhs(uint gidx) const
{
    assert(gidx < pStatedef->countSpecs());
    return pSpec_RHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int ssolver::Reacdef::upd(uint gidx) const
{
	assert(pSetupdone == true);
    assert(gidx < pStatedef->countSpecs());
    return pSpec_UPD[gidx];
}

////////////////////////////////////////////////////////////////////////////////

bool ssolver::Reacdef::reqspec(uint gidx) const
{
	assert(pSetupdone == true);
    assert(gidx < pStatedef->countSpecs());
    if (pSpec_DEP[gidx] != DEP_NONE) return true;
    if (pSpec_RHS[gidx] != 0) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

// END
