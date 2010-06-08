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


/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

// Standard library & STL headers.
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <sstream>
#include <algorithm>

// STEPS headers.
#include "../common.h"
#include "tetexact.hpp"
#include "kproc.hpp"
#include "tet.hpp"
#include "tri.hpp"
#include "reac.hpp"
#include "sreac.hpp"
#include "diff.hpp"
#include "comp.hpp"
#include "patch.hpp"
#include "../math/constants.hpp"
#include "../error.hpp"
#include "../solver/statedef.hpp"
#include "../solver/compdef.hpp"
#include "../solver/patchdef.hpp"
#include "../solver/reacdef.hpp"
#include "../solver/sreacdef.hpp"
#include "../solver/diffdef.hpp"
#include "../solver/types.hpp"
#include "../geom/tetmesh.hpp"
#include "../geom/tet.hpp"
#include "../geom/tri.hpp"

////////////////////////////////////////////////////////////////////////////////

#define SCHEDULEWIDTH 16
//#define MAXLEVELS 10

////////////////////////////////////////////////////////////////////////////////


NAMESPACE_ALIAS(steps::tetexact, stex);
NAMESPACE_ALIAS(steps::solver, ssolver);
NAMESPACE_ALIAS(steps::math, smath);

////////////////////////////////////////////////////////////////////////////////

/// Unary function that calls the array delete[] operator on pointers. Easy
/// to use with STL/Boost (see steps::tools::DeletePointer).
///
struct DeleteArray
{
    template <typename Type> void operator() (Type * pointer) const
    {
        delete[] pointer;
    }
};

////////////////////////////////////////////////////////////////////////////////

void stex::schedIDXSet_To_Vec(stex::SchedIDXSet const & s, stex::SchedIDXVec & v)
{
    v.resize(s.size());
    std::copy(s.begin(), s.end(), v.begin());
}

////////////////////////////////////////////////////////////////////////////////

stex::Tetexact::Tetexact(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r)
: API(m, g, r)
, pMesh(0)
, pKProcs()
, pComps()
, pCompMap()
, pPatches()
, pTets()
, pTris()
, pA0(0.0)
, pLevelSizes()
, pLevels()
, pBuilt(false)
, pIndices(0)
, pMaxUpSize(0)
, pRannum(0)
{
	// Perform upcast.
	pMesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom());
	// TODO: deal with situation that it's not a mesh in a proper way -> throw exception

	// First initialise the pTets, pTris vector, because
	// want tets and tris to maintain indexing from Geometry
	uint ntets = mesh()->countTets();
	uint ntris = mesh()->countTris();
	pTets = std::vector<steps::tetexact::Tet *>(ntets, NULL);
	pTris = std::vector<steps::tetexact::Tri *>(ntris, NULL);

	/*
	// Now check that every Tetrahedron has been added to a compartment.
	// This is a requirement for the forseeable future
	uint ntmcomps = mesh()->_countComps();
	uint tetcount = 0;

	for (uint c = 0; c < ntmcomps; ++c)
	{
		steps::wm::Comp * wmcomp = mesh()->_getComp(c);
		steps::tetmesh::TmComp * tmcomp = dynamic_cast<steps::tetmesh::TmComp*>(wmcomp);
		tetcount += tmcomp->countTets();
	}
    if (tetcount != ntets)
    {
    	std::ostringstream os;
    	os << "Cannot create solver object with this ";
    	os << "steps.geom.Tetmesh geometry description object";
    	os << " All tetrahedrons in Tetmesh must be assigned to a";
    	os << " steps.geom.TmComp compartment object.";
    	throw steps::ArgErr(os.str());
    }
	*/

	// Now create the actual compartments.
	ssolver::CompDefPVecCI c_end = statedef()->endComp();
    for (ssolver::CompDefPVecCI c = statedef()->bgnComp(); c != c_end; ++c)
    {
        uint compdef_gidx = (*c)->gidx();
        uint comp_idx = _addComp(*c);
        assert(compdef_gidx == comp_idx);
    }
    // Create the actual patches.
    ssolver::PatchDefPVecCI p_end = statedef()->endPatch();
    for (ssolver::PatchDefPVecCI p = statedef()->bgnPatch(); p != p_end; ++p)
    {
        uint patchdef_gidx = (*p)->gidx();
        uint patch_idx = _addPatch(*p);
        assert(patchdef_gidx == patch_idx);
    }

    uint ncomps = pComps.size();
    assert (mesh()->_countComps() == ncomps);
    for (uint c = 0; c < ncomps; ++c)
    {
        // Now add the tets for this comp
    	// We have checked the indexing- c is the global index
        steps::wm::Comp * wmcomp = mesh()->_getComp(c);
        // Perform upcast
        steps::tetmesh::TmComp * tmcomp = dynamic_cast<steps::tetmesh::TmComp*>(wmcomp);
        steps::tetexact::Comp * localcomp = pComps[c];
        std::vector<uint> const & tetindcs = tmcomp->_getAllTetIndices();
        std::vector<uint>::const_iterator t_end = tetindcs.end();
        for (std::vector<uint>::const_iterator t = tetindcs.begin();
             t != t_end; ++t)
        {
        	steps::tetmesh::Tet * tet = new steps::tetmesh::Tet(mesh(), (*t));
        	assert (tet->getComp() == tmcomp);
        	double vol = tet->getVol();
        	double a0 = tet->getTri0Area();
        	double a1 = tet->getTri1Area();
        	double a2 = tet->getTri2Area();
        	double a3 = tet->getTri3Area();
        	double d0 = tet->getTet0Dist();
        	double d1 = tet->getTet1Dist();
        	double d2 = tet->getTet2Dist();
        	double d3 = tet->getTet3Dist();
        	// At this point fetch the indices of neighbouring tets too
        	int tet0 = tet->getTet0Idx();
        	int tet1 = tet->getTet1Idx();
        	int tet2 = tet->getTet2Idx();
        	int tet3 = tet->getTet3Idx();

        	_addTet((*t), localcomp, vol, a0, a1, a2, a3, d0, d1, d2, d3,
        			tet0, tet1, tet2, tet3);
        	delete tet;
        }
    }

    uint npatches = pPatches.size();
    assert (mesh()->_countPatches() == npatches);
    for (uint p = 0; p < npatches; ++p)
    {
    	// Now add the tris for this patch
    	// We have checked the indexing- p is the global index
    	steps::wm::Patch * wmpatch = mesh()->_getPatch(p);
        // Perform upcast
        steps::tetmesh::TmPatch * tmpatch = dynamic_cast<steps::tetmesh::TmPatch*>(wmpatch);
    	steps::tetexact::Patch * localpatch = pPatches[p];
    	std::vector<uint> const & triindcs = tmpatch->_getAllTriIndices();
    	std::vector<uint>::const_iterator t_end = triindcs.end();
    	for (std::vector<uint>::const_iterator t = triindcs.begin();
			 t != t_end; ++t)
    	{
    		steps::tetmesh::Tri * tri = new steps::tetmesh::Tri(mesh(), (*t));
    		assert (tri->getPatch() == tmpatch);
    		double area = tri->getArea();
    		//// TO DO: For 2D diffusion find length and distance information here
    		// e.g.
    		// double l0 = tri->getBar0Length();
    		// ..
    		// double d2 = tri->getTri2Dist();
    		int tetinner = tri->getTet0Idx();
    		int tetouter = tri->getTet1Idx();

    		_addTri((*t), localpatch, area, tetinner, tetouter);
    		delete tri;
    	}
    }

    // All tets and tris that belong to some comp or patch have been created
    // locally- now we can connect them locally
    // NOTE: currently if a tetrahedron's neighbour belongs to a different
    // comp they do not talk to each other (see stex::Tet::setNextTet())
    //

    assert (ntets == pTets.size());
    // pTets member size of all tets in geometry, but may not be filled with
    // local tets if they have not been added to a compartment
    /*
	 std::vector<steps::tetexact::Tet *>::const_iterator t_end = pTets.end();
	 for (std::vector<steps::tetexact::Tet *>::const_iterator t = pTets.begin();
	 t < nlocaltets; ++t)
	 */
    for (uint t = 0; t < ntets; ++t)
    {
    	if (pTets[t] == 0) continue;
    	int tet0 = pTets[t]->tet(0);
    	int tet1 = pTets[t]->tet(1);
    	int tet2 = pTets[t]->tet(2);
    	int tet3 = pTets[t]->tet(3);
    	// DEBUG 19/03/09 : again, incorrectly didn't allow for tet index == 0
    	if (tet0 >= 0)
    	{
			if (pTets[tet0] != 0) pTets[t]->setNextTet(0, pTets[tet0]);
    	}
    	if (tet1 >= 0)
    	{
			if (pTets[tet1] != 0) pTets[t]->setNextTet(1, pTets[tet1]);
    	}
    	if (tet2 >= 0)
    	{
    		if (pTets[tet2] != 0) pTets[t]->setNextTet(2, pTets[tet2]);
    	}
    	if (tet3 >= 0)
    	{
    		if (pTets[tet3] != 0) pTets[t]->setNextTet(3, pTets[tet3]);
    	}
    	// Not setting Tet triangles at this point- only want to set
    	// for surface triangles
    }
    assert (ntris == pTris.size());

    for (uint t = 0; t < ntris; ++t)
    {
    	// Looping over all possible tris, but only some have been added to a patch
    	if (pTris[t] == 0) continue;

    	// By convention, triangles in a patch should have an inner tetrahedron defined
    	// (neighbouring tets 'flipped' if necessary in Tetmesh)
    	// but not necessarily an outer tet
    	int tetinner = pTris[t]->tet(0);
    	int tetouter = pTris[t]->tet(1);

    	// DEBUG 18/03/09:
    	// Now correct check, previously didn't allow for tet index == 0
    	assert(tetinner >= 0);
    	assert(pTets[tetinner] != 0 );

		pTris[t]->setInnerTet(pTets[tetinner]);
    	// Now add this triangle to inner tet's list of neighbours
    	for (uint i=0; i <= 4; ++i)
    	{
    		// include assert for debugging purposes and remove
    		// once this is tested
    		assert (i < 4);														//////////
    		// check if there is already a neighbouring tet or tri
    		// In theory if there is a tri to add, the tet should
    		// have less than 4 neighbouring tets added because
    		// a neighbouring tet(s) is in a different compartment
    		// Also check tris because in some cases a surface tet
    		// may have more than 1 neighbouring tri
    		// NOTE: The order here will end up being different to the
    		// neighbour order at the Tetmesh level
    		if (pTets[tetinner]->nextTet(i) != 0) continue;
    		if (pTets[tetinner]->nextTri(i) != 0) continue;
    		pTets[tetinner]->setNextTri(i, pTris[t]);
    		break;
    	}



    	// DEBUG 18/03/09:
    	// Now correct check, previously didn't allow for tet index == 0
    	if (tetouter >= 0)
    	{
    		if (pTets[tetouter] != 0)
    		{
    			pTris[t]->setOuterTet(pTets[tetouter]);
    			// Add this triangle to outer tet's list of neighbours
    			for (uint i=0; i <= 4; ++i)
    			{
    				assert (i < 4);
    				if (pTets[tetouter]->nextTet(i) != 0) continue;
    				if (pTets[tetouter]->nextTri(i) != 0) continue;
    				pTets[tetouter]->setNextTri(i, pTris[t]);
    				break;
				}
    		}
    	}
    }

    _setup();
}

////////////////////////////////////////////////////////////////////////////////

stex::Tetexact::~Tetexact(void)
{
    CompPVecCI comp_e = pComps.end();
    for (CompPVecCI c = pComps.begin(); c != comp_e; ++c) delete *c;
    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI p = pPatches.begin(); p != patch_e; ++p) delete *p;

    TetPVecCI tet_e = pTets.end();
    for (TetPVecCI t = pTets.begin(); t != tet_e; ++t) delete *t;
    TriPVecCI tri_e = pTris.end();
    for (TriPVecCI t = pTris.end(); t != tri_e; ++t) delete *t;

    std::for_each(pLevels.begin(), pLevels.end(), DeleteArray());

    delete[] pIndices;
    delete[] pRannum;
}

////////////////////////////////////////////////////////////////////////////////

std::string stex::Tetexact::getSolverName(void) const
{
	return "tetexact";
}

////////////////////////////////////////////////////////////////////////////////

std::string stex::Tetexact::getSolverDesc(void) const
{
	return "SSA Direct Method in tetrahedral mesh";
}

////////////////////////////////////////////////////////////////////////////////

std::string stex::Tetexact::getSolverAuthors(void) const
{
	return "Stefan Wils and Iain Hepburn";
}

////////////////////////////////////////////////////////////////////////////////

std::string stex::Tetexact::getSolverEmail(void) const
{
	return "stefan@tnb.ua.ac.be, ihepburn@oist.jp";
}


////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setup(void)
{
	TetPVecCI tet_end = pTets.end();
	for (TetPVecCI t = pTets.begin(); t != tet_end; ++t)
	{
		// DEBUG: vector holds all possible tetrahedrons,
		// but they have not necessarily been added to a compartment.
		if ((*t) == 0) continue;

		(*t)->setupKProcs(this);
	}

	TriPVecCI tri_end = pTris.end();
	for (TriPVecCI t = pTris.begin(); t != tri_end; ++t)
	{
		// DEBUG: vector holds all possible triangles, but
		// only patch triangles are filled
		if ((*t) == 0) continue;

		(*t)->setupKProcs(this);
	}
	// Resolve all dependencies
	for (TetPVecCI t = pTets.begin(); t != tet_end; ++t)
	{
		// DEBUG: vector holds all possible tetrahedrons,
		// but they have not necessarily been added to a compartment.
		if ((*t) == 0) continue;

		KProcPVecCI kprocend = (*t)->kprocEnd();
		for (KProcPVecCI k = (*t)->kprocBegin(); k != kprocend; ++k)
		{
		    (*k)->setupDeps();
		}
	}
	for (TriPVecCI t = pTris.begin(); t != tri_end; ++t)
	{
		// DEBUG: vector holds all possible triangles, but
		// only patch triangles are filled
		if ((*t) == 0) continue;

	    KProcPVecCI kprocend = (*t)->kprocEnd();
	    for (KProcPVecCI k = (*t)->kprocBegin(); k != kprocend; ++k)
	    {
	        (*k)->setupDeps();
	    }
	}

	_build();
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_addComp(steps::solver::Compdef * cdef)
{
    stex::Comp * comp = new Comp(cdef);
    assert(comp != 0);
    uint compidx = pComps.size();
    pComps.push_back(comp);
    pCompMap[cdef] = comp;
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_addPatch(steps::solver::Patchdef * pdef)
{
    /* Comp * icomp = 0;
	 Comp * ocomp = 0;
	 if (pdef->icompdef()) icomp = pCompMap[pdef->icompdef()];
	 if (pdef->ocompdef()) ocomp = pCompMap[pdef->ocompdef()];
	 */
    stex::Patch * patch = new Patch(pdef);
    assert(patch != 0);
    uint patchidx = pPatches.size();
    pPatches.push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_addTet(uint tetidx,
							 steps::tetexact::Comp * comp, double vol,
							 double a1, double a2, double a3, double a4,
							 double d1, double d2, double d3, double d4,
							 int tet0, int tet1, int tet2, int tet3)
{
	steps::solver::Compdef * compdef  = comp->def();
    stex::Tet * localtet = new stex::Tet(compdef, vol, a1, a2, a3, a4, d1, d2, d3, d4,
									     tet0, tet1, tet2, tet3);
    assert(localtet != 0);
    assert(tetidx < pTets.size());
    assert(pTets[tetidx] == 0);
    pTets[tetidx] = localtet;
    comp->addTet(localtet);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_addTri(uint triidx, steps::tetexact::Patch * patch, double area,
							 int tinner, int touter)
{
    steps::solver::Patchdef * patchdef = patch->def();
    stex::Tri * tri = new stex::Tri(patchdef, area, tinner, touter);
    assert(tri != 0);
    assert (triidx < pTris.size());
    assert (pTris[triidx] == 0);
    pTris[triidx] = tri;
    patch->addTri(tri);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::reset(void)
{
	std::for_each(pComps.begin(), pComps.end(), std::mem_fun(&Comp::reset));
	std::for_each(pPatches.begin(), pPatches.end(), std::mem_fun(&Patch::reset));

    // std::for_each(pTets.begin(), pTets.end(), std::mem_fun(&Tet::reset));

    TetPVecCI tet_end = pTets.end();
    for (TetPVecCI t = pTets.begin(); t != tet_end; ++t)
    {
    	if ((*t) == 0) continue;
    	(*t)->reset();
    }

	TriPVecCI tri_end = pTris.end();
	for (TriPVecCI t = pTris.begin(); t != tri_end; ++t)
	{
		if ((*t) == 0) continue;
		(*t)->reset();
	}

	statedef()->resetTime();
	statedef()->resetNSteps();

	_reset();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::run(double endtime)
{
	if (endtime < statedef()->time())
	{
		std::ostringstream os;
		os << "Endtime is before current simulation time";
	    throw steps::ArgErr(os.str());
	}
	while (statedef()->time() < endtime)
	{
		stex::KProc * kp = _getNext();
		if (kp == 0) break;
		double a0 = getA0();
		if (a0 == 0.0) break;
		double dt = rng()->getExp(a0);
		if ((statedef()->time() + dt) > endtime) break;
		_executeStep(kp, dt);
	}
	statedef()->setTime(endtime);
}

////////////////////////////////////////////////////////////////////////

void stex::Tetexact::advance(double adv)
{
	if (adv < 0.0)
	{
		std::ostringstream os;
		os << "Time to advance cannot be negative";
	    throw steps::ArgErr(os.str());
	}

	double endtime = statedef()->time() + adv;
	run(endtime);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::step(void)
{
	stex::KProc * kp = _getNext();
	if (kp == 0) return;
	double a0 = getA0();
	if (a0 == 0.0) return;
	double dt = rng()->getExp(a0);
	_executeStep(kp, dt);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::getTime(void) const
{
	return statedef()->time();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setTime(double time)
{
	statedef()->setTime(time);
}

////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::getNSteps(void) const
{
    return statedef()->nsteps();
}

////////////////////////////////////////////////////////////////////////

void stex::Tetexact::setNSteps(uint nsteps)
{
    statedef()->setNSteps(nsteps);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompVol(uint cidx) const
{
	assert(cidx < statedef()->countComps());
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	return comp->vol();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompCount(uint cidx, uint sidx) const
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint slidx = comp->def()->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return 0.0;

	uint count = 0;
	TetPVecCI t_end = comp->endTet();
	for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		count += (*t)->pools()[slidx];
	}

	return count;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompCount(uint cidx, uint sidx, double n)
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	assert (n >= 0.0);
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	uint slidx = comp->def()->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return;

	double totalvol = comp->def()->vol();
	TetPVecCI t_end = comp->endTet();

	double n_int = std::floor(n);
	double n_frc = n - n_int;
	uint c = static_cast<uint>(n_int);
	if (n_frc > 0.0)
	{
		double rand01 = rng()->getUnfIE();
		if (rand01 < n_frc) c++;
	}
	if (c >= comp->countTets())
	{
		uint nremoved = 0;
		for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
		{
			Tet * tet = *t;
			double fract = static_cast<double>(c) * (tet->vol() / totalvol);
			uint n3 = static_cast<uint>(std::floor(fract));
			tet->setCount(slidx, n3);
            // BUGFIX 18/11/09 IH. By reducing c here we were not giving all
            // tets an equal share. Tets with low index would have a
            // greater share than those with higher index.
			//c -= n3;
			nremoved += n3;
		}
		assert(nremoved <= c);
		c -= nremoved;
	}

	while (c != 0)
	{
		Tet * tet = comp->pickTetByVol(rng()->getUnfIE());
        assert (tet != 0);
        tet->setCount(slidx, (tet->pools()[slidx] + 1.0));
        c--;
	}

	for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		_updateSpec(*t, slidx);
	}

	// Rates have changed
	_reset();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompAmount(uint cidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getCompCount(cidx, sidx);
	return (count / smath::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompAmount(uint cidx, uint sidx, double a)
{
	// convert amount in mols to number of molecules
	double a2 = a * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setCompCount(cidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompConc(uint cidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getCompCount(cidx, sidx);
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	double vol = comp->vol();
	return count/ (1.0e3 * vol * steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompConc(uint cidx, uint sidx, double c)
{
	assert(c >= 0.0);
	assert (cidx < statedef()->countComps());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	double count = c * (1.0e3 * comp->vol() * steps::math::AVOGADRO);
	// the following method does all the necessary argument checking
	_setCompCount(cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getCompClamped(uint cidx, uint sidx) const
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	assert(statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lsidx = comp->def()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return false;

    TetPVecCI t_end = comp->endTet();
    for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
    	if ((*t)->clamped(lsidx) == false) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompClamped(uint cidx, uint sidx, bool b)
{
	assert(cidx < statedef()->countComps());
	assert(sidx < statedef()->countSpecs());
	assert(statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lsidx = comp->def()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED) return;

	// Set the flag in def object, though this may not be necessary
	comp->def()->setClamped(lsidx, b);

	TetPVecCI t_end = comp->endTet();
	for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		(*t)->setClamped(lsidx, b);
	}
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacK(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	assert(statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lridx = comp->def()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// We're just returning the default value for this comp, individual
	// tets may have different Kcsts set individually
	return (comp->def()->kcst(lridx));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompReacK(uint cidx, uint ridx, double kf)
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	assert(statedef()->countComps() == pComps.size());
	assert(kf >= 0.0);
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lridx = comp->def()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	// First set the default value for the comp
	comp->def()->setKcst(lridx, kf);

	// Now update all tetrahedra in this comp
	TetPVecCI t_end = comp->endTet();
	for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		(*t)->reac(lridx)->setKcst(kf);
	}

	// Rates have changed
	_reset();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getCompReacActive(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	assert(statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lridx = comp->def()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return false;

	TetPVecCI t_end = comp->endTet();
	for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		if ((*t)->reac(lridx)->inactive() == true) return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompReacActive(uint cidx, uint ridx, bool a)
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	assert(statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert(comp != 0);
	uint lridx = comp->def()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	// Set the default value for the comp, though this is not entirely
	// necessary
	comp->def()->setActive(lridx, a);

	TetPVecCI t_end = comp->endTet();
	for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		(*t)->reac(lridx)->setActive(a);
	}
    // It's cheaper to just recompute everything.
    _reset();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompDiffD(uint cidx, uint didx) const
{
	assert (cidx < statedef()->countComps());
	assert (didx < statedef()->countDiffs());
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert (comp != 0);
	uint ldidx = comp->def()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED) return 0.0;

	// We're just returning the default value for this comp, individual
	// tets may have different Dcsts set individually
	return (comp->def()->dcst(ldidx));

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompDiffD(uint cidx, uint didx, double dk)
{
	assert (cidx < statedef()->countComps());
	assert (didx < statedef()->countDiffs());
	assert (statedef()->countComps() == pComps.size());
	assert(dk >= 0.0);
	stex::Comp * comp = _comp(cidx);
	assert (comp != 0);
	uint ldidx = comp->def()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED) return;

	// First set the default value for the comp
	comp->def()->setDcst(ldidx, dk);

	// Now update all tets in this comp
	TetPVecCI t_end = comp->endTet();
	for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		(*t)->diff(ldidx)->setDcst(dk);
	}

	// Rates have changed
	_reset();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getCompDiffActive(uint cidx, uint didx) const
{
	assert (cidx < statedef()->countComps());
	assert (didx < statedef()->countDiffs());
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert (comp != 0);
	uint ldidx = comp->def()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED) return false;

	TetPVecCI t_end = comp->endTet();
	for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		if ((*t)->diff(ldidx)->inactive() == true) return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setCompDiffActive(uint cidx, uint didx, bool act)
{
	assert (cidx < statedef()->countComps());
	assert (didx < statedef()->countDiffs());
	assert (statedef()->countComps() == pComps.size());
	stex::Comp * comp = _comp(cidx);
	assert (comp != 0);
	uint ldidx = comp->def()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED) return;

	TetPVecCI t_end = comp->endTet();
	for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
	{
		(*t)->diff(ldidx)->setActive(act);
	}
    // It's cheaper to just recompute everything.
    _reset();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchArea(uint pidx) const
{
	assert (pidx < statedef()->countPatches());
	assert (statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert (patch != 0);
	return patch->area();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchCount(uint pidx, uint sidx) const
{
	assert (pidx < statedef()->countPatches());
	assert (sidx < statedef()->countSpecs());
	assert (statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert (patch != 0);
	uint slidx = patch->def()->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return 0.0;

	uint count = 0;
	TriPVecCI t_end = patch->endTri();
	for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
	{
		count += (*t)->pools()[slidx];
	}

	return count;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchCount(uint pidx, uint sidx, double n)
{
	assert (pidx < statedef()->countPatches());
	assert (sidx < statedef()->countSpecs());
	assert (statedef()->countPatches() == pPatches.size());
	assert (n >= 0.0);
	stex::Patch * patch = _patch(pidx);
	assert (patch != 0);
	uint slidx = patch->def()->specG2L(sidx);
	if (slidx == ssolver::LIDX_UNDEFINED) return;

	double totalarea = patch->def()->area();
	TriPVecCI t_end = patch->endTri();

	double n_int = std::floor(n);
	double n_frc = n - n_int;
	uint c = static_cast<uint>(n_int);
	if (n_frc > 0.0)
	{
		double rand01 = rng()->getUnfIE();
		if (rand01 < n_frc) c++;
	}

	if (c >= patch->countTris())
	{
		uint nremoved = 0;
		for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
		{
			Tri * tri = *t;
			double fract = static_cast<double>(c) * (tri->area()/totalarea);
            uint n3 = static_cast<uint>(std::floor(fract));
            tri->setCount(slidx, n3);
            // BUGFIX 18/11/09 IH. By reducing c here we were not giving all
            // triangles an equal share. Triangles with low index would have a
            // greater share than those with higher index.
            // c -= n3;
            nremoved += n3;
		}
		assert(nremoved <= c);
		c -= nremoved;
	}

	while(c != 0)
	{
		Tri * tri = patch->pickTriByArea(rng()->getUnfIE());
		assert (tri != 0);
		tri->setCount(slidx, (tri->pools()[slidx] + 1.0));
		c--;
	}
	for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
	{
		_updateSpec(*t, slidx);
	}

	// Rates have changed
	_reset();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchAmount(uint pidx, uint sidx) const
{
	// the following method does all the necessary argument checking
	double count = _getPatchCount(pidx, sidx);
	return (count / steps::math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchAmount(uint pidx, uint sidx, double a)
{
	assert(a > 0.0);
	// convert amount in mols to number of molecules
	double a2 = a * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setPatchCount(pidx, sidx, a2);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getPatchClamped(uint pidx, uint sidx) const
{
	assert (pidx < statedef()->countPatches());
	assert (sidx < statedef()->countSpecs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsidx = patch->def()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return false;

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if ((*t)->clamped(lsidx) == false) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchClamped(uint pidx, uint sidx, bool buf)
{
	assert (pidx < statedef()->countPatches());
	assert (sidx < statedef()->countSpecs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsidx = patch->def()->specG2L(sidx);
    if (lsidx == ssolver::LIDX_UNDEFINED) return;

    // Set the flag in def object for consistency, though this is not
    // entirely necessary
    patch->def()->setClamped(lsidx, buf);

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        (*t)->setClamped(lsidx, buf);
    }

}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacK(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED) return 0.0;

    // We're just returning the default value for this patch, individual
    // triangles may have different Kcsts set
    return (patch->def()->kcst(lsridx));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchSReacK(uint pidx, uint ridx, double kf)
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	assert(statedef()->countPatches() == pPatches.size());
	assert(kf >= 0.0);
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED) return;

    // First set the default values for this patch
    patch->def()->setKcst(lsridx, kf);

    // Now update all triangles in this patch
    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        (*t)->sreac(lsridx)->setKcst(kf);
    }

	// Rates have changed
	_reset();
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getPatchSReacActive(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED) return false;

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        if ((*t)->sreac(lsridx)->inactive() == true) return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setPatchSReacActive(uint pidx, uint ridx, bool a)
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	assert(statedef()->countPatches() == pPatches.size());
	stex::Patch * patch = _patch(pidx);
	assert(patch != 0);
	uint lsridx = patch->def()->sreacG2L(ridx);
    if (lsridx == ssolver::LIDX_UNDEFINED) return;

    // First set the flags in def object for consistency, though this is
    // not entirely necessary for this solver
    patch->def()->setActive(lsridx, a);

    TriPVecCI t_end = patch->endTri();
    for (TriPVecCI t = patch->bgnTri(); t != t_end; ++t)
    {
        (*t)->sreac(lsridx)->setActive(a);
    }
    // It's cheaper to just recompute everything.
    _reset();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::addKProc(steps::tetexact::KProc * kp)
{
	assert (kp != 0);

	SchedIDX nidx = pKProcs.size();
	pKProcs.push_back(kp);
	kp->setSchedIDX(nidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_build(void)
{
	assert (pBuilt == false);

    // Setup level.
    uint clsize = pKProcs.size();
    if (clsize == 0) return;

    // Work up.
    uint clevel = 0;
    do
    {
        // Make sure the new size is a multiple of SCHEDULEWIDTH.
        uint extra = clsize % SCHEDULEWIDTH;
        if (extra != 0) clsize += SCHEDULEWIDTH - extra;

        // Create the level and add it.
        double * level = new double[clsize];
        std::fill_n(level, clsize, 0.0);
        pLevelSizes.push_back(clsize);
        pLevels.push_back(level);

        // Prepare for next level.
        clevel++;
        clsize = clsize / SCHEDULEWIDTH;
    }
    while (clsize > 1);

    // Set top level.
    pA0 = 0.0;

    // Time to create ONE indices table to hold the run's present reaction of
    // choice's update vector. This will be re-used and replace old-version's
    // hard-coded table in _update. Size is the maximum possible, found by looping
    // over all KProcs. This little bit of computational time is well worth all
    // that needless memory allocation
	uint maxupdvecsize = 0;
	TetPVecCI tet_end = pTets.end();
	for (TetPVecCI tet = pTets.begin(); tet != tet_end; ++tet)
	{
		// DEBUG: pTets element may be zero pointer if tet not added to a compartment
		if ((*tet) == 0) continue;

		// To find the maximum size of an individual kproc's upd vec size
		uint kpuv_size = 0;
		SchedIDXSet updset;
		KProcPVecCI kproc_end = (*tet)->kprocEnd();
		for (KProcPVecCI k = (*tet)->kprocBegin(); k != kproc_end; ++k)
		{
			updset.insert((*k)->schedIDX());
			if ((*k)->updVecSize() > kpuv_size) kpuv_size = (*k)->updVecSize();
		}

	    // Loop over neighbouring triangles.
		for (uint i = 0; i < 4; ++i)
		{
			Tri * tri = (*tet)->nextTri(i);
			if (tri == 0) continue;
			kproc_end = tri->kprocEnd();
			for (KProcPVecCI k = tri->kprocBegin(); k != kproc_end; ++k)
			{
				updset.insert((*k)->schedIDX());
			}
		}
		SchedIDXVec updvec;
		schedIDXSet_To_Vec(updset, updvec);
		if (maxupdvecsize < updvec.size()) maxupdvecsize = updvec.size();
		if (maxupdvecsize < kpuv_size) maxupdvecsize = kpuv_size;
	}
	TriPVecCI tri_end = pTris.end();
	for (TriPVecCI tri = pTris.begin(); tri != tri_end; ++tri)
	{
		uint kpuv_size = 0;
		//DEBUG: Didn't allow for a zero pointer (not added to geometry i.e. not
		// belonging to a patch
		if ((*tri) == 0) continue;

		SchedIDXSet updset;
	    KProcPVecCI kproc_end = (*tri)->kprocEnd();
		for (KProcPVecCI k = (*tri)->kprocBegin(); k != kproc_end; ++k)
		{
			updset.insert((*k)->schedIDX());
			if ((*k)->updVecSize() > kpuv_size) kpuv_size = (*k)->updVecSize();
		}
		SchedIDXVec updvec;
		schedIDXSet_To_Vec(updset, updvec);
		if (maxupdvecsize < updvec.size()) maxupdvecsize = updvec.size();
		if (maxupdvecsize < kpuv_size) maxupdvecsize = kpuv_size;
	}

	pMaxUpSize = maxupdvecsize;
	pIndices = new uint[pMaxUpSize];

    // Also let's create a random number holder-table,
    // size of number of KProcs % SCHEDULEWIDTH or pLevels.size()
    // This will be re-used in _getNext as apposed to hard-coded (again maximum
    // limit).
    uint lsize = pLevels.size();
    pRannum = new double[lsize];

    pBuilt = true;
}

////////////////////////////////////////////////////////////////////////////////

steps::tetexact::KProc * stex::Tetexact::_getNext(void) const
{
    assert(pA0 >= 0.0);
    // Quick check to see whether nothing is there.
    if (pA0 == 0.0) return 0;

    // Start at top level.
    uint clevel = pLevels.size();
    // And start at the first node of that level.
    uint cur_node = 0;

    // Prepare random numbers.
    for (uint i = 0; i < clevel; ++i)
    {
        pRannum[i] = rng()->getUnfIE();
    }

    // Run until top level.
    double a0 = pA0;
    double * level = 0;
    while (clevel != 0)
    {
        // Decrease the current level.
        clevel--;
        // and start looking in the right place.
        cur_node *= SCHEDULEWIDTH;
        uint max_node = cur_node + SCHEDULEWIDTH;

        // Fetch the level.
        level = pLevels[clevel];

        // Compute local selector.
        double selector = pRannum[clevel] * a0;

        // Compare.
        double accum = 0.0;
        // 27/10/09 I.H. 'old' removed from for loop because not used.
        // double old = 0.0;
        double curval = 0.0;
        for (uint i = 0; i < SCHEDULEWIDTH; ++i)
        {
            curval = level[cur_node];
            if (selector < curval + accum) break;
            accum += curval;
            // old = accum;
            cur_node++;
        }

        // Checks.
        assert(cur_node < max_node);
        assert(curval > 0.0);
        a0 = curval;
    }

    // Check.
    assert(cur_node < pKProcs.size());
    return pKProcs[cur_node];
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_reset(void)
{
    if (pKProcs.size() == 0) return;

    // Reset the basic level: compute rates.
    double * oldlevel = pLevels[0];
    uint cur_node = 0;
    std::vector<KProc*>::iterator kp_end = pKProcs.end();
    for (std::vector<KProc*>::iterator kp = pKProcs.begin(); kp != kp_end; ++kp)
    {
        oldlevel[cur_node++] = (*kp)->rate();
    }

    // Work up.
    for (uint cur_level = 1; cur_level < pLevels.size(); ++cur_level)
    {
        // Compute the number of nodes to reset on this level.
        uint numnodes = pLevelSizes[cur_level - 1] / SCHEDULEWIDTH;

        // Fetch a pointer to this level.
        double * level = pLevels[cur_level];

        // Recompute them.
        uint child_node = 0;
        for (cur_node = 0; cur_node < numnodes; ++cur_node)
        {
            double val = 0.0;
            for (uint i = 0; i < SCHEDULEWIDTH; ++i)
            {
                val += oldlevel[child_node++];
            }
            level[cur_node] = val;
        }

        // Copy the level.
        oldlevel = level;
    }

    // Compute zero propensity.
    pA0 = 0.0;
    for (uint i = 0; i < SCHEDULEWIDTH; ++i)
    {
        pA0 += oldlevel[i];
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_executeStep(steps::tetexact::KProc * kp, double dt)
{
	SchedIDXVec const & upd = kp->apply(rng());
	_update(upd);
	statedef()->incTime(dt);
	statedef()->incNSteps(1);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_update(SchedIDXVec const & entries)
{
    if (countKProcs() == 0) return;

    // Prefetch zero level.
    double * level0 = pLevels[0];
    // Number of entries.
    assert(entries.size() <= pMaxUpSize);								/////////

    // Recompute rates.
    SchedIDXVecCI sidx_end = entries.end();
    uint prev_e = 0xFFFFFFFF;
    uint cur_e = 0;
    for (SchedIDXVecCI sidx = entries.begin(); sidx != sidx_end; ++sidx)
    {
        // Fetch index.
        uint idx = *sidx;
        // Recompute rate, get difference, and store.
        double newrate = pKProcs[idx]->rate();
        level0[idx] = newrate;

        // Store and collapse if possible.
        idx /= SCHEDULEWIDTH;
        if (prev_e == 0xFFFFFFFF)
        {
            prev_e = 0;
            pIndices[cur_e++] = idx;
        }
        else if (pIndices[prev_e] != idx)
        {
            prev_e = cur_e;
            pIndices[cur_e++] = idx;
        }
    }
    uint nentries = cur_e;

    // Update upper levels.
    uint nlevels = pLevels.size();
    double * prevlevel = pLevels[0];
    for (uint l = 1; l < nlevels; ++l)
    {
        // Update the first entry.
        cur_e = 0;
        prev_e = 0xFFFFFFFF;

        // Fetch a pointer to the current level.
        double * currlevel = pLevels[l];

        // Recompute the entries.
        for (uint e = 0; e < nentries; ++e)
        {
            // Fetch index.
            uint idx = pIndices[e];

            // Recompute.
            double val = 0.0;
            uint idx2 = idx * SCHEDULEWIDTH;
            for (uint i = 0; i < SCHEDULEWIDTH; ++i)
            {
                val += prevlevel[idx2++];
            }
            currlevel[idx] = val;

            // Store and collapse if possible.
            idx /= SCHEDULEWIDTH;
            if (prev_e == 0xFFFFFFFF)
            {
                prev_e = 0;
                pIndices[cur_e++] = idx;
            }
            else if (pIndices[prev_e] != idx)
            {
                prev_e = cur_e;
                pIndices[cur_e++] = idx;
            }
        }

        // Update the pointer to the previous level.
        prevlevel = currlevel;

        // cur_e now is the new number of entries to handle.
        nentries = cur_e;
    }

    // Update zero propensity.
    double * toplevel = pLevels[pLevels.size() - 1];
    pA0 = 0.0;
    for (uint i = 0; i < SCHEDULEWIDTH; ++i)
    {
        pA0 += toplevel[i];
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_updateSpec(steps::tetexact::Tet * tet, uint spec_lidx)
{
    SchedIDXSet updset;

    // Loop over tet.
    KProcPVecCI kproc_end = tet->kprocEnd();
    for (KProcPVecCI k = tet->kprocBegin(); k != kproc_end; ++k)
    {
        updset.insert((*k)->schedIDX());
    }

    // Loop over neighbouring triangles.
    for (uint i = 0; i < 4; ++i)
    {
        Tri * tri = tet->nextTri(i);
        if (tri == 0) continue;
        kproc_end = tri->kprocEnd();
        for (KProcPVecCI k = tri->kprocBegin(); k != kproc_end; ++k)
        {
            updset.insert((*k)->schedIDX());
        }
    }

    // Send the list of kprocs that need to be updated to the schedule.
    if (updset.empty()) return;
    SchedIDXVec updvec;
    schedIDXSet_To_Vec(updset, updvec);

    _update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_updateSpec(steps::tetexact::Tri * tri, uint spec_lidx)
{
    SchedIDXSet updset;

    KProcPVecCI kproc_end = tri->kprocEnd();
    for (KProcPVecCI k = tri->kprocBegin(); k != kproc_end; ++k)
    {
        updset.insert((*k)->schedIDX());
    }

    // Send the list of kprocs that need to be updated to the schedule.
    if (updset.empty()) return;
    SchedIDXVec updvec;
    schedIDXSet_To_Vec(updset, updvec);

    _update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacH(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Comp object has same index as solver::Compdef object
	stex::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);

	TetPVecCI t_bgn = lcomp->bgnTet();
	TetPVecCI t_end = lcomp->endTet();
	if (t_bgn == t_end) return 0.0;

	double h = 0.0;
	for (TetPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::Reac * reac = (*t)->reac(lridx);
		h += reac->h();
	}

	return h;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacC(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Comp object has same index as solver::Compdef object
	stex::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);

	TetPVecCI t_bgn = lcomp->bgnTet();
	TetPVecCI t_end = lcomp->endTet();
	if (t_bgn == t_end) return 0.0;
	double c2 = 0.0;
	double v = 0.0;
	for (TetPVecCI t = t_bgn; t != t_end; ++t)
	{
		double v2 = (*t)->vol();
		stex::Reac * reac = (*t)->reac(lridx);
		c2 += (reac->c() * v2);
		v += v2;
	}
	assert(v > 0.0);
	return c2/v;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getCompReacA(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Comp object has same index as solver::Compdef object
	stex::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);

	TetPVecCI t_bgn = lcomp->bgnTet();
	TetPVecCI t_end = lcomp->endTet();
	if (t_bgn == t_end) return 0.0;

	double a = 0.0;
	for (TetPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::Reac * reac = (*t)->reac(lridx);
		a += reac->rate();
	}

	return a;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_getCompReacExtent(uint cidx, uint ridx) const
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Comp object has same index as solver::Compdef object
	stex::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);

	TetPVecCI t_bgn = lcomp->bgnTet();
	TetPVecCI t_end = lcomp->endTet();
	if (t_bgn == t_end) return 0.0;

	uint x = 0.0;
	for (TetPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::Reac * reac = (*t)->reac(lridx);
		x += reac->getExtent();
	}

	return x;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_resetCompReacExtent(uint cidx, uint ridx)
{
	assert(cidx < statedef()->countComps());
	assert(ridx < statedef()->countReacs());
	ssolver::Compdef * comp = statedef()->compdef(cidx);
	assert(comp != 0);
	uint lridx = comp->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	// The 'local' Comp object has same index as solver::Compdef object
	stex::Comp * lcomp = pComps[cidx];
	assert (lcomp->def() == comp);

	TetPVecCI t_bgn = lcomp->bgnTet();
	TetPVecCI t_end = lcomp->endTet();
	if (t_bgn == t_end) return;

	for (TetPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::Reac * reac = (*t)->reac(lridx);
		reac->resetExtent();
	}
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacH(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	uint lsridx = patch->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Patch object has same index as solver::Patchdef object
	stex::Patch * lpatch = pPatches[pidx];
	assert (lpatch->def() == patch);

	TriPVecCI t_bgn = lpatch->bgnTri();
	TriPVecCI t_end = lpatch->endTri();
	if (t_bgn == t_end) return 0.0;

	double h = 0.0;
	for (TriPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::SReac * sreac = (*t)->sreac(lsridx);
		h += sreac->h();
	}

	return h;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacC(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	uint lsridx = patch->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Patch object has same index as solver::Patchdef object
	stex::Patch * lpatch = pPatches[pidx];
	assert (lpatch->def() == patch);

	TriPVecCI t_bgn = lpatch->bgnTri();
	TriPVecCI t_end = lpatch->endTri();
	if (t_bgn == t_end) return 0.0;

	double c = 0.0;
	double a = 0.0;
	for (TriPVecCI t = t_bgn; t != t_end; ++t)
	{
		double a2 = (*t)->area();
		stex::SReac * sreac = (*t)->sreac(lsridx);
		c += (sreac->c() * a2);
		a += a2;
	}
	assert (a > 0.0);
	return c/a;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getPatchSReacA(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	uint lsridx = patch->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Patch object has same index as solver::Patchdef object
	stex::Patch * lpatch = pPatches[pidx];
	assert (lpatch->def() == patch);

	TriPVecCI t_bgn = lpatch->bgnTri();
	TriPVecCI t_end = lpatch->endTri();
	if (t_bgn == t_end) return 0.0;

	double a = 0.0;
	for (TriPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::SReac * sreac = (*t)->sreac(lsridx);
		a += sreac->rate();
	}
	return a;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Tetexact::_getPatchSReacExtent(uint pidx, uint ridx) const
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	uint lsridx = patch->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return 0.0;

	// The 'local' Patch object has same index as solver::Patchdef object
	stex::Patch * lpatch = pPatches[pidx];
	assert (lpatch->def() == patch);

	TriPVecCI t_bgn = lpatch->bgnTri();
	TriPVecCI t_end = lpatch->endTri();
	if (t_bgn == t_end) return 0.0;

	double x = 0.0;
	for (TriPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::SReac * sreac = (*t)->sreac(lsridx);
		x += sreac->getExtent();
	}
	return x;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_resetPatchSReacExtent(uint pidx, uint ridx)
{
	assert (pidx < statedef()->countPatches());
	assert (ridx < statedef()->countSReacs());
	ssolver::Patchdef * patch = statedef()->patchdef(pidx);
	uint lsridx = patch->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return;

	// The 'local' Patch object has same index as solver::Patchdef object
	stex::Patch * lpatch = pPatches[pidx];
	assert (lpatch->def() == patch);

	TriPVecCI t_bgn = lpatch->bgnTri();
	TriPVecCI t_end = lpatch->endTri();
	if (t_bgn == t_end) return;

	for (TriPVecCI t = t_bgn; t != t_end; ++t)
	{
		stex::SReac * sreac = (*t)->sreac(lsridx);
		sreac->resetExtent();
	}
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetVol(uint tidx) const
{
	assert (tidx < pTets.size());
	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}
	return pTets[tidx]->vol();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetVol(uint tidx, double vol)
{
    throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetCount(uint tidx, uint sidx) const
{
	assert (tidx < pTets.size());
	assert (sidx < statedef()->countSpecs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];
	uint lsidx = tet->compdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED) return 0.0;
	return tet->pools()[lsidx];
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetCount(uint tidx, uint sidx, double n)
{
	assert (tidx < pTets.size());
	assert (sidx < statedef()->countSpecs());
	assert (n >= 0.0);

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lsidx = tet->compdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED) return;

	double n_int = std::floor(n);
	double n_frc = n - n_int;
	uint c = static_cast<uint>(n_int);
	if (n_frc > 0.0)
	{
		double rand01 = rng()->getUnfIE();
		if (rand01 < n_frc) c++;
	}

	// Tet object updates def level Comp object counts
	tet->setCount(lsidx, c);
	_updateSpec(tet, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetAmount(uint tidx, uint sidx) const
{
	// following method does all necessary argument checking
	double count = _getTetCount(tidx, sidx);
	return count/steps::math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetAmount(uint tidx, uint sidx, double m)
{
	// convert amount in mols to number of molecules
	double m2 = m * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setTetCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetConc(uint tidx, uint sidx) const
{
	// following method does all necessary argument checking
	double count = _getTetCount(tidx, sidx);
	stex::Tet * tet = pTets[tidx];
	double vol = tet->vol();
	return (count/(1.0e3 * vol * steps::math::AVOGADRO));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetConc(uint tidx, uint sidx, double c)
{
	assert (c >= 0.0);
	assert (tidx < pTets.size());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];
	double count = c * (1.0e3 * tet->vol() * steps::math::AVOGADRO);
	// the following method does all the necessary argument checking
	_setTetCount(tidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetClamped(uint tidx, uint sidx) const
{
	assert (tidx < pTets.size());
	assert (sidx < statedef()->countSpecs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lsidx = tet->compdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED) return false;
	return tet->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetClamped(uint tidx, uint sidx, bool buf)
{
	assert (tidx < pTets.size());
	assert (sidx < statedef()->countSpecs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lsidx = tet->compdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED) return;

	tet->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacK(uint tidx, uint ridx) const
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return (tet->reac(lridx)->kcst());
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetReacK(uint tidx, uint ridx, double kf)
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());
	assert (kf >= 0.0);

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	tet->reac(lridx)->setKcst(kf);

	// Update the rates for all KProcs affected by this change
	SchedIDXVec updvec;
	updvec.push_back(tet->reac(lridx)->schedIDX());
	_update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetReacActive(uint tidx, uint ridx) const
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return false;

	if (tet->reac(lridx)->inactive() == true) return false;
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetReacActive(uint tidx, uint ridx, bool act)
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return;

	tet->reac(lridx)->setActive(act);

	SchedIDXVec updvec;
	updvec.push_back(tet->reac(lridx)->schedIDX());
	_update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetDiffD(uint tidx, uint didx) const
{
	assert (tidx < pTets.size());
	assert (didx < statedef()->countDiffs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint ldidx = tet->compdef()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED) return 0.0;

	return (tet->diff(ldidx)->dcst());

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetDiffD(uint tidx, uint didx, double dk)
{
	assert (tidx < pTets.size());
	assert (didx < statedef()->countDiffs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint ldidx = tet->compdef()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED) return;

	tet->diff(ldidx)->setDcst(dk);

	// Update the rates for all KProcs affected by this change
	SchedIDXVec updvec;
	updvec.push_back(tet->diff(ldidx)->schedIDX());
	_update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTetDiffActive(uint tidx, uint didx) const
{
	assert (tidx < pTets.size());
	assert (didx < statedef()->countDiffs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint ldidx = tet->compdef()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED) return false;

	if (tet->diff(ldidx)->inactive() == true) return false;
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTetDiffActive(uint tidx, uint didx, bool act)
{
	assert (tidx < pTets.size());
	assert (didx < statedef()->countDiffs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint ldidx = tet->compdef()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED) return;

	tet->diff(ldidx)->setActive(act);

	SchedIDXVec updvec;
	updvec.push_back(tet->diff(ldidx)->schedIDX());
	_update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacH(uint tidx, uint ridx) const
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return tet->reac(lridx)->h();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacC(uint tidx, uint ridx) const
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return tet->reac(lridx)->c();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetReacA(uint tidx, uint ridx) const
{
	assert (tidx < pTets.size());
	assert (ridx < statedef()->countReacs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint lridx = tet->compdef()->reacG2L(ridx);
	if (lridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return tet->reac(lridx)->rate();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTetDiffA(uint tidx, uint didx) const
{
	assert (tidx < pTets.size());
	assert (didx < statedef()->countDiffs());

	if (pTets[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Tetrahedron " << tidx << " has not been assigned to a compartment.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tet * tet = pTets[tidx];

	uint ldidx = tet->compdef()->diffG2L(didx);
	if (ldidx == ssolver::LIDX_UNDEFINED) return 0.0;

	return tet->diff(ldidx)->rate();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriArea(uint tidx) const
{
	assert (tidx < pTris.size());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	return pTris[tidx]->area();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriArea(uint tidx, double area)
{
	throw steps::NotImplErr();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriCount(uint tidx, uint sidx) const
{
	assert (tidx < pTris.size());
	assert (sidx < statedef()->countSpecs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];
	uint lsidx = tri->patchdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED) return 0.0;

	return tri->pools()[lsidx];
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriCount(uint tidx, uint sidx, double n)
{
	assert (tidx < pTris.size());
	assert (sidx < statedef()->countSpecs());
	assert (n >= 0.0);

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];
	uint lsidx = tri->patchdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED) return;

	double n_int = std::floor(n);
	double n_frc = n - n_int;
	uint c = static_cast<uint>(n_int);
	if (n_frc > 0.0)
	{
		double rand01 = rng()->getUnfIE();
		if (rand01 < n_frc) c++;
	}

	// Tri object updates counts in def level Comp object
	tri->setCount(lsidx, c);
	_updateSpec(tri, lsidx);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriAmount(uint tidx, uint sidx) const
{
	// following method does all necessary argument checking
	double count = _getTriCount(tidx, sidx);
	return count/steps::math::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriAmount(uint tidx, uint sidx, double m)
{
	// convert amount in mols to number of molecules
	double m2 = m * steps::math::AVOGADRO;
	// the following method does all the necessary argument checking
	_setTriCount(tidx, sidx, m2);
}

////////////////////////////////////////////////////////////////////////////////


bool stex::Tetexact::_getTriClamped(uint tidx, uint sidx) const
{
	assert (tidx < pTris.size());
	assert (sidx < statedef()->countSpecs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsidx = tri->patchdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED) return false;

	return tri->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriClamped(uint tidx, uint sidx, bool buf)
{
	assert (tidx < pTris.size());
	assert (sidx < statedef()->countSpecs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsidx = tri->patchdef()->specG2L(sidx);
	if (lsidx == ssolver::LIDX_UNDEFINED) return;

	tri->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacK(uint tidx, uint ridx) const
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return (tri->sreac(lsridx)->kcst());
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriSReacK(uint tidx, uint ridx, double kf)
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return;

	tri->sreac(lsridx)->setKcst(kf);

	// Update the rates for all KProcs affected by this change
	SchedIDXVec updvec;
	updvec.push_back(tri->sreac(lsridx)->schedIDX());
	_update(updvec);

}

////////////////////////////////////////////////////////////////////////////////

bool stex::Tetexact::_getTriSReacActive(uint tidx, uint ridx) const
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return false;

	if (tri->sreac(lsridx)->inactive() == true) return false;
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tetexact::_setTriSReacActive(uint tidx, uint ridx, bool act)
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return;

	tri->sreac(lsridx)->setActive(act);

	SchedIDXVec updvec;
	updvec.push_back(tri->sreac(lsridx)->schedIDX());
	_update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacH(uint tidx, uint ridx) const
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return tri->sreac(lsridx)->h();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacC(uint tidx, uint ridx) const
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return tri->sreac(lsridx)->c();
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tetexact::_getTriSReacA(uint tidx, uint ridx) const
{
	assert (tidx < pTris.size());
	assert (ridx < statedef()->countSReacs());

	if (pTris[tidx] == 0)
	{
    	std::ostringstream os;
    	os << "Triangle " << tidx << " has not been assigned to a patch.";
    	throw steps::ArgErr(os.str());
	}

	stex::Tri * tri = pTris[tidx];

	uint lsridx = tri->patchdef()->sreacG2L(ridx);
	if (lsridx == ssolver::LIDX_UNDEFINED) return 0.0;

	return tri->sreac(lsridx)->rate();
}

////////////////////////////////////////////////////////////////////////////////

// END
