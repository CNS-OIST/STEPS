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


// STL headers.
#include <cassert>
#include <algorithm>
#include <vector>
#include <sstream>

// STEPS headers.
#include "../common.h"
#include "tmcomp.hpp"
#include "../error.hpp"

NAMESPACE_ALIAS(steps::tetmesh, stetmesh);

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmComp::TmComp(std::string const & id, Tetmesh * container,
		                 std::vector<uint> const & tets)
: steps::wm::Comp(id, container, 0.0)
, pTetmesh(container)
, pTetsN(0)
, pTet_indices(0)
, pXmin(0.0)
, pXmax(0.0)
, pYmin(0.0)
, pYmax(0.0)
, pZmin(0.0)
, pZmax(0.0)
{
    if (pTetmesh == 0)
    {
    std::ostringstream os;
    os << "No mesh provided to Comp initializer function";
    throw steps::ArgErr(os.str());
    }

	// The maximum tetrahedron index in tetrahedral mesh
	uint maxidx = (pTetmesh-> countTets())-1;

	for(uint i=0; i< tets.size(); ++i)
	{
		// perform some checks on this tet
		bool included = false;
		for (uint j=0; j < pTetsN; ++j)
		{
			// check if tet has already occurred in this list (duplicate)
			if (tets[i] == pTet_indices[j])
			{
				included = true;
				break;
			}
		}
		if (included == true) continue;

		if (tets[i] > maxidx)
		{
			std::ostringstream os;
			os << "Invalid index supplied for tetrahedron #" << i << " in list.";
			throw steps::ArgErr(os.str());
		}
		if (pTetmesh->getTetComp(tets[i]) != 0)
		{
			std::ostringstream os;
			os << "Cannot add tetrahedron with index " << tets[i] << "(# " << i;
			os << " in list) to compartment; ";
			os << "tetrahedron belongs to a different compartment.";
			throw steps::ArgErr(os.str());
		}
		// add the tetrahedron to this compartment
		pTet_indices.push_back(tets[i]);
		// add this tetrahedron volume to the total
		pVol += pTetmesh->getTetVol(tets[i]);
		// perform annotation of tetrahedron
		pTetmesh->setTetComp(tets[i], this);
		// increment the number of tets in this compartment
		++ pTetsN;
	}
	assert (pTetsN == pTet_indices.size());

	// Compute the bounds of this compartment
	// first fetch vector of first tetrahedron's vertices
	uint * tet = pTetmesh->_getTet(pTet_indices[0]);
	// initialise min and max x coordinates with first x coordinate of first tet
	pXmin = pTetmesh->_getVertex(tet[0])[0];
	pXmax = pTetmesh->_getVertex(tet[0])[0];
	// initialise min and max y coordinates with first y coordinate of first tet
	pYmax = pTetmesh->_getVertex(tet[0])[1];
	pYmin = pTetmesh->_getVertex(tet[0])[1];
	// initialise min and max z coordinates with first z coordinate of first tet
	pZmin = pTetmesh->_getVertex(tet[0])[2];
	pZmax = pTetmesh->_getVertex(tet[0])[2];
	for (uint i=0; i< pTetsN; ++i)
	{
		// fetch the 4 vertices of the ith tet
		uint * tet = pTetmesh->_getTet(pTet_indices[i]);
		// compare each vertex to current values
		for (uint j=0; j<4; ++j)
		{
			double xtemp = pTetmesh->_getVertex(tet[j])[0];
			if (xtemp < pXmin) pXmin = xtemp;
			if (xtemp > pXmax) pXmax = xtemp;
			double ytemp = pTetmesh->_getVertex(tet[j])[1];
			if (ytemp < pYmin) pYmin = ytemp;
			if (ytemp > pYmax) pYmax  = ytemp;
			double ztemp = pTetmesh->_getVertex(tet[j])[2];
			if (ztemp < pZmin) pZmin = ztemp;
			if (ztemp > pZmax) pZmax = ztemp;
		}
	}
	// check values with some simple asserts
	assert (pZmin < pZmax);
	assert (pYmin < pYmax);
	assert (pXmin < pXmax);

}

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmComp::~TmComp(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::TmComp::setVol(double vol)
{
	std::ostringstream os;
	os << "Cannot set volume of Tetmesh comp object; vol calculated internally";
	throw steps::NotImplErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::TmComp::getBoundMin(void) const
{
    std::vector<double> b_min(3);
    b_min[0] = pXmin;
    b_min[1] = pYmin;
    b_min[2] = pZmin;
    return b_min;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::TmComp::getBoundMax(void) const
{
    std::vector<double> b_max(3);
    b_max[0] = pXmax;
    b_max[1] = pYmax;
    b_max[2] = pZmax;
    return b_max;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> stetmesh::TmComp::isTetInside(std::vector<uint> tet) const
{
	uint notets = tet.size();
	std::vector<bool> inside(notets);
	for (uint i=0; i < notets; ++i)
	{
		bool tetinside = false;
		for (uint j=0; j< pTet_indices.size() ; ++j)
		{
			if (tet[i] == pTet_indices[j])
			{
				tetinside = true;
				break;
			}
		}
		if (tetinside == true) inside[i] = true;
		else inside[i] = false;
	}
	return inside;
}

////////////////////////////////////////////////////////////////////////////////

// END


