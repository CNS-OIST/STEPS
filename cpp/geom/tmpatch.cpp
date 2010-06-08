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
#include <string>

// STEPS headers.
#include "../common.h"
#include "tmpatch.hpp"
#include "tmcomp.hpp"
#include "../error.hpp"
#include "tri.hpp"
#include "tet.hpp"

NAMESPACE_ALIAS(steps::tetmesh, stetmesh);

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmPatch::TmPatch(std::string const & id, Tetmesh * container,
       	 std::vector<uint> const & tris, steps::tetmesh::TmComp* icomp,
       	 steps::tetmesh::TmComp* ocomp)
: steps::wm::Patch(id, container, icomp, ocomp,  0.0)
, pTetmesh(container)
, pTri_indices()
, pTrisN(0)

{
    if (pTetmesh == 0)
    {
    std::ostringstream os;
    os << "No mesh provided to Patch initializer function";
    throw steps::ArgErr(os.str());
    }

    // The maximum triangle index in tetrahedral mesh
    uint maxidx = (pTetmesh->countTris() -1);
    // vector to store all triangles to be 'flipped'
    std::vector<uint> fliplist = std::vector<uint>();
    for (uint i=0; i <tris.size(); ++i)
    {
    	// perform some checks on this triangle
    	bool included = false;
    	for(uint j = 0; j <pTrisN; ++j)
    	{
    		// check if tri has already occured in this list (duplicate)
    		if (tris[i] == pTri_indices[j])
    		{
    			included = true;
    			break;
    		}
    	}
    	if (included == true) continue;

    	if (tris[i] > maxidx)
    	{
    		std::ostringstream os;
    		os << "Invalid index supplied for triangle #" << i << " in list.";
    		throw steps::ArgErr(os.str());
    	}
    	if (pTetmesh->getTriPatch(tris[i]) != 0)
    	{
    		std::ostringstream os;
    		os << "Cannot add triangle with index " << tris[i] << "(#" << i;
    		os << " in list) to patch; ";
    		os << "triangle belongs to a different patch.";
    		throw steps::ArgErr(os.str());
    	}
    	// Check if the triangle neighbours belong to the correct compartments
    	// (bearing in mind their relative position might need to be 'flipped';
    	// add to 'flip list' if necessary
    	//
    	// fetch the tet indices of the inner and outer comps (may be -1, no comp)
    	int icmpidx = pTetmesh->_getTriTetNeighb(tris[i])[0];
    	int ocmpidx = pTetmesh->_getTriTetNeighb(tris[i])[1];
    	stetmesh::TmComp * icmp = 0;
    	stetmesh::TmComp * ocmp = 0;
    	if (icmpidx != -1) icmp = pTetmesh->getTetComp(icmpidx);
    	if (ocmpidx != -1) ocmp = pTetmesh->getTetComp(ocmpidx);
    	if (icmp == icomp && ocmp == ocomp)
    	{
    		// add the triangle to this patch
    		pTri_indices.push_back(tris[i]);
    		// add the triangle area to the total
    		pArea += pTetmesh->getTriArea(tris[i]);
    		// increment the number of triangles in this patch
    		pTetmesh->setTriPatch(tris[i], this);
    		++pTrisN;
    	}
    	else if (ocmp == icomp && icmp == ocomp)
    	{
    		// add the triangle to this patch
    		pTri_indices.push_back(tris[i]);
    		// add the triangle area to the total
    		pArea += pTetmesh->getTriArea(tris[i]);
    		// increment the number of triangles in this patch
    		pTetmesh->setTriPatch(tris[i], this);
    		// add this to the list to be 'flipped'
    		fliplist.push_back(tris[i]);
    		++pTrisN;
    	}
    	else
    	{
    		std::ostringstream os;
    		os << "Triangle (index " << tris[i] << ") cannot belong to this patch; ";
    		os << "inner and outer compartments don't match;";
    		throw steps::ArgErr(os.str());
    	}

    } // end of loop over all tris (argument to constructor)

    for (uint i=0; i< fliplist.size(); ++i)
    {
    	pTetmesh->_flipTriTetNeighb(fliplist[i]);
    }

    // Now loop over ALL triangles and check if normal points away from inner
    // tetrahedron. If it doesn't, flip the vertices
    for (uint t = 0; t < pTri_indices.size(); ++t)
    {
    	// first create the triangle object
    	steps::tetmesh::Tri tri(pTetmesh, pTri_indices[t]);
    	// now fetch the barycentre of the triangle
    	double * baryctri = tri._getBarycenter();
    	// now create the tetrahedron object of the INNER tet
    	steps::tetmesh::Tet tet = tri.getInnerTet();
    	// now fetch the barycentre of the inner tet
    	double * baryctet = tet._getBarycenter();
    	// now find the vector from the tet barycentre to the
    	// tri barycentre. If dot prod with tri normal is positive, tri verts don't have to be flipped
    	double vec1[3];
		std::vector<double> vec2(3);
    	vec1[0] = baryctri[0] - baryctet[0];
    	vec1[1] = baryctri[1] - baryctet[1];
    	vec1[2] = baryctri[2] - baryctet[2];
    	vec2 = tri.getNorm();
    	// calculate the dot product
    	double dp = 0.0;
    	dp += vec1[0]*vec2[0];
    	dp += vec1[1]*vec2[1];
    	dp += vec1[2]*vec2[2];
    	// flip the vertices if dot product is is negative
    	if (dp < 0.0) pTetmesh->_flipTriVerts(pTri_indices[t]);
    	// now check the dot product is positive
    	vec2 = tri.getNorm();
    	dp = 0.0;
    	dp += vec1[0]*vec2[0];
    	dp += vec1[1]*vec2[1];
    	dp += vec1[2]*vec2[2];
    	assert(dp > 0.0);
    }
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmPatch::TmPatch(std::string const & id, Tetmesh * container,
						   std::vector<uint> const & tris, steps::wm::Comp* wmicomp,
						   steps::wm::Comp* wmocomp)
: steps::wm::Patch(id, container, wmicomp, wmocomp, 0.0)
, pTetmesh(container)
, pTri_indices()
, pTrisN(0)

{
    if (pTetmesh == 0)
    {
		std::ostringstream os;
		os << "No mesh provided to Patch initializer function";
		throw steps::ArgErr(os.str());
    }

	// upcast the compartment pointers for this overloaded constructor
	steps::tetmesh::TmComp * icomp = 0;
	icomp = dynamic_cast<steps::tetmesh::TmComp *>(wmicomp);
	steps::tetmesh::TmComp * ocomp = 0;
	if (wmocomp != 0) ocomp = dynamic_cast<steps::tetmesh::TmComp *> (wmocomp);

    // The maximum triangle index in tetrahedral mesh
    uint maxidx = (pTetmesh->countTris() -1);
    // vector to store all triangles to be 'flipped'
    std::vector<uint> fliplist = std::vector<uint>();
    for (uint i=0; i <tris.size(); ++i)
    {
    	// perform some checks on this triangle
    	bool included = false;
    	for(uint j = 0; j <pTrisN; ++j)
    	{
    		// check if tri has already occured in this list (duplicate)
    		if (tris[i] == pTri_indices[j])
    		{
    			included = true;
    			break;
    		}
    	}
    	if (included == true) continue;

    	if (tris[i] > maxidx)
    	{
    		std::ostringstream os;
    		os << "Invalid index supplied for triangle #" << i << " in list.";
    		throw steps::ArgErr(os.str());
    	}
    	if (pTetmesh->getTriPatch(tris[i]) != 0)
    	{
    		std::ostringstream os;
    		os << "Cannot add triangle with index " << tris[i] << "(#" << i;
    		os << " in list) to patch; ";
    		os << "triangle belongs to a different patch.";
    		throw steps::ArgErr(os.str());
    	}
    	// Check if the triangle neighbours belong to the correct compartments
    	// (bearing in mind their relative position might need to be 'flipped';
    	// add to 'flip list' if necessary
    	//
    	// fetch the tet indices of the inner and outer comps (may be -1, no comp)
    	int icmpidx = pTetmesh->_getTriTetNeighb(tris[i])[0];
    	int ocmpidx = pTetmesh->_getTriTetNeighb(tris[i])[1];
    	stetmesh::TmComp * icmp = 0;
    	stetmesh::TmComp * ocmp = 0;
    	if (icmpidx != -1) icmp = pTetmesh->getTetComp(icmpidx);
    	if (ocmpidx != -1) ocmp = pTetmesh->getTetComp(ocmpidx);
    	if (icmp == icomp && ocmp == ocomp)
    	{
    		// add the triangle to this patch
    		pTri_indices.push_back(tris[i]);
    		// add the triangle area to the total
    		pArea += pTetmesh->getTriArea(tris[i]);
    		// increment the number of triangles in this patch
    		pTetmesh->setTriPatch(tris[i], this);
    		++pTrisN;
    	}
    	else if (ocmp == icomp && icmp == ocomp)
    	{
    		// add the triangle to this patch
    		pTri_indices.push_back(tris[i]);
    		// add the triangle area to the total
    		pArea += pTetmesh->getTriArea(tris[i]);
    		// increment the number of triangles in this patch
    		pTetmesh->setTriPatch(tris[i], this);
    		// add this to the list to be 'flipped'
    		fliplist.push_back(tris[i]);
    		++pTrisN;
    	}
    	else
    	{
    		std::ostringstream os;
    		os << "Triangle (index " << tris[i] << ") cannot belong to this patch; ";
    		os << "inner and outer compartments don't match;";
    		throw steps::ArgErr(os.str());
    	}

    } // end of loop over all tris (argument to constructor)

    for (uint i=0; i< fliplist.size(); ++i)
    {
    	pTetmesh->_flipTriTetNeighb(fliplist[i]);
    }

    // Now loop over ALL triangles and check if normal points away from inner
    // tetrahedron. If it doesn't, flip the vertices
    for (uint t = 0; t < pTri_indices.size(); ++t)
    {
    	// first create the triangle object
    	steps::tetmesh::Tri tri(pTetmesh, pTri_indices[t]);
    	// now fetch the barycentre of the triangle
    	double * baryctri = tri._getBarycenter();
    	// now create the tetrahedron object of the INNER tet
    	steps::tetmesh::Tet tet = tri.getInnerTet();
    	// now fetch the barycentre of the inner tet
    	double * baryctet = tet._getBarycenter();
    	// now find the vector from the tet barycentre to the
    	// tri barycentre. If dot prod with tri normal is positive, tri verts don't have to be flipped
    	double vec1[3];
		std::vector<double> vec2(3);
    	vec1[0] = baryctri[0] - baryctet[0];
    	vec1[1] = baryctri[1] - baryctet[1];
    	vec1[2] = baryctri[2] - baryctet[2];
    	vec2 = tri.getNorm();
    	// calculate the dot product
    	double dp = 0.0;
    	dp += vec1[0]*vec2[0];
    	dp += vec1[1]*vec2[1];
    	dp += vec1[2]*vec2[2];
    	// flip the vertices if dot product is is negative
    	if (dp < 0.0) pTetmesh->_flipTriVerts(pTri_indices[t]);
    	// now check the dot product is positive
    	vec2 = tri.getNorm();
    	dp = 0.0;
    	dp += vec1[0]*vec2[0];
    	dp += vec1[1]*vec2[1];
    	dp += vec1[2]*vec2[2];
    	assert(dp > 0.0);
    }
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmPatch::~TmPatch(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::TmPatch::setArea(double area)
{
	std::ostringstream os;
	os << "Cannot set area of Tetmesh patch object; area calculated internally";
	throw steps::NotImplErr(os.str());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> stetmesh::TmPatch::isTriInside(std::vector<uint> tri) const
{
	uint notris = tri.size();
	std::vector<bool> inside(notris);
	for (uint i=0; i < notris; ++i)
	{
		bool triinside = false;
		for (uint j=0; j< pTri_indices.size() ; ++j)
		{
			if (tri[i] == pTri_indices[j])
			{
				triinside = true;
				break;
			}
		}
		if (triinside == true) inside[i] = true;
		else inside[i] = false;
	}
	return inside;
}

////////////////////////////////////////////////////////////////////////////////

// END
