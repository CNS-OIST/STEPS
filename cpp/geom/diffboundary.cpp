////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2011ÊOkinawa Institute of Science and Technology, Japan.
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


// STL headers.
#include <cassert>
#include <algorithm>
#include <vector>
#include <sstream>
#include <string>

// STEPS headers.
#include "../common.h"
#include "diffboundary.hpp"
#include "tmcomp.hpp"
#include "../error.hpp"
#include "tri.hpp"
#include "tet.hpp"

NAMESPACE_ALIAS(steps::tetmesh, stetmesh);

////////////////////////////////////////////////////////////////////////////////

stetmesh::DiffBoundary::DiffBoundary(std::string const & id, Tetmesh * container,
       	 std::vector<uint> const & tris)
: pID(id)
, pTetmesh(container)
, pTri_indices()
, pComps()
, pTrisN(0)

{
    if (pTetmesh == 0)
    {
    	std::ostringstream os;
    	os << "No mesh provided to Diffusion Boundary initializer function";
    	throw steps::ArgErr(os.str());
    }

    if (tris.size() == 0)
    {
    	std::ostringstream os;
    	os << "No triangles provided to Diffusion Boundary initializer function";
    	throw steps::ArgErr(os.str());
    }

    int * tri_tet_neighb = pTetmesh->_getTriTetNeighb(tris[0]);

	if (tri_tet_neighb[0] == -1 or tri_tet_neighb[1] == -1)
	{
    	std::ostringstream os;
		os << "Cannot add triangle with index " << tris[0] << "(#" << 0;
		os << " in list) to diffusion boundary; ";
		os << "triangle is on the mesh surface.";
		throw steps::ArgErr(os.str());
    }

	steps::tetmesh::TmComp * icmp = pTetmesh->getTetComp(tri_tet_neighb[0]);
	steps::tetmesh::TmComp * ocmp = pTetmesh->getTetComp(tri_tet_neighb[1]);

	if (icmp == 0 or ocmp == 0)
	{
    	std::ostringstream os;
		os << "Cannot add triangle with index " << tris[0] << "(#" << 0;
		os << " in list) to diffusion boundary; ";
		os << "triangle does not have an inner and outer compartment.";
		throw steps::ArgErr(os.str());
    }

	// Fill the vector now, but do the checks and throw exceptions in the loop
	pComps.push_back(icmp);
	pComps.push_back(ocmp);

    // The maximum triangle index in tetrahedral mesh
    uint maxidx = (pTetmesh->countTris() -1);
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

    	if (pTetmesh->getTriDiffBoundary(tris[i]) != 0)
    	{
    		std::ostringstream os;
    		os << "Cannot add triangle with index " << tris[i] << "(#" << i;
    		os << " in list) to diffusion boundary; ";
    		os << "triangle belongs to a different diffusion boundary.";
    		throw steps::ArgErr(os.str());
    	}

    	if (pTetmesh->getTriPatch(tris[i]) != 0)
    	{
    		std::ostringstream os;
    		os << "Cannot add triangle with index " << tris[i] << "(#" << i;
    		os << " in list) to diffusion boundary; ";
    		os << "triangle belongs to a patch.";
    		throw steps::ArgErr(os.str());
    	}


        int * tri_tet_neighb_temp = pTetmesh->_getTriTetNeighb(tris[i]);

    	if (tri_tet_neighb_temp[0] == -1 or tri_tet_neighb_temp[1] == -1)
    	{
        	std::ostringstream os;
    		os << "Cannot add triangle with index " << tris[i] << "(#" << i;
    		os << " in list) to diffusion boundary; ";
    		os << "triangle is on the mesh surface.";
    		throw steps::ArgErr(os.str());
        }

    	steps::tetmesh::TmComp * icmp_temp = pTetmesh->getTetComp(tri_tet_neighb_temp[0]);
    	steps::tetmesh::TmComp * ocmp_temp = pTetmesh->getTetComp(tri_tet_neighb_temp[1]);

    	if (icmp_temp == 0 or ocmp_temp == 0)
    	{
        	std::ostringstream os;
    		os << "Cannot add triangle with index " << tris[i] << "(#" << i;
    		os << " in list) to diffusion boundary; ";
    		os << "triangle does not have an inner and outer compartment.";
    		throw steps::ArgErr(os.str());
        }

    	if (icmp_temp != icmp)
    	{
    		// Try the other way around
    		steps::tetmesh::TmComp * holder = icmp_temp;
    		icmp_temp = ocmp_temp;;
    		ocmp_temp = holder;
    	}

    	if (icmp_temp == icmp && ocmp_temp == ocmp)
    	{
    		// add the triangle to this diffusion boundary
    		pTri_indices.push_back(tris[i]);
    		++pTrisN;
    		pTetmesh->setTriDiffBoundary(tris[i], this);

    	}
    	else
    	{
        	std::ostringstream os;
    		os << "Cannot add triangle with index " << tris[i] << "(#" << i;
    		os << " in list) to diffusion boundary; ";
    		os << "triangle does not have an inner and outer compartment.";
    		throw steps::ArgErr(os.str());
    	}

    } // end of loop over all tris (argument to constructor)

    /*
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
    }\*/

    pTetmesh->_handleDiffBoundaryAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::DiffBoundary::~DiffBoundary(void)
{
	if (pTetmesh == 0) return;
	_handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::DiffBoundary::setID(std::string const & id)
{
	assert(pTetmesh != 0);
	if (id == pID) return;
	// The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pTetmesh->_handleDiffBoundaryIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> stetmesh::DiffBoundary::isTriInside(std::vector<uint> tri) const
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

void stetmesh::DiffBoundary::_handleSelfDelete(void)
{
	pTetmesh->_handleDiffBoundaryDel(this);
	pComps.clear();
	pTri_indices.clear();
	pTetmesh = 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
