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
#include <iostream>
#include <vector>
#include <sstream>

// STEPS headers.
#include "../common.h"
#include "../math/tetrahedron.hpp"
#include "tetmesh.hpp"
#include "../math/triangle.hpp"
#include "../error.hpp"

NAMESPACE_ALIAS(steps::tetmesh, stetmesh);

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tetmesh::Tetmesh(uint nverts, uint ntets, uint ntris)
: Geom()
, pSetupDone(false)
, pVertsN(nverts)
, pVerts(0)
, pTrisN(ntris) // initialise this member with user supplied number, but this will be modified later
, pTris(0)
, pTri_areas(0)
, pTri_barycs(0)
, pTri_norms(0)
, pTri_patches(0)
, pTri_tet_neighbours(0)
, pTetsN(ntets)
, pTets(0)
, pTet_vols(0)
, pTet_barycentres(0)
, pTet_comps(0)
, pTet_tri_neighbours(0)
, pTet_tet_neighbours(0)
, pTris_user(0)
, pXmin(0.0)
, pXmax(0.0)
, pYmin(0.0)
, pYmax(0.0)
, pZmin(0.0)
, pZmax(0.0)
{	/*
    assert(pVertsN > 0);
    assert(pTrisN > 0);
    assert(pTetsN > 0);
    pVerts = new double[pVertsN * 3];
    pTris = new uint[pTrisN * 3];
    pTri_areas = new double[pTrisN];
    pTri_norms = new double[pTrisN * 3];
    pTri_patches = new stetmesh::TmPatch*[pTrisN];
    // Initialise patch pointers to zero (fill_n doesn't work for zero pointers)
    for (uint i=0; i<pTrisN; ++i) pTri_patches[i] = 0;
    pTri_tet_neighbours = new int[pTrisN * 2];
    // Initialise each triangle's tetrahedron neighbours to -1, which implies no tetrahedron
    // After setup each triangle should have 2 neighbouring tetrahedra,
    // or 1 (if triangle is a surface triangle)
    // if both neighbouring tetrahedra are -1 this is an error
    std::fill_n(pTri_tet_neighbours, (pTrisN *2), -1);
    pTets = new uint[pTetsN * 4];
    pTet_vols = new double[pTetsN];
    pTet_comps = new stetmesh::TmComp*[pTetsN];
    // Initialise comp pointers to zero (fill_n doesn't work for zero pointers)
    for (uint i=0; i<pTetsN; ++i) pTet_comps[i] = 0;
    pTet_tri_neighbours = new uint[pTetsN * 4];
    pTet_tet_neighbours = new int[pTetsN * 4];
    // Initialise the Tet_tet_neighbour array with -1, which implies no neighbour
    // (tetrahedron is on the surface)
    // This is currently the only method of setting the -1, the main algorithm
    // does not explicitly check if it is a surface tet. Normally only 1 face could be on the surface,
    // but in some cases 2 or even 3 faces could be on the surface
    std::fill_n(pTet_tet_neighbours, (pTetsN * 4), -1);
    */

	// DEBUG 7/4/09 Constructor didn't have enough margin for user error, particularly the
	// user had to supply the total number of surface triangles in the mesh (i.e. not just
	// surface ones)
	// Now this constructor creates tetrahedron and vertex arrays, but makes pTris_user
	// available to user to supply triangle information. Total number of triangles will then
	// be found by setup() after user has supplied all tetrahedron (and vertex, and possible
	// some triangle) data
	//
	if (nverts < 4)
	{
	    std::ostringstream os;
	    os << "Number of vertices must be positive and greater than 3.";
	    throw steps::ArgErr(os.str());
	}
	if (ntris < 0)
	{
	    std::ostringstream os;
	    os << "Number of triangles must be positive or zero.";
	    throw steps::ArgErr(os.str());
	}
	if (ntets < 1)
	{
	    std::ostringstream os;
	    os << "Number of tetrahedra must be positive and greater than zero.";
	    throw steps::ArgErr(os.str());
	}

    pVerts = new double[pVertsN * 3];
    pTets = new uint[pTetsN * 4];
    pTet_vols = new double[pTetsN];
    pTet_barycentres = new double[pTetsN*3];
    pTet_comps = new stetmesh::TmComp*[pTetsN];
    // Initialise comp pointers to zero (fill_n doesn't work for zero pointers)
    for (uint i=0; i<pTetsN; ++i) pTet_comps[i] = 0;
    pTet_tri_neighbours = new uint[pTetsN * 4];
    pTet_tet_neighbours = new int[pTetsN * 4];
    // Initialise the Tet_tet_neighbour array with -1, which implies no neighbour
    // (tetrahedron is on the surface)
    // This is currently the only method of setting the -1, the main algorithm
    // does not explicitly check if it is a surface tet. Normally only 1 face could be on the surface,
    // but in some cases 2 or even 3 faces could be on the surface
    std::fill_n(pTet_tet_neighbours, (pTetsN * 4), -1);

    // Make the pTris_user array available for user entry. Actual number of ALL triangles
    // will be found later
    pTris_user = new uint[ntris * 3];

}

////////////////////////////////////////////////////////////////////////////////////////////

stetmesh::Tetmesh::Tetmesh(std::vector<double> const & verts,
						   std::vector<uint> const & tets,
		                   std::vector<uint> const & tris)
: Geom()
, pSetupDone(false)
, pVertsN(0)
, pVerts(0)
, pTrisN(0)
, pTris(0)
, pTri_areas(0)
, pTri_norms(0)
, pTri_barycs(0)
, pTri_patches(0)
, pTri_tet_neighbours(0)
, pTris_user(0) // not used by this contructor
, pTetsN(0)
, pTets(0)
, pTet_vols(0)
, pTet_comps(0)
, pTet_tri_neighbours(0)
, pTet_tet_neighbours(0)
, pXmin(0.0)
, pXmax(0.0)
, pYmin(0.0)
, pYmax(0.0)
, pZmin(0.0)
, pZmax(0.0)
{
	// check the vectors are of the expected size
	if ((verts.size() % 3) != 0 || (tris.size() % 3) != 0
		|| (tets.size() % 4) != 0)
	{
	    std::ostringstream os;
	    os << "Tables supplied to Tet mesh initialiser function ";
	    os << " are not of the expected dimensions";
	    throw steps::ArgErr(os.str());
	}
	pVertsN = verts.size() / 3;
	pTetsN = tets.size() / 4;
	if (pVertsN == 0 || pTetsN == 0)
	{
		std::ostringstream os;
	    os << "Vertex table or Tet table not supplied to Tet mesh initialiser function.";
	    throw steps::ArgErr(os.str());
	}

	pVerts = new double[pVertsN * 3];
	// copy the supplied vertices information to pVerts member
	for (uint i = 0; i < pVertsN*3; ++i) pVerts[i] = verts[i];

	pTets = new uint[pTetsN * 4];
	// copy the supplied tetrahedron information to the pTets member
	for (uint i = 0; i< pTetsN*4; ++i) pTets[i] = tets[i];
	pTet_vols = new double[pTetsN];
    pTet_barycentres = new double[pTetsN*3];
	pTet_comps = new stetmesh::TmComp*[pTetsN];
	// Initialise comp pointers to zero (fill_n doesn't appear to work for zero pointers)
	for (uint i=0; i<pTetsN; ++i) pTet_comps[i] = 0;
	pTet_tri_neighbours = new uint[pTetsN * 4];
	pTet_tet_neighbours = new int[pTetsN * 4];
	// Initialise the Tet_tet_neighbour array with -1, which implies no neighbour
	// (tetrahedron is on the surface)
	// this is currently the only method of setting the -1, the main algorithm
	// does not explicitly check if it is a surface tet. Normally only 1 face could be on the surface,
	// but in some cases 2 or even 3 faces could be on the surface
	std::fill_n(pTet_tet_neighbours, (pTetsN * 4), -1);

	// Create some triangle data arrays initially with maximum theoretical size. Each
	// tetrahedron must share at least 2 triangles with a neighbour -> maximum
	// number of triangles is (TetsN x 3) + 1
	uint trisn_max = (pTetsN*3) + 1;

	uint * tris_temp = new uint[trisn_max * 3];
	int * tri_tet_neighbours_temp = new int[trisn_max * 2];
	std::fill_n(tri_tet_neighbours_temp, (trisn_max *2), -1);

	// Add any triangles supplied first; maintaining indices
	uint tris_added = 0;
	for (uint i=0; i<tris.size(); i+=3)
	{
		int triele[3] = {tris[i+0], tris[i+1], tris[i+2]};
		std::sort(triele, triele + 3);
		tris_temp[i+0] = triele[0];
		tris_temp[i+1] = triele[1];
		tris_temp[i+2] = triele[2];
		tris_added++;
	}

	uint tettetadded = 0;
	// Loop over all tetrahedra and fill tris_temp, pTet_tri_neighbours,
	// pTet_tet_neighbours, tri_tet_neighbours_temp
	for (uint tet=0; tet < pTetsN; ++tet)
	{
		/// set this tetrahedron's volume
		///
		// find pointers to this tetrahedron's 4 vertices
		double * vert0 = pVerts + (3 * pTets[tet*4]);
		double * vert1 = pVerts + (3 * pTets[(tet*4)+1]);
		double * vert2 = pVerts + (3 * pTets[(tet*4)+2]);
		double * vert3 = pVerts + (3 * pTets[(tet*4)+3]);
		pTet_vols[tet] = steps::math::tet_vol(vert0, vert1, vert2, vert3);
		steps::math::tet_barycenter(vert0, vert1, vert2, vert3,pTet_barycentres + tet*3);

		// create array for triangle formed by edges (0,1,2)
		// will also sort for ease of comparison
		uint tri_vert0[3] = {pTets[tet*4], pTets[(tet*4)+1], pTets[(tet*4)+2]};
		std::sort(tri_vert0, tri_vert0 + 3);
		// create array for triangle formed by edges (0,1,3)
		uint tri_vert1[3] = {pTets[tet*4], pTets[(tet*4)+1], pTets[(tet*4)+3]};
		std::sort(tri_vert1, tri_vert1 + 3);
		// create array for triangle formed by edges (0,2,3)
		uint tri_vert2[3] = {pTets[tet*4], pTets[(tet*4)+2], pTets[(tet*4)+3]};
		std::sort(tri_vert2, tri_vert2 + 3);
		// create array for triangle formed by edges (1,2,3)
		uint tri_vert3[3] = {pTets[(tet*4)+1], pTets[(tet*4)+2], pTets[(tet*4)+3]};
		std::sort(tri_vert3, tri_vert3 + 3);

		int tri0idx = -1;
		int tri1idx = -1;
		int tri2idx = -1;
		int tri3idx = -1;
		uint trisfound = 0;	/// unused at the moment
		for (uint tri = 0; tri < tris_added; ++tri)
		{
			uint thistri[3] = {tris_temp[tri*3], tris_temp[(tri*3)+1], tris_temp[(tri*3)+2]};
			// should be no need to sort them now // std::sort(thistri, thistri + 3);
			if (tri0idx == -1 && array_srt_cmp(tri_vert0, thistri, 3)) tri0idx = tri;
			if (tri1idx == -1 && array_srt_cmp(tri_vert1, thistri, 3)) tri1idx = tri;
			if (tri2idx == -1 && array_srt_cmp(tri_vert2, thistri, 3)) tri2idx = tri;
			if (tri3idx == -1 && array_srt_cmp(tri_vert3, thistri, 3)) tri3idx = tri;
		}
		// First add the triangles that are not already included
		if (tri0idx == -1)
		{
			tris_temp[tris_added*3] = tri_vert0[0];
			tris_temp[(tris_added*3)+1] = tri_vert0[1];
			tris_temp[(tris_added*3)+2] = tri_vert0[2];
			// label this triangle and increment triangle counter
			tri0idx = tris_added++;
		}
		if (tri1idx == -1)
		{
			tris_temp[tris_added*3] = tri_vert1[0];
			tris_temp[(tris_added*3)+1] = tri_vert1[1];
			tris_temp[(tris_added*3)+2] = tri_vert1[2];
			// label this triangle and increment triangle counter
			tri1idx = tris_added++;
		}
		if (tri2idx == -1)
		{
			tris_temp[tris_added*3] = tri_vert2[0];
			tris_temp[(tris_added*3)+1] = tri_vert2[1];
			tris_temp[(tris_added*3)+2] = tri_vert2[2];
			// label this triangle and increment triangle counter
			tri2idx = tris_added++;
		}
		if (tri3idx == -1)
		{
			tris_temp[tris_added*3] = tri_vert3[0];
			tris_temp[(tris_added*3)+1] = tri_vert3[1];
			tris_temp[(tris_added*3)+2] = tri_vert3[2];
			// label this triangle and increment triangle counter
			tri3idx = tris_added++;
		}

		// Use this information to fill neighbours information
		pTet_tri_neighbours[tet*4] = tri0idx;
		pTet_tri_neighbours[(tet*4)+1] = tri1idx;
		pTet_tri_neighbours[(tet*4)+2] = tri2idx;
		pTet_tri_neighbours[(tet*4)+3] = tri3idx;


		/* BUGFIX 13/08/09 In following if else section:
		 * Indexing of triangle and tet neighbours were not aligned.
		 * This was causing descrepancies in the setting of the scaled dcsts in the simulation
		 * layer algorithm, causing concentration differences at steady-state.
		 * It was decided to fix the indexing here rather than alter the sim algorithm.
		 */

		// Add this tet to neighbours of first triangle
		if (tri_tet_neighbours_temp[tri0idx*2] == -1) tri_tet_neighbours_temp[tri0idx*2] = tet;
		else if (tri_tet_neighbours_temp[(tri0idx*2)+1] == -1)
		{

			tri_tet_neighbours_temp[(tri0idx*2)+1] = tet;

			// This triangle has a 2 tet neighbours now- time to tell those tets they are neighbours
			uint tet1 = tri_tet_neighbours_temp[tri0idx*2];

			assert(pTet_tet_neighbours[(tet*4)] == -1);
			pTet_tet_neighbours[(tet*4)] = tet1;

			// Find this triangle's index from tet1's triangle neighbours
			for (uint idx = 0; idx <=4; ++idx)
			{
				assert(idx != 4);
				int tet1tri0idx = pTet_tri_neighbours[(tet1*4)+idx];
				assert(tet1tri0idx >= 0);
				if (tet1tri0idx == tri0idx)
				{
					assert(pTet_tet_neighbours[(tet1*4)+idx] == -1);
					pTet_tet_neighbours[(tet1*4)+idx] = tet;
					break;
				}
			}

			tettetadded+=2;
		}
		else assert(false);
		if (tri_tet_neighbours_temp[tri1idx*2] == -1) tri_tet_neighbours_temp[tri1idx*2] = tet;
		else if (tri_tet_neighbours_temp[(tri1idx*2)+1] == -1)
		{
			tri_tet_neighbours_temp[(tri1idx*2)+1] = tet;
			// This triangle has a 2 tet neighbours now- time to tell those tets they are neighbours
			uint tet1 = tri_tet_neighbours_temp[tri1idx*2];

			assert(pTet_tet_neighbours[(tet*4)+1] == -1);
			pTet_tet_neighbours[(tet*4)+1] = tet1;

			// Find this triangle's index from tet1's triangle neighbours
			for (uint idx = 0; idx <= 4; ++idx)
			{
				assert(idx != 4);
				int tet1tri1idx = pTet_tri_neighbours[(tet1*4)+idx];
				assert(tet1tri1idx >= 0);
				if (tet1tri1idx == tri1idx)
				{
					assert(pTet_tet_neighbours[(tet1*4)+idx] == -1);
					pTet_tet_neighbours[(tet1*4)+idx] = tet;
					break;
				}
			}

			tettetadded+=2;
		}
		else assert(false);
		if (tri_tet_neighbours_temp[tri2idx*2] == -1) tri_tet_neighbours_temp[tri2idx*2] = tet;
		else if (tri_tet_neighbours_temp[(tri2idx*2)+1] == -1)
		{
			tri_tet_neighbours_temp[(tri2idx*2)+1] = tet;
			// This triangle has a 2 tet neighbours now- time to tell those tets they are neighbours
			uint tet1 = tri_tet_neighbours_temp[tri2idx*2];

			assert(pTet_tet_neighbours[(tet*4)+2] == -1);
			pTet_tet_neighbours[(tet*4)+2] = tet1;

			// Find this triangle's index from tet1's triangle neighbours
			for (uint idx = 0; idx <= 4; ++idx)
			{
				assert(idx != 4);
				int tet1tri2idx = pTet_tri_neighbours[(tet1*4)+idx];
				assert(tet1tri2idx >= 0);
				if (tet1tri2idx == tri2idx)
				{
					assert(pTet_tet_neighbours[(tet1*4)+idx] == -1);
					pTet_tet_neighbours[(tet1*4)+idx] = tet;
					break;
				}
			}

			tettetadded+=2;
		}
		else assert(false);
		if (tri_tet_neighbours_temp[tri3idx*2] == -1) tri_tet_neighbours_temp[tri3idx*2] = tet;
		else if (tri_tet_neighbours_temp[(tri3idx*2)+1] == -1)
		{
			tri_tet_neighbours_temp[(tri3idx*2)+1] = tet;
			// This triangle has a 2 tet neighbours now- time to tell those tets they are neighbours
			uint tet1 = tri_tet_neighbours_temp[tri3idx*2];

			assert(pTet_tet_neighbours[(tet*4)+3] == -1);
			pTet_tet_neighbours[(tet*4)+3] = tet1;

			// Find this triangle's index from tet1's triangle neighbours
			for (uint idx = 0; idx <= 4; ++idx)
			{
				assert(idx != 4);
				int tet1tri3idx = pTet_tri_neighbours[(tet1*4)+idx];
				assert(tet1tri3idx >= 0);
				if (tet1tri3idx == tri3idx)
				{
					assert(pTet_tet_neighbours[(tet1*4)+idx] == -1);
					pTet_tet_neighbours[(tet1*4)+idx] = tet;
					break;
				}
			}

			tettetadded+=2;
		}
		else assert(false);
	} // end of loop over all tetrahedra

	// Most tetrahedra should have had 4 tetneighbours added, but fewer for surface tets
	assert (tettetadded < pTetsN*4);

    ////////////////////////////////////////////////////////////////////////

	// We know know how many triangles are in the mesh
	pTrisN = tris_added;

	pTris = new uint[pTrisN * 3];
	// copy the supplied triangles information to pTris member
	for (uint i = 0; i < pTrisN*3; ++i) pTris[i] = tris_temp[i];
	pTri_areas = new double[pTrisN];
	pTri_barycs = new double[pTrisN * 3];
	pTri_norms = new double[pTrisN * 3];
	pTri_patches = new stetmesh::TmPatch*[pTrisN];
	// Initialise patch pointers to zero (fill_n doesn't work for zero pointers)
	for (uint i=0; i<pTrisN; ++i) pTri_patches[i] = 0;
	pTri_tet_neighbours = new int[pTrisN * 2];
	for (uint i=0; i < pTrisN*2; ++i) pTri_tet_neighbours[i] = tri_tet_neighbours_temp[i];

	// Free up memory from temporary tables
	delete[] tris_temp;
	delete[] tri_tet_neighbours_temp;

	/// loop over all triangles and set pTri_areas and pTri_norms
	for (uint tri = 0; tri < pTrisN; ++tri)
	{
		/// set this triangle's area
		///
		// find pointers to this triangle's 3 vertices
		double * vert0 = pVerts + (3 * pTris[tri*3]);
		double * vert1 = pVerts + (3 * pTris[(tri*3)+1]);
		double * vert2 = pVerts + (3 * pTris[(tri*3)+2]);
		// call steps::math method to set the area
		pTri_areas[tri] = steps::math::triArea(vert0, vert1, vert2);

		// create arrays to store this triangle's barycenter and normal
		double baryc[3];
		double norm[3];
		// call math methods to find the barycenter, normal and store in arrays
		steps::math::triBarycenter(vert0, vert1, vert2, baryc);
		steps::math::triNormal(vert0, vert1, vert2, norm);

		// set the barycenter
		pTri_barycs[tri*3] = baryc[0];
		pTri_barycs[(tri*3)+1] = baryc[1];
		pTri_barycs[(tri*3)+2] = baryc[2];
		// set the normal
		pTri_norms[tri*3] = norm[0];
		pTri_norms[(tri*3)+1] = norm[1];
		pTri_norms[(tri*3)+2] = norm[2];
	}

    ////////////////////////////////////////////////////////////////////////

	/// Find the minimal and maximal boundary values
	double xmin = pVerts[0];
	double xmax = pVerts[0];
	double ymin = pVerts[1];
	double ymax = pVerts[1];
	double zmin = pVerts[2];
	double zmax = pVerts[2];

	for (uint i=1; i<pVertsN; ++i)
	{
		if (pVerts[i*3] < xmin) xmin = pVerts[i*3];
		if (pVerts[i*3] > xmax) xmax = pVerts[i*3];
		if (pVerts[(i*3)+1] < ymin) ymin = pVerts[(i*3)+1];
		if (pVerts[(i*3)+1] > ymax) ymax = pVerts[(i*3)+1];
		if (pVerts[(i*3)+2] < zmin) zmin = pVerts[(i*3)+2];
		if (pVerts[(i*3)+2] > zmax) zmax = pVerts[(i*3)+2];
	}
	pXmin = xmin;
	pXmax = xmax;
	pYmin = ymin;
	pYmax = ymax;
	pZmin = zmin;
	pZmax = zmax;

	// Constructor has completed all necessary setting-up: Set flag.
	pSetupDone = true;
}


////////////////////////////////////////////////////////////////////////////////

stetmesh::Tetmesh::Tetmesh(std::vector<double> const & verts,
		std::vector<uint> const & tris,
   		std::vector<double> const & tri_areas,
   		std::vector<double> const & tri_norms,
   		std::vector<int> const & tri_tet_neighbs,
   		std::vector<uint> const & tets,
   		std::vector<double> const & tet_vols,
   		std::vector<double> const & tet_barycs,
   		std::vector<uint> const & tet_tri_neighbs,
   		std::vector<int> const & tet_tet_neighbs)
: Geom()
, pSetupDone(false)
, pVertsN(0)
, pVerts(0)
, pTrisN(0)
, pTris(0)
, pTri_areas(0)
, pTri_norms(0)
, pTri_barycs(0)
, pTri_patches(0)
, pTri_tet_neighbours(0)
, pTris_user(0) // not used by this contructor
, pTetsN(0)
, pTets(0)
, pTet_vols(0)
, pTet_comps(0)
, pTet_tri_neighbours(0)
, pTet_tet_neighbours(0)
, pXmin(0.0)
, pXmax(0.0)
, pYmin(0.0)
, pYmax(0.0)
, pZmin(0.0)
, pZmax(0.0)
{
	// Check all vectors are of expected size
   	if ((verts.size() % 3) != 0 )
   	{
   	    std::ostringstream os;
   	    os << "Vertex coordinate data supplied to Tet mesh initialiser function ";
   	    os << " is not of the expected dimensions (3)";
   	    throw steps::ArgErr(os.str());
   	}
   	if ((tris.size() % 3) != 0 )
   	{
   	    std::ostringstream os;
   	    os << "Triangle node data supplied to Tet mesh initialiser function ";
   	    os << " is not of the expected dimensions (3)";
   	    throw steps::ArgErr(os.str());
   	}
   	if ((tets.size() % 4) != 0 )
   	{
   	    std::ostringstream os;
   	    os << "Tetrahedron node data supplied to Tet mesh initialiser function ";
   	    os << " is not of the expected dimensions (4)";
   	    throw steps::ArgErr(os.str());
   	}

	pVertsN = verts.size() / 3;
	pTrisN = tris.size() /3;
	pTetsN = tets.size() / 4;
	if (pVertsN == 0 || pTrisN == 0 || pTetsN == 0)
	{
		std::ostringstream os;
	    os << "Vertex, Triangle or Tet table not supplied to Tet mesh initialiser function.";
	    throw steps::ArgErr(os.str());
	}

	if (tri_areas.size() != pTrisN)
	{
		std::ostringstream os;
	    os << "Triangle areas table is not of expected size.";
	    throw steps::ArgErr(os.str());
	}
	if (tri_norms.size() != pTrisN*3)
	{
		std::ostringstream os;
	    os << "Triangle normals table is not of expected size.";
	    throw steps::ArgErr(os.str());
	}
	if (tri_tet_neighbs.size() != pTrisN*2)
	{
		std::ostringstream os;
	    os << "Triangle tet neighbours table is not of expected size.";
	    throw steps::ArgErr(os.str());
	}
	if (tet_vols.size() != pTetsN)
	{
		std::ostringstream os;
	    os << "Tetrahedron volumes table is not of expected size.";
	    throw steps::ArgErr(os.str());
	}
	if (tet_barycs.size() != pTetsN*3)
	{
		std::ostringstream os;
	    os << "Tetrahedron barycenters table is not of expected size.";
	    throw steps::ArgErr(os.str());
	}
	if (tet_tri_neighbs.size() != pTetsN*4)
	{
		std::ostringstream os;
	    os << "Tetrahedron tri neighbours table is not of expected size.";
	    throw steps::ArgErr(os.str());
	}
	if (tet_tet_neighbs.size() != pTetsN*4)
	{
		std::ostringstream os;
	    os << "Tetrahedron tet neighbours table is not of expected size.";
	    throw steps::ArgErr(os.str());
	}

	//

	pVerts = new double[pVertsN * 3];
	pTris = new uint[pTrisN * 3];
	pTri_areas = new double[pTrisN];
	pTri_norms = new double[pTrisN*3];
	pTri_barycs = new double[pTrisN*3];
	pTri_tet_neighbours = new int[pTrisN*2];
	pTets = new uint[pTetsN * 4];
	pTet_vols = new double[pTetsN];
	pTet_barycentres = new double[pTetsN*3];
	pTet_tri_neighbours = new uint[pTetsN*4];
	pTet_tet_neighbours = new int[pTetsN*4];

	pTet_comps = new stetmesh::TmComp*[pTetsN];
	// Initialise comp pointers to zero (fill_n doesn't appear to work for zero pointers)
	for (uint i=0; i<pTetsN; ++i) pTet_comps[i] = 0;
	pTri_patches = new stetmesh::TmPatch*[pTrisN];
	// Initialise patch pointers to zero (fill_n doesn't appear to work for zero pointers)
	for (uint i=0; i<pTrisN; ++i) pTri_patches[i] = 0;

	// copy the supplied vertices information to pVerts member
	for (uint i = 0; i < pVertsN*3; ++i) pVerts[i] = verts[i];

	// Loop once over triangle indexes and copy the information to the local structures
	for (uint t = 0; t < pTrisN; ++t)
	{
		pTris[t*3] = tris[t*3];
		pTris[(t*3)+1] = tris[(t*3)+1];
		pTris[(t*3)+2] = tris[(t*3)+2];

		pTri_areas[t] = tri_areas[t];

		pTri_norms[t*3] = tri_norms[t*3];
		pTri_norms[(t*3)+1] = tri_norms[(t*3)+1];
		pTri_norms[(t*3)+2] = tri_norms[(t*3)+2];

		pTri_tet_neighbours[t*2] = tri_tet_neighbs[t*2];
		pTri_tet_neighbours[(t*2)+1] = tri_tet_neighbs[(t*2)+1];

		// The barycenter of the triangle is not store in the text file right now
		/// set this triangle's barycenter
		///
		// find pointers to this triangle's 3 vertices
		double * vert0 = pVerts + (3 * tris[t*3]);
		double * vert1 = pVerts + (3 * tris[(t*3)+1]);
		double * vert2 = pVerts + (3 * tris[(t*3)+2]);
		// create array to store this triangle's barycenter
		double baryc[3];
		steps::math::triBarycenter(vert0, vert1, vert2, baryc);
		// set the barycenter
		pTri_barycs[t*3] = baryc[0];
		pTri_barycs[(t*3)+1] = baryc[1];
		pTri_barycs[(t*3)+2] = baryc[2];
	}

	for (uint t = 0; t < pTetsN; ++t)
	{
		pTets[t*4] = tets[t*4];
		pTets[(t*4)+1] = tets[(t*4)+1];
		pTets[(t*4)+2] = tets[(t*4)+2];
		pTets[(t*4)+3] = tets[(t*4)+3];

		pTet_vols[t] = tet_vols[t];

		pTet_barycentres[t*3] = tet_barycs[t*3];
		pTet_barycentres[(t*3)+1] = tet_barycs[(t*3)+1];
		pTet_barycentres[(t*3)+2] = tet_barycs[(t*3)+2];

		pTet_tri_neighbours[t*4] = tet_tri_neighbs[t*4];
		pTet_tri_neighbours[(t*4)+1] = tet_tri_neighbs[(t*4)+1];
		pTet_tri_neighbours[(t*4)+2] = tet_tri_neighbs[(t*4)+2];
		pTet_tri_neighbours[(t*4)+3] = tet_tri_neighbs[(t*4)+3];

		pTet_tet_neighbours[t*4] = tet_tet_neighbs[t*4];
		pTet_tet_neighbours[(t*4)+1] = tet_tet_neighbs[(t*4)+1];
		pTet_tet_neighbours[(t*4)+2] = tet_tet_neighbs[(t*4)+2];
		pTet_tet_neighbours[(t*4)+3] = tet_tet_neighbs[(t*4)+3];

	}


	/// Find the minimal and maximal boundary values
	double xmin = pVerts[0];
	double xmax = pVerts[0];
	double ymin = pVerts[1];
	double ymax = pVerts[1];
	double zmin = pVerts[2];
	double zmax = pVerts[2];

	for (uint i=1; i<pVertsN; ++i)
	{
		if (pVerts[i*3] < xmin) xmin = pVerts[i*3];
		if (pVerts[i*3] > xmax) xmax = pVerts[i*3];
		if (pVerts[(i*3)+1] < ymin) ymin = pVerts[(i*3)+1];
		if (pVerts[(i*3)+1] > ymax) ymax = pVerts[(i*3)+1];
		if (pVerts[(i*3)+2] < zmin) zmin = pVerts[(i*3)+2];
		if (pVerts[(i*3)+2] > zmax) zmax = pVerts[(i*3)+2];
	}
	pXmin = xmin;
	pXmax = xmax;
	pYmin = ymin;
	pYmax = ymax;
	pZmin = zmin;
	pZmax = zmax;

	// Constructor has completed all necessary setting-up: Set flag.
	pSetupDone = true;
}


////////////////////////////////////////////////////////////////////////////////

stetmesh::Tetmesh::~Tetmesh(void)
{
	// Memory created in 1st constructor must be freed if setup hasn't been called
	if (pSetupDone == false) delete[] pTris_user;

	delete[] pVerts;
	delete[] pTris;
	delete[] pTri_areas;
	delete[] pTri_norms;
	delete[] pTri_patches;
	delete[] pTri_tet_neighbours;
	delete[] pTets;
	delete[] pTet_vols;
	delete[] pTet_comps;
	delete[] pTet_tri_neighbours;
	delete[] pTet_tet_neighbours;
	delete[] pTet_barycentres;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::setVertex(uint vidx, double x, double y, double z)
{
	if (vidx >= pVertsN)
	{
		std::ostringstream os;
		os << "Vertex index is out of range.";
		throw steps::ArgErr(os.str());
	}
	if (pSetupDone == true)
	{
		std::ostringstream os;
		os << "Cannot set vertex after mesh has been setup.";
		throw steps::ArgErr(os.str());
	}

    uint vidx2 = vidx * 3;
    pVerts[vidx2++] = x;
    pVerts[vidx2++] = y;
    pVerts[vidx2] = z;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::setTri(uint tidx, uint vidx0, uint vidx1, uint vidx2)
{
	if (tidx >= pTrisN)
	{
		std::ostringstream os;
		os << "Triangle index out of user-supplied range.";
		throw steps::ArgErr(os.str());
	}
	if (vidx0 >= pVertsN || vidx1 >= pVertsN || vidx2 >= pVertsN)
	{
		std::ostringstream os;
		os << "Vertex index is out of range.";
		throw steps::ArgErr(os.str());
	}
	if (pSetupDone == true)
	{
		std::ostringstream os;
		os << "Cannot set triangle after mesh has been setup.";
		throw steps::ArgErr(os.str());
	}

    uint tidx2 = tidx * 3;
    pTris_user[tidx2++] = vidx0;
    pTris_user[tidx2++] = vidx1;
    pTris_user[tidx2] = vidx2;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::setTet
(
 uint tidx,
 uint vidx0, uint vidx1,
 uint vidx2, uint vidx3
)
{
	if (tidx >= pTetsN)
	{
		std::ostringstream os;
		os << "Tetrahedron index out of range.";
		throw steps::ArgErr(os.str());
	}
	if (vidx0 >= pVertsN || vidx1 >= pVertsN || vidx2 >= pVertsN || vidx3 >= pVertsN)
	{
		std::ostringstream os;
		os << "Vertex index out of range.";
		throw steps::ArgErr(os.str());
	}
	if (pSetupDone == true)
	{
		std::ostringstream os;
		os << "Cannot set tetrahedron after mesh has been setup.";
		throw steps::ArgErr(os.str());
	}

    uint tidx2 = tidx * 4;
    pTets[tidx2++] = vidx0;
    pTets[tidx2++] = vidx1;
    pTets[tidx2++] = vidx2;
    pTets[tidx2] = vidx3;
}

////////////////////////////////////////////////////////////////////////////////

/// This method only makes sense if the first constructor is used and user has
/// supplied all mesh information prior to call,
/// i.e. all set***() methods should have filled pVerts and pTets
//  pTris may not be complete, or even empty, so setup initialises pTris.
/// Set all members that can be determined from the mesh information:
/// bounding box values, areas, volumes, and neighbouring triangle and
/// tetrahedral information
//
// DEBUG: User may have only supplied surface triangle information. Ok up to now,
// but this code must find the 'real' number of triangles, similar to the
// second constructor
//
void stetmesh::Tetmesh::setup(void)
{
	if (pSetupDone == true)
	{
		std::ostringstream os;
		os << "Setup completed.";
		throw steps::ArgErr(os.str());
	}
	/* replaced the below code with code similar to second constructor. See notes
	// in first constructor

	/// Find the minimal and maximal boundary values
	double xmin = pVerts[0];
	double xmax = pVerts[0];
	double ymin = pVerts[1];
	double ymax = pVerts[1];
	double zmin = pVerts[2];
	double zmax = pVerts[2];

	for (uint i=1; i < pVertsN; ++i)
	{
		if (pVerts[i*3] < xmin) xmin = pVerts[i*3];
		if (pVerts[i*3] > xmax) xmax = pVerts[i*3];
		if (pVerts[(i*3)+1] < ymin) ymin = pVerts[(i*3)+1];
		if (pVerts[(i*3)+1] > ymax) ymax = pVerts[(i*3)+1];
		if (pVerts[(i*3)+2] < zmin) zmin = pVerts[(i*3)+2];
		if (pVerts[(i*3)+2] > zmax) zmax = pVerts[(i*3)+2];
	}
	pXmin = xmin;
	pXmax = xmax;
	pYmin = ymin;
	pYmax = ymax;
	pZmin = zmin;
	pZmax = zmax;

	/// loop over all triangles and set pTri_areas and pTri_norms
	for (uint tri = 0; tri < pTrisN; ++tri)
	{
		/// set this triangle's area
		///
		// find pointers to this triangle's 3 vertices
		double * vert0 = pVerts + (3 * pTris[tri*3]);
		double * vert1 = pVerts + (3 * pTris[(tri*3)+1]);
		double * vert2 = pVerts + (3 * pTris[(tri*3)+2]);
		// call steps::math method to set the area
		pTri_areas[tri] = steps::math::triArea(vert0, vert1, vert2);

		// create an array to store this triangle's normal
		double norm[3];
		// call steps::math method to find the normal and store in norm array
		steps::math::triNormal(vert0, vert1, vert2, norm);
		// set the normal
		pTri_norms[tri*3] = norm[0];
		pTri_norms[(tri*3)+1] = norm[1];
		pTri_norms[(tri*3)+2] = norm[2];
	}


	/// loop over all tetrahedra and set pTet_vols, pTet_tri_neighbours,
	/// pTri_tet_neighbours and pTet_tet_neighbours
	for (uint tet = 0; tet < pTetsN; ++tet)
	{
		/// set this tetrahedron's volume
		///
		// find pointers to this tetrahedron's 4 vertices
		double * vert0 = pVerts + (3 * pTets[tet*4]);
		double * vert1 = pVerts + (3 * pTets[(tet*4)+1]);
		double * vert2 = pVerts + (3 * pTets[(tet*4)+2]);
		double * vert3 = pVerts + (3 * pTets[(tet*4)+3]);
		pTet_vols[tet] = steps::math::tet_vol(vert0, vert1, vert2, vert3);


		/// find the 4 triangle faces of this tetrahedron, by vertex indices
		/// and sort into ascending order by applying algorithm sort()
		/// function for ease of comparing
		///
		// create array for triangle formed by edges (0,1,2)
		uint tri_vert0[3] = {pTets[tet*4], pTets[(tet*4)+1], pTets[(tet*4)+2]};
		std::sort(tri_vert0, tri_vert0 + 3);
		// create array for triangle formed by edges (0,1,3)
		uint tri_vert1[3] = {pTets[tet*4], pTets[(tet*4)+1], pTets[(tet*4)+3]};
		std::sort(tri_vert1, tri_vert1 + 3);
		// create array for triangle formed by edges (0,2,3)
		uint tri_vert2[3] = {pTets[tet*4], pTets[(tet*4)+2], pTets[(tet*4)+3]};
		std::sort(tri_vert2, tri_vert2 + 3);
		// create array for triangle formed by edges (1,2,3)
		uint tri_vert3[3] = {pTets[(tet*4)+1], pTets[(tet*4)+2], pTets[(tet*4)+3]};
		std::sort(tri_vert3, tri_vert3 + 3);

		// loop over all tetrahedra to find this tetrahedron's 4 neighbours
		uint tets_added = 0;
		for(uint tet_neighb = 0; tet_neighb < pTetsN; ++tet_neighb)
		{
			// pass if this is the same tetrahedron
			if (tet_neighb == tet) continue;

			/// similar to above, store arrays of this tetrahedron's triangle faces and sort
			// create array for triangle formed by edges (0,1,2)
			uint tri_vert0_neighb[3] = {pTets[tet_neighb*4], pTets[(tet_neighb*4)+1],
			pTets[(tet_neighb*4)+2]};
			std::sort(tri_vert0_neighb, tri_vert0_neighb + 3);
			// create array for triangle formed by edges (0,1,3)
			uint tri_vert1_neighb[3] = {pTets[tet_neighb*4], pTets[(tet_neighb*4)+1],
			pTets[(tet_neighb*4)+3]};
			std::sort(tri_vert1_neighb, tri_vert1_neighb + 3);
			// create array for triangle formed by edges (0,2,3)
			uint tri_vert2_neighb[3] = {pTets[tet_neighb*4], pTets[(tet_neighb*4)+2],
			pTets[(tet_neighb*4)+3]};
			std::sort(tri_vert2_neighb, tri_vert2_neighb + 3);
			// create array for triangle formed by edges (1,2,3)
			uint tri_vert3_neighb[3] = {pTets[(tet_neighb*4)+1], pTets[(tet_neighb*4)+2],
			pTets[(tet_neighb*4)+3]};
			std::sort(tri_vert3_neighb, tri_vert3_neighb + 3);

			// compare the faces of the 2 tetrahedra to each other.
			// If one matches they are neighbours
			// first face of tetrahedron of interest (the one we are filling information for)
			if (array_srt_cmp(tri_vert0, tri_vert0_neighb, 3)
				|| array_srt_cmp(tri_vert0, tri_vert1_neighb, 3)
				|| array_srt_cmp(tri_vert0, tri_vert2_neighb, 3)
				|| array_srt_cmp(tri_vert0, tri_vert3_neighb, 3))
			{
				// fill Tet_tet_neighbour with this information.
				// NOT also filling the information for the neighbouring tet
				// (index tet_neighb) at this point
				pTet_tet_neighbours[tet*4] = tet_neighb;
				++ tets_added;
				// break if all 4 neighbouring tets have been found
				if (tets_added == 4) break;
			}
			// now compare second face
			else if (array_srt_cmp(tri_vert1, tri_vert0_neighb, 3)
					 || array_srt_cmp(tri_vert1, tri_vert1_neighb, 3)
					 || array_srt_cmp(tri_vert1, tri_vert2_neighb, 3)
					 || array_srt_cmp(tri_vert1, tri_vert3_neighb, 3))
			{
				pTet_tet_neighbours[(tet*4)+1] = tet_neighb;
				++ tets_added;
				if (tets_added == 4) break;
			}
			else if (array_srt_cmp(tri_vert2, tri_vert0_neighb, 3)
					 || array_srt_cmp(tri_vert2, tri_vert1_neighb, 3)
					 || array_srt_cmp(tri_vert2, tri_vert2_neighb, 3)
					 || array_srt_cmp(tri_vert2, tri_vert3_neighb, 3))
			{
				pTet_tet_neighbours[(tet*4)+2] = tet_neighb;
				++ tets_added;
				if (tets_added == 4) break;
			}
			else if(array_srt_cmp(tri_vert3, tri_vert0_neighb, 3)
					|| array_srt_cmp(tri_vert3, tri_vert1_neighb, 3)
					|| array_srt_cmp(tri_vert3, tri_vert2_neighb, 3)
					|| array_srt_cmp(tri_vert3, tri_vert3_neighb, 3))
			{
				pTet_tet_neighbours[(tet*4)+3] = tet_neighb;
				++ tets_added;
				if(tets_added == 4) break;
			}
		}
		// A surface tet will usually have 3 neighbours, and in some cases a tet may
		// have only 2 or even 1 neighbour; but no neighbours is an error.
		assert(tets_added > 0);

		// now compare to each triangle
		uint tris_added = 0;
		for (uint tri = 0; tri < pTrisN; ++tri)
		{
			// create array to sort triangle vertices into ascending order
			uint tri_sorted[3] = {pTris[tri*3], pTris[(tri*3)+1], pTris[(tri*3)+2]};
			std::sort(tri_sorted, tri_sorted + 3);
			// compare triangles
			if(array_srt_cmp(tri_sorted, tri_vert0, 3))
			{
				pTet_tri_neighbours[(tet*4) + 0] = tri;
				// also add this tetrahedron to this triangle's neighbours
				// no test if tetrahedron is inner or outer at this point;
				// this may have to be 'flipped' when triangle is added to a patch
				if(pTri_tet_neighbours[tri*2] == -1) pTri_tet_neighbours[tri*2] = tet;
				else if (pTri_tet_neighbours[(tri*2)+1] == -1) pTri_tet_neighbours[(tri*2)+1] = tet;
				else assert(false);

				++ tris_added;
				if (tris_added == 4) break;
			}
			else if(array_srt_cmp(tri_sorted, tri_vert1, 3))
			{
				pTet_tri_neighbours[(tet*4) + 1] = tri;
				if(pTri_tet_neighbours[tri*2] == -1) pTri_tet_neighbours[tri*2] = tet;
				else if (pTri_tet_neighbours[(tri*2)+1] == -1) pTri_tet_neighbours[(tri*2)+1] = tet;
				else assert(false);

				++ tris_added;
				if (tris_added == 4) break;
			}
			else if(array_srt_cmp(tri_sorted, tri_vert2, 3))
			{
				pTet_tri_neighbours[(tet*4) + 2] = tri;
				if(pTri_tet_neighbours[tri*2] == -1) pTri_tet_neighbours[tri*2] = tet;
				else if (pTri_tet_neighbours[(tri*2)+1] == -1) pTri_tet_neighbours[(tri*2)+1] = tet;
				else assert(false);

				++ tris_added;
				if (tris_added == 4) break;
			}
			else if(array_srt_cmp(tri_sorted, tri_vert3, 3))
			{
				pTet_tri_neighbours[(tet*4) + 3] = tri;
				if(pTri_tet_neighbours[tri*2] == -1) pTri_tet_neighbours[tri*2] = tet;
				else if (pTri_tet_neighbours[(tri*2)+1] == -1) pTri_tet_neighbours[(tri*2)+1] = tet;
				else assert(false);

				++ tris_added;
				if (tris_added == 4) break;
			}
		}
		// All tets should have 4 neighbouring triangles
		assert (tris_added == 4);

	} // end of loop over all tetrahedra
	*/
	// End of replaced code


	// Create some triangle data arrays initially with maximum theoretical size. Each
	// tetrahedron must share at least 2 triangles with a neighbour -> maximum
	// number of triangles is (TetsN x 3) + 1
	uint trisn_max = (pTetsN*3) + 1;

	uint * tris_temp = new uint[trisn_max * 3];
	int * tri_tet_neighbours_temp = new int[trisn_max * 2];
	std::fill_n(tri_tet_neighbours_temp, (trisn_max *2), -1);

	// Add any triangles supplied first; maintaining indices
	uint tris_added = 0;
	for (uint i=0; i<pTrisN*3; i+=3)
	{
		int triele[3] = {pTris_user[i+0], pTris_user[i+1], pTris_user[i+2]};
		std::sort(triele, triele + 3);
		tris_temp[i+0] = triele[0];
		tris_temp[i+1] = triele[1];
		tris_temp[i+2] = triele[2];
		tris_added++;
	}
	// Now can free memory for array of user-supplied triangle information
	delete[] pTris_user;

	uint tettetadded = 0;
	// Loop over all tetrahedra and fill tris_temp, pTet_tri_neighbours,
	// pTet_tet_neighbours, tri_tet_neighbours_temp
	for (uint tet=0; tet < pTetsN; ++tet)
	{
		/// set this tetrahedron's volume
		///
		// find pointers to this tetrahedron's 4 vertices
		double * vert0 = pVerts + (3 * pTets[tet*4]);
		double * vert1 = pVerts + (3 * pTets[(tet*4)+1]);
		double * vert2 = pVerts + (3 * pTets[(tet*4)+2]);
		double * vert3 = pVerts + (3 * pTets[(tet*4)+3]);
		pTet_vols[tet] = steps::math::tet_vol(vert0, vert1, vert2, vert3);
		steps::math::tet_barycenter(vert0, vert1, vert2, vert3,pTet_barycentres + tet*3);

		// create array for triangle formed by edges (0,1,2)
		// will also sort for ease of comparison
		uint tri_vert0[3] = {pTets[tet*4], pTets[(tet*4)+1], pTets[(tet*4)+2]};
		std::sort(tri_vert0, tri_vert0 + 3);
		// create array for triangle formed by edges (0,1,3)
		uint tri_vert1[3] = {pTets[tet*4], pTets[(tet*4)+1], pTets[(tet*4)+3]};
		std::sort(tri_vert1, tri_vert1 + 3);
		// create array for triangle formed by edges (0,2,3)
		uint tri_vert2[3] = {pTets[tet*4], pTets[(tet*4)+2], pTets[(tet*4)+3]};
		std::sort(tri_vert2, tri_vert2 + 3);
		// create array for triangle formed by edges (1,2,3)
		uint tri_vert3[3] = {pTets[(tet*4)+1], pTets[(tet*4)+2], pTets[(tet*4)+3]};
		std::sort(tri_vert3, tri_vert3 + 3);

		int tri0idx = -1;
		int tri1idx = -1;
		int tri2idx = -1;
		int tri3idx = -1;
		uint trisfound = 0;	/// unused at the moment
		for (uint tri = 0; tri < tris_added; ++tri)
		{
			uint thistri[3] = {tris_temp[tri*3], tris_temp[(tri*3)+1], tris_temp[(tri*3)+2]};
			// should be no need to sort them now // std::sort(thistri, thistri + 3);
			if (tri0idx == -1 && array_srt_cmp(tri_vert0, thistri, 3)) tri0idx = tri;
			if (tri1idx == -1 && array_srt_cmp(tri_vert1, thistri, 3)) tri1idx = tri;
			if (tri2idx == -1 && array_srt_cmp(tri_vert2, thistri, 3)) tri2idx = tri;
			if (tri3idx == -1 && array_srt_cmp(tri_vert3, thistri, 3)) tri3idx = tri;
		}
		// First add the triangles that are not already included
		if (tri0idx == -1)
		{
			tris_temp[tris_added*3] = tri_vert0[0];
			tris_temp[(tris_added*3)+1] = tri_vert0[1];
			tris_temp[(tris_added*3)+2] = tri_vert0[2];
			// label this triangle and increment triangle counter
			tri0idx = tris_added++;
		}
		if (tri1idx == -1)
		{
			tris_temp[tris_added*3] = tri_vert1[0];
			tris_temp[(tris_added*3)+1] = tri_vert1[1];
			tris_temp[(tris_added*3)+2] = tri_vert1[2];
			// label this triangle and increment triangle counter
			tri1idx = tris_added++;
		}
		if (tri2idx == -1)
		{
			tris_temp[tris_added*3] = tri_vert2[0];
			tris_temp[(tris_added*3)+1] = tri_vert2[1];
			tris_temp[(tris_added*3)+2] = tri_vert2[2];
			// label this triangle and increment triangle counter
			tri2idx = tris_added++;
		}
		if (tri3idx == -1)
		{
			tris_temp[tris_added*3] = tri_vert3[0];
			tris_temp[(tris_added*3)+1] = tri_vert3[1];
			tris_temp[(tris_added*3)+2] = tri_vert3[2];
			// label this triangle and increment triangle counter
			tri3idx = tris_added++;
		}

		// Use this information to fill neighbours information
		pTet_tri_neighbours[tet*4] = tri0idx;
		pTet_tri_neighbours[(tet*4)+1] = tri1idx;
		pTet_tri_neighbours[(tet*4)+2] = tri2idx;
		pTet_tri_neighbours[(tet*4)+3] = tri3idx;

		// Add this tet to neighbours of first triangle
		if (tri_tet_neighbours_temp[tri0idx*2] == -1) tri_tet_neighbours_temp[tri0idx*2] = tet;
		else if (tri_tet_neighbours_temp[(tri0idx*2)+1] == -1)
		{
			tri_tet_neighbours_temp[(tri0idx*2)+1] = tet;

			// This triangle has a 2 tet neighbours now- time to tell those tets they are neighbours
			uint tet1 = tri_tet_neighbours_temp[tri0idx*2];

			// Add the tet in the main loop to it's neighbour's neighbours
			if (pTet_tet_neighbours[tet1*4] == -1) pTet_tet_neighbours[tet1*4] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+1] == -1) pTet_tet_neighbours[(tet1*4)+1] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+2] == -1) pTet_tet_neighbours[(tet1*4)+2] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+3] == -1) pTet_tet_neighbours[(tet1*4)+3] = tet;
			else assert(false);
			tettetadded++;

			// Add tet1 to tet's neighbours
			if (pTet_tet_neighbours[tet*4] == -1) pTet_tet_neighbours[tet*4] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+1] == -1) pTet_tet_neighbours[(tet*4)+1] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+2] == -1) pTet_tet_neighbours[(tet*4)+2] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+3] == -1) pTet_tet_neighbours[(tet*4)+3] = tet1;
			else assert(false);
			tettetadded++;
		}
		else assert(false);
		if (tri_tet_neighbours_temp[tri1idx*2] == -1) tri_tet_neighbours_temp[tri1idx*2] = tet;
		else if (tri_tet_neighbours_temp[(tri1idx*2)+1] == -1)
		{
			tri_tet_neighbours_temp[(tri1idx*2)+1] = tet;
			// This triangle has a 2 tet neighbours now- time to tell those tets they are neighbours
			uint tet1 = tri_tet_neighbours_temp[tri1idx*2];

			// Add the tet in the main loop to it's neighbour's neighbours
			if (pTet_tet_neighbours[tet1*4] == -1) pTet_tet_neighbours[tet1*4] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+1] == -1) pTet_tet_neighbours[(tet1*4)+1] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+2] == -1) pTet_tet_neighbours[(tet1*4)+2] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+3] == -1) pTet_tet_neighbours[(tet1*4)+3] = tet;
			else assert(false);
			tettetadded++;

			// Add tet1 to tet's neighbours
			if (pTet_tet_neighbours[tet*4] == -1) pTet_tet_neighbours[tet*4] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+1] == -1) pTet_tet_neighbours[(tet*4)+1] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+2] == -1) pTet_tet_neighbours[(tet*4)+2] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+3] == -1) pTet_tet_neighbours[(tet*4)+3] = tet1;
			else assert(false);
			tettetadded++;
		}
		else assert(false);
		if (tri_tet_neighbours_temp[tri2idx*2] == -1) tri_tet_neighbours_temp[tri2idx*2] = tet;
		else if (tri_tet_neighbours_temp[(tri2idx*2)+1] == -1)
		{
			tri_tet_neighbours_temp[(tri2idx*2)+1] = tet;
			// This triangle has a 2 tet neighbours now- time to tell those tets they are neighbours
			uint tet1 = tri_tet_neighbours_temp[tri2idx*2];

			// Add the tet in the main loop to it's neighbour's neighbours
			if (pTet_tet_neighbours[tet1*4] == -1) pTet_tet_neighbours[tet1*4] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+1] == -1) pTet_tet_neighbours[(tet1*4)+1] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+2] == -1) pTet_tet_neighbours[(tet1*4)+2] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+3] == -1) pTet_tet_neighbours[(tet1*4)+3] = tet;
			else assert(false);
			tettetadded++;

			// Add tet1 to tet's neighbours
			if (pTet_tet_neighbours[tet*4] == -1) pTet_tet_neighbours[tet*4] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+1] == -1) pTet_tet_neighbours[(tet*4)+1] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+2] == -1) pTet_tet_neighbours[(tet*4)+2] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+3] == -1) pTet_tet_neighbours[(tet*4)+3] = tet1;
			tettetadded++;
		}
		else assert(false);
		if (tri_tet_neighbours_temp[tri3idx*2] == -1) tri_tet_neighbours_temp[tri3idx*2] = tet;
		else if (tri_tet_neighbours_temp[(tri3idx*2)+1] == -1)
		{
			tri_tet_neighbours_temp[(tri3idx*2)+1] = tet;
			// This triangle has a 2 tet neighbours now- time to tell those tets they are neighbours
			uint tet1 = tri_tet_neighbours_temp[tri3idx*2];

			// Add the tet in the main loop to it's neighbour's neighbours
			if (pTet_tet_neighbours[tet1*4] == -1) pTet_tet_neighbours[tet1*4] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+1] == -1) pTet_tet_neighbours[(tet1*4)+1] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+2] == -1) pTet_tet_neighbours[(tet1*4)+2] = tet;
			else if (pTet_tet_neighbours[(tet1*4)+3] == -1) pTet_tet_neighbours[(tet1*4)+3] = tet;
			else assert(false);
			tettetadded++;

			// Add tet1 to tet's neighbours
			if (pTet_tet_neighbours[tet*4] == -1) pTet_tet_neighbours[tet*4] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+1] == -1) pTet_tet_neighbours[(tet*4)+1] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+2] == -1) pTet_tet_neighbours[(tet*4)+2] = tet1;
			else if (pTet_tet_neighbours[(tet*4)+3] == -1) pTet_tet_neighbours[(tet*4)+3] = tet1;
			tettetadded++;
		}
		else assert(false);
	} // end of loop over all tetrahedra

	// Most tetrahedra should have had 4 tetneighbours added, but fewer for surface tets
	assert (tettetadded < pTetsN*4);

    ////////////////////////////////////////////////////////////////////////

	// We know know how many triangles are in the mesh
	// Replace previous user-supplied number with total number
	pTrisN = tris_added;

	pTris = new uint[pTrisN * 3];
	// copy the supplied triangles information to pTris member
	for (uint i = 0; i < pTrisN*3; ++i) pTris[i] = tris_temp[i];
	pTri_areas = new double[pTrisN];
	pTri_norms = new double[pTrisN * 3];
	pTri_patches = new stetmesh::TmPatch*[pTrisN];
	// Initialise patch pointers to zero (fill_n doesn't work for zero pointers)
	for (uint i=0; i<pTrisN; ++i) pTri_patches[i] = 0;
	pTri_tet_neighbours = new int[pTrisN * 2];
	for (uint i=0; i < pTrisN*2; ++i) pTri_tet_neighbours[i] = tri_tet_neighbours_temp[i];

	// Free up memory from temporary tables
	delete[] tris_temp;
	delete[] tri_tet_neighbours_temp;

	/// loop over all triangles and set pTri_areas and pTri_norms
	for (uint tri = 0; tri < pTrisN; ++tri)
	{
		/// set this triangle's area
		///
		// find pointers to this triangle's 3 vertices
		double * vert0 = pVerts + (3 * pTris[tri*3]);
		double * vert1 = pVerts + (3 * pTris[(tri*3)+1]);
		double * vert2 = pVerts + (3 * pTris[(tri*3)+2]);
		// call steps::math method to set the area
		pTri_areas[tri] = steps::math::triArea(vert0, vert1, vert2);

		// create an array to store this triangle's normal
		double norm[3];
		// call steps::math method to find the normal and store in norm array
		steps::math::triNormal(vert0, vert1, vert2, norm);
		// set the normal
		pTri_norms[tri*3] = norm[0];
		pTri_norms[(tri*3)+1] = norm[1];
		pTri_norms[(tri*3)+2] = norm[2];
	}

    ////////////////////////////////////////////////////////////////////////

	/// Find the minimal and maximal boundary values
	double xmin = pVerts[0];
	double xmax = pVerts[0];
	double ymin = pVerts[1];
	double ymax = pVerts[1];
	double zmin = pVerts[2];
	double zmax = pVerts[2];

	for (uint i=1; i<pVertsN; ++i)
	{
		if (pVerts[i*3] < xmin) xmin = pVerts[i*3];
		if (pVerts[i*3] > xmax) xmax = pVerts[i*3];
		if (pVerts[(i*3)+1] < ymin) ymin = pVerts[(i*3)+1];
		if (pVerts[(i*3)+1] > ymax) ymax = pVerts[(i*3)+1];
		if (pVerts[(i*3)+2] < zmin) zmin = pVerts[(i*3)+2];
		if (pVerts[(i*3)+2] > zmax) zmax = pVerts[(i*3)+2];
	}
	pXmin = xmin;
	pXmax = xmax;
	pYmin = ymin;
	pYmax = ymax;
	pZmin = zmin;
	pZmax = zmax;


	pSetupDone = true;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double>  stetmesh::Tetmesh::getVertex(uint vidx) const
{
	if (vidx >= pVertsN)
	{
		std::ostringstream os;
		os << "Vertex index is out of range.";
		throw steps::ArgErr(os.str());
	}
	std::vector<double> vert(3);
	vert[0] = pVerts[vidx*3];
	vert[1] = pVerts[(vidx*3)+1];
	vert[2] = pVerts[(vidx*3)+2];
	return vert;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> stetmesh::Tetmesh::getTri(uint tidx) const
{
	if (tidx >= pTrisN)
	{
		std::ostringstream os;
		os << "Triangle index is out of range.";
		throw steps::ArgErr(os.str());
	}
	std::vector<uint> tri(3);
	tri[0] = pTris[tidx*3];
	tri[1] = pTris[(tidx*3)+1];
	tri[2] = pTris[(tidx*3)+2];
	return tri;
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tetmesh::getTriArea(uint tidx) const
{
	if (tidx >= pTrisN)
	{
		std::ostringstream os;
		os << "Triangle index is out of range.";
		throw steps::ArgErr(os.str());
	}
	return (pTri_areas[tidx]);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getTriBarycenter(uint tidx) const
{
	if (tidx >= pTrisN)
	{
		std::ostringstream os;
		os << "Triangle index is out of range.";
		throw steps::ArgErr(os.str());
	}
	std::vector<double> baryc(3);
	baryc[0] = pTri_barycs[(tidx*3)];
	baryc[1] = pTri_barycs[(tidx*3)+1];
	baryc[2] = pTri_barycs[(tidx*3)+2];

	return baryc;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getTriNorm(uint tidx) const
{
    assert(pSetupDone == true);
	if (tidx >= pTrisN)
	{
		std::ostringstream os;
		os << "Triangle index is out of range.";
		throw steps::ArgErr(os.str());
	}
	std::vector<double> norm(3);
	norm[0] = pTri_norms[(tidx*3)];
	norm[1] = pTri_norms[(tidx*3)+1];
	norm[2] = pTri_norms[(tidx*3)+2];

	return norm;
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmPatch * stetmesh::Tetmesh::getTriPatch(uint tidx) const
{
	if (tidx >= pTrisN)
	{
		std::ostringstream os;
		os << "Triangle index is out of range.";
		throw steps::ArgErr(os.str());
	}
    return pTri_patches[tidx];
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::setTriPatch(uint tidx, stetmesh::TmPatch * patch)
{
	if (tidx >= pTrisN)
	{
		std::ostringstream os;
		os << "Triangle index is out of range.";
		throw steps::ArgErr(os.str());
	}
	pTri_patches[tidx] = patch;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> stetmesh::Tetmesh::getTriTetNeighb(uint tidx) const
{
    assert(pSetupDone == true);
	if (tidx >= pTrisN)
	{
		std::ostringstream os;
		os << "Triangle index is out of range.";
		throw steps::ArgErr(os.str());
	}
	std::vector<int> tettemp(2);
	tettemp[0] = pTri_tet_neighbours[tidx*2];
	tettemp[1] = pTri_tet_neighbours[(tidx*2)+1];
	return tettemp;
}

////////////////////////////////////////////////////////////////////////////////

// Created by weiliang 2010.02.02
std::vector<int> stetmesh::Tetmesh::getTriBoundary(void) const 
{
    assert(pSetupDone == true);
    std::vector<int> tribounds;
    for (int t = 0; t < pTrisN; t++) {
        std::vector<int> trineighbor = getTriTetNeighb(t);
        if (trineighbor[0] == -1 || trineighbor[1] == -1) {
            tribounds.push_back(t);
        }
    }
    return tribounds;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_flipTriTetNeighb(uint tidx)
{
    assert(pSetupDone == true);
	assert (tidx < pTrisN);

	int tettemp = pTri_tet_neighbours[tidx*2];
	pTri_tet_neighbours[tidx*2] = pTri_tet_neighbours[(tidx*2)+1];
	pTri_tet_neighbours[(tidx*2)+1] = tettemp;
}

void stetmesh::Tetmesh::_flipTriVerts(uint tidx)
{
    assert(pSetupDone == true);
	assert(tidx < pTrisN);

	uint tmp = pTris[tidx*3];
	pTris[tidx*3] = pTris[(tidx*3)+1];
	pTris[(tidx*3)+1] = tmp;
	// also recalculate the triangle's normal; the long way (though should be
	// negative of previous normal
	double * vert0 = pVerts + (3 * pTris[tidx*3]);
	double * vert1 = pVerts + (3 * pTris[(tidx*3)+1]);
	double * vert2 = pVerts + (3 * pTris[(tidx*3)+2]);
	double norm[3];
	// call steps::math method to find the normal and store in norm array
	steps::math::triNormal(vert0, vert1, vert2, norm);
	// set the normal
	pTri_norms[tidx*3] = norm[0];
	pTri_norms[(tidx*3)+1] = norm[1];
	pTri_norms[(tidx*3)+2] = norm[2];


}
////////////////////////////////////////////////////////////////////////////////
/*
 uint stetmesh::Tetmesh::addTet(uint vidx0, uint vidx1, uint vidx2, uint vidx3)
 {

 }

 ////////////////////////////////////////////////////////////////////////////////

 void stetmesh::Tetmesh::delTet(uint tidx)
 {

 }

 ////////////////////////////////////////////////////////////////////////////////

 void stetmesh::Tetmesh::copyTet(uint tidx, uint * tet)
 {

 }
 */
////////////////////////////////////////////////////////////////////////////////

std::vector<uint>  stetmesh::Tetmesh::getTet(uint tidx) const
{
	if (tidx >= pTetsN)
	{
		std::ostringstream os;
		os << "Tetrahedron index out of range.";
		throw steps::ArgErr(os.str());
	}
	std::vector<uint> tet(4);
	tet[0] = pTets[tidx*4];
	tet[1] = pTets[(tidx*4)+1];
	tet[2] = pTets[(tidx*4)+2];
	tet[3] = pTets[(tidx*4)+3];
	return tet;
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tetmesh::getMeshVolume(void) const
{
    assert(pSetupDone == true);

	double totvol = 0.0;
	for (uint i=0; i<pTetsN; ++i) totvol += pTet_vols[i];
	return totvol;
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tetmesh::getTetVol(uint tidx) const
{
    assert(pSetupDone == true);

	if (tidx >= pTetsN)
	{
		std::ostringstream os;
		os << "Tetrahedron index out of range.";
		throw steps::ArgErr(os.str());
	}
	return (pTet_vols[tidx]);
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tetmesh::getTetQualityRER(uint tidx) const
{
	assert(pSetupDone == true);
	if (tidx >= pTetsN)
	{
		std::ostringstream os;
		os << "Tetrahedron index out of range.";
		throw steps::ArgErr(os.str());
	}

	uint * verts = _getTet(tidx);

    double * v0 = _getVertex(verts[0]);
    double * v1 = _getVertex(verts[1]);
    double * v2 = _getVertex(verts[2]);
    double * v3 = _getVertex(verts[3]);

    return (steps::math::tet_circumrad(v0, v1, v2, v3) /
            steps::math::tet_shortestedge(v0, v1, v2, v3));
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getTetBarycenter(uint tidx) const
{
	assert(pSetupDone == true);
	if (tidx >= pTetsN)
	{
		std::ostringstream os;
		os << "Tetrahedron index out of range.";
		throw steps::ArgErr(os.str());
	}
	std::vector<double> baryctemp(3);
	baryctemp[0] = pTet_barycentres[tidx*3];
	baryctemp[1] = pTet_barycentres[(tidx*3)+1];
	baryctemp[2] = pTet_barycentres[(tidx*3)+2];
	return baryctemp;
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmComp * stetmesh::Tetmesh::getTetComp(uint tidx) const
{
	if (tidx >= pTetsN)
	{
		std::ostringstream os;
		os << "Tetrahedron index out of range.";
		throw steps::ArgErr(os.str());
	}
	return pTet_comps[tidx];
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::setTetComp(uint tidx, stetmesh::TmComp * comp)
{
	if (tidx >= pTetsN)
	{
		std::ostringstream os;
		os << "Tetrahedron index out of range.";
		throw steps::ArgErr(os.str());
	}
	pTet_comps[tidx] = comp;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> stetmesh::Tetmesh::getTetTriNeighb(uint tidx) const
{
	assert(pSetupDone == true);
	if (tidx >= pTetsN)
	{
		std::ostringstream os;
		os << "Tetrahedron index out of range.";
		throw steps::ArgErr(os.str());
	}
	std::vector<uint> tritemp(4);
	tritemp[0] = pTet_tri_neighbours[tidx*4];
	tritemp[1] = pTet_tri_neighbours[(tidx*4)+1];
	tritemp[2] = pTet_tri_neighbours[(tidx*4)+2];
	tritemp[3] = pTet_tri_neighbours[(tidx*4)+3];
	return tritemp;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> stetmesh::Tetmesh::getTetTetNeighb(uint tidx) const
{
    assert(pSetupDone == true);
	if (tidx >= pTetsN)
	{
		std::ostringstream os;
		os << "Tetrahedron index out of range.";
		throw steps::ArgErr(os.str());
	}
	std::vector<int> tettemp(4);
	tettemp[0] = pTet_tet_neighbours[tidx*4];
	tettemp[1] = pTet_tet_neighbours[(tidx*4)+1];
	tettemp[2] = pTet_tet_neighbours[(tidx*4)+2];
	tettemp[3] = pTet_tet_neighbours[(tidx*4)+3];
	return tettemp;
}

////////////////////////////////////////////////////////////////////////////////

int stetmesh::Tetmesh::findTetByPoint(std::vector<double> p) const
{
    assert(pSetupDone == true);
    int tetidx = -1;
    // initial check to see if point is outside boundary box
    if (p[0] < pXmin || p[1] < pYmin || p[2] < pZmin
    	|| p[0] > pXmax || p[1] > pYmax || p[2] > pZmax)
    {
		return tetidx;
    }
    double pnt[3] = {p[0], p[1], p[2]};
    for (uint tidx = 0; tidx < pTetsN; ++tidx)
    {
    	double * vert0 = pVerts + (3 * pTets[tidx*4]);
    	double * vert1 = pVerts + (3 * pTets[(tidx*4)+1]);
    	double * vert2 = pVerts + (3 * pTets[(tidx*4)+2]);
    	double * vert3 = pVerts + (3 * pTets[(tidx*4)+3]);
    	if (steps::math::tet_inside(vert0, vert1, vert2, vert3, pnt))
    	{
    		tetidx = tidx;
    		break;
    	}
    }

	return tetidx;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getBoundMin(void) const
{
	assert(pSetupDone == true);

	std::vector<double> b_min(3);
	b_min[0] = pXmin;
	b_min[1] = pYmin;
	b_min[2] = pZmin;
	return b_min;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getBoundMax(void) const
{
	assert(pSetupDone == true);

	std::vector<double> b_max(3);
	b_max[0] = pXmax;
	b_max[1] = pYmax;
	b_max[2] = pZmax;
	return b_max;
}

////////////////////////////////////////////////////////////////////////////////

double * stetmesh::Tetmesh::_getVertex(uint vidx) const
{
    return pVerts + (vidx * 3);
}

////////////////////////////////////////////////////////////////////////////////

uint * stetmesh::Tetmesh::_getTri(uint tidx) const
{
    return pTris + (tidx * 3);
}

////////////////////////////////////////////////////////////////////////////////

uint * stetmesh::Tetmesh::_getTet(uint tidx) const
{
    return pTets + (tidx * 4);
}

////////////////////////////////////////////////////////////////////////////////

int * stetmesh::Tetmesh::_getTriTetNeighb(uint tidx) const
{
	return pTri_tet_neighbours + (tidx * 2);
}

////////////////////////////////////////////////////////////////////////////////

uint * stetmesh::Tetmesh::_getTetTriNeighb(uint tidx) const
{
	return pTet_tri_neighbours + (tidx * 4);
}

////////////////////////////////////////////////////////////////////////////////

int * stetmesh::Tetmesh::_getTetTetNeighb(uint tidx) const
{
	return pTet_tet_neighbours + (tidx * 4);
}

////////////////////////////////////////////////////////////////////////////////

double * stetmesh::Tetmesh::_getTriNorm(uint tidx) const
{
	return pTri_norms + (tidx * 3);
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
bool stetmesh::array_srt_cmp(T ar1[], T ar2[], uint ar_size)
{
	for(unsigned int i = 0; i < ar_size; ++i)
	{
		if (ar1[i] != ar2[i]) return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////

// END
