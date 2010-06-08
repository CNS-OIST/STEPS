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
#include <cmath>
#include <sstream>
#include <iostream>

// STEPS headers.
#include "../common.h"
#include "../math/tetrahedron.hpp"
#include "../math/tools.hpp"
#include "tet.hpp"
#include "tetmesh.hpp"
#include "../rng/rng.hpp"
#include "../error.hpp"

NAMESPACE_ALIAS(steps::rng, srng);
NAMESPACE_ALIAS(steps::tetmesh, stetmesh);
USING(stetmesh, Tetmesh);

//////////////////////////////////////////////////////////////////////////////// vector

stetmesh::Tet::Tet(Tetmesh * mesh, uint tidx)
: pTetmesh(mesh)
, pTidx(tidx)
, pVerts()
, pBaryc()
{
	if (pTetmesh == 0)
    {
		std::ostringstream os;
        os << "No mesh provided to Tet initializer function";
        throw steps::ArgErr(os.str());
    }

	uint * tet_temp = pTetmesh->_getTet(tidx);
	pVerts[0] = tet_temp[0];
	pVerts[1] = tet_temp[1];
	pVerts[2] = tet_temp[2];
	pVerts[3] = tet_temp[3];

	pBaryc = new double[3];
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tet::~Tet(void)
{
	delete pBaryc;
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tet::getVol(void) const
{
	assert(pTetmesh != 0);
	return (pTetmesh->getTetVol(pTidx));
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tet::getBarycenter(void) const
{
	double * v0 = pTetmesh->_getVertex(pVerts[0]);
	double * v1 = pTetmesh->_getVertex(pVerts[1]);
	double * v2 = pTetmesh->_getVertex(pVerts[2]);
	double * v3 = pTetmesh->_getVertex(pVerts[3]);
	/*double v0[3], v1[3], v2[3], v3[3];		// Defunct code- now use internal
	for (uint i=0; i < 3; ++i)					// method which returns pointer
	{
		v0[i] = v0vec[i];
		v1[i] = v1vec[i];
		v2[i] = v2vec[i];
		v3[i] = v3vec[i];
	}*/
	steps::math::tet_barycenter(v0, v1, v2, v3, pBaryc);
	std::vector<double> barycentre(3);
	barycentre[0] = pBaryc[0];
	barycentre[1] = pBaryc[1];
	barycentre[2] = pBaryc[2];
	return barycentre;
 }

////////////////////////////////////////////////////////////////////////////////
/*
double stetmesh::Tet::getQualityAR(void) const
{

}
*/
////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tet::getQualityRER(void) const
{
	//std::vector<double> pnt0 = getVertex(0);
	//std::vector<double> pnt1 = getVertex(1);
	//std::vector<double> pnt2 = getVertex(2);
	//std::vector<double> pnt3 = getVertex(3);
	//double p0[3] = {pnt0[0], pnt0[1], pnt0[2]};
	//double p1[3] = {pnt1[0], pnt1[2], pnt1[2]};
	//double p2[3] = {pnt2[0], pnt2[1], pnt2[2]};
	//double p3[3] = {pnt3[0], pnt3[1], pnt3[2]};
	//return (steps::math::tet_circumrad(p0, p1, p2, p3)/
	//		steps::math::tet_shortestedge(p0, p1, p2, p3));
    // DEBUG 5-Apr-2009
    double * v0 = pTetmesh->_getVertex(pVerts[0]);
    double * v1 = pTetmesh->_getVertex(pVerts[1]);
    double * v2 = pTetmesh->_getVertex(pVerts[2]);
    double * v3 = pTetmesh->_getVertex(pVerts[3]);
    return (steps::math::tet_circumrad(v0, v1, v2, v3) /
            steps::math::tet_shortestedge(v0, v1, v2, v3));
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmComp * stetmesh::Tet::getComp(void) const
{
	assert(pTetmesh != 0);
	return (pTetmesh->getTetComp(pTidx));
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tet stetmesh::Tet::getTet(uint i) const
{
	assert (i <= 3);
	int tetidx = pTetmesh->_getTetTetNeighb(pTidx)[i];
	assert(tetidx != -1);
	return (Tet(pTetmesh, tetidx));
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tet::getTetDist(uint i) const
{
	assert(i <= 3);
	int tetidx = getTetIdx(i);
	if (tetidx == -1)
	{
		return 0.0;
	}
	Tet * tettemp = new Tet(pTetmesh, tetidx);
	double * bary1 = _getBarycenter();
	double * bary2 = tettemp->_getBarycenter();
	double xdist = bary1[0] - bary2[0];
	double ydist = bary1[1] - bary2[1];
	double zdist = bary1[2] - bary2[2];
	delete tettemp;
	return (sqrt((xdist*xdist) + (ydist*ydist) + (zdist*zdist)));
}

////////////////////////////////////////////////////////////////////////////////

int stetmesh::Tet::getTetIdx(uint i) const
{
	assert(i <= 3);
	return (pTetmesh->_getTetTetNeighb(pTidx)[i]);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tri stetmesh::Tet::getTri(uint i) const
{
	assert (i<=3);
	uint triidx = pTetmesh->_getTetTriNeighb(pTidx)[i];
	return (Tri(pTetmesh, triidx));
}

////////////////////////////////////////////////////////////////////////////////

uint stetmesh::Tet::getTriIdx(uint i) const
{
	assert(i<=3);
	return (pTetmesh->_getTetTriNeighb(pTidx)[i]);
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tet::getTriDist(uint i) const
{
	assert(i <= 3);
	uint triidx = getTriIdx(i);
	Tri * tritemp = new Tri(pTetmesh, triidx);
	double * bary1 = _getBarycenter();
	double * bary2 = tritemp->_getBarycenter();
	double xdist = bary1[0] - bary2[0];
	double ydist = bary1[1] - bary2[1];
	double zdist = bary1[2] - bary2[2];
	delete tritemp;
	return (sqrt((xdist*xdist) + (ydist*ydist) + (zdist*zdist)));
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tet::getTriArea(uint i) const
{
	assert(i <= 3);
	uint triidx = pTetmesh->_getTetTriNeighb(pTidx)[i];
	return (pTetmesh->getTriArea(triidx));
}

////////////////////////////////////////////////////////////////////////////////

uint stetmesh::Tet::getVertexIdx(uint i) const
{
	assert(i <=3);
	return pVerts[i];
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tet::getVertex(uint i) const
{
	assert (i <= 3);
	return (pTetmesh->getVertex(getVertexIdx(i)));
}

////////////////////////////////////////////////////////////////////////////////

bool stetmesh::Tet::isInside(std::vector<double> p) const
{
	double * p0 = _getVertex(0);
	double * p1 = _getVertex(1);
	double * p2 = _getVertex(2);
	double * p3 = _getVertex(3);

	double pnt[3] = {p[0], p[1], p[2]};
	/*double p0[3] = {pnt0[0], pnt0[1], pnt0[2]};	// Defunct code
	double p1[3] = {pnt1[0], pnt1[2], pnt1[2]};
	double p2[3] = {pnt2[0], pnt2[1], pnt2[2]};
	double p3[3] = {pnt3[0], pnt3[1], pnt3[2]}; */
	return (steps::math::tet_inside(p0, p1, p2, p3, pnt));
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tet::getRanPnt(srng::RNG * r, uint n) const
{
	double * p = new double[n*3];

	double * p0 = _getVertex(0);
	double * p1 = _getVertex(1);
	double * p2 = _getVertex(2);
	double * p3 = _getVertex(3);
	/*double p0[3] = {pnt0[0], pnt0[1], pnt0[2]}; 	// Defunct code
	double p1[3] = {pnt1[0], pnt1[2], pnt1[2]};
	double p2[3] = {pnt2[0], pnt2[1], pnt2[2]};
	double p3[3] = {pnt3[0], pnt3[1], pnt3[2]}; */
	for (uint i=0; i< n; ++i)
	{
		double rn1 = r->getUnfIE();
		double rn2 = r->getUnfIE();
		double rn3 = r->getUnfIE();
		steps::math::tet_ranpnt(p0, p1, p2, p3, rn1, rn2, rn3, (p+(i*3)));
	}

	std::vector<double> pnt(n*3);
	for (uint i=0; i< (n*3); ++i)
	{
		pnt[i] = p[i];
	}
	delete[] p;

	return pnt;
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tri stetmesh::Tet::getTri0(void) const
{
	return getTri(0);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tri stetmesh::Tet::getTri1(void) const
{
	return getTri(1);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tri stetmesh::Tet::getTri2(void) const
{
	return getTri(2);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tri stetmesh::Tet::getTri3(void) const
{
	return getTri(3);
}


////////////////////////////////////////////////////////////////////////////////

double * stetmesh::Tet::_getBarycenter(void) const
{
	double * v0 = pTetmesh->_getVertex(pVerts[0]);
	double * v1 = pTetmesh->_getVertex(pVerts[1]);
	double * v2 = pTetmesh->_getVertex(pVerts[2]);
	double * v3 = pTetmesh->_getVertex(pVerts[3]);
	/*double v0[3], v1[3], v2[3], v3[3];		// Defunct code- now use internal
	for (uint i=0; i < 3; ++i)					// method which returns pointer
	{
		v0[i] = v0vec[i];
		v1[i] = v1vec[i];
		v2[i] = v2vec[i];
		v3[i] = v3vec[i];
	}*/
	steps::math::tet_barycenter(v0, v1, v2, v3, pBaryc);

	return pBaryc;
}

////////////////////////////////////////////////////////////////////////////////

double * stetmesh::Tet::_getVertex(uint i) const
{
	assert (i <= 3);
	return (pTetmesh->_getVertex(getVertexIdx(i)));
}

////////////////////////////////////////////////////////////////////////////////

// END
