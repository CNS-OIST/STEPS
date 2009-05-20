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
#include <cassert>
#include <cmath>
#include <sstream>
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/tetrahedron.hpp>
#include <steps/math/tools.hpp>
#include <steps/geom/tet.hpp>
#include <steps/geom/tetmesh.hpp>
#include <steps/rng/rng.hpp>
#include <steps/error.hpp>

NAMESPACE_ALIAS(steps::rng, srng);
NAMESPACE_ALIAS(steps::tetmesh, stetmesh);
USING(stetmesh, Tetmesh);

//////////////////////////////////////////////////////////////////////////////// vector

stetmesh::Tet::Tet(Tetmesh * mesh, uint tidx)
: pTetmesh(mesh)
, pTidx(tidx)
, pVerts()
{
	if (pTetmesh == 0)
    {
		std::ostringstream os;
        os << "No mesh provided to Tet initializer function";
        throw steps::ArgErr(os.str());
    }

	std::vector<uint> tet_temp = pTetmesh->getTet(tidx);
	pVerts[0] = tet_temp[0];
	pVerts[1] = tet_temp[1];
	pVerts[2] = tet_temp[2];
	pVerts[3] = tet_temp[3];
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tet::~Tet(void)
{

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
	std::vector<double> v0vec = pTetmesh->getVertex(pVerts[0]);
	std::vector<double> v1vec = pTetmesh->getVertex(pVerts[1]);
	std::vector<double> v2vec = pTetmesh->getVertex(pVerts[2]);
	std::vector<double> v3vec = pTetmesh->getVertex(pVerts[3]);
	double v0[3], v1[3], v2[3], v3[3];
	for (uint i=0; i < 3; ++i)
	{
		v0[i] = v0vec[i];
		v1[i] = v1vec[i];
		v2[i] = v2vec[i];
		v3[i] = v3vec[i];
	}
	double baryc[3];
	steps::math::tet_barycenter(v0, v1, v2, v3, baryc);
	std::vector<double> barycentre(3);
	barycentre[0] = baryc[0];
	barycentre[1] = baryc[1];
	barycentre[2] = baryc[2];
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
	int tetidx = pTetmesh->getTetTetNeighb(pTidx)[i];
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
	std::vector<double> bary1 = getBarycenter();
	std::vector<double> bary2 = tettemp->getBarycenter();
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
	return (pTetmesh->getTetTetNeighb(pTidx)[i]);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tri stetmesh::Tet::getTri(uint i) const
{
	assert (i<=3);
	uint triidx = pTetmesh->getTetTriNeighb(pTidx)[i];
	return (Tri(pTetmesh, triidx));
}

////////////////////////////////////////////////////////////////////////////////

uint stetmesh::Tet::getTriIdx(uint i) const
{
	assert(i<=3);
	return (pTetmesh->getTetTriNeighb(pTidx)[i]);
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tet::getTriDist(uint i) const
{
	assert(i <= 3);
	uint triidx = getTriIdx(i);
	Tri * tritemp = new Tri(pTetmesh, triidx);
	std::vector<double> bary1 = getBarycenter();
	std::vector<double> bary2 = tritemp->getBarycenter();
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
	uint triidx = pTetmesh->getTetTriNeighb(pTidx)[i];
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

bool stetmesh::Tet::isInside(double * p) const
{
	std::vector<double> pnt0 = getVertex(0);
	std::vector<double> pnt1 = getVertex(1);
	std::vector<double> pnt2 = getVertex(2);
	std::vector<double> pnt3 = getVertex(3);
	double p0[3] = {pnt0[0], pnt0[1], pnt0[2]};
	double p1[3] = {pnt1[0], pnt1[2], pnt1[2]};
	double p2[3] = {pnt2[0], pnt2[1], pnt2[2]};
	double p3[3] = {pnt3[0], pnt3[1], pnt3[2]};
	return (steps::math::tet_inside(p0, p1, p2, p3, p));
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tet::getRanPnt(srng::RNG * r, uint n) const
{
	double * p = new double[n*3];

	std::vector<double> pnt0 = getVertex(0);
	std::vector<double> pnt1 = getVertex(1);
	std::vector<double> pnt2 = getVertex(2);
	std::vector<double> pnt3 = getVertex(3);
	double p0[3] = {pnt0[0], pnt0[1], pnt0[2]};
	double p1[3] = {pnt1[0], pnt1[2], pnt1[2]};
	double p2[3] = {pnt2[0], pnt2[1], pnt2[2]};
	double p3[3] = {pnt3[0], pnt3[1], pnt3[2]};
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

// END
