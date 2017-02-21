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
#include <iostream>
#include <cassert>
#include <sstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/efield/efield.hpp"
#include "steps/solver/efield/tetmesh.hpp"
#include "steps/solver/efield/tetcoupler.hpp"

#include "third_party/easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace sefield = steps::solver::efield;

////////////////////////////////////////////////////////////////////////////////

// All parameters are supplied in base s.i. units. EField object converts to
// different units for matrix calculation:
// Specific membrane capacitance: microfarad per square centimetre
// Current: picoAmpere
// Electric potential: millivolt
// Coordinates of vertices are stored in microns in solver object to avoid
// copying here.

////////////////////////////////////////////////////////////////////////////////

sefield::EField::EField(std::unique_ptr<EFieldSolver> impl):
    pVProp(std::move(impl)),
    pMesh(0),
    pNVerts(0), pNTris(0), pNTets(0),
    pCPerm()
{}
    
void sefield::EField::initMesh(uint nverts, double * verts,
                               uint ntris, uint * tris,
                               uint ntets, uint * tets,
                               uint opt_method,
                               std::string const & opt_file_name,
                               double search_percent)
{
    pNVerts = nverts;
    pNTris = ntris;
    pNTets = ntets;

    // First, the mesh is constructed -- VertexElements are created
    // and triangle and tetrahedron arrays are copied
    // TODO: (copying tris and tets is maybe not necessary?).
    pMesh = new sefield::TetMesh(pNVerts, verts, pNTris, tris, pNTets, tets);
    assert(pMesh != 0);

    // Extract all unique connections, by looping over vertices and
    // analyzing their edges.
    pMesh->extractConnections();

    // For each triangle, compute its surface area and add one third
    // of this value to each vertex associated with that triangle. The
    // surface area of each vertex therefore receives contributions from
    // each triangle that it is part of.
    pMesh->allocateSurface();

    // "Couple the mesh": this means that the coupling constant between
    // each vertex-vertex connection gets computed.
    TetCoupler tc(pMesh);
    tc.coupleMesh();

    pMesh->axisOrderElements(opt_method, opt_file_name, search_percent);
    pCPerm = pMesh->getVertexPermutation();

    // Geometry is in microns, calculation uses pF, so we need to supply
    // specific capacitance in pF/um2. Default 1 uF/cm^2 = 0.01 pF/um^2
    pMesh->applySurfaceCapacitance(0.01);

    // Geometry is in microns, times in ms and voltages in mV, so here we
    // need a conductivity in units of nS/um. or resistivity in Gohm_micron:
    // 1 ohm_m = 10-9 Gohm_m = 10-3 Gohm_micron.
    // Default is 100 ohm.cm  = 1 ohm.m
    pMesh->applyConductance(1.0/1.0e-3);

    // Default value for the membrane potential is -65mV but may be changed with
    // solver method setPotential
    assert(static_cast<bool>(pVProp));
    pVProp->initMesh(pMesh);
    pVProp->setPotential(-65);

    //pTritoVert = new uint[pNTris*3];
    pTritoVert.resize(pNTris*3);
    for (uint i=0, j=0; i< pNTris; ++i, j+=3) 
        for (uint v=0; v<3; ++v) 
            pTritoVert[j+v] = pMesh->getTriangleVertex(i,v);
}

////////////////////////////////////////////////////////////////////////////////

sefield::EField::~EField(void)
{
    delete pMesh;
}


////////////////////////////////////////////////////////////////////////////////

void sefield::EField::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&pNVerts, sizeof(uint));
    cp_file.write((char*)&pNTris, sizeof(uint));
    cp_file.write((char*)&pNTets, sizeof(uint));

    uint nCPerm = pCPerm.size();
    cp_file.write((char*)&nCPerm, sizeof(uint));
    cp_file.write((char*)&pCPerm.front(), sizeof(uint) * nCPerm);

    pMesh->checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&pNVerts, sizeof(uint));
    cp_file.read((char*)&pNTris, sizeof(uint));
    cp_file.read((char*)&pNTets, sizeof(uint));

    uint nCPerm;
    cp_file.read((char*)&nCPerm, sizeof(uint));
    pCPerm.resize(nCPerm);
    cp_file.read((char*)&pCPerm.front(), sizeof(uint) * nCPerm);

    pMesh->restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::setMembCapac(uint midx, double cm)
{
    // Currently midx should be zero until multiple membranes are supported
    assert (midx == 0);

    assert (cm >= 0.0);
    // Geometry is in microns, calculation uses pF, so we need to supply
    // specific capacitance in pF/um2.
    // Argument is in F/m^2: 1 F/m^2 = 1 pF / um^2 so no conversion needed!
    pMesh->applySurfaceCapacitance(cm);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::setMembPotential(uint midx, double v)
{
    // Currently midx should be zero until multiple membranes are supported
    assert (midx == 0);

    // We require mV
    pVProp->setPotential(v*1.0e3);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::setMembVolRes(uint midx, double ro)
{
    assert(ro >= 0.0);
    pMesh->applyConductance(1.0/(ro*1.0e-3));
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::setSurfaceResistivity(uint midx,double rspec, double vrev)
{
    // arguments in ohm_m^2 and V
    // want Gohm_um^2, so divide by 10^9 for Gohm, multiply by 10^12 for um^2

    double rs = 1.0e3*rspec;
    double vext = vrev*1.0e3;

    pVProp->setSurfaceConductance(1.0/rs, vext);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::setVertIClamp(uint vidx, double cur)
{
    // vidx argument converted to local index in Tetexact.
    assert(vidx < pNVerts);

    //Convert the index to one that makes sense to vprop
    vidx = pCPerm[vidx];

    // let VProp do any further necessary argument checking (should be less than number of vertices)
    // convert current to picoamps and maintain convention that positive applied current
    // is depolarising
    pVProp->setVertIClamp(vidx, -cur*1.0e12);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::advance(double dt)
{
    assert(dt >= 0.0);

    //Convert to ms
    pVProp->advance(dt*1.0e3);
}

////////////////////////////////////////////////////////////////////////////////

double sefield::EField::getVertV(uint vidx)
{
    // vidx argument converted to local index in Tetexact.
    assert(vidx < pNVerts);
    // let VProp do any further necessary argument checking (should be less than numer of vertices)

    //Convert the index to one that makes sense to vprop
    vidx = pCPerm[vidx];

    // convert from mV to V
    return pVProp->getV(vidx)*1.0e-3;
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::setVertV(uint vidx, double v)
{
    // vidx argument converted to local index in Tetexact.
    assert(vidx < pNVerts);

    //Convert the index to one that makes sense to vprop
    vidx = pCPerm[vidx];

    // convert to millivolts
    pVProp->setV(vidx, v*1.0e3);
}

////////////////////////////////////////////////////////////////////////////////

bool sefield::EField::getVertVClamped(uint vidx)
{
    // vidx argument converted to local index in Tetexact.
    assert(vidx < pNVerts);

    //Convert the index to one that makes sense to vprop
    vidx = pCPerm[vidx];

    return pVProp->getClamped(vidx);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::setVertVClamped(uint vidx, bool cl)
{
    // vidx argument converted to local index in Tetexact.
    assert(vidx < pNVerts);

    //Convert the index to one that makes sense to vprop
    vidx = pCPerm[vidx];

    pVProp->setClamped(vidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

double  sefield::EField::getTriV(uint tidx)
{
    // tidx argument converted to local index in Tetexact
    assert(tidx < pNTris);

    double pot = 0.0;

    pot += pVProp->getV(pTritoVert[tidx*3]);
    pot += pVProp->getV(pTritoVert[(tidx*3)+1]);
    pot += pVProp->getV(pTritoVert[(tidx*3)+2]);

    // getV returns in milliVolts
    return ((pot*1.0e-3)/3.0);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::setTriV(uint tidx, double v)
{
    // tidx argument converted to local index in Tetexact
    assert(tidx < pNTris);

    // Lets convert to millivolts
    v *= 1.0e3;

    // Directly use VProp since we'll use mesh indices
    pVProp->setV(pMesh->getTriangleVertex(tidx, 0), v);
    pVProp->setV(pMesh->getTriangleVertex(tidx, 1), v);
    pVProp->setV(pMesh->getTriangleVertex(tidx, 2), v);
}

////////////////////////////////////////////////////////////////////////////////

bool sefield::EField::getTriVClamped(uint tidx)
{
    // tidx argument converted to local index in Tetexact
    assert(tidx < pNTris);

    // Directly use VProp since we'll use mesh indices
    if (pVProp->getClamped(pMesh->getTriangleVertex(tidx, 0)) == false) return false;
    if (pVProp->getClamped(pMesh->getTriangleVertex(tidx, 1)) == false) return false;
    if (pVProp->getClamped(pMesh->getTriangleVertex(tidx, 2)) == false) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::setTriVClamped(uint tidx, bool cl)
{
    // tidx argument converted to local index in Tetexact
    assert(tidx < pNTris);

    // Directly use VProp since we'll use mesh indices
    pVProp->setClamped(pMesh->getTriangleVertex(tidx, 0), cl);
    pVProp->setClamped(pMesh->getTriangleVertex(tidx, 1), cl);
    pVProp->setClamped(pMesh->getTriangleVertex(tidx, 2), cl);
}
////////////////////////////////////////////////////////////////////////////////

double sefield::EField::getTriI(uint tidx)
{
    // tidx argument converted to local index in Tetexact
    assert(tidx < pNTris);

    // convert from picoamp to amp
    return (pVProp->getTriI(tidx)*1.0e-12);
}

////////////////////////////////////////////////////////////////////////////////

void    sefield::EField::setTriI(uint tidx, double cur)
{
    // tidx argument converted to local index in Tetexact
    assert(tidx < pNTris);

    pVProp->setTriI(tidx, cur*1.0e12);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::EField::setTriIClamp(uint tidx, double cur)
{
    // tidx argument converted to local index in Tetexact
    assert(tidx < pNTris);
    // maintain convention that positive applied current is depolarising
    pVProp->setTriIClamp(tidx, -cur*1.0e12);
}

////////////////////////////////////////////////////////////////////////////////
/*
void    sefield::EField::setTriI(double * cur)
{
    for (uint i = 0; i < pNTris; ++i)
    {
        // convert to picoamp
        pVProp->setTriI(i, cur[i]*1.0e12);
    }
}
*/
////////////////////////////////////////////////////////////////////////////////

double    sefield::EField::getTetV(uint tidx)
{
    // TODO: improve on this function, just returning the mean of the vertices at the moment

    // tidx argument converted to local index in Tetexact
    assert(tidx < pNTets);

    double pot = 0.0;

    uint v0 = pMesh->getTetrahedronVertex(tidx, 0);
    uint v1 = pMesh->getTetrahedronVertex(tidx, 1);
    uint v2 = pMesh->getTetrahedronVertex(tidx, 2);
    uint v3 = pMesh->getTetrahedronVertex(tidx, 3);

    // Directly use VProp since we have mesh indices
    pot += pVProp->getV(v0);
    pot += pVProp->getV(v1);
    pot += pVProp->getV(v2);
    pot += pVProp->getV(v3);

    // getV returns mV
    return ((pot*1.0e-3)/4.0);
}

////////////////////////////////////////////////////////////////////////////////

void     sefield::EField::setTetV(uint tidx, double v)
{
    // tidx argument converted to local index in Tetexact
    assert(tidx < pNTets);

    // Lets convert to millivolts here
    v *= 1.0e3;

    // Directly use VProp since we'll use mesh indices
    pVProp->setV(pMesh->getTetrahedronVertex(tidx, 0), v);
    pVProp->setV(pMesh->getTetrahedronVertex(tidx, 1), v);
    pVProp->setV(pMesh->getTetrahedronVertex(tidx, 2), v);
    pVProp->setV(pMesh->getTetrahedronVertex(tidx, 3), v);
}

////////////////////////////////////////////////////////////////////////////////

bool    sefield::EField::getTetVClamped(uint tidx)
{
    // tidx argument converted to local index in Tetexact
    assert(tidx < pNTets);

    // Directly use VProp since we'll use mesh indices
    if (pVProp->getClamped(pMesh->getTetrahedronVertex(tidx, 0)) == false) return false;
    if (pVProp->getClamped(pMesh->getTetrahedronVertex(tidx, 1)) == false) return false;
    if (pVProp->getClamped(pMesh->getTetrahedronVertex(tidx, 2)) == false) return false;
    if (pVProp->getClamped(pMesh->getTetrahedronVertex(tidx, 3)) == false) return false;

    return true;
}

////////////////////////////////////////////////////////////////////////////////

void    sefield::EField::setTetVClamped(uint tidx, bool cl)
{
    // tidx argument converted to local index in Tetexact
    assert(tidx < pNTets);

    pVProp->setClamped(pMesh->getTetrahedronVertex(tidx, 0), cl);
    pVProp->setClamped(pMesh->getTetrahedronVertex(tidx, 1), cl);
    pVProp->setClamped(pMesh->getTetrahedronVertex(tidx, 2), cl);
    pVProp->setClamped(pMesh->getTetrahedronVertex(tidx, 3), cl);

}

void sefield::EField::saveOptimal(std::string const & opt_file_name)
{
    pMesh->saveOptimal(opt_file_name);
}

////////////////////////////////////////////////////////////////////////////////

// END

