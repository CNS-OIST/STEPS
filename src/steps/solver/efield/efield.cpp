/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#include "efield.hpp"

#include "tetcoupler.hpp"
#include "util/checkpointing.hpp"
#include "util/error.hpp"

namespace steps::solver::efield {

// All parameters are supplied in base s.i. units. EField object converts to
// different units for matrix calculation:
// Specific membrane capacitance: microfarad per square centimetre
// Current: picoAmpere
// Electric potential: millivolt
// Coordinates of vertices are stored in microns in solver object to avoid
// copying here.

////////////////////////////////////////////////////////////////////////////////

EField::EField(std::unique_ptr<EFieldSolver> impl)
    : pVProp(std::move(impl)) {}

void EField::initMesh(const std::vector<double>& verts,
                      const std::vector<vertex_id_t>& tris,
                      const std::vector<vertex_id_t>& tets,
                      uint opt_method,
                      std::string const& opt_file_name,
                      double search_percent) {
    pNVerts = verts.size() / 3;
    pNTris = tris.size() / 3;
    pNTets = tets.size() / 4;

    // First, the mesh is constructed -- VertexElements are created
    // and triangle and tetrahedron arrays are copied
    // TODO: (copying tris and tets is maybe not necessary?).
    pMesh = new TetMesh(pNVerts, verts.data(), pNTris, tris.data(), pNTets, tets.data());
    AssertLog(pMesh != nullptr);

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
    pMesh->applyConductance(1.0 / 1.0e-3);

    // Default value for the membrane potential is -65mV but may be changed with
    // solver method setPotential
    AssertLog(static_cast<bool>(pVProp));
    pVProp->initMesh(pMesh);
    setMembPotential(solver::membrane_global_id(0), DEFAULT_MEMB_POT);

    // pTritoVert = new uint[pNTris*3];
    pTritoVert.resize(pNTris * 3);
    for (uint i = 0, j = 0; i < pNTris; ++i, j += 3) {
        for (uint v = 0; v < 3; ++v) {
            pTritoVert[j + v] = pMesh->getTriangleVertex(triangle_local_id(i), v);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

EField::~EField() {
    delete pMesh;
}

////////////////////////////////////////////////////////////////////////////////

void EField::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pNVerts);
    util::checkpoint(cp_file, pNTris);
    util::checkpoint(cp_file, pNTets);
    util::checkpoint(cp_file, pCPerm);
    pMesh->checkpoint(cp_file);
    pVProp->checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void EField::restore(std::fstream& cp_file) {
    util::compare(cp_file, pNVerts, "Mismatched EField pNVerts restore value.");
    util::compare(cp_file, pNTris, "Mismatched EField pNTris restore value.");
    util::compare(cp_file, pNTets, "Mismatched EField pNTets restore value.");
    util::restore(cp_file, pCPerm);
    pMesh->restore(cp_file);
    pVProp->restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void EField::setMembCapac(solver::membrane_global_id midx, double cm) {
    // Currently midx should be zero until multiple membranes are supported
    AssertLog(midx.get() == 0);

    AssertLog(cm >= 0.0);
    // Geometry is in microns, calculation uses pF, so we need to supply
    // specific capacitance in pF/um2.
    // Argument is in F/m^2: 1 F/m^2 = 1 pF / um^2 so no conversion needed!
    pMesh->applySurfaceCapacitance(cm);
}

void EField::setTriCapac(triangle_local_id tidx, double cm) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTris);

    AssertLog(cm >= 0.0);

    // Geometry is in microns, calculation uses pF, so we need to supply
    // specific capacitance in pF/um2.
    // Argument is in F/m^2: 1 F/m^2 = 1 pF / um^2 so no conversion needed!

    pMesh->applyTriCapacitance(tidx, cm);
}

////////////////////////////////////////////////////////////////////////////////

void EField::setMembPotential(solver::membrane_global_id midx, double v) {
    // Currently midx should be zero until multiple membranes are supported
    AssertLog(midx.get() == 0);

    // We require mV
    pVProp->setPotential(v * 1.0e3);
}

////////////////////////////////////////////////////////////////////////////////

void EField::setMembVolRes(solver::membrane_global_id /*midx*/, double ro) {
    AssertLog(ro >= 0.0);
    pMesh->applyConductance(1.0 / (ro * 1.0e-3));
}

////////////////////////////////////////////////////////////////////////////////

void EField::setSurfaceResistivity(solver::membrane_global_id /*midx*/, double rspec, double vrev) {
    // arguments in ohm.m^2 and V
    // want Gohm.um^2, so divide by 10^9 for Gohm, multiply by 10^12 for um^2

    double rs = 1.0e3 * rspec;
    double vext = vrev * 1.0e3;

    pVProp->setSurfaceConductance(1.0 / rs, vext);
}

std::pair<double, double> EField::getSurfaceResistivity() const {
    // resistivity and reversal potential are stored in Gohm.um^2 and mV
    // return them in ohm.m^2 and V

    const auto& cond = pVProp->getSurfaceConductance();
    return {1 / (cond.first * 1.0e3), cond.second / 1.0e3};
}

////////////////////////////////////////////////////////////////////////////////

double EField::getVertIClamp(vertex_id_t vidx) {
    // vidx argument converted to local index in Tetexact.
    AssertLog(vidx < pNVerts);

    // Convert the index to one that makes sense to vprop
    vidx = pCPerm[vidx.get()];

    // let VProp do any further necessary argument checking (should be less than
    // number of vertices)
    return pVProp->getVertIClamp(vidx) / 1.0e12;
}

////////////////////////////////////////////////////////////////////////////////

void EField::setVertIClamp(vertex_id_t vidx, double cur) {
    // vidx argument converted to local index in Tetexact.
    AssertLog(vidx < pNVerts);

    // Convert the index to one that makes sense to vprop
    vidx = pCPerm[vidx.get()];

    // let VProp do any further necessary argument checking (should be less than
    // number of vertices) convert current to picoamps and maintain convention
    // that positive applied current is depolarising
    pVProp->setVertIClamp(vidx, -cur * 1.0e12);
}

////////////////////////////////////////////////////////////////////////////////

void EField::advance(double dt) {
    AssertLog(dt >= 0.0);

    // Convert to ms
    pVProp->advance(dt * 1.0e3);
}

////////////////////////////////////////////////////////////////////////////////

double EField::getVertV(vertex_id_t vidx) {
    // vidx argument converted to local index in Tetexact.
    AssertLog(vidx < pNVerts);
    // let VProp do any further necessary argument checking (should be less than
    // numer of vertices)

    // Convert the index to one that makes sense to vprop
    vidx = pCPerm[vidx.get()];

    // convert from mV to V
    return pVProp->getV(vidx) * 1.0e-3;
}

////////////////////////////////////////////////////////////////////////////////

void EField::setVertV(vertex_id_t vidx, double v) {
    // vidx argument converted to local index in Tetexact.
    AssertLog(vidx < pNVerts);

    // Convert the index to one that makes sense to vprop
    vidx = pCPerm[vidx.get()];

    // convert to millivolts
    pVProp->setV(vidx, v * 1.0e3);
}

////////////////////////////////////////////////////////////////////////////////

bool EField::getVertVClamped(vertex_id_t vidx) {
    // vidx argument converted to local index in Tetexact.
    AssertLog(vidx < pNVerts);

    // Convert the index to one that makes sense to vprop
    vidx = pCPerm[vidx.get()];

    return pVProp->getClamped(vidx);
}

////////////////////////////////////////////////////////////////////////////////

void EField::setVertVClamped(vertex_id_t vidx, bool cl) {
    // vidx argument converted to local index in Tetexact.
    AssertLog(vidx < pNVerts);

    // Convert the index to one that makes sense to vprop
    vidx = pCPerm[vidx.get()];

    pVProp->setClamped(vidx, cl);
}

////////////////////////////////////////////////////////////////////////////////

double EField::getTriV(triangle_local_id tidx) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTris);

    double pot = 0.0;

    pot += pVProp->getV(pTritoVert[tidx.get() * 3]);
    pot += pVProp->getV(pTritoVert[(tidx.get() * 3) + 1]);
    pot += pVProp->getV(pTritoVert[(tidx.get() * 3) + 2]);

    // getV returns in milliVolts
    return (pot * 1.0e-3) / 3.0;
}

////////////////////////////////////////////////////////////////////////////////

void EField::setTriV(triangle_local_id tidx, double v) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTris);

    // Lets convert to millivolts
    v *= 1.0e3;

    // Directly use VProp since we'll use mesh indices
    pVProp->setV(pMesh->getTriangleVertex(tidx, 0), v);
    pVProp->setV(pMesh->getTriangleVertex(tidx, 1), v);
    pVProp->setV(pMesh->getTriangleVertex(tidx, 2), v);
}

////////////////////////////////////////////////////////////////////////////////

bool EField::getTriVClamped(triangle_local_id tidx) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTris);

    // Directly use VProp since we'll use mesh indices
    if (!pVProp->getClamped(pMesh->getTriangleVertex(tidx, 0))) {
        return false;
    }
    if (!pVProp->getClamped(pMesh->getTriangleVertex(tidx, 1))) {
        return false;
    }
    return pVProp->getClamped(pMesh->getTriangleVertex(tidx, 2));
}

////////////////////////////////////////////////////////////////////////////////

void EField::setTriVClamped(triangle_local_id tidx, bool cl) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTris);

    // Directly use VProp since we'll use mesh indices
    pVProp->setClamped(pMesh->getTriangleVertex(tidx, 0), cl);
    pVProp->setClamped(pMesh->getTriangleVertex(tidx, 1), cl);
    pVProp->setClamped(pMesh->getTriangleVertex(tidx, 2), cl);
}
////////////////////////////////////////////////////////////////////////////////

double EField::getTriI(triangle_local_id tidx) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTris);

    // convert from picoamp to amp
    return pVProp->getTriI(tidx) * 1.0e-12;
}

////////////////////////////////////////////////////////////////////////////////

void EField::setTriI(triangle_local_id tidx, double cur) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTris);

    pVProp->setTriI(tidx, cur * 1.0e12);
}

////////////////////////////////////////////////////////////////////////////////

double EField::getTriIClamp(triangle_local_id tidx) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTris);
    // maintain convention that positive applied current is depolarising
    return pVProp->getTriIClamp(tidx) / 1.0e12;
}

////////////////////////////////////////////////////////////////////////////////

void EField::setTriIClamp(triangle_local_id tidx, double cur) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTris);
    // maintain convention that positive applied current is depolarising
    pVProp->setTriIClamp(tidx, -cur * 1.0e12);
}

////////////////////////////////////////////////////////////////////////////////
/*
void    EField::setTriI(double * cur)
{
    for (uint i = 0; i < pNTris; ++i)
    {
        // convert to picoamp
        pVProp->setTriI(i, cur[i]*1.0e12);
    }
}
*/
////////////////////////////////////////////////////////////////////////////////

double EField::getTetV(tetrahedron_local_id tidx) {
    // TODO: improve on this function, just returning the mean of the vertices at
    // the moment

    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTets);

    double pot = 0.0;

    const auto v0 = pMesh->getTetrahedronVertex(tidx, 0);
    const auto v1 = pMesh->getTetrahedronVertex(tidx, 1);
    const auto v2 = pMesh->getTetrahedronVertex(tidx, 2);
    const auto v3 = pMesh->getTetrahedronVertex(tidx, 3);

    // Directly use VProp since we have mesh indices
    pot += pVProp->getV(v0);
    pot += pVProp->getV(v1);
    pot += pVProp->getV(v2);
    pot += pVProp->getV(v3);

    // getV returns mV
    return (pot * 1.0e-3) / 4.0;
}

////////////////////////////////////////////////////////////////////////////////

void EField::setTetV(tetrahedron_local_id tidx, double v) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTets);

    // Lets convert to millivolts here
    v *= 1.0e3;

    // Directly use VProp since we'll use mesh indices
    pVProp->setV(pMesh->getTetrahedronVertex(tidx, 0), v);
    pVProp->setV(pMesh->getTetrahedronVertex(tidx, 1), v);
    pVProp->setV(pMesh->getTetrahedronVertex(tidx, 2), v);
    pVProp->setV(pMesh->getTetrahedronVertex(tidx, 3), v);
}

////////////////////////////////////////////////////////////////////////////////

bool EField::getTetVClamped(tetrahedron_local_id tidx) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTets);

    // Directly use VProp since we'll use mesh indices
    if (!pVProp->getClamped(pMesh->getTetrahedronVertex(tidx, 0))) {
        return false;
    }
    if (!pVProp->getClamped(pMesh->getTetrahedronVertex(tidx, 1))) {
        return false;
    }
    if (!pVProp->getClamped(pMesh->getTetrahedronVertex(tidx, 2))) {
        return false;
    }
    return pVProp->getClamped(pMesh->getTetrahedronVertex(tidx, 3));
}

////////////////////////////////////////////////////////////////////////////////

void EField::setTetVClamped(tetrahedron_local_id tidx, bool cl) {
    // tidx argument converted to local index in Tetexact
    AssertLog(tidx < pNTets);

    pVProp->setClamped(pMesh->getTetrahedronVertex(tidx, 0), cl);
    pVProp->setClamped(pMesh->getTetrahedronVertex(tidx, 1), cl);
    pVProp->setClamped(pMesh->getTetrahedronVertex(tidx, 2), cl);
    pVProp->setClamped(pMesh->getTetrahedronVertex(tidx, 3), cl);
}

void EField::saveOptimal(std::string const& opt_file_name) {
    pMesh->saveOptimal(opt_file_name);
}

}  // namespace steps::solver::efield
