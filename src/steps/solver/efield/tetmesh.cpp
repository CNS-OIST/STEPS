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

// STEPS headers.
#include "tetmesh.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"
#include "util/error.hpp"
#include "vertexconnection.hpp"
#include "vertexelement.hpp"

// logging
#include <easylogging++.h>

namespace steps::solver::efield {

static std::array<vertex_id_t, 4> sort_vertices(vertex_id_t v1,
                                                vertex_id_t v2,
                                                vertex_id_t v3,
                                                vertex_id_t v4) {
    std::array<vertex_id_t, 4> eax{{v1, v2, v3, v4}};
    std::sort(eax.begin(), eax.end());
    return eax;
}

TetStub::TetStub(vertex_id_t v1, vertex_id_t v2, vertex_id_t v3, vertex_id_t v4)
    : pSortedVerts(sort_vertices(v1, v2, v3, v4)) {}

////////////////////////////////////////////////////////////////////////////////

TetStub::TetStub(vertex_id_t* v)
    : TetStub(v[0], v[1], v[2], v[3]) {}

////////////////////////////////////////////////////////////////////////////////

bool TetStub::operator<(TetStub const& t) const {
    if (pSortedVerts[0] < t.pSortedVerts[0]) {
        return true;
    }
    if (pSortedVerts[0] > t.pSortedVerts[0]) {
        return false;
    }
    if (pSortedVerts[1] < t.pSortedVerts[1]) {
        return true;
    }
    if (pSortedVerts[1] > t.pSortedVerts[1]) {
        return false;
    }
    if (pSortedVerts[2] < t.pSortedVerts[2]) {
        return true;
    }
    if (pSortedVerts[2] > t.pSortedVerts[2]) {
        return false;
    }
    if (pSortedVerts[3] < t.pSortedVerts[3]) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

TetMesh::TetMesh(uint nv,
                 const double* vpos,
                 uint ntr,
                 const vertex_id_t* trivi,
                 uint ntet,
                 const vertex_id_t* tetvi)
    : pElements(nv)
    , pNTri(ntr)
    , pNTet(ntet)
    , pTetrahedrons(nullptr)
    , pTriangles(nullptr)

{
    // Copy vertices.
    const double* vpos2 = vpos;
    for (uint i = 0; i < nv; ++i, vpos2 += 3) {
        auto* velt = new VertexElement(i, vpos2);
        pElements[i] = velt;
    }

    // Copy triangles.
    AssertLog(pNTri > 0);
    pTriangles = new vertex_id_t[pNTri * 3];
    memcpy(pTriangles, trivi, 3 * pNTri * sizeof(vertex_id_t));

    // Copy tetrahedrons.
    AssertLog(pNTet > 0);
    pTetrahedrons = new vertex_id_t[pNTet * 4];
    memcpy(pTetrahedrons, tetvi, 4 * pNTet * sizeof(vertex_id_t));

    // Make a look up table for checking if a tetrahedron exists.
    for (uint i = 0; i < pNTet; ++i) {
        TetStub t(pTetrahedrons + (i * 4));
        // We don't look kindly on duplicate tetrahedrons.
        AssertLog(pTetLUT.find(t) == pTetLUT.end());
        pTetLUT.insert(t);
    }
}

////////////////////////////////////////////////////////////////////////////////

TetMesh::~TetMesh() {
    delete[] pTriangles;
    delete[] pTetrahedrons;

    // Free memory in maps and sets.
    // delete pTetHS;
    /*
    uint n_hms = pNbrHMS->size();
    for (uint i = 0; i < n_hms; ++i)
    {
        delete (*pNbrHMS)[i];
    }
    delete pNbrHMS;
    */

    //
    // Copied from Mesh::~Mesh.
    //
    //              |
    //              |
    //             \|/
    //              '
    for (auto const& pElement: pElements) {
        delete pElement;
    }

    for (auto const& pConnection: pConnections) {
        delete pConnection;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetMesh::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pElements.size());
    for (auto const& e: pElements) {
        e->checkpoint(cp_file);
    }
    util::checkpoint(cp_file, pConnections.size());
    for (auto const& c: pConnections) {
        c->checkpoint(cp_file);
    }
    util::checkpoint(cp_file, pTriAreas);
    util::checkpoint(cp_file, pTriPrevCapac);
    util::checkpoint(cp_file, pVertexPerm);
}

////////////////////////////////////////////////////////////////////////////////

void TetMesh::restore(std::fstream& cp_file) {
    util::compare(cp_file, pElements.size(), "Mismatched pElements restore size.");
    for (auto const& e: pElements) {
        e->restore(cp_file);
    }
    util::compare(cp_file, pConnections.size(), "Mismatched pConnections restore size.");

    for (auto const& c: pConnections) {
        c->restore(cp_file);
    }
    util::restore(cp_file, pTriAreas);
    util::restore(cp_file, pTriPrevCapac);
    util::restore(cp_file, pVertexPerm);
}

////////////////////////////////////////////////////////////////////////////////

void TetMesh::extractConnections() {
    typedef std::set<ConnStub> ConnSet;

    // Keeps track of connections already added, so we don't do them twice.
    ConnSet conns;

    // Loop over all tetrahedrons and identify all edges.
    for (uint itet = 0; itet < pNTet; ++itet) {
        // Loop over all elements
        for (uint j = 0; j < 3; ++j) {
            // And loop over all possible connections
            // I.H. 24/11/09 Changes here- I didn't really get what was intended
            // before and I think vertices not belonging to same tetrahedron were
            // being connected
            for (uint k = j + 1; k < 4; ++k) {
                AssertLog(k != j);
                VertexElement* vj = getVertex(pTetrahedrons[(itet * 4) + j]);
                VertexElement* vk = getVertex(pTetrahedrons[(itet * 4) + k]);
                conns.insert(ConnStub(vj, vk));
            }
        }
    }
    // Build connection objects for all edges.
    for (auto const& conn: conns) {
        newConnection(conn.fVertex1, conn.fVertex2);
    }

    // Do some additional set up stuff.
    // Why here?
    auto nelems = pElements.size();
    pVertexPerm.resize(nelems);
    for (auto i = 0u; i < nelems; ++i) {
        // This allows the vertex to set up some additional data structures
        // depending on the number of other vertices it is connected to
        // through an edge.
        pElements[i]->fix();
        pVertexPerm[i] = vertex_id_t(i);
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetMesh::allocateSurface() {
    AssertLog(pTriangles != nullptr);

    pTriAreas.resize(pNTri);
    pTriPrevCapac.resize(pNTri);
    auto* itr = pTriangles;
    for (uint iitr = 0; iitr < pNTri; ++iitr, itr += 3) {
        VertexElement* v0 = getVertex(itr[0]);
        VertexElement* v1 = getVertex(itr[1]);
        VertexElement* v2 = getVertex(itr[2]);

        double a0 = v1->getX() - v0->getX();
        double a1 = v1->getY() - v0->getY();
        double a2 = v1->getZ() - v0->getZ();

        double b0 = v2->getX() - v0->getX();
        double b1 = v2->getY() - v0->getY();
        double b2 = v2->getZ() - v0->getZ();

        double cx = a1 * b2 - a2 * b1;
        double cy = a2 * b0 - a0 * b2;
        double cz = a0 * b1 - a1 * b0;

        double area = 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
        double area3 = area / 3.0;
        pTriAreas[iitr] = area;
        pTriPrevCapac[iitr] = 0;

        v0->incrementSurfaceArea(area3);
        v1->incrementSurfaceArea(area3);
        v2->incrementSurfaceArea(area3);
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetMesh::reindexElements() {
    uint i = 0;
    for (auto const& e: pElements) {
        e->setIDX(i++);
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::array<uint, 3>> TetMesh::getNeighboringTetrahedra(VertexElement* ve) const {
    std::vector<std::array<uint, 3>> outputvec;
    // Check input.
    AssertLog(ve != nullptr);

    vertex_id_t v0x(ve->getIDX());
    uint nn = ve->getNCon();
    for (uint i = 0; i < nn; ++i) {
        for (uint j = i + 1; j < nn; ++j) {
            for (uint k = j + 1; k < nn; ++k) {
                vertex_id_t v1x(ve->nbrIdx(i));
                vertex_id_t v2x(ve->nbrIdx(j));
                vertex_id_t v3x(ve->nbrIdx(k));
                if (pTetLUT.find(TetStub(v0x, v1x, v2x, v3x)) != pTetLUT.end()) {
                    outputvec.push_back(std::array<uint, 3>{{i, j, k}});
                }
            }
        }
    }

    return outputvec;
}

////////////////////////////////////////////////////////////////////////////////

// Applies a surface capacitance to a built mesh. This sets the capacitance for
// each vertex from its associated area. Assuming the dimension units are
// microns, then the capacitance should be supplied in pF per square micron. So,
// for example, to set a capacitance of 1 uF/cm2, you should call this with
// argument 0.01
//
void TetMesh::applySurfaceCapacitance(double d) {
    for (auto& pElement: pElements) {
        pElement->applySurfaceCapacitance(d);
    }
    for (auto& capac: pTriPrevCapac) {
        capac = d;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetMesh::applyTriCapacitance(triangle_local_id tidx, double cm) {
    AssertLog(tidx < pNTri);

    double cap = (cm - pTriPrevCapac[tidx.get()]) * pTriAreas[tidx.get()] / 3.0;
    pTriPrevCapac[tidx.get()] = cm;

    pElements[getTriangleVertex(tidx, 0).get()]->updateCapacitance(cap);
    pElements[getTriangleVertex(tidx, 1).get()]->updateCapacitance(cap);
    pElements[getTriangleVertex(tidx, 2).get()]->updateCapacitance(cap);
}

////////////////////////////////////////////////////////////////////////////////

void TetMesh::applyConductance(double d) {
    for (auto& pElement: pElements) {
        pElement->applyConductance(d);
    }
}

////////////////////////////////////////////////////////////////////////////////

double TetMesh::getTotalCapacitance() {
    double ret = 0.;
    for (auto& pElement: pElements) {
        ret += pElement->getCapacitance();
    }
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

double TetMesh::getTotalArea() {
    double ret = 0.;
    for (auto& pElement: pElements) {
        ret += pElement->getSurfaceArea();
    }
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

std::array<double, 3> TetMesh::centerOfMass() {
    std::array<double, 3> ret{};
    const auto nelt = pElements.size();

    double f = 1. / nelt;
    for (auto i = 0u; i < nelt; i++) {
        ret[0] += f * pElements[i]->getX();
        ret[1] += f * pElements[i]->getY();
        ret[2] += f * pElements[i]->getZ();
    }
    //      E.info("cof g " + ret[0] + " " + ret[1] + " " + ret[2]);
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

void TetMesh::axisOrderElements(uint opt_method,
                                std::string const& opt_file_name,
                                double search_percent) {
    // Now this method provides a choice between Stefan and Robert's method
    // and the new method by Iain. The original method is fast and suffices for
    // simple geometries, Iain's method is superior and important for complex
    // geometries, but slow.

    if (!opt_file_name.empty()) {
        std::fstream opt_file;

        opt_file.open(opt_file_name.c_str(), std::fstream::in | std::fstream::binary);
        opt_file.seekg(0);

        uint nelems = 0;
        opt_file.read(reinterpret_cast<char*>(&nelems), sizeof(uint));
        if (pElements.size() != nelems) {
            std::ostringstream os;
            os << "optimal data mismatch with simulator parameters: "
                  "Tetmesh::nelems, ";
            os << nelems << ":" << pElements.size();
            ArgErrLog(os.str());
        }

        util::restore(opt_file, nelems, pVertexPerm);

        VertexElementPVec elements_temp = pElements;

        for (uint vidx = 0; vidx < nelems; ++vidx) {
            VertexElementP vep = elements_temp[vidx];
            // sanity check
            AssertLog(vep->getIDX() == vidx);
            const auto new_idx = pVertexPerm[vidx];
            pElements[new_idx.get()] = vep;
        }

        reindexElements();
        reordered();
        opt_file.close();
        return;
    }

    if (opt_method == 2) {
        // The breadth first search with Cuthill-McKee improvement

        auto pNVerts = pElements.size();

        // the best vertex(narrowest width)
        uint bestone = 0;
        auto bestwidth = static_cast<int>(pNVerts);

        std::vector<VertexElement*> orig_indices = pElements;

        std::stringstream ss;
        ss << "\n\n-- " << search_percent << "%  Samples of breadth first search --" << std::endl;
        CLOG(INFO, "general_log") << ss.str() << std::endl;

        // time_t btime;
        // time_t etime;
        // time(&btime);

        for (auto vidx = 0u; vidx < pNVerts; vidx += static_cast<uint>(100 / search_percent)) {
            std::set<VertexElement*> verteleset;
            std::vector<VertexElement*> vertelevec;
            std::queue<VertexElement*> vertelqueue;

            verteleset.insert(orig_indices[vidx]);
            vertelevec.push_back(orig_indices[vidx]);
            uint ve0ncons = orig_indices[vidx]->getNCon();
            VertexElement** ve0neighbs = orig_indices[vidx]->getNeighbours();

            fill_ve_vec(verteleset, vertelevec, vertelqueue, ve0ncons, ve0neighbs);
            pElements.clear();

            uint ielt = 0;

            for (auto const& vertele: vertelevec) {
                pElements.push_back(vertele);
                pVertexPerm[vertele->getIDX()] = vertex_id_t(ielt);
                ielt++;
            }

            // Note: this will reorder the vertex indices as we go- not a problem in
            // this loop but they must be reset for the final index setting to work
            reindexElements();

            int maxdi = 0;
            for (auto iv: vertex_id_t::range(pNVerts)) {
                VertexElement* ve = getVertex(iv);
                int ind = ve->getIDX();

                for (auto i = 0u; i < ve->getNCon(); ++i) {
                    auto inbr = static_cast<int>(ve->nbrIdx(i));
                    auto di = ind - inbr;
                    if (di < 0) {
                        di = -di;
                    }
                    if (di > maxdi) {
                        maxdi = di;
                    }
                }
            }
            if (maxdi < bestwidth) {
                bestone = vidx;
                bestwidth = maxdi;
            }
        }

        // reset the vertex indices to their original value: important
        for (uint v = 0; v < pNVerts; ++v) {
            orig_indices[v]->setIDX(v);
        }

        std::set<VertexElement*> verteleset;
        std::vector<VertexElement*> vertelevec;
        std::queue<VertexElement*> vertelqueue;

        verteleset.insert(orig_indices[bestone]);
        vertelevec.push_back(orig_indices[bestone]);
        uint ve0ncons = orig_indices[bestone]->getNCon();
        VertexElement** ve0neighbs = orig_indices[bestone]->getNeighbours();

        fill_ve_vec(verteleset, vertelevec, vertelqueue, ve0ncons, ve0neighbs);
        pElements.clear();

        uint ielt = 0;
        for (auto const& vertele: vertelevec) {
            pElements.push_back(vertele);
            pVertexPerm[vertele->getIDX()] = vertex_id_t(ielt);
            ielt++;
        }

    }
    // / / / / / / / / / / / /  / / / / / / / / / / / / / / / / / / / / / / //

    else if (opt_method == 1) {
        // time_t btime;
        // time_t etime;
        // time(&btime);

        auto const& com = centerOfMass();
        auto** m = new double*[3];
        for (uint i = 0; i < 3; i++) {
            m[i] = new double[3];
            m[i][0] = 0.0;
            m[i][1] = 0.0;
            m[i][2] = 0.0;
        }

        for (auto ivert: vertex_id_t::range(countVertices())) {
            VertexElement* ve = getVertex(ivert);
            double x = ve->getX() - com[0];
            double y = ve->getY() - com[1];
            double z = ve->getZ() - com[2];

            m[0][0] += x * x;
            m[0][1] += x * y;
            m[0][2] += x * z;

            m[1][0] += y * x;
            m[1][1] += y * y;
            m[1][2] += y * z;

            m[2][0] += z * x;
            m[2][1] += z * y;
            m[2][2] += z * z;
        }

        uint iret;
        double gamma;
        double evec[3];
        mainEvec(3, m, &iret, &gamma, evec);
        if (iret != 1) {
            CLOG(INFO, "general_log") << "\nWarning - eigenvalue faliure " << std::endl;
        }

        for (uint i = 0; i < 3; i++) {
            delete[] m[i];
        }
        delete[] m;

        double vx = evec[0];
        double vy = evec[1];
        double vz = evec[2];
        vx += 1.e-8;
        vy += 1.e-9;
        vz += 1.e-10;

        std::map<double, VertexElement*> hm;
        // vector<double> da;
        for (auto ielt: vertex_id_t::range(countVertices())) {
            VertexElement* ve = getVertex(ielt);
            double d = vx * ve->getX() + vy * ve->getY() + vz * ve->getZ();

            hm[d] = ve;
            // da.push_back(d);
        }

        // sort(da.begin(), da.end());

        pElements.clear();
        uint ielt = 0;
        auto iter = hm.begin();

        while (iter != hm.end()) {
            double d = iter->first;
            pElements.push_back(hm[d]);
            ++iter;

            // pVertexPerm[i] contains the new index in the elements array
            // for the vertex that was originally at index i
            pVertexPerm[hm[d]->getIDX()] = vertex_id_t(ielt);
            ielt++;
        }
    } else {
        std::ostringstream os;
        os << "Unknown optimization method.\n";
        ArgErrLog(os.str());
    }

    reindexElements();
    reordered();
}

///////////////////////////////////////////////////////////////////////////////

void TetMesh::saveOptimal(std::string const& opt_file_name) {
    std::fstream opt_file;

    opt_file.open(opt_file_name.c_str(),
                  std::fstream::out | std::fstream::binary | std::fstream::trunc);
    uint nelems = static_cast<uint>(pElements.size());
    opt_file.write(reinterpret_cast<char*>(&nelems), sizeof(uint));

    util::checkpoint(opt_file, pVertexPerm, false /* with_size */);
    opt_file.close();
}

////////////////////////////////////////////////////////////////////////////////

void TetMesh::fill_ve_vec(std::set<VertexElement*>& veset,
                          std::vector<VertexElement*>& vevec,
                          std::queue<VertexElement*>& vequeue,
                          uint ncons,
                          VertexElement** nbrs) {
    while (true) {
        std::multimap<uint, uint> neighb_connections;
        for (uint i = 0; i < ncons; ++i) {
            bool inserted = veset.insert(nbrs[i]).second;
            if (inserted) {
                neighb_connections.emplace(nbrs[i]->getNCon(), i);
            }
        }

        const auto it_end = neighb_connections.end();
        for (auto it = neighb_connections.begin(); it != it_end; ++it) {
            vevec.push_back(nbrs[it->second]);
            vequeue.push(nbrs[it->second]);
        }

        if (vequeue.empty()) {
            return;
        } else {
            VertexElement* nextve = vequeue.front();
            vequeue.pop();

            nbrs = nextve->getNeighbours();
            ncons = nextve->getNCon();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetMesh::mainEvec(uint n, double** A, uint* it, double* gamma, double* x1) {
    double eps = 1.e-10;
    double dta = 1.e-10;
    int maxit = 100;

    std::vector<double> x0(n);
    for (auto i = 0u; i < n; i++) {
        x0[i] = 1.0 / std::sqrt(i + 1);
    }
    *it = UINT_MAX;
    int l = 1;
    while (*it == UINT_MAX && l < maxit) {
        *gamma = 0.;
        for (auto i = 0u; i < n; i++) {
            x1[i] = 0.;
            for (auto j = 0u; j < n; j++) {
                x1[i] += A[i][j] * x0[j];
            }

            if (std::abs(x1[i]) > std::abs(*gamma)) {
                *gamma = x1[i];
            }
        }
        if (std::abs(*gamma) < eps) {
            *it = 0;
        } else {
            for (auto i = 0u; i < n; i++) {
                x1[i] /= *gamma;
            }
            double phi = 0.0;
            for (auto i = 0u; i < n; i++) {
                double s = std::abs(x1[i] - x0[i]);
                if (s > phi) {
                    phi = s;
                }
            }
            if (phi < dta) {
                *it = 1;
            } else {
                for (auto i = 0u; i < n; i++) {
                    x0[i] = x1[i];
                }
                l++;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetMesh::reordered() {
    auto* tris = pTriangles;
    for (uint i = 0; i < pNTri; ++i, tris += 3) {
        tris[0] = pVertexPerm[tris[0].get()];
        tris[1] = pVertexPerm[tris[1].get()];
        tris[2] = pVertexPerm[tris[2].get()];
    }

    auto* tets = pTetrahedrons;
    for (uint i = 0; i < pNTet; ++i, tets += 4) {
        tets[0] = pVertexPerm[tets[0].get()];
        tets[1] = pVertexPerm[tets[1].get()];
        tets[2] = pVertexPerm[tets[2].get()];
        tets[3] = pVertexPerm[tets[3].get()];
    }
}

////////////////////////////////////////////////////////////////////////////////

VertexConnection* TetMesh::newConnection(VertexElement* v1, VertexElement* v2) {
    auto con = new VertexConnection(v1, v2);
    pConnections.push_back(con);
    return con;
}

}  // namespace steps::solver::efield
