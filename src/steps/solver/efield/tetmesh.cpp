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
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <cassert>
#include <queue>
#include <fstream>
#include <time.h>       /* time_t, struct tm, difftime, time, mktime */

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/efield/tetmesh.hpp"
#include "steps/solver/efield/vertexconnection.hpp"
#include "steps/solver/efield/vertexelement.hpp"


////////////////////////////////////////////////////////////////////////////////

namespace sefield = steps::solver::efield;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

sefield::TetStub::TetStub(uint v1, uint v2, uint v3, uint v4)
{
    pSortedVerts[0] = v1;
    pSortedVerts[1] = v2;
    pSortedVerts[2] = v3;
    pSortedVerts[3] = v4;
    sort(pSortedVerts, pSortedVerts + 4);
}

////////////////////////////////////////////////////////////////////////////////

sefield::TetStub::TetStub(uint * v)
{
    pSortedVerts[0] = v[0];
    pSortedVerts[1] = v[1];
    pSortedVerts[2] = v[2];
    pSortedVerts[3] = v[3];
    sort(pSortedVerts, pSortedVerts + 4);
}

////////////////////////////////////////////////////////////////////////////////

bool sefield::TetStub::operator< (sefield::TetStub const & t) const
{
    if (pSortedVerts[0] < t.pSortedVerts[0]) return true;
    if (pSortedVerts[0] > t.pSortedVerts[0]) return false;
    if (pSortedVerts[1] < t.pSortedVerts[1]) return true;
    if (pSortedVerts[1] > t.pSortedVerts[1]) return false;
    if (pSortedVerts[2] < t.pSortedVerts[2]) return true;
    if (pSortedVerts[2] > t.pSortedVerts[2]) return false;
    if (pSortedVerts[3] < t.pSortedVerts[3]) return true;
    //if (pSortedVerts[3] > t.pSortedVerts[3]) return false;
    return false;
}

////////////////////////////////////////////////////////////////////////////////
/*
ostream & operator<< (ostream & os, sefield::TetStub const & ts)
{
    os << "(";
    os << ts.pSortedVerts[0] << ", ";
    os << ts.pSortedVerts[1] << ", ";
    os << ts.pSortedVerts[2] << ", ";
    os << ts.pSortedVerts[3] << ")";
    return os;
}
*/
////////////////////////////////////////////////////////////////////////////////

sefield::TetMesh::TetMesh
(
    uint nv, double * vpos,
    uint ntr, uint * trivi,
    uint ntet, uint * tetvi
)
: pElements(nv)
, pConnections()
, pVertexPerm(0)
, pNTri(ntr)
, pNTet(ntet)
, pTetrahedrons(0)
, pTriangles(0)
, pTetLUT()
//, pTetHS(0)
//, pNbrHMS(0)
{
    // Copy vertices.
    double * vpos2 = vpos;
    for (uint i = 0; i < nv; ++i, vpos2 += 3)
    {
        VertexElement * velt = new VertexElement(i, vpos2);
        pElements[i] = velt;
    }

    // Copy triangles.
    assert(pNTri > 0);
    pTriangles = new uint[pNTri * 3];
    memcpy(pTriangles, trivi, 3 * pNTri * sizeof(uint));

    // Copy tetrahedrons.
    assert(pNTet > 0);
    pTetrahedrons = new uint[pNTet * 4];
    memcpy(pTetrahedrons, tetvi, 4 * pNTet * sizeof(uint));

    // Make a look up table for checking if a tetrahedron exists.
    for (uint i = 0; i < pNTet; ++i)
    {
        TetStub t(pTetrahedrons + (i * 4));
        // We don't look kindly on duplicate tetrahedrons.
        assert(pTetLUT.find(t) == pTetLUT.end());
        pTetLUT.insert(t);
    }

}

////////////////////////////////////////////////////////////////////////////////

sefield::TetMesh::~TetMesh(void)
{
    delete[] pTriangles;
    delete[] pTetrahedrons;

    // Free memory in maps and sets.
    //delete pTetHS;
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

    uint nelems = pElements.size();
    for (unsigned i = 0; i < nelems; ++i)
    {
        delete pElements[i];
    }

    uint nconns = pConnections.size();
    for (unsigned i = 0; i < nconns; ++i)
    {
        delete pConnections[i];
    }

    delete[] pVertexPerm;
}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetMesh::checkpoint(std::fstream & cp_file)
{
    uint nelems = pElements.size();
    cp_file.write((char*)&nelems, sizeof(uint));
    for (uint e = 0; e < nelems; e++) {
        pElements[e]->checkpoint(cp_file);
    }
    uint nconns = pConnections.size();
    cp_file.write((char*)&nconns, sizeof(uint));
    for (uint c = 0; c < nconns; c++) {
        pConnections[c]->checkpoint(cp_file);
    }
    cp_file.write((char*)pVertexPerm, sizeof(uint) * nelems);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetMesh::restore(std::fstream & cp_file)
{
    uint nelems = 0;
    cp_file.read((char*)&nelems, sizeof(uint));
    if (pElements.size() != nelems) {
        std::ostringstream os;
        os << "checkpoint data mismatch with simulator parameters: sefield::Tetmesh::nelems, ";
        os << nelems << ":" << pElements.size();
        throw steps::ArgErr(os.str());
    }
    for (uint e = 0; e < nelems; e++) {
        pElements[e]->restore(cp_file);
    }

    uint nconns = 0;
    cp_file.read((char*)&nconns, sizeof(uint));
    if (pConnections.size() != nconns) {
        std::ostringstream os;
        os << "checkpoint data mismatch with simulator parameters: sefield::Tetmesh::nconns.";
        os << nconns << ":" << pConnections.size();
        throw steps::ArgErr(os.str());
    }
    for (uint c = 0; c < nconns; c++) {
        pConnections[c]->restore(cp_file);
    }
    cp_file.read((char*)pVertexPerm, sizeof(uint) * nelems);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetMesh::extractConnections(void)
{
    typedef set<ConnStub> ConnSet;
    typedef ConnSet::iterator ConnSetI;
    typedef ConnSet::const_iterator ConnSetCI;


    // Keeps track of connections already added, so we don't do them twice.
    ConnSet conns;

    // Loop over all tetrahedrons and identify all edges.
    for (uint itet = 0; itet < pNTet; ++itet)
    {
        // Loop over all elements
        for (uint j = 0; j < 3; ++j)
        {
            // And loop over all possible connections
            // I.H. 24/11/09 Changes here- I didn't really get what was intended before
            // and I think vertices not belonging to same tetrahedron were being connected
            for (uint k = j+1; k < 4; ++k)
            {
                assert (k != j);
                VertexElement * vj = getVertex(pTetrahedrons[(itet * 4) + j]);
                VertexElement * vk = getVertex(pTetrahedrons[(itet * 4) + k]);
                conns.insert(ConnStub(vj, vk));
            }
        }

    }
    // Build connection objects for all edges.
    ConnSetI conns_end = conns.end();

    for (ConnSetI conn = conns.begin(); conn != conns_end; ++conn)
    {
        newConnection(conn->fVertex1, conn->fVertex2);
    }


    // Do some additional set up stuff.
    // Why here?
    uint nelems = pElements.size();
    pVertexPerm = new uint[nelems];
    for (uint i = 0; i < nelems; ++i)
    {
        // This allows the vertex to set up some additional data structures
        // depending on the number of other vertices it is connected to
        // through an edge.
        pElements[i]->fix();
        pVertexPerm[i] = i;
    }
}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetMesh::allocateSurface(void)
{
    assert(pTriangles != 0);
    /*
    double * areas = new double[pElements.size()];

    std::cout << "\nbefore...\n[";
    for (uint i=0; i< pElements.size(); ++i)
    {
        std::cout << areas[i] << ",";
    }
    std::cout << "]";
    std::cin.get();
    */
    uint * itr = pTriangles;
    for (uint iitr = 0; iitr < pNTri; ++iitr, itr += 3)
    {
        VertexElement * v0 = getVertex(itr[0]);
        VertexElement * v1 = getVertex(itr[1]);
        VertexElement * v2 = getVertex(itr[2]);

        double a0 = v1->getX() - v0->getX();
        double a1 = v1->getY() - v0->getY();
        double a2 = v1->getZ() - v0->getZ();

        double b0 = v2->getX() - v0->getX();
        double b1 = v2->getY() - v0->getY();
        double b2 = v2->getZ() - v0->getZ();

        double cx = a1 * b2 - a2 * b1;
        double cy = a2 * b0 - a0 * b2;
        double cz = a0 * b1 - a1 * b0;

        double area = 0.5 * sqrt(cx * cx + cy * cy + cz * cz);
        double area3 = area / 3.0;

        v0->incrementSurfaceArea(area3);
        v1->incrementSurfaceArea(area3);
        v2->incrementSurfaceArea(area3);
        /*
        areas[v0->getIDX()]+=area3;
        areas[v1->getIDX()]+=area3;
        areas[v2->getIDX()]+=area3;
         */
    }
    /*
    std::cout << "\nafter...\n[";
    for (uint i=0; i< pElements.size(); ++i)
    {
        std::cout << areas[i] << ",";
    }
    std::cout << "]";
    std::cin.get();
     delete[] areas;
     */
}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetMesh::reindexElements(void)
{
    uint i = 0;
    VertexElementPVecI e_end = pElements.end();
    for (VertexElementPVecI e = pElements.begin(); e != e_end; ++e, ++i)
    {
        (*e)->setIDX(i);
    }
}

////////////////////////////////////////////////////////////////////////////////

vector<vector<uint> > sefield::TetMesh::getNeighboringTetrahedra(VertexElement * ve)
{
    // Check input.
    assert(ve != 0);

    uint v0x = ve->getIDX();
    uint nn = ve->getNCon();
    vector<vector<uint> >  outputvec = vector<vector<uint> >();
    for (uint i = 0; i < nn; ++i)
    {
        for (uint j = i + 1; j < nn; ++j)
        {
            for (uint k = j + 1; k < nn; ++k)
            {
                uint v1x = ve->nbrIdx(i);
                uint v2x = ve->nbrIdx(j);
                uint v3x = ve->nbrIdx(k);
                if (pTetLUT.find(TetStub(v0x, v1x, v2x, v3x)) != pTetLUT.end())
                {
                    std::vector<uint> tri = std::vector<uint>(3);
                    tri[0] = i;
                    tri[1] = j;
                    tri[2] = k;
                    outputvec.push_back(tri);
                }
            }
        }
    }

    return outputvec;
}

////////////////////////////////////////////////////////////////////////////////

// Applies a surface capacitance to a built mesh. This sets the capacitance for
// each vertex from its associated area. Assuming the dimension units are microns,
// then the capacitance should be supplied in pF per square micron. So, for example,
// to set a capacitance of 1 uF/cm2, you should call this with argument 0.01
//
void sefield::TetMesh::applySurfaceCapacitance(double d)
{
    for (unsigned i = 0; i < pElements.size(); i++)
    {
        pElements[i]->applySurfaceCapacitance(d);
    }
}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetMesh::applyConductance(double d)
{
    for (unsigned i = 0; i < pElements.size(); i++)
    {
        pElements[i]->applyConductance(d);
    }
}

////////////////////////////////////////////////////////////////////////////////

double sefield::TetMesh::getTotalCapacitance(void)
{
    double ret = 0.;
    for (unsigned i = 0; i < pElements.size(); i++)
    {
        ret += pElements[i]->getCapacitance();
    }
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

double sefield::TetMesh::getTotalArea(void)
{
    double ret = 0.;
    for (unsigned i = 0; i < pElements.size(); i++)
    {
        ret += pElements[i]->getSurfaceArea();
    }
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> sefield::TetMesh::centerOfMass(void)
{
    int nelt = pElements.size();
    std::vector<double> ret = std::vector<double>(3);
    ret[0] = 0.0;
    ret[1] = 0.0;
    ret[2] = 0.0;

    double f = 1. / nelt;
    for (int i = 0; i < nelt; i++)
    {
        ret[0] += f * pElements[i]->getX();
        ret[1] += f * pElements[i]->getY();
        ret[2] += f * pElements[i]->getZ();
    }
    //      E.info("cof g " + ret[0] + " " + ret[1] + " " + ret[2]);
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetMesh::axisOrderElements(uint opt_method, std::string const & opt_file_name, double search_percent)
{

    // Now this method provides a choice between Stefan and Robert's method
    // and the new method by Iain. The original method is fast and suffices for
    // simple geometries, Iain's method is superior and important for complex
    // geometries, but slow.

    if (opt_file_name != "")
    {
        std::fstream opt_file;

        opt_file.open(opt_file_name.c_str(),
                    std::fstream::in | std::fstream::binary);
        opt_file.seekg(0);

        uint nelems = 0;
        opt_file.read((char*)&nelems, sizeof(uint));
        if (pElements.size() != nelems) {
            std::ostringstream os;
            os << "optimal data mismatch with simulator parameters: sefield::Tetmesh::nelems, ";
            os << nelems << ":" << pElements.size();
            throw steps::ArgErr(os.str());
        }
        opt_file.read((char*)pVertexPerm, sizeof(uint) * nelems);

        VertexElementPVec elements_temp = pElements;

        for (uint vidx = 0; vidx < nelems; ++vidx)
        {
            VertexElementP vep = elements_temp[vidx];
            // sanity check
            assert(vep->getIDX() == vidx);
            uint new_idx = pVertexPerm[vidx];
            pElements[new_idx] = vep;
        }

        reindexElements();
        reordered();
        opt_file.close();
        return;

    }


    if (opt_method == 2)
    {
        // The breadth first search with Cuthill-McKee improvement

        uint pNVerts = pElements.size();


        // the best vertex(narrowest width)
        uint bestone = 0;
        uint bestwidth = pNVerts;

        std::vector<VertexElement*> orig_indices = pElements;

        stringstream ss;
        ss << "\n\n-- " << search_percent << "%  Samples of breadth first search --" << std::endl;
        cout << ss.str() << endl;

        //time_t btime;
        //time_t etime;
        //time(&btime);

        for (uint vidx = 0; vidx < pNVerts; vidx += 100/search_percent)
        {
            set<VertexElement*> verteleset = set<VertexElement*>();
            vector<VertexElement*> vertelevec = vector<VertexElement*>();
            queue<VertexElement*> vertelqueue = queue<VertexElement*>();

            verteleset.insert(orig_indices[vidx]);
            vertelevec.push_back(orig_indices[vidx]);
            uint ve0ncons = orig_indices[vidx]->getNCon();
            VertexElement ** ve0neighbs = orig_indices[vidx]->getNeighbours();

            fill_ve_vec(verteleset, vertelevec, vertelqueue, ve0ncons, ve0neighbs);
            pElements.clear();

            vector<VertexElement*>::iterator vertele_end = vertelevec.end();
            uint ielt = 0;

            for (vector<VertexElement*>::iterator vertele = vertelevec.begin(); vertele != vertele_end; ++vertele)
            {
                pElements.push_back(*vertele);
                pVertexPerm[(*vertele)->getIDX()]=ielt;
                ielt++;
            }

            // Note: this will reorder the vertex indices as we go- not a problem in this loop
            // but they must be reset for the final index setting to work
            reindexElements();

            uint maxdi = 0;
            for (int iv = 0; iv < pNVerts; ++iv)
            {
                VertexElement* ve = getVertex(iv);
                int ind = ve->getIDX();

                for (int i = 0; i < ve->getNCon(); ++i)
                {
                    int inbr = ve->nbrIdx(i);
                    int di = ind - inbr;
                    if (di < 0)
                    {
                        di = -di;
                    }
                    if (di > maxdi)
                    {
                        maxdi = di;
                    }
                }
            }
            if (maxdi < bestwidth)
            {
                bestone = vidx;
                bestwidth = maxdi;
                //cout << "\nOriginal vertex "<< vidx << " gives banded matrix half-width " << maxdi;
            }
        }

        // reset the vertex indices to their original value: important
        for (uint v = 0; v < pNVerts; ++v)
        {
            orig_indices[v]->setIDX(v);
        }

        set<VertexElement*> verteleset = set<VertexElement*>();
        vector<VertexElement*> vertelevec = vector<VertexElement*>();
        queue<VertexElement*> vertelqueue = queue<VertexElement*>();

        verteleset.insert(orig_indices[bestone]);
        vertelevec.push_back(orig_indices[bestone]);
        uint ve0ncons = orig_indices[bestone]->getNCon();
        VertexElement ** ve0neighbs = orig_indices[bestone]->getNeighbours();

        fill_ve_vec(verteleset, vertelevec, vertelqueue, ve0ncons, ve0neighbs);
        pElements.clear();

        vector<VertexElement*>::iterator vertele_end = vertelevec.end();
        uint ielt = 0;
        for (vector<VertexElement*>::iterator vertele = vertelevec.begin(); vertele != vertele_end; ++vertele)
        {
            pElements.push_back(*vertele);
            pVertexPerm[(*vertele)->getIDX()]=ielt;
            ielt++;
        }

        //time(&etime);
        //std::cout << "\nTime taken: "<< difftime(etime, btime) << std::endl;
    }
    // / / / / / / / / / / / /  / / / / / / / / / / / / / / / / / / / / / / //

    else if (opt_method == 1)
    {
        //time_t btime;
        //time_t etime;
        //time(&btime);

        std::vector<double> com = centerOfMass();
        double ** m = new double*[3];
        for (uint i = 0; i < 3; i++)
        {
            m[i] = new double[3];
            m[i][0] = 0.0;
            m[i][1] = 0.0;
            m[i][2] = 0.0;
        }

        for (uint ivert = 0; ivert < countVertices(); ivert++)
        {
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
        if (iret != 1)
        {
            cout << "\nWarning - eigenvalue faliure " << endl;
        }

        for (uint i = 0; i < 3; i++)
        {
            delete[] m[i];
        }
        delete[] m;

        double vx = evec[0];
        double vy = evec[1];
        double vz = evec[2];
        vx += 1.e-8;
        vy += 1.e-9;
        vz += 1.e-10;
        stringstream ss;
        ss << "aligning to axis " << vx << " " << vy << " " << vz << endl;
        //cout << ss.str();

        map<double, VertexElement*> hm;
        //vector<double> da;
        for (uint ielt = 0; ielt < countVertices(); ielt++)
        {
            VertexElement * ve = getVertex(ielt);
            double d = vx * ve->getX() + vy * ve->getY() + vz * ve->getZ();

            hm[d] = ve;
            //da.push_back(d);
        }

        //sort(da.begin(), da.end());

        pElements.clear();
        uint ielt = 0;
        map<double, VertexElement*>::const_iterator iter = hm.begin();

        while (iter != hm.end())
        {
            double d = iter->first;
            pElements.push_back(hm[d]);
            ++iter;

            // pVertexPerm[i] contains the new index in the elements array
            // for the vertex that was originally at index i
            pVertexPerm[hm[d]->getIDX()] = ielt;
            ielt++;
        }
        //time(&etime);
        //std::cout << "\nTime taken: "<< difftime(etime, btime) << std::endl;
    }
    else
    {
        std::ostringstream os;
        os << "Unknown optimization method.\n";
        throw steps::ArgErr(os.str());
    }


    reindexElements();
    reordered();

}

///////////////////////////////////////////////////////////////////////////////

void sefield::TetMesh::saveOptimal(std::string const & opt_file_name)
{

    std::fstream opt_file;

    opt_file.open(opt_file_name.c_str(),
                std::fstream::out | std::fstream::binary | std::fstream::trunc);
    uint nelems = pElements.size();
    opt_file.write((char*)&nelems, sizeof(uint));

    opt_file.write((char*)pVertexPerm, sizeof(uint) * nelems);
    opt_file.close();

}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetMesh::fill_ve_vec(set<VertexElement*> & veset, vector<VertexElement*> & vevec, queue<VertexElement*> & vequeue, uint ncons, VertexElement ** nbrs)
{

    /*

    vector<uint>newindcs=std::vector<uint>();
    for (uint i=0; i < ncons; ++i)
    {
        bool inserted = veset.insert(nbrs[i]).second;
        if (inserted)
        {
            newindcs.push_back(i);
            // This is the important one- in a vector to keep the ordering
            vevec.push_back(nbrs[i]);
            // And add to the queue too, to find their neighbours:
            vequeue.push(nbrs[i]);
        }
    }

    // Go to the next element in the queue
    if (vequeue.empty() == true) return;
    else
    {
        VertexElement * nextve = vequeue.front();
        vequeue.pop();
        VertexElement ** newneighbs = nextve->getNeighbours();
        uint nextcons = nextve->getNCon();
        fill_ve_vec(veset, vevec, vequeue, nextcons, newneighbs);
    }
    */

    while(1)
    {
        multimap<uint, uint> neighb_connections = multimap<uint, uint>();
        for (uint i=0; i < ncons; ++i)
        {
            bool inserted = veset.insert(nbrs[i]).second;
            if (inserted)
            {
                neighb_connections.insert(std::pair<uint, uint>(nbrs[i]->getNCon(), i));
            }
        }


        multimap<uint,uint>::const_iterator it_end=neighb_connections.end();
        for (multimap<uint,uint>::const_iterator it=neighb_connections.begin(); it!=it_end; ++it)
        {
            vevec.push_back(nbrs[(*it).second]);
            vequeue.push(nbrs[(*it).second]);

        }

        if (vequeue.empty() == true) return;
        else
        {


            VertexElement * nextve = vequeue.front();
            vequeue.pop();

            nbrs = nextve->getNeighbours();
            ncons = nextve->getNCon();

        }

    }

}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetMesh::mainEvec
(
    uint n,
    double ** A,
    uint * it,
    double * gamma,
    double * x1
)
{
    double eps = 1.e-10;
    double dta = 1.e-10;
    int maxit = 100;

    double x0[n];
    for (int i = 0; i < n; i++)
    {
        x0[i] = 1.0 / sqrt(i+1);
    }
    *it=-1;
    int l = 1;
    while (*it == -1 && l < maxit)
    {
        *gamma=0.;
        for (int i = 0; i < n; i++)
        {
            x1[i] = 0.;
            for (int j = 0; j < n; j++)
            {
                x1[i] += A[i][j] * x0[j];
            }

            if (fabs(x1[i]) > fabs(*gamma))
            {
                *gamma=x1[i];
            }
        }
        if (fabs(*gamma) < eps)
        {
            *it=0;
        }
        else
        {
            for (int i=0; i < n; i++)
            {
                x1[i] /= *gamma;
            }
            double phi = 0.0;
            for (int i = 0; i < n; i++)
            {
                double s = fabs(x1[i] - x0[i]);
                if (s > phi)
                {
                    phi = s;
                }
            }
            if (phi < dta)
            {
                *it=1;
            }
            else
            {
                for (int i=0; i < n; i++)
                {
                    x0[i] = x1[i];
                }
                l++;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetMesh::reordered(void)
{
    uint * tris = pTriangles;
    for (uint i = 0; i < pNTri; ++i, tris += 3)
    {
        tris[0] = pVertexPerm[tris[0]];
        tris[1] = pVertexPerm[tris[1]];
        tris[2] = pVertexPerm[tris[2]];
    }

    uint * tets = pTetrahedrons;
    for (uint i = 0; i < pNTet; ++i, tets += 4)
    {
        tets[0] = pVertexPerm[tets[0]];
        tets[1] = pVertexPerm[tets[1]];
        tets[2] = pVertexPerm[tets[2]];
        tets[3] = pVertexPerm[tets[3]];
    }
}

////////////////////////////////////////////////////////////////////////////////

sefield::VertexConnection * sefield::TetMesh::newConnection
(
    VertexElement * v1,
    VertexElement * v2
)
{
    VertexConnection * con = new VertexConnection(v1, v2);
    pConnections.push_back(con);
    return con;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> sefield::TetMesh::getVertexPermutation(void)
{
    uint nverts = countVertices();

    // Return a copy of this vector
    std::vector<uint> vert_ret = std::vector<uint>(nverts);

    for (uint i = 0; i < nverts; ++i)
    {
        vert_ret[i] = pVertexPerm[i];
    }

    return vert_ret;
}

////////////////////////////////////////////////////////////////////////////////
/*
void sefield::TetMesh::savematrix(void)
{
    ofstream fout;
    fout.open("matrixdata.txt");

    uint nverts = countVertices();
    uint ncons =    pConnections.size();

    std::vector<std::vector<int> > cmatrix = std::vector<std::vector<int> >(nverts);
    for (uint i=0; i < nverts; ++i)
    {
        cmatrix[i] = std::vector<int>(nverts);
    }

    for (uint v = 0; v < nverts; ++v)
    {
        fout << v << " ";
        VertexElement * nextve = pElements[v];
        VertexElement ** newneighbs = nextve->getNeighbours();
        uint nextcons = nextve->getNCon();
        for (uint c = 0; c < nextcons; ++c)
        {
            fout << newneighbs[c]->getIDX();
            fout << " ";
        }
        fout << "\n";
    }

    fout.close();

}
*/
////////////////////////////////////////////////////////////////////////////////

// END
