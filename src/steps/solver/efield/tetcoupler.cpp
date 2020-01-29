/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/efield/matrix.hpp"
#include "steps/solver/efield/tetcoupler.hpp"
#include "steps/solver/efield/tetmesh.hpp"
#include "steps/solver/efield/vertexconnection.hpp"
#include "steps/solver/efield/vertexelement.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace sefield = steps::solver::efield;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

sefield::TetCoupler::TetCoupler(TetMesh * mesh)
: pMesh(mesh)
{
}

////////////////////////////////////////////////////////////////////////////////

sefield::TetCoupler::~TetCoupler()
= default;

////////////////////////////////////////////////////////////////////////////////

void sefield::TetCoupler::coupleMesh()
{

    typedef vector<double*> doublePVec;

    // For each vertex allocate space to accumulate coupling coefficients
    // to neighbors.
    uint nvertices = pMesh->countVertices();
    doublePVec vccs(nvertices);
//    auto vccs_bgn = vccs.begin();
//    auto vccs_end = vccs.end();
#pragma omp parallel for
    for (uint i = 0; i < nvertices; ++i)
    {
        VertexElement * vertex = pMesh->getVertex(i);
        AssertLog(vertex->getIDX() == i);
        uint ncons = vertex->getNCon();
        //double * new_vccs = new double[ncons]();
        // initialize the array.
        // std::fill_n(new_vccs, ncons, 0.0);
        //vccs[i] = new_vccs;
        vccs[i] = new double[ncons]();
    }

    // Loop over all vertices in the mesh. For each vertex:
    //
    //   * Fetch a list of neighbouring tetrahedra ('neighbouring'
    //     meaning a tetrahedron that includes the current vertex).
    //
    //   * For each neighbouring tetrahedron:
    //
    //       * Transform its vertex indices into pointers to
    //         corresponding VertexElement objects.
    //
    //       * Compute flux coefficients and distribute them over the
    //         'neighbouring' vertex nodes.
    //
#pragma omp parallel for
    for (uint ivert = 0; ivert < nvertices; ++ivert)
    {
        VertexElement* ve = pMesh->getVertex(ivert);

        // Now we need the following information about vertex ve:
        // a list of tetrahedrons of which ve is a part. This list
        // is returned as indices in the ve's neighbour's list.
        const vector<std::array<uint, 3>>  vti = pMesh->getNeighboringTetrahedra(ve);

        // for (int[] tetinds : pMesh.neighboringTetrahedra(ve)) {

        // Loop over neighoubring tetrahedra
        for (auto tetinds : vti) {
            // Get the tetrahedron by vertex neighbour indices (e.g if the vertex has 6 nighbours one
            // tetrahedron may be 2,4,5 or something

          // Get an array of the neighbouring vertices.
            // (As said before, the ints in tetinds are indices in
            // the neighbours list of vertex ve, not global indices.
            auto ves = new VertexElement*[3];
            for (uint i = 0; i < 3; ++i)
            {
                ves[i] = ve->getNeighbor(tetinds[i]);
            }
            // Compute the flux into the polyhedron around the vertex in
            // terms of the potential difference to each of the corners.
            double facs[3] = {0., 0., 0.};
            //facs[0]=0.0;
            //facs[1]=0.0;
            //facs[2]=0.0;


            fluxCoeficients(ve, ves, facs);
            // Destroy the ves array
            delete[] ves;

            // Accumulate coefficients for neighbours in vccs.
            for (uint i = 0; i < 3; ++i)
            {
                //AssertLog(facs[i] > 0.0);
                vccs[ivert][tetinds[i]] += facs[i];
            }

            // Delete the return vector.
            //delete[] facs;
        }

        /*
        // Destroy the list of neighbouring tetrahedrons.
        vector<vector<uint> >::iterator vti_end = vti.end();
        for (vector<vector<uint> >::iterator vti_r = vti.begin();
            vti_r != vti_end; ++vti_r)
        {
            delete[] (*vti_r);
        }
        */
    }

    // If all has gone according to plan, then the fluxes are symmetric
    // and fall into the form flux_j = sum (w_i,j (v_i - v_j))
    // with w_i,j = w_j,i
    // the w_i,j is then the coupling constant for the connection from i to j
    uint ntot = 0;
    uint ndif = 0;

#pragma omp parallel for
    for (uint icon = 0; icon < pMesh->ncon(); ++icon)
    {
        VertexConnection * vc = pMesh->getConnection(icon);
        VertexElement * va = vc->getA();
        VertexElement * vb = vc->getB();

        uint va_idx = va->getIDX();
        uint vb_idx = vb->getIDX();

        // vccs only contains elements for the adjacent points - need to
        // look through ti find which one corresponds to the other vertex
        // in this connection
        double wab = 0.0;
        for (uint i = 0; i < va->getNCon(); ++i)
        {
            if (va->getNeighbor(i)->getIDX() == vb_idx)
            {
                wab = vccs[va_idx][i];
            }
        }

        // do the same the other way round, just to check
        double wba = 0.0;
        for (uint i = 0; i < vb->getNCon(); ++i)
        {
            if (vb->getNeighbor(i)->getIDX() == va_idx)
            {
                wba = vccs[vb_idx][i];
            }
        }

        if (dblsDiffer(wab, wba))
        {
#pragma omp atomic
            ndif += 1;
#ifdef _OPENMP            
            auto tid = omp_get_thread_num();
            if (!tid) CLOG_N_TIMES(100, DEBUG, "general_log") << "symmetry miscount " << wab << " " << wba; 
#endif 
}
        else
        {
            vc->setGeomCouplingConstant(wab);
        }
#pragma omp atomic
         ntot += 1;
    }

    // should ndif > 0 throw an exception?
    if (ndif > 0) {
        std::ostringstream os;
        os << ndif << " out of " << ntot << " failed sym test. Nvert=" << pMesh->countVertices();
        ProgErrLog(os.str());
    }

    // Deallocate vccs.
#pragma omp parallel for
    for (auto i = 0u; i < vccs.size(); ++i) {
      delete[] vccs[i];
    }
    //for (doublePVecI vccs_i = vccs_bgn; vccs_i != vccs_end; ++vccs_i)
    //{
    //    delete[] (*vccs_i);
    //}
}

////////////////////////////////////////////////////////////////////////////////

bool sefield::TetCoupler::dblsDiffer(double a, double b)
{
// Old code:
//
//    bool ret = false;
//    if (a == 0 && b == 0)
//    {
//        // ok, same
//    }
//    else if (a + b == 0)
//    {
//        ret = true;
//    }
//    else if (abs((a - b) / (a + b)) > 1.e-12)
//    {
//        ret = true;
//    }
//    return ret;

    bool ret = false;
    if ((a == 0.0) && (b == 0.0))
    {
        // ok, same
    }
    else if ((a + b) == 0.0)
    {
        ret = true;
    }
    else if (fabs((a - b) / (a + b)) > 1.e-7)
    {
        ret = true;
    }
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetCoupler::cross_product(double * a, double * b, double * c)
{

    c[0] = (a[1]*b[2]) - (a[2]*b[1]);
    c[1] = (a[2]*b[0]) - (a[0]*b[2]);
    c[2] = (a[0]*b[1]) - (a[1]*b[0]);

}

////////////////////////////////////////////////////////////////////////////////

void sefield::TetCoupler::fluxCoeficients(sefield::VertexElement * ve, sefield::VertexElement ** ves, double * ret)
{
    // Assume that the potential varies linearly across the element between
    // the values at the vertices. Then the procedure is:
    //
    // First, find an expression for the field in terms of the vertex values.
    //
    // Second, integrate this component of the field that is normal to the
    // boundary over the boundary of the polyhedron centered on ve.
    //
    // The integral is simple because the field is constant - it just
    // needs a cross product with the triangles forming the part of the
    // polyhedron boundary that is in this tetrahedron.
    //
    //
    // The result will be a linear expression in terms of the vertex
    // values. The return value is an array containing the coefficients
    // in this expression.

    // Matrix ordering is [row][colunmn], so m[2][0] is the third row,
    // first column of matrix with the vectors to the three adjacent


    // vertices in its rows.
    auto** m = new double*[3];
    for (uint i = 0; i < 3; ++i)
    {
        m[i] = new double[3];
    }

    // This matrix is already the transpose of matrix in equation 13
    for (uint iv = 0; iv < 3; ++iv)
    {
        m[iv][0] = ves[iv]->getX() - ve->getX();
        m[iv][1] = ves[iv]->getY() - ve->getY();
        m[iv][2] = ves[iv]->getZ() - ve->getZ();
    }
    // Create a matrix with m as initialization.
    auto* vm = new Matrix(3, m);

    // Need to be consistent about the orientation of tetrahedra. One way
    // is to take the dot product of one vector with the cross product of
    // the other two and make sure it always has the same sign by swapping
    // rows two and three if neessary. The dot-cross product is the same
    // as the determinant up to a scaling factor


    bool swap = false;
    if (vm->det() < 0.)
    {
        // switch second and third rows of m;
        for (uint i = 0; i < 3; ++i)
        {
            double w = m[1][i];
            m[1][i] = m[2][i];
            m[2][i] = w;
            swap = true;
        }
        delete vm;
        vm = new Matrix(3, m);
    }


    // CORRECT METHOD
    Matrix * vmi = vm->inverse();
    delete vm;

    // consider the other vertices in turn
    for (int ivert = 0; ivert < 3; ivert++)
    {
        int ip = ivert;
        int iq = ivert + 1;
        int ir = ivert + 2;

        if (iq >= 3) {
            iq -= 3;
        }

        if (ir >= 3) {
            ir -= 3;
        }

        double c1[3];
        double c2[3];
        double c3[3];

        cross_product(m[ip], m[iq], c1);
        cross_product(m[iq], m[ir], c2);
        cross_product(m[ir], m[ip], c3);

        double vec[3];
        double f = 1.0/12.0;
        double g = 2.0/12.0;

        vec[0] = (f * c1[0]) + (g * c2[0]) + (f * c3[0]);
        vec[1] = (f * c1[1]) + (g * c2[1]) + (f * c3[1]);
        vec[2] = (f * c1[2]) + (g * c2[2]) + (f * c3[2]);

        double * wk = vmi->lvprod(vec);

        // accumulate the contribution from these two triangles into the final return array
        // the 0.5 is because the area of the triangle is half the cross product
        for (int i = 0; i < 3; i++)
        {
            ret[i] += (0.5 * wk[i]);
        }
        delete wk;
    }

    // if we swapped the sense of the tetrahedron above, swap the results
    // back so they go in the right slots when returned
    if (swap)
    {
        double w = ret[1];
        ret[1] = ret[2];
        ret[2] = w;
    }

    // memory allocated for m not cleaned up. Lets do that now
    for (uint i = 0; i < 3; ++i)
    {
        delete[] m[i];
    }
    delete[] m;

    delete vmi;

}

////////////////////////////////////////////////////////////////////////////////

// END
