
/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/efield/matrix.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace sefield = steps::solver::efield;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

sefield::Matrix::Matrix(uint n0)
: pA(0)
, pWS(0)
, pN(n0)
, pPerm(0)
, pSign(0)
{
    // Check size argument.
    AssertLog(pN != 0);

    pA = new double*[pN];
    for (uint i = 0; i < pN; ++i)
    {
        pA[i] = new double[pN];
    }
    pPerm = new int[pN];
    pWS = new double[pN];
}

////////////////////////////////////////////////////////////////////////////////

sefield::Matrix::Matrix(uint nn, double ** da)
: pA(0)
, pWS(0)
, pN(nn)
, pPerm(0)
, pSign(0)
{
    // Check input arguments.
    AssertLog(pN != 0);
    AssertLog(da != 0);

    pA = new double*[pN];
    for (uint i = 0; i < pN; ++i)
    {
        pA[i] = new double[pN];
        for (uint j = 0; j < pN; ++j)
        {
            pA[i][j] = da[i][j];
        }
    }

    pPerm = new int[pN];
    pWS = new double[pN];
}

////////////////////////////////////////////////////////////////////////////////

sefield::Matrix::~Matrix(void)
{
    delete[] pPerm;
    delete[] pWS;

    for (uint i = 0; i < pN; ++i)
    {
        delete[] pA[i];
    }
    delete[] pA;
}

////////////////////////////////////////////////////////////////////////////////

void sefield::Matrix::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&pN, sizeof(uint));
    cp_file.write((char*)&pSign, sizeof(int));
    cp_file.write((char*)pA, sizeof(double) * pN * pN);
    cp_file.write((char*)pWS, sizeof(double) * pN);
    cp_file.write((char*)pPerm, sizeof(int) * pN);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::Matrix::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&pN, sizeof(uint));
    cp_file.read((char*)&pSign, sizeof(int));
    cp_file.read((char*)pA, sizeof(double) * pN * pN);
    cp_file.read((char*)pWS, sizeof(double) * pN);
    cp_file.read((char*)pPerm, sizeof(int) * pN);
}

////////////////////////////////////////////////////////////////////////////////

sefield::Matrix * sefield::Matrix::copy(void)
{
    // NOTE: The memory is cleaned up in det() and inverse()
    Matrix * m = new Matrix(pN, pA);
    for (uint i = 0; i < pN; ++i)
    {
        m->pPerm[i] = pPerm[i];
    }
    m->pSign = pSign;
    return m;
}

////////////////////////////////////////////////////////////////////////////////

double * sefield::Matrix::lvprod(double * v)
{
    // NOTE: THe memory allocated for this vector is eventually cleaned up in
    // Tetcoupler::coupleMesh()
    double * r = new double[pN];
    fill_n(r, pN, 0.0);

    for (uint i = 0; i < pN; ++i)
    {
        for (uint j = 0; j < pN; ++j)
        {
            r[j] += v[i] * pA[i][j];
        }
    }

    return r;
}

////////////////////////////////////////////////////////////////////////////////

double sefield::Matrix::det(void)
{
    Matrix * t = copy();
    t->LU();
    double d = 1.0 * t->pSign;
    for (uint i = 0; i < pN; ++i)
    {
        d *= t->pA[i][i];
    }
    delete t;
    return d;
}

////////////////////////////////////////////////////////////////////////////////

void sefield::Matrix::LU(void)
{
    int i, imax, j, k;
    double big, dum, sum, temp;
    double * vv = new double[pN];
    double TINY = 1.0e-20;

    pSign = 1;

    imax = -1;
    for (i = 0; i < pN; ++i)
    {
        big = 0.0;
        for (j = 0; j < pN; ++j)
        {
            if ((temp = fabs(pA[i][j])) > big)
            {
                big = temp;
            }
        }
        if (big == 0.0)
        {
            //  Sp("Singular Matrix in routine LUDCMP");
        }
        vv[i] = 1.0 / big;
    }

    for (j = 0; j < pN; ++j)
    {
        for (i = 0; i < j; ++i)
        {
            sum = pA[i][j];
            for (k = 0; k < i; ++k)
            {
                sum -= pA[i][k] * pA[k][j];
            }
            pA[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i < pN; ++i)
        {
            sum = pA[i][j];
            for (k = 0; k < j; ++k)
            {
                sum -= pA[i][k] * pA[k][j];
            }
            pA[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big)
            {
                big = dum;
                imax = i;
            }
        }
        if (j != imax)
        {
            for (k = 0; k < pN; ++k)
            {
                dum = pA[imax][k];
                pA[imax][k] = pA[j][k];
                pA[j][k] = dum;
            }
            pSign = -pSign;
            vv[imax] = vv[j];
        }
        pPerm[j] = imax;
        if (pA[j][j] == 0.0)
        {
            pA[j][j] = TINY;
        }
        if (j != pN)
        {
            dum = 1.0 / (pA[j][j]);
            for (i = j + 1; i < pN; ++i)
            {
                pA[i][j] *= dum;
            }
        }
    }
    delete[] vv;
}

////////////////////////////////////////////////////////////////////////////////

sefield::Matrix * sefield::Matrix::inverse(void)
{
    Matrix * t = copy();
    Matrix * r = copy();
    t->LU();

    double * c = new double[pN];
    for (uint j = 0; j < pN; ++j)
    {
        for (uint i = 0; i < pN; ++i)
        {
            c[i] = 0.0;
        }
        c[j] = 1.0;
        t->lubksb(c);
        for (uint i = 0; i < pN; ++i)
        {
            r->pA[i][j] = c[i];
        }
    }
    delete t;
    delete[] c;

    // NOTE: memory allocated for r is cleaned in TetCoupler::fluxCoeficients
    return r;
}

////////////////////////////////////////////////////////////////////////////////

double * sefield::Matrix::lubksb(double * b)
{
    int ip;
    int ii = -1;
    double sum;

    for (int i = 0; i < pN; ++i)
    {
        ip = pPerm[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii >= 0)
        {
            for (int j = ii; j < i; ++j)
            {
                sum -= pA[i][j] * b[j];
            }
        } else if (sum != 0.0)
        {
            ii = i;
        }
        b[i] = sum;
    }

    for (int i = pN - 1; i >= 0; --i)
    {
        sum = b[i];
        for (int j = i + 1; j < pN; ++j)
        {
            sum -= pA[i][j] * b[j];
        }
        b[i] = sum / pA[i][i];
    }

    return b;
}

////////////////////////////////////////////////////////////////////////////////

// END
