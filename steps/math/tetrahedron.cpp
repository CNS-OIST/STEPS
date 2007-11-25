////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
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
// $Id$
////////////////////////////////////////////////////////////////////////////////

// Standard headers.
#include <Python.h>
#include <numpy/arrayobject.h>
#include <cmath>
#include <cstdio>
#include <cstring>

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/rng/rng.hpp>

NAMESPACE_ALIAS(steps::rng, srng);

// A SERVICE ANNOUNCEMENT FROM THE TECH DEPT.
// 
// This module was written directly in C for a number of reasons:
//
// 1. We found that these routines become a bottle neck for many 
//    mesh-related operations, prompting an implementation in C.
// 2. We do not want to involve SWIG here, because this introduces
//    two extra layers of function calls, one in Python and one
//    in C. Whether this is a real problem remains open for debate,
//    but we play on safe. Furthermore, we still need direct access
//    to the NumPy C API, so why not use the Python C/API directly
//    as well...? 
// 3. We are currently reluctant to introduce dependencies on other 
//    wrapper or assistance tools, such as Pyrex or Weave. While such
//    tools greatly simplify the glue code that needs to be written
//    to make the module behave well from the point of view of Python
//    (e.g. adding docstrings, argument parsing and exception handling),
//    they would have to be installed whenever, wherever somebody wants 
//    to work on STEPS code. We want to avoid making too many demands
//    on the build process or on the UNIX environment required for working 
//    on STEPS, because this must be possible even in the most adverse 
//    conditions, where time is limited and computing facilities might be 
//    basic (e.g. during crash visits to other labs, conferences, 
//    workshops, ...). 
//     
// For optimal understanding of this code, please refer to these works:
//
// 1. Python/C API Reference Manual
//    http://docs.python.org/api/api.html
//
// 2. Travis E. Oliphant, Guide to Numpy, E-book (2006)
//    http://www.tramy.us/guidetoscipy.html
//
// Thank you.

////////////////////////////////////////////////////////////////////////////////

static int linsolve(int n, int rhs_num, double a[])
{
    // Precompute n+ rhs_num
    int n_plus_rhs_num = n + rhs_num;

    // Loop over all rows.
    for (int j = 0; j < n; ++j)
    {
        // Choose a pivot row: first we select j.
        int ipivot = j;
        double apivot = a[j + j * n];
        // But we really want the largest.
        for (int i = j + 1; i < n; ++i)
        {
            if (fabs(apivot) < fabs(a[i + j * n]))
            {
                apivot = a[i + j * n];
                ipivot = i;
            }
        }

        // Singular system: report!
        if (apivot == 0.0)
        {
            return j;
        }

        // Swap.
        for (int i = 0; i < n_plus_rhs_num; ++i)
        {
            double temp          = a[ipivot + i * n];
            a[ipivot + i * n]    = a[j + i * n];
            a[j + i * n]         = temp;
        }

        // a[j,j] becomes 1.
        // a[j + j * n] = 1.0;
        for (int k = j; k < n_plus_rhs_num; ++k)
        {
            a[j + k * n] = a[j + k * n] / apivot;
        }

        // a[i,j] becomes 0.
        for (int i = 0; i < n; ++i)
        {
            if (i != j)
            {
                double factor = a[i + j * n];
                // a[i + j * n] = 0.0;
                for (int k = j; k < n_plus_rhs_num; ++k)
                {
                    a[i + k * n] = a[i + k * n] - factor * a[j + k * n];
                }
            }
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

// CODE NOTE:
// 
// I'm not 100% sure whether the following code, in which computation of 
// a 4x4 determinant has been written out explicitly, is absolutely the 
// fastest. Modern optimizing compilers might in fact do better when given 
// the algorithm expressed in loops. This might give them the necessary
// information to perform automatic loop vectorization and to avoid stalls
// in the pipeline. 
// 
// If the code in this module remains a bottleneck, this should be 
// investigated first.

// Directly computes the determinant for the following 4x4 matrix:
// 
//      (0,0)     (0,1)     (0,2)     1.0
//      =m[0]     =m[1]     =m[2]
// 
//      (1,0)     (1,1)     (1,2)     1.0
//      =m[3]     =m[4]     =m[5]
// 
//      (1,0)     (1,1)     (1,2)     1.0
//      =m[6]     =m[7]     =m[8]
// 
//      (1,0)     (1,1)     (1,2)     1.0
//      =m[9]     =m[10]    =m[11]

static inline npy_double DET4X3_1(npy_double * m)
{
    return m[5]*m[7]*m[9]-m[2]*m[7]*m[9]-1.0*m[4]*m[8]*m[9]+m[1]*
           m[8]*m[9]+m[2]*m[4]*m[9]-m[1]*m[5]*m[9]-1.0*m[5]*m[6]*
           m[10]+m[2]*m[6]*m[10]+1.0*m[3]*m[8]*m[10]-m[0]*m[8]*
           m[10]-m[2]*m[3]*m[10]+m[0]*m[5]*m[10]+1.0*m[4]*m[6]*
           m[11]-m[1]*m[6]*m[11]-1.0*m[3]*m[7]*m[11]+m[0]*m[7]*
           m[11]+m[1]*m[3]*m[11]-m[0]*m[4]*m[11]-m[2]*m[4]*m[6]*
           1.0+m[1]*m[5]*m[6]*1.0+m[2]*m[3]*m[7]*1.0-m[0]*m[5]*
           m[7]*1.0-m[1]*m[3]*m[8]*1.0+m[0]*m[4]*m[8];
}

// Directly computes the determinant for the following 4x4 matrix:
// 
//      (0,0)     (0,1)     (0,2)     1.0
//      =m0[0]    =m0[1]    =m0[2]
// 
//      (1,0)     (1,1)     (1,2)     1.0
//      =m1[0]    =m1[1]    =m1[2]
// 
//      (1,0)     (1,1)     (1,2)     1.0
//      =m2[0]    =m2[1]    =m2[2]
// 
//      (1,0)     (1,1)     (1,2)     1.0
//      =m3[0]    =m3[1]    =m3[2]

static inline npy_double DET4X3_2
(
    npy_double * m0, npy_double * m1,
    npy_double * m2, npy_double * m3
)
{
    return m1[2]*m2[1]*m3[0]-m0[2]*m2[1]*m3[0]-1.0*m1[1]*m2[2]*m3[0]+m0[1]*
           m2[2]*m3[0]+m0[2]*m1[1]*m3[0]-m0[1]*m1[2]*m3[0]-1.0*m1[2]*m2[0]*
           m3[1]+m0[2]*m2[0]*m3[1]+1.0*m1[0]*m2[2]*m3[1]-m0[0]*m2[2]*
           m3[1]-m0[2]*m1[0]*m3[1]+m0[0]*m1[2]*m3[1]+1.0*m1[1]*m2[0]*
           m3[2]-m0[1]*m2[0]*m3[2]-1.0*m1[0]*m2[1]*m3[2]+m0[0]*m2[1]*
           m3[2]+m0[1]*m1[0]*m3[2]-m0[0]*m1[1]*m3[2]-m0[2]*m1[1]*m2[0]*
           1.0+m0[1]*m1[2]*m2[0]*1.0+m0[2]*m1[0]*m2[1]*1.0-m0[0]*m1[2]*
           m2[1]*1.0-m0[1]*m1[0]*m2[2]*1.0+m0[0]*m1[1]*m2[2];
}

////////////////////////////////////////////////////////////////////////////////

// Provides access to the internals of a SWIG-exposed extension module.
// Source: the wizzards at Boost-Python.
//
// When you want to use a C++ object that has been exposed to Python by
// using SWIG as an intermediary, the only way to get access to that object
// from within a non-SWIG C++ extension, is to mimick some SWIG internals
// locally. This makes up dependent on knowing certain details of SWIG's
// internal kitchen. Practically speaking, this internal kitchen seems to 
// have stabilized. From the point of view of engineering principles, this 
// is clearly not an ideal situation.
//
// Maybe we should never have used SWIG in the first place... Als ik ooit
// eens vijf minuten tijd had...

struct PySwigObject 
{
    PyObject_HEAD 
    void * ptr;
    const char * desc;
};

void * getSWIGptr(PyObject* obj)
{
    char thisStr[] = "this";
    // First we need to get the this attribute from the Python object.
    if (!PyObject_HasAttrString(obj, thisStr))
        return 0; 
    PyObject* thisAttr = PyObject_GetAttrString(obj, thisStr);
    if (thisAttr == 0)
        return 0;
    // This Python Object is a SWIG wrapper and contains our pointer.
    return (((PySwigObject*)thisAttr)->ptr);
}

////////////////////////////////////////////////////////////////////////////////

// Since many functions in this module expect the same input arguments, we 
// decided to put the parsing and checking code in a separate function.
// 
// OVERVIEW:
//
// Case 1: only p
//     Check size (4N, 3) and dtype (floating)
//     Set case1flag, ntets
// Case 2: p ant t
//     Check size (N, 3; M, 4); and dtype (floating, integer)
//     Set case2flag, ntets
//
// RETURNS:
//      0       Error occured, e.g. parameters do not conform.
//      1       Only p was specified
//      2       Both p and t were specified.

enum CallTypePT 
{ 
    CTPT_ERR = 0, 
    CTPT_P_NOT_T = 1, 
    CTPT_P_AND_T = 2 
};

static CallTypePT parse_input_CallTypePT
(
    PyObject * args,
    PyObject ** p,
    PyObject ** t,
    npy_intp & npnts,
    npy_intp & ntets,
    npy_double ** pptr,
    npy_double ** pptr_max,
    npy_intp ** tptr,
    npy_intp ** tptr_max
)
{
    PyObject * in_p = 0;
    PyObject * in_t = 0;
    
    // Unpack arguments & tentatively check whether we have "only-p" or
    // "p-and-t" style input.
    if (!PyArg_ParseTuple(args, "O|O", &in_p, &in_t)) return CTPT_ERR;
    CallTypePT rescode = (in_t == 0 ? CTPT_P_NOT_T : CTPT_P_AND_T);
    
    // Deal with p.
    *p = PyArray_ContiguousFromAny(in_p, PyArray_DOUBLE, 1, 2);
    if (*p == 0) return CTPT_ERR;
    *pptr = (npy_double*)PyArray_DATA(*p);
    if (*pptr == 0)
    {
        PyErr_Format(PyExc_RuntimeError, "Cannot access ndarray for arg 'p'");
        return CTPT_ERR;
    }
    if (PyArray_NDIM(*p) == 1)
    {
        int ncols = PyArray_DIM(*p, 0);
        if ((ncols % 3) != 0)
        {
            PyErr_Format(PyExc_ValueError, 
                "Size of arg 'p' must be a multiple of 3");
            return CTPT_ERR;
        }
        npnts = ncols / 3;
    }
    else
    {
        if (PyArray_DIM(*p, 1) != 3)
        {
            PyErr_Format(PyExc_ValueError, "Array 'p' must have 3 columns");
            return CTPT_ERR;
        }
        npnts = PyArray_DIM(*p, 0);
    }
    *pptr_max = *pptr + (3 * npnts);
    
    // If only p was specified, we can exit now.
    if (rescode == CTPT_P_NOT_T)
    {
        if ((npnts % 4) != 0)
        {
            PyErr_Format(PyExc_ValueError, 
                "Size of arg 'p' must be a multiple of 4");
            return CTPT_ERR;
        }
        ntets = npnts / 4;
        return rescode;
    }
    
    // Else, continue to deal with t.
    *t = PyArray_ContiguousFromAny(in_t, PyArray_INTP, 1, 2);
    if (*t == 0) return CTPT_ERR;
    *tptr = (npy_intp*)PyArray_DATA(*t);
    if (*tptr == 0)
    {
        PyErr_Format(PyExc_RuntimeError, "Cannot access ndarray for arg 't'");
        return CTPT_ERR;
    }
    if (PyArray_NDIM(*t) == 1)
    {
        int ncols = PyArray_DIM(*t, 0);
        if ((ncols % 4) != 0)
        {
            PyErr_Format(PyExc_ValueError,
                "Size of arg 't' must be a multiple of 4");
            return CTPT_ERR;
        }
        ntets = ncols / 4;
    }
    else
    {
        if (PyArray_DIM(*t, 1) != 4)
        {
            PyErr_Format(PyExc_ValueError, "Array 't' must have 4 columns");
            return CTPT_ERR;
        }
        ntets = PyArray_DIM(*t, 0);
    }
    *tptr_max = *tptr + (ntets * 4);
    return rescode;
}

// Code for creating output matrix of doubles.
#define CREATE_OUTPUT_NDARRAY_DOUBLE(ARG_DIMS, ARG_NROWS) \
    out_dims[0] = ARG_NROWS; \
    out = PyArray_SimpleNew(ARG_DIMS, out_dims, PyArray_DOUBLE); \
    if (out == 0) goto ret_err; \
    outptr = (npy_double*)PyArray_DATA(out); \
    if (outptr == 0) \
    { \
        PyErr_Format(PyExc_RuntimeError, "Cannot access result ndarray"); \
        goto ret_err; \
    }

// This macro transforms a strip of 4 indices in a t(et) matrix into 
// 4 pointers into the p(oint) matrix. Performs out-of-bounds checking.
//
// Requires function definition of:
//     npy_intp npnts;
//     {Goto label} ret_err
// And local scope definition of:
//     npy_intp tidx;
//     npy_double * p0;
//     npy_double * p1;
//     npy_double * p2;
//     npy_double * p3;
#define TPTR_TO_4PPTR \
    tidx = *(tptr++); \
    if (tidx >= npnts) \
    { \
        PyErr_Format(PyExc_IndexError, \
            "Point index out of range (%d)", tidx); \
        goto ret_err; \
    } \
    p0 = pptr + (3 * tidx); \
    tidx = *(tptr++); \
    if (tidx >= npnts) \
    { \
        PyErr_Format(PyExc_IndexError, \
            "Point index out of range (%d)", tidx); \
        goto ret_err; \
    } \
    p1 = pptr + (3 * tidx); \
    tidx = *(tptr++); \
    if (tidx >= npnts) \
    { \
        PyErr_Format(PyExc_IndexError, \
            "Point index out of range (%d)", tidx); \
        goto ret_err; \
    } \
    p2 = pptr + (3 * tidx); \
    tidx = *(tptr++); \
    if (tidx >= npnts) \
    { \
        PyErr_Format(PyExc_IndexError, \
            "Point index out of range (%d)", tidx); \
        goto ret_err; \
    } \
    p3 = pptr + (3 * tidx);

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_vol__doc__,
"Compute the volume for one or more tetrahedrons.\n\
\n\
PARAMETERS:\n\
    p\n\
        An N_pnt * 3 array of points, or similar sized sequence\n\
        object that can be converted to an array. Only 2-dimensional\n\
        arrays/sequences are allowed.\n\
    t\n\
        An optional argument. \n\
        If specified, it is interpreted as an N_tet * 4 array of \n\
        integer indices into p (or any nested sequence object that\n\
        can be converted to such an array). Only 1- or 2-dimensional\n\
        arrays_sequences are allowed.\n\
        If not specified, N_pnt must be a multiple of four, meaning\n\
        that each consecutive block of four rows in p defines one\n\
        tetrahedron.\n\
\n\
RETURNS:\n\
    A 1-dimensional array of size N_tet with volumes.\n");

static PyObject * tet_vol(PyObject * self, PyObject * args)
{
    // OVERVIEW:
    // Parse input arguments
    // Prepare output array (2D ntets * 3, double)
    // Perform volume computation
    //   If case1flag:
    //     Loop over all elements of p in blocks of 4*3=12 doubles
    //       Use DET4X3_1 for each block
    //   If case2flag:
    //     Loop over each row of t
    //       Check index of each column in each strip
    //       Use DET4X3_2 for 4 indexed rows in p 
    // Return output
    
    PyObject * out = 0;
    npy_intp out_dims[1] = { 0 };
    npy_double * outptr;
    
    PyObject * p = 0;
    npy_intp npnts = 0;
    npy_double * pptr;
    npy_double * pptr_max;
    PyObject * t = 0;
    npy_intp ntets = 0;
    npy_intp * tptr;
    npy_intp * tptr_max;
    CallTypePT presult = 
        parse_input_CallTypePT(args, &p, &t, npnts, ntets, 
        &pptr, &pptr_max, &tptr, &tptr_max);
    if (presult == CTPT_ERR) goto ret_err;
    
    CREATE_OUTPUT_NDARRAY_DOUBLE(1, ntets);
    // Perform volume computation.
    if (presult == CTPT_P_NOT_T)
    {
        for (; pptr < pptr_max; pptr += 12)
            *(outptr++) = fabs(DET4X3_1(pptr) / 6.0);
    }
    else if (presult == CTPT_P_AND_T)
    {
        npy_intp tidx;
        npy_double * p0;
        npy_double * p1;
        npy_double * p2;
        npy_double * p3;
        while (tptr < tptr_max)
        {
            TPTR_TO_4PPTR
            *(outptr++) = fabs(DET4X3_2(p0, p1, p2, p3) / 6.0);
        }
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError, "Bug: code shouldn't reach here (0)");
        goto ret_err;
    }
    
    // Return output.
    Py_XDECREF(p);
    Py_XDECREF(t);
    return out;
    
ret_err:
    // Clean up code after error.
    Py_XDECREF(p);
    Py_XDECREF(t);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_barycenter__doc__,
"Compute the barycenter for one or more tetrahedrons.\n\
\n\
PARAMETERS:\n\
    p\n\
        An N_pnt * 3 array of points, or similar sized sequence\n\
        object that can be converted to an array. Only 2-dimensional\n\
        arrays/sequences are allowed.\n\
    t\n\
        An optional argument. \n\
        If specified, it is interpreted as an N_tet * 4 array of \n\
        integer indices into p (or any nested sequence object that\n\
        can be converted to such an array). Only 1- or 2-dimensional\n\
        arrays/sequences are allowed.\n\
        If not specified, N_pnt must be a multiple of four, meaning\n\
        that each consecutive block of four rows in p defines one\n\
        tetrahedron.\n\
\n\
RETURNS:\n\
    A 2-dimensional array of size N_tet * 3 with barycenters.\n");

static PyObject * tet_barycenter(PyObject * self, PyObject * args)
{
    // OVERVIEW:
    // Parse input arguments
    // Prepare output array (2D ntets// 3, double)
    // Perform barycenter computation
    //   If case1flag:
    //     Loop over all elements of p in blocks of 4*3=12 doubles
    //       Compute barycenter
    //   If case2flag:
    //     Loop over each row of t
    //       Check index of each column in each strip
    //       Compute barycenter 
    // Return output
    
    npy_intp out_dims[2] = { 0, 3 };
    PyObject * out = 0;
    npy_double * outptr;
    
    PyObject * p = 0;
    npy_intp npnts = 0;
    npy_double * pptr;
    npy_double * pptr_max;
    PyObject * t = 0;
    npy_intp ntets = 0;
    npy_intp * tptr;
    npy_intp * tptr_max;
    CallTypePT presult = 
        parse_input_CallTypePT(args, &p, &t, npnts, ntets, 
        &pptr, &pptr_max, &tptr, &tptr_max);
    if (presult == CTPT_ERR) goto ret_err;
    
    CREATE_OUTPUT_NDARRAY_DOUBLE(2, ntets)
    
    // Compute barycenters.
    if (presult == CTPT_P_NOT_T)
    {
        for (; pptr < pptr_max; pptr += 12)
        {
            *(outptr++) = (pptr[ 0] + pptr[ 3] + pptr[ 6] + pptr[ 9]) / 4.0;
            *(outptr++) = (pptr[ 1] + pptr[ 4] + pptr[ 7] + pptr[10]) / 4.0;
            *(outptr++) = (pptr[ 2] + pptr[ 5] + pptr[ 8] + pptr[11]) / 4.0;
        }
    }
    else if (presult == CTPT_P_AND_T)
    {
        npy_intp tidx;
        npy_double * p0;
        npy_double * p1;
        npy_double * p2;
        npy_double * p3;
        while (tptr < tptr_max)
        {
            TPTR_TO_4PPTR
            *(outptr++) = (p0[0] + p1[0] + p2[0] + p3[0]) / 4.0;
            *(outptr++) = (p0[1] + p1[1] + p2[1] + p3[1]) / 4.0;
            *(outptr++) = (p0[2] + p1[2] + p2[2] + p3[2]) / 4.0;
        }
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError, "Bug: code shouldn't reach here (0)");
        goto ret_err;
    }
    
    // Return output.
    Py_XDECREF(p);
    Py_XDECREF(t);
    return out;
    
ret_err:
    // Clean up code after error.
    Py_XDECREF(p);
    Py_XDECREF(t);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_toBarycentric__doc__,
"Transform one or more 3D points into their barycentric coordinates,\n\
determined by the corner points of a 3D simplex (i.e., a tetrahedron).\n\
\n\
PARAMETERS:\n\
    tcp\n\
        Tetrahedron corner points: a 4*3 array, or a sequence object\n\
        that can be transformed into such an array. Only 2-dimensional\n\
        sequences are allowed.\n\
    p\n\
        An N_pnt * 3 array object of points that must be transformed,\n\
        or a sequence object that can be transformed into such an\n\
        array. Only 1- and 2-dimensional sequences are allowed.\n\
\n\
RETURNS:\n\
    A 2-dimensional array of size N_pnt * 4 giving the barycentric \n\
    coordinates for each point in p.\n\
\n\
TODO:\n\
    The precision of this method appears to be lacking a bit at this \n\
    moment, at least when comparing to results obtained with Matlab.\n\
    (In tetrahedron_test.py, we can only expect tolerance to 3 digits\n\
    when comparing the ratio of STEPS- and Matlab computed values\n\
    to 1.0...!) The sort() doesn't seem to be doing much about this, \n\
    so we'll have to look into what linalg.solve() does, and how to \n\
    improve on that.\n");

static PyObject * tet_toBarycentric(PyObject * self, PyObject * args)
{
    // OVERVIEW:
    // Check input arguments
    //   tcp must have 4 rows and 3 columns, dtype must be double
    //   p must have 3 columns and > 1 row, dtype must be double
    // Copy p0 = [ x0, y0, z0 ] into local variables
    // Prepare 3*3 'basis' matrix [p1-p0; p2-p0; p3 - p0]
    // Prepare output matrix
    // For each point pt:
    //   Set up linear system
    //     Copy 'basis' matrix to first three rows of a
    //     Set fourth row to: pt[0]-x0, pt[1]-y0, pt[2]-z0 
    //     Solve system a -> this gives us (b1,b2,b3) in fourth row
    //     Find fourth barycentric coordinate b0
    //     Copy (b0,b1,b2,b3) to output matrix
    // Return output matrix.

    PyObject * in_tcp = 0;
    PyObject * tcp = 0;
    npy_double * tcpptr;
    PyObject * in_p = 0;
    PyObject * p = 0;
    int npnts;
    npy_double * pptr;
    npy_double * pptr_max;
    double x0, y0, z0;
    double base[9];
    size_t base_nbytes = 9 * sizeof(double);
    npy_intp out_dims[2] = { 0, 4 };
    PyObject * out = 0;
    npy_double * outptr;
    double a[12];
    int info;
    double b1, b2, b3, temp;

    // Parse arguments.
    if (!PyArg_ParseTuple(args, "OO", &in_tcp, &in_p)) goto ret_err;
    
    tcp = PyArray_ContiguousFromAny(in_tcp, PyArray_DOUBLE, 2, 2);
    if (tcp == 0) goto ret_err;
    if (PyArray_DIM(tcp, 0) != 4)
    {
        PyErr_Format(PyExc_ValueError, "Array tcp must have 4 rows");
        goto ret_err;
    }
    if (PyArray_DIM(tcp, 1) != 3)
    {
        PyErr_Format(PyExc_ValueError, "Array tcp must have 3 columns");
        goto ret_err;
    }
    tcpptr = (npy_double*)PyArray_DATA(tcp);
    if (tcpptr == 0)
    {
        PyErr_Format(PyExc_RuntimeError, "Cannot access ndarray for tcp");
        goto ret_err;
    }
    
    p = PyArray_ContiguousFromAny(in_p, PyArray_DOUBLE, 1, 2);
    if (p == 0) goto ret_err;
    if (PyArray_NDIM(p) == 1)
    {
        if ((PyArray_DIM(p,0) % 3) != 0)
        {
            PyErr_Format(PyExc_ValueError, 
               "1-dimensional array p must be multiple of 3 in size");
            goto ret_err;
        }
        npnts = PyArray_DIM(p, 0) / 3;
    }
    else
    {
        npnts = PyArray_DIM(p, 0);
        if (PyArray_DIM(p, 1) != 3)
        {
            PyErr_Format(PyExc_ValueError, "Array p must have 3 columns");
                goto ret_err;
        }
    }
    if (npnts < 1)
    {
        PyErr_Format(PyExc_ValueError, "Error in number of points specified");
        goto ret_err;
    }
    pptr = (npy_double*)PyArray_DATA(p);
    if (pptr == 0)
    {
        PyErr_Format(PyExc_RuntimeError, "Cannot access ndarray for p");
        goto ret_err;
    }
    pptr_max = pptr + (3 * npnts);
    
    // Copy first tet corner point in local variables.
    x0 = tcpptr[0];
    y0 = tcpptr[1];
    z0 = tcpptr[2];
    // Prepase base matrix.
    base[0] = tcpptr[ 3] - x0;
    base[1] = tcpptr[ 4] - y0;
    base[2] = tcpptr[ 5] - z0;
    base[3] = tcpptr[ 6] - x0;
    base[4] = tcpptr[ 7] - y0;
    base[5] = tcpptr[ 8] - z0;
    base[6] = tcpptr[ 9] - x0;
    base[7] = tcpptr[10] - y0;
    base[8] = tcpptr[11] - z0;
    
    // Prepare output matrix. 
    CREATE_OUTPUT_NDARRAY_DOUBLE(2, npnts)
    
    // Loop over points, compute barycentric coords.
    while (pptr < pptr_max)
    {
        memcpy(a, base, base_nbytes);
        a[ 9] = *(pptr++) - x0;
        a[10] = *(pptr++) - y0;
        a[11] = *(pptr++) - z0;
        info = linsolve(3, 1, a);
        if (info != 0)
        {
            PyErr_Format(PyExc_ZeroDivisionError, "Tetrahedron is degenerate");
            goto ret_err;
        }
        b1 = a[ 9];
        b2 = a[10];
        b3 = a[11];
        if (b2 < b1)
        {
            temp = b2;
            b2 = b1;
            b1 = temp;
        }
        if (b3 < b2)
        {
            temp = b3;
            b3 = b2;
            b2 = temp;
            if (b2 < b1)
            {
                temp = b2;
                b2 = b1;
                b1 = temp;
            }
        }
        *(outptr++) = 1.0 - (b1 + b2 + b3);
        *(outptr++) = b1;
        *(outptr++) = b2;
        *(outptr++) = b3;
    }
    
    // Return. 
    Py_XDECREF(tcp);
    Py_XDECREF(p);
    return out;    
    
ret_err:
    // Return with error.
    Py_XDECREF(tcp);
    Py_XDECREF(p);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_inside__doc__,
"Test whether one or more 3D points are inside a tetrahedron.\n\
\n\
PARAMETERS:\n\
    tcp\n\
        Tetrahedron corner points: a 4*3 array, or a sequence object\n\
        that can be transformed into such an array.\n\
    p\n\
        An N_pnt * 3 array object of points that must be transformed,\n\
        or a sequence object that can be transformed into such an\n\
        array.\n\
        \n\
RETURNS:\n\
    A 1-dimensional N_pnt-sized array of boolean values.\n");

static PyObject * tet_inside(PyObject * self, PyObject * args)
{
    // OVERVIEW:
    // Call tet_toBarycentric to get an array of barycentric coordinates.
    // If not successful, return with error.
    // Else, create a boolean 1D output array of size N_pnts.
    // For each point:
    //   if all 4 barycentric coordinates are positive: inside
    //   else, outside
    // Clean up barycentric coordinates.
    // Return output array.

    PyObject * bary = 0;
    npy_double * baryptr;
    npy_intp npnts;
    npy_intp out_dims[1];
    PyObject * out = 0;
    npy_bool * outptr;
    npy_bool * outptr_max;
    
    // Fetch barycentric coords.
    bary = tet_toBarycentric(0, args);
    if (bary == 0) goto ret_err;
    baryptr = (npy_double*)PyArray_DATA(bary);
    if (baryptr == 0)
    {
        PyErr_Format(PyExc_RuntimeError, "Cannot access barycenter array");
        goto ret_err;
    }
    npnts = PyArray_DIM(bary, 0);
    
    // Prepare output array.
    out_dims[0] = npnts;
    out = PyArray_SimpleNew(1, out_dims, PyArray_BOOL);
    if (out == 0) goto ret_err;
    outptr = (npy_bool*)PyArray_DATA(out);
    if (outptr == 0)
    {
        PyErr_Format(PyExc_RuntimeError, "Cannot access result ndarray");
        goto ret_err;
    }
    outptr_max = outptr + npnts;
    
    // Check each point for being inside.
    while (outptr < outptr_max)
    {
        if (*(baryptr++) < 0.0)
        {
            *(outptr++) = NPY_FALSE;
            continue;
        }
        if (*(baryptr++) < 0.0)
        {
            *(outptr++) = NPY_FALSE;
            continue;
        }
        if (*(baryptr++) < 0.0)
        {
            *(outptr++) = NPY_FALSE;
            continue;
        }
        if (*(baryptr++) < 0.0)
        {
            *(outptr++) = NPY_FALSE;
            continue;
        }
        *(outptr++) = NPY_TRUE;
    }

    // Return output array.
    Py_XDECREF(bary);
    return out;
    
ret_err:
    // Return on error.
    Py_XDECREF(bary);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_ranpnt__doc__,
"Generate a number of random points in a tetrahedron.\n\
\n\
Currently not implemented!!!\n\
\n\
The default number of random points generated is 1.\n\
\n\
PARAMETERS:\n\
    r\n\
        A random number generator (from the steps.rng\n\
        interface).\n\
    tcp\n\
        Tetrahedron corner points (a 4*3 array).\n\
    num\n\
        The number of random points to generate (default = 1).\n\
\n\
RETURNS:\n\
    A num*3 Numpy array of tetrahedron-bounded random points.\n\
\n\
RAISES:\n\
    ---\n");

// SOURCE:
// http://vcg.isti.cnr.it/activities/geometryegraphics/pointintetraedro.html

static PyObject * tet_ranpnt(PyObject * self, PyObject * args)
{   
    PyObject * in_rng = 0;
    void * rng1;
    srng::RNG * rng;
    PyObject * in_tcp = 0;
    PyObject * tcp = 0;
    npy_double * tcpptr;
    int npnts = 1;
    npy_intp out_dims[2] = { 0, 3 };
    PyObject * out = 0;
    npy_double * outptr;
    npy_double * outptr_max;
    
    // Parse input.
    if (!PyArg_ParseTuple(args, "OO|i", &in_rng, &in_tcp, &npnts)) 
        goto ret_err;
    rng1 = getSWIGptr(in_rng);
    if (rng1 == 0) 
    {
        PyErr_Format(PyExc_ValueError, 
            "Can't extract STEPS rng from parameter");
        goto ret_err;
    }
    rng = static_cast<srng::RNG*>(rng1);
    tcp = PyArray_ContiguousFromAny(in_tcp, PyArray_DOUBLE, 2, 2);
    if (tcp == 0) goto ret_err;
    if (PyArray_DIM(tcp, 0) != 4)
    {
        PyErr_Format(PyExc_ValueError, "Array tcp must have 4 rows");
        goto ret_err;
    }
    if (PyArray_DIM(tcp, 1) != 3)
    {
        PyErr_Format(PyExc_ValueError, "Array tcp must have 3 columns");
        goto ret_err;
    }
    tcpptr = (npy_double*)PyArray_DATA(tcp);
    if (tcpptr == 0)
    {
        PyErr_Format(PyExc_RuntimeError, "Cannot access ndarray for tcp");
        goto ret_err;
    }
    if (npnts <= 0)
    {
        PyErr_Format(PyExc_ValueError, 
            "Invalid number of points (%d) requested", npnts);
        goto ret_err;
    }
    
    // Prepare output array.
    CREATE_OUTPUT_NDARRAY_DOUBLE(2, npnts)
    outptr_max = outptr + (npnts * 3);

    // Generate points.
    while (outptr < outptr_max)
    {
        // double alpha = pow(rng->getUnfEE(), 1.0/3.0);
        // double alpha1 = 1.0 - alpha;
        // double p0x = alpha1 * tcpptr[ 0];
        // double p0y = alpha1 * tcpptr[ 1];
        // double p0z = alpha1 * tcpptr[ 2];
        // double p01x = p0x + (alpha * tcpptr[ 3]);
        // double p01y = p0y + (alpha * tcpptr[ 4]);
        // double p01z = p0z + (alpha * tcpptr[ 5]);
        // double p02x = p0x + (alpha * tcpptr[ 6]);
        // double p02y = p0y + (alpha * tcpptr[ 7]);
        // double p02z = p0z + (alpha * tcpptr[ 8]);
        // double p03x = p0x + (alpha * tcpptr[ 9]);
        // double p03y = p0y + (alpha * tcpptr[10]);
        // double p03z = p0z + (alpha * tcpptr[11]);
        // alpha = sqrt(rng->getUnfEE()); 
        // alpha1 = 1.0 - alpha;
        // p01x *= alpha1;
        // p01y *= alpha1;
        // p01z *= alpha1;
        // double p12x = p01x + (alpha * p02x); 
        // double p12y = p01y + (alpha * p02y);
        // double p12z = p01z + (alpha * p02z);
        // double p13x = p01x + (alpha * p03x); 
        // double p13y = p01y + (alpha * p03y);
        // double p13z = p01z + (alpha * p03z);
        // alpha = rng->getUnfEE();
        // alpha1 = 1.0 - alpha;
        // *(outptr++) = (alpha * p12x) + (alpha1 * p13x);
        // *(outptr++) = (alpha * p12y) + (alpha1 * p13y);
        // *(outptr++) = (alpha * p12z) + (alpha1 * p13z);
        
        double s = rng->getUnfEE();
        double t = rng->getUnfEE();
        double u = rng->getUnfEE();
        if (s+t>1.0) 
        {
            // cut'n fold the cube into a prism
            s = 1.0 - s;
            t = 1.0 - t;
        }
        if (t+u>1.0) 
        { 
            // cut'n fold the prism into a tetrahedron
            double tmp = u;
            u = 1.0 - s - t;
            t = 1.0 - tmp;
        } 
        else if (s+t+u>1.0)
        {
            double tmp = u;
            u = s + t + u - 1.0;
            s = 1 - t - tmp;
        }
        // a,s,t,u are the barycentric coordinates of the random point.
        double a = 1 - s - t - u;
        *(outptr++) = (a * tcpptr[ 0]) + (s * tcpptr[ 3]) + 
                      (t * tcpptr[ 6]) + (u * tcpptr[ 9]);
        *(outptr++) = (a * tcpptr[ 1]) + (s * tcpptr[ 4]) + 
                      (t * tcpptr[ 7]) + (u * tcpptr[10]);
        *(outptr++) = (a * tcpptr[ 2]) + (s * tcpptr[ 5]) + 
                      (t * tcpptr[ 8]) + (u * tcpptr[11]);
    }
    
    // Return normally.
    Py_XDECREF(tcp);
    return out;

ret_err:
    // Clean up and report error.
    Py_XDECREF(tcp);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_circumsphere__doc__,
"Compute the circumsphere for one or more tetrahedrons.\n\
\n\
The circumsphere of a tetrahedron passes through all 4 corner points.\n\
It can be found by solving the linear system A.X = B, with A given by:\n\
\n\
    |  x2 - x1   x3 - x1   x4 - x1  |\n\
    |  y2 - y1   y3 - y1   y4 - y1  |\n\
    |  z2 - z1   z3 - z1   z4 - z1  |\n\
\n\
and B:\n\
\n\
    | (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 |\n\
    | (x3-x1)^2 + (y3-y1)^2 + (z3-z1)^2 |\n\
    | (x4-x1)^2 + (y4-y1)^2 + (z4-z1)^2 |\n\
\n\
The radius is computed using the solution (x,y,z) to the linear system:\n\
\n\
    r = 1/2 * sqrt(x^2 + y^2 + z^2)\n\
\n\
The coordinates of the sphere's center:\n\
\n\
        | x1 |         | x |\n\
    C = | y1 | + 1/2 * | y |\n\
        | z1 |         | z |\n\
\n\
PARAMETERS:\n\
    p\n\
        An N_pnt * 3 array of points, or a sequence object\n\
        that can be transformed into such an array.\n\
    t\n\
        An optional argument.\n\
        If specified, it is interpreted as an N_tet * 4 array of\n\
        integer indices into p (or any nested sequence object that\n\
        can be converted to such an array).\n\
        If not specified, N_pnt must be a multiple of four, meaning\n\
        that each consecutive block of four rows in p defines a\n\
        tetrahedron.\n\
\n\
RETURNS:\n\
    A N_tet * 4 array of floating point values. Columns 0 to 2\n\
    are the coordinates of the circumsphere's centers, column 3\n\
    are the radii.\n\
\n\
NOTES:\n\
    Even though this is the standard way of computing the properties\n\
    of a circumsphere, imo this has one main problem: it involves a \n\
    lot of subtraction, at least in the part where the linear system \n\
    is being set up. This might have an impact on the precision of the \n\
    result. Maybe there are better ways of doing it? (Since this method\n\
    will be mostly used in quality measures, i.e. outside of inner \n\
    loops or responsive GUI code, efficiency is of less importance.)\n\
\n\
    In general: results appear to become worse with increased \n\
    coplanarity.\n\
\n\
SEE ALSO:\n\
    steps.math.tetrahedron.circumradius()\n\
    steps.math.tetrahedron.circumradius2()\n");

static PyObject * tet_circumsphere(PyObject * self, PyObject * args)
{
    // OVERVIEW:
    // Parse input arguments
    // Prepare output array (2D ntets// 4, double)
    // Perform volume computation
    //   If case1flag:
    //     Loop over all elements of p in blocks of 4*3=12 doubles
    //       Compute circumradius & center
    //   If case2flag:
    //     Loop over each row of t
    //       Check index of each column in each strip
    //       Compute circumradius & center 
    // Return output
    
    PyObject * out = 0;
    npy_intp out_dims[2] = { 0, 4 };
    npy_double * outptr;
    
    PyObject * p = 0;
    npy_intp npnts = 0;
    npy_double * pptr;
    npy_double * pptr_max;
    PyObject * t = 0;
    npy_intp ntets = 0;
    npy_intp * tptr;
    npy_intp * tptr_max;
    CallTypePT presult = 
        parse_input_CallTypePT(args, &p, &t, npnts, ntets, 
        &pptr, &pptr_max, &tptr, &tptr_max);
    if (presult == CTPT_ERR) goto ret_err;
    
    CREATE_OUTPUT_NDARRAY_DOUBLE(2, ntets)
    
    // Perform circumradius computation.
    if (presult == CTPT_P_NOT_T)
    {
        while (pptr < pptr_max)
        {
            npy_double x0 = *(pptr++);
            npy_double y0 = *(pptr++);
            npy_double z0 = *(pptr++);
            
            double a[12];
            
            npy_double pp = a[ 0] = *(pptr++) - x0;
            npy_double qq = a[ 3] = *(pptr++) - y0;
            npy_double ss = a[ 6] = *(pptr++) - z0;
            a[ 9] = (pp * pp) + (qq * qq) + (ss * ss);
            
            pp = a[ 1] = *(pptr++) - x0;
            qq = a[ 4] = *(pptr++) - y0;
            ss = a[ 7] = *(pptr++) - z0;
            a[10] = (pp * pp) + (qq * qq) + (ss * ss);
            
            pp = a[ 2] = *(pptr++) - x0;
            qq = a[ 5] = *(pptr++) - y0;
            ss = a[ 8] = *(pptr++) - z0;
            a[11] = (pp * pp) + (qq * qq) + (ss * ss);
            
            int info = linsolve(3, 1, a);
            if (info != 0)
            {
                PyErr_Format(PyExc_ZeroDivisionError, 
                    "Tetrahedron is degenerate");
                goto ret_err;
            }
            
            npy_double xx0 = a[ 9];
            npy_double yy0 = a[10];
            npy_double zz0 = a[11];
            
            *(outptr++) = x0 + (0.5 * xx0);
            *(outptr++) = y0 + (0.5 * yy0);
            *(outptr++) = z0 + (0.5 * zz0);
            *(outptr++) = 0.5 * sqrt((xx0 * xx0) + (yy0 * yy0) + (zz0 * zz0));
        }
    }
    else if (presult == CTPT_P_AND_T)
    {
        while (tptr < tptr_max)
        {
            npy_double a[12];
            
            npy_intp tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            npy_double * pptr2 = pptr + (3 * tidx);
            npy_double x0 = *(pptr2++);
            npy_double y0 = *(pptr2++);
            npy_double z0 = *pptr2;
            
            tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            pptr2 = pptr + (3 * tidx);
            npy_double pp = a[ 0] = *(pptr2++) - x0;
            npy_double qq = a[ 3] = *(pptr2++) - y0;
            npy_double ss = a[ 6] = *pptr2 - z0;
            a[ 9] = (pp * pp) + (qq * qq) + (ss * ss);
            
            tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            pptr2 = pptr + (3 * tidx);
            pp = a[ 1] = *(pptr2++) - x0;
            qq = a[ 4] = *(pptr2++) - y0;
            ss = a[ 7] = *pptr2 - z0;
            a[10] = (pp * pp) + (qq * qq) + (ss * ss);
                        
            tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            pptr2 = pptr + (3 * tidx);
            pp = a[ 2] = *(pptr2++) - x0;
            qq = a[ 5] = *(pptr2++) - y0;
            ss = a[ 8] = *pptr2 - z0;
            a[11] = (pp * pp) + (qq * qq) + (ss * ss);
            
            int info = linsolve(3, 1, a);
            if (info != 0)
            {
                PyErr_Format(PyExc_ZeroDivisionError, 
                    "Tetrahedron is degenerate");
                goto ret_err;
            }
            
            npy_double xx0 = a[ 9];
            npy_double yy0 = a[10];
            npy_double zz0 = a[11];
            
            *(outptr++) = x0 + (0.5 * xx0);
            *(outptr++) = y0 + (0.5 * yy0);
            *(outptr++) = z0 + (0.5 * zz0);
            *(outptr++) = 0.5 * sqrt((xx0 * xx0) + (yy0 * yy0) + (zz0 * zz0));
        }
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError, "Bug: code shouldn't reach here (0)");
        goto ret_err;
    }

    // Return output.
    Py_XDECREF(p);
    Py_XDECREF(t);
    return out;
    
ret_err:
    // Signal error.
    Py_XDECREF(p);
    Py_XDECREF(t);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_circumradius__doc__,
"Compute the radius of the circumsphere for one or more tetrahedrons.\n\
\n\
This value can be used for computing tetmesh quality measures.\n\
\n\
PARAMETERS:\n\
    p\n\
        An N_pnt * 3 array of points, or a sequence object\n\
        that can be transformed into such an array.\n\
    t\n\
        An optional argument.\n\
        If specified, it is interpreted as an N_tet * 4 array of\n\
        integer indices into p (or any nested sequence object that\n\
        can be converted to such an array).\n\
        If not specified, N_pnt must be a multiple of four, meaning\n\
        that each consecutive block of four rows in p defines a\n\
        tetrahedron.\n\
\n\
RETURNS:\n\
    A 1D array of N_tet radii.\n\
\n\
SEE ALSO:\n\
    steps.math.tetrahedron.circumsphere()\n\
    steps.math.tetrahedron.circumradius2()\n");

static PyObject * tet_circumradius(PyObject * self, PyObject * args)
{
    // OVERVIEW:
    // Parse input arguments
    // Prepare output array (2D ntets// 4, double)
    // Perform volume computation
    //   If case1flag:
    //     Loop over all elements of p in blocks of 4*3=12 doubles
    //       Compute circumradius
    //   If case2flag:
    //     Loop over each row of t
    //       Check index of each column in each strip
    //       Compute circumradius 
    // Return output
    
    npy_intp tidx;
    npy_double * pptr2;
    PyObject * out = 0;
    npy_intp out_dims[1] = { 0 };
    npy_double * outptr;
    npy_double x0, y0, z0;
    npy_double pp, qq, ss;
    npy_double a[12];
    int info;
    
    PyObject * p = 0;
    npy_intp npnts = 0;
    npy_double * pptr;
    npy_double * pptr_max;
    PyObject * t = 0;
    npy_intp ntets = 0;
    npy_intp * tptr;
    npy_intp * tptr_max;
    CallTypePT presult = 
        parse_input_CallTypePT(args, &p, &t, npnts, ntets, 
        &pptr, &pptr_max, &tptr, &tptr_max);
    if (presult == CTPT_ERR) goto ret_err;
        
    CREATE_OUTPUT_NDARRAY_DOUBLE(1, ntets)
    
    // Perform circumradius computation.
    if (presult == CTPT_P_NOT_T)
    {
        while (pptr < pptr_max)
        {
            x0 = *(pptr++);
            y0 = *(pptr++);
            z0 = *(pptr++);
            
            pp = a[ 0] = *(pptr++) - x0;
            qq = a[ 3] = *(pptr++) - y0;
            ss = a[ 6] = *(pptr++) - z0;
            a[ 9] = (pp * pp) + (qq * qq) + (ss * ss);
            
            pp = a[ 1] = *(pptr++) - x0;
            qq = a[ 4] = *(pptr++) - y0;
            ss = a[ 7] = *(pptr++) - z0;
            a[10] = (pp * pp) + (qq * qq) + (ss * ss);
            
            pp = a[ 2] = *(pptr++) - x0;
            qq = a[ 5] = *(pptr++) - y0;
            ss = a[ 8] = *(pptr++) - z0;
            a[11] = (pp * pp) + (qq * qq) + (ss * ss);
            
            info = linsolve(3, 1, a);
            if (info != 0)
            {
                PyErr_Format(PyExc_ZeroDivisionError, 
                    "Tetrahedron is degenerate");
                goto ret_err;
            }

            *(outptr++) = 
                0.5 * sqrt((a[9]*a[9]) + (a[10]*a[10]) + (a[11]*a[11]));
        }
    }
    else if (presult == CTPT_P_AND_T)
    {
        while (tptr < tptr_max)
        {
            tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            pptr2 = pptr + (3 * tidx);
            x0 = *(pptr2++);
            y0 = *(pptr2++);
            z0 = *pptr2;
            
            tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            pptr2 = pptr + (3 * tidx);
            pp = a[ 0] = *(pptr2++) - x0;
            qq = a[ 3] = *(pptr2++) - y0;
            ss = a[ 6] = *pptr2 - z0;
            a[ 9] = (pp * pp) + (qq * qq) + (ss * ss);
            
            tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            pptr2 = pptr + (3 * tidx);
            pp = a[ 1] = *(pptr2++) - x0;
            qq = a[ 4] = *(pptr2++) - y0;
            ss = a[ 7] = *pptr2 - z0;
            a[10] = (pp * pp) + (qq * qq) + (ss * ss);
                        
            tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            pptr2 = pptr + (3 * tidx);
            pp = a[ 2] = *(pptr2++) - x0;
            qq = a[ 5] = *(pptr2++) - y0;
            ss = a[ 8] = *pptr2 - z0;
            a[11] = (pp * pp) + (qq * qq) + (ss * ss);
            
            info = linsolve(3, 1, a);
            if (info != 0)
            {
                PyErr_Format(PyExc_ZeroDivisionError, 
                    "Tetrahedron is degenerate");
                goto ret_err;
            }
            
            *(outptr++) = 
                0.5 * sqrt((a[9]*a[9]) + (a[10]*a[10]) + (a[11]*a[11]));
        }
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError, "Bug: code shouldn't reach here (0)");
        goto ret_err;
    }

    // Return output.
    Py_XDECREF(p);
    Py_XDECREF(t);
    return out;
    
ret_err:
    // Signal error.
    Py_XDECREF(p);
    Py_XDECREF(t);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_circumradius2__doc__,
"Compute the square radius of the circumsphere for one or more\n\
tetrahedrons.\n\
\n\
PARAMETERS:\n\
    p\n\
        An N_pnt * 3 array of points, or a sequence object\n\
        that can be transformed into such an array.\n\
    t\n\
        An optional argument.\n\
        If specified, it is interpreted as an N_tet * 4 array of\n\
        integer indices into p (or any nested sequence object that\n\
        can be converted to such an array).\n\
        If not specified, N_pnt must be a multiple of four, meaning\n\
        that each consecutive block of four rows in p defines a\n\
        tetrahedron.\n\
\n\
RETURNS:\n\
    A 1D array of N_tet square radii.\n\
\n\
SEE ALSO:\n\
    steps.math.tetrahedron.circumsphere()\n\
    steps.math.tetrahedron.circumradius()\n");

static PyObject * tet_circumradius2(PyObject * self, PyObject * args)
{
    // OVERVIEW:
    // Parse input arguments
    // Prepare output array (2D ntets// 4, double)
    // Perform volume computation
    //   If case1flag:
    //     Loop over all elements of p in blocks of 4*3=12 doubles
    //       Compute circumradius
    //   If case2flag:
    //     Loop over each row of t
    //       Check index of each column in each strip
    //       Compute circumradius 
    // Return output
    
    npy_intp tidx;
    npy_double * pptr2;
    PyObject * out = 0;
    npy_intp out_dims[1] = { 0 };
    npy_double * outptr;
    npy_double x0, y0, z0;
    double pp, qq, ss;
    double a[12];
    int info;
    
    PyObject * p = 0;
    npy_intp npnts = 0;
    npy_double * pptr;
    npy_double * pptr_max;
    PyObject * t = 0;
    npy_intp ntets = 0;
    npy_intp * tptr;
    npy_intp * tptr_max;
    CallTypePT presult = 
        parse_input_CallTypePT(args, &p, &t, npnts, ntets, 
        &pptr, &pptr_max, &tptr, &tptr_max);
    if (presult == CTPT_ERR) goto ret_err;
    
    CREATE_OUTPUT_NDARRAY_DOUBLE(1, ntets)
    
    // Perform circumradius computation.
    if (presult == CTPT_P_NOT_T)
    {
        while (pptr < pptr_max)
        {
            x0 = *(pptr++);
            y0 = *(pptr++);
            z0 = *(pptr++);
            
            pp = a[ 0] = *(pptr++) - x0;
            qq = a[ 3] = *(pptr++) - y0;
            ss = a[ 6] = *(pptr++) - z0;
            a[ 9] = (pp * pp) + (qq * qq) + (ss * ss);
            
            pp = a[ 1] = *(pptr++) - x0;
            qq = a[ 4] = *(pptr++) - y0;
            ss = a[ 7] = *(pptr++) - z0;
            a[10] = (pp * pp) + (qq * qq) + (ss * ss);
            
            pp = a[ 2] = *(pptr++) - x0;
            qq = a[ 5] = *(pptr++) - y0;
            ss = a[ 8] = *(pptr++) - z0;
            a[11] = (pp * pp) + (qq * qq) + (ss * ss);
            
            info = linsolve(3, 1, a);
            if (info != 0)
            {
                PyErr_Format(PyExc_ZeroDivisionError, 
                    "Tetrahedron is degenerate");
                goto ret_err;
            }

            *(outptr++) = 0.25 * ((a[9]*a[9]) + (a[10]*a[10]) + (a[11]*a[11]));
        }
    }
    else if (presult == CTPT_P_AND_T)
    {
        while (tptr < tptr_max)
        {
            tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            pptr2 = pptr + (3 * tidx);
            x0 = *(pptr2++);
            y0 = *(pptr2++);
            z0 = *pptr2;
            
            tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            pptr2 = pptr + (3 * tidx);
            pp = a[ 0] = *(pptr2++) - x0;
            qq = a[ 3] = *(pptr2++) - y0;
            ss = a[ 6] = *pptr2 - z0;
            a[ 9] = (pp * pp) + (qq * qq) + (ss * ss);
            
            tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            pptr2 = pptr + (3 * tidx);
            pp = a[ 1] = *(pptr2++) - x0;
            qq = a[ 4] = *(pptr2++) - y0;
            ss = a[ 7] = *pptr2 - z0;
            a[10] = (pp * pp) + (qq * qq) + (ss * ss);
                        
            tidx = *(tptr++);
            if (tidx >= npnts)
            {
                PyErr_Format(PyExc_IndexError, 
                    "Point index out of range (%d)", tidx);
                goto ret_err;
            }
            pptr2 = pptr + (3 * tidx);
            pp = a[ 2] = *(pptr2++) - x0;
            qq = a[ 5] = *(pptr2++) - y0;
            ss = a[ 8] = *pptr2 - z0;
            a[11] = (pp * pp) + (qq * qq) + (ss * ss);
            
            info = linsolve(3, 1, a);
            if (info != 0)
            {
                PyErr_Format(PyExc_ZeroDivisionError, 
                    "Tetrahedron is degenerate");
                goto ret_err;
            }
            
            *(outptr++) = 0.25 * ((a[9]*a[9]) + (a[10]*a[10]) + (a[11]*a[11]));
        }
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError, "Bug: code shouldn't reach here (0)");
        goto ret_err;
    }

    // Return output.
    Py_XDECREF(p);
    Py_XDECREF(t);
    return out;
    
ret_err:
    // Signal error.
    Py_XDECREF(p);
    Py_XDECREF(t);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_shortestEdge__doc__,
"Compute the length of the shortest edge of one or more tetrahedrons.\n\
\n\
PARAMETERS:\n\
    p\n\
        An N_pnt * 3 array of points, or a sequence object\n\
        that can be transformed into such an array.\n\
    t\n\
        An optional argument.\n\
        If specified, it is interpreted as an N_tet * 4 array of\n\
        integer indices into p (or any nested sequence object that\n\
        can be converted to such an array).\n\
        If not specified, N_pnt must be a multiple of four, meaning\n\
        that each consecutive block of four rows in p defines a\n\
        tetrahedron.\n\
\n\
RETURNS:\n\
    A 1D array of N_tet shortest edge lengths.\n\
\n\
SEE ALSO:\n\
    steps.math.tetrahedron.shortestEdge2()\n\
    steps.math.tetrahedron.longestEdge()\n\
    steps.math.tetrahedron.longestEdge2()\n");

static PyObject * tet_shortestEdge(PyObject * self, PyObject * args)
{
    // OVERVIEW:
    // Parse input arguments
    // Prepare output array (1D ntets, double)
    // Find shortest edge
    //   If case1flag:
    //     Loop over all elements of p in blocks of 4*3=12 doubles
    //      Find shortest (r01, r02, r03, r12, r13, r23) 
    //   If case2flag:
    //     Loop over each row of t
    //       Check index of each column in each strip
    //       Find shortest 
    // Return output
    
    PyObject * out = 0;
    npy_intp out_dims[1] = { 0 };
    npy_double * outptr;
    double dx, dy, dz;

    double ls, ln;
    
    PyObject * p = 0;
    npy_intp npnts = 0;
    npy_double * pptr;
    npy_double * pptr_max;
    PyObject * t = 0;
    npy_intp ntets = 0;
    npy_intp * tptr;
    npy_intp * tptr_max;
    CallTypePT presult = 
        parse_input_CallTypePT(args, &p, &t, npnts, ntets, 
        &pptr, &pptr_max, &tptr, &tptr_max);
    if (presult == CTPT_ERR) goto ret_err;

    CREATE_OUTPUT_NDARRAY_DOUBLE(1, ntets)
    
    if (presult == CTPT_P_NOT_T)
    {
        for (;pptr < pptr_max; pptr += 12)
        {
            // Edge 0: p1 - p0
            dx = pptr[ 3] - pptr[ 0];
            dy = pptr[ 4] - pptr[ 1];
            dz = pptr[ 5] - pptr[ 2];
            ls = (dx * dx) + (dy * dy) + (dz * dz);
            
            // Edge 1: p2 - p0
            dx = pptr[ 6] - pptr[ 0];
            dy = pptr[ 7] - pptr[ 1];
            dz = pptr[ 8] - pptr[ 2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 2: p3 - p0
            dx = pptr[ 9] - pptr[ 0];
            dy = pptr[10] - pptr[ 1];
            dz = pptr[11] - pptr[ 2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 3: p2 - p1
            dx = pptr[ 6] - pptr[ 3];
            dy = pptr[ 7] - pptr[ 4];
            dz = pptr[ 8] - pptr[ 5];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 4: p3 - p1
            dx = pptr[ 9] - pptr[ 3];
            dy = pptr[10] - pptr[ 4];
            dz = pptr[11] - pptr[ 5];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 5: p3 - p2
            dx = pptr[ 9] - pptr[ 6];
            dy = pptr[10] - pptr[ 7];
            dz = pptr[11] - pptr[ 8];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            *(outptr++) = sqrt(ls);
        }
    }
    else if (presult == CTPT_P_AND_T)
    {
        npy_intp tidx;
        npy_double * p0; 
        npy_double * p1;
        npy_double * p2;
        npy_double * p3;
        while (tptr < tptr_max)
        {
            TPTR_TO_4PPTR
            
            // Edge 0: p1 - p0
            dx = p1[0] - p0[0];
            dy = p1[1] - p0[1];
            dz = p1[2] - p0[2];
            ls = (dx * dx) + (dy * dy) + (dz * dz);
            
            // Edge 1: p2 - p0
            dx = p2[0] - p0[0];
            dy = p2[1] - p0[1];
            dz = p2[2] - p0[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 2: p3 - p0
            dx = p3[0] - p0[0];
            dy = p3[1] - p0[1];
            dz = p3[2] - p0[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 3: p2 - p1
            dx = p2[0] - p1[0];
            dy = p2[1] - p1[1];
            dz = p2[2] - p1[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 4: p3 - p1
            dx = p3[0] - p1[0];
            dy = p3[1] - p1[1];
            dz = p3[2] - p1[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 5: p3 - p2
            dx = p3[0] - p2[0];
            dy = p3[1] - p2[1];
            dz = p3[2] - p2[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            *(outptr++) = sqrt(ls);
        }
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError, "Bug: code shouldn't reach here (0)");
        goto ret_err;
    }

    // Return output.
    Py_XDECREF(p);
    Py_XDECREF(t);
    return out;
    
ret_err:
    // Signal error.
    Py_XDECREF(p);
    Py_XDECREF(t);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_shortestEdge2__doc__,
"Compute the square length of the shortest edge of one or more \n\
tetrahedrons.\n\
\n\
PARAMETERS:\n\
    p\n\
        An N_pnt * 3 array of points, or a sequence object\n\
        that can be transformed into such an array.\n\
    t\n\
        An optional argument.\n\
        If specified, it is interpreted as an N_tet * 4 array of\n\
        integer indices into p (or any nested sequence object that\n\
        can be converted to such an array).\n\
        If not specified, N_pnt must be a multiple of four, meaning\n\
        that each consecutive block of four rows in p defines a\n\
        tetrahedron.\n\
\n\
RETURNS:\n\
    A 1D array of N_tet shortest square edge lengths.\n\
\n\
SEE ALSO:\n\
    steps.math.tetrahedron.shortestEdge()\n\
    steps.math.tetrahedron.longestEdge()\n\
    steps.math.tetrahedron.longestEdge2()\n");

static PyObject * tet_shortestEdge2(PyObject * self, PyObject * args)
{
    // OVERVIEW:
    // Parse input arguments
    // Prepare output array (1D ntets, double)
    // Find shortest edge
    //   If case1flag:
    //     Loop over all elements of p in blocks of 4*3=12 doubles
    //      Find square shortest (r01, r02, r03, r12, r13, r23) 
    //   If case2flag:
    //     Loop over each row of t
    //       Check index of each column in each strip
    //       Find square shortest 
    // Return output
    
    PyObject * out = 0;
    npy_intp out_dims[1] = { 0 };
    npy_double * outptr;
    double dx, dy, dz;
    double ls, ln;
    
    PyObject * p = 0;
    npy_intp npnts = 0;
    npy_double * pptr;
    npy_double * pptr_max;
    PyObject * t = 0;
    npy_intp ntets = 0;
    npy_intp * tptr;
    npy_intp * tptr_max;
    CallTypePT presult = 
        parse_input_CallTypePT(args, &p, &t, npnts, ntets, 
        &pptr, &pptr_max, &tptr, &tptr_max);
    if (presult == CTPT_ERR) goto ret_err;
    
    CREATE_OUTPUT_NDARRAY_DOUBLE(1, ntets)
    
    if (presult == CTPT_P_NOT_T)
    {
        for (;pptr < pptr_max; pptr += 12)
        {
            // Edge 0: p1 - p0
            dx = pptr[ 3] - pptr[ 0];
            dy = pptr[ 4] - pptr[ 1];
            dz = pptr[ 5] - pptr[ 2];
            ls = (dx * dx) + (dy * dy) + (dz * dz);
            
            // Edge 1: p2 - p0
            dx = pptr[ 6] - pptr[ 0];
            dy = pptr[ 7] - pptr[ 1];
            dz = pptr[ 8] - pptr[ 2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 2: p3 - p0
            dx = pptr[ 9] - pptr[ 0];
            dy = pptr[10] - pptr[ 1];
            dz = pptr[11] - pptr[ 2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 3: p2 - p1
            dx = pptr[ 6] - pptr[ 3];
            dy = pptr[ 7] - pptr[ 4];
            dz = pptr[ 8] - pptr[ 5];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 4: p3 - p1
            dx = pptr[ 9] - pptr[ 3];
            dy = pptr[10] - pptr[ 4];
            dz = pptr[11] - pptr[ 5];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 5: p3 - p2
            dx = pptr[ 9] - pptr[ 6];
            dy = pptr[10] - pptr[ 7];
            dz = pptr[11] - pptr[ 8];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            *(outptr++) = ls;
        }
    }
    else if (presult == CTPT_P_AND_T)
    {
        npy_intp tidx;
        npy_double * p0; 
        npy_double * p1;
        npy_double * p2;
        npy_double * p3;
        while (tptr < tptr_max)
        {
            TPTR_TO_4PPTR
            
            // Edge 0: p1 - p0
            dx = p1[0] - p0[0];
            dy = p1[1] - p0[1];
            dz = p1[2] - p0[2];
            ls = (dx * dx) + (dy * dy) + (dz * dz);
            
            // Edge 1: p2 - p0
            dx = p2[0] - p0[0];
            dy = p2[1] - p0[1];
            dz = p2[2] - p0[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 2: p3 - p0
            dx = p3[0] - p0[0];
            dy = p3[1] - p0[1];
            dz = p3[2] - p0[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 3: p2 - p1
            dx = p2[0] - p1[0];
            dy = p2[1] - p1[1];
            dz = p2[2] - p1[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 4: p3 - p1
            dx = p3[0] - p1[0];
            dy = p3[1] - p1[1];
            dz = p3[2] - p1[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            // Edge 5: p3 - p2
            dx = p3[0] - p2[0];
            dy = p3[1] - p2[1];
            dz = p3[2] - p2[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln < ls) ls = ln;
            
            *(outptr++) = ls;
        }
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError, "Bug: code shouldn't reach here (0)");
        goto ret_err;
    }

    // Return output.
    Py_XDECREF(p);
    Py_XDECREF(t);
    return out;
    
ret_err:
    // Signal error.
    Py_XDECREF(p);
    Py_XDECREF(t);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_longestEdge__doc__,
"Compute the length of the longest edge of one or more tetrahedrons.\n\
\n\
PARAMETERS:\n\
    p\n\
        An N_pnt * 3 array of points, or a sequence object\n\
        that can be transformed into such an array.\n\
    t\n\
        An optional argument.\n\
        If specified, it is interpreted as an N_tet * 4 array of\n\
        integer indices into p (or any nested sequence object that\n\
        can be converted to such an array).\n\
        If not specified, N_pnt must be a multiple of four, meaning\n\
        that each consecutive block of four rows in p defines a\n\
        tetrahedron.\n\
\n\
RETURNS:\n\
    A 1D array of N_tet longest edge lengths.\n\
\n\
SEE ALSO:\n\
    steps.math.tetrahedron.shortestEdge()\n\
    steps.math.tetrahedron.shortestEdge2()\n\
    steps.math.tetrahedron.longestEdge2()\n");

static PyObject * tet_longestEdge(PyObject * self, PyObject * args)
{
    // OVERVIEW:
    // Parse input arguments
    // Prepare output array (1D ntets, double)
    // Find longest edge
    //   If case1flag:
    //     Loop over all elements of p in blocks of 4*3=12 doubles
    //      Find longest (r01, r02, r03, r12, r13, r23) 
    //   If case2flag:
    //     Loop over each row of t
    //       Check index of each column in each strip
    //       Find longest 
    // Return output
    
    PyObject * out = 0;
    npy_intp out_dims[1] = { 0 };
    npy_double * outptr;
    double dx, dy, dz;
    double ls, ln;
    
    PyObject * p = 0;
    npy_intp npnts = 0;
    npy_double * pptr;
    npy_double * pptr_max;
    PyObject * t = 0;
    npy_intp ntets = 0;
    npy_intp * tptr;
    npy_intp * tptr_max;
    CallTypePT presult = 
        parse_input_CallTypePT(args, &p, &t, npnts, ntets, 
        &pptr, &pptr_max, &tptr, &tptr_max);
    if (presult == CTPT_ERR) goto ret_err;
    
    CREATE_OUTPUT_NDARRAY_DOUBLE(1, ntets)
    
    if (presult == CTPT_P_NOT_T)
    {
        for (;pptr < pptr_max; pptr += 12)
        {
            // Edge 0: p1 - p0
            dx = pptr[ 3] - pptr[ 0];
            dy = pptr[ 4] - pptr[ 1];
            dz = pptr[ 5] - pptr[ 2];
            ls = (dx * dx) + (dy * dy) + (dz * dz);
            
            // Edge 1: p2 - p0
            dx = pptr[ 6] - pptr[ 0];
            dy = pptr[ 7] - pptr[ 1];
            dz = pptr[ 8] - pptr[ 2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 2: p3 - p0
            dx = pptr[ 9] - pptr[ 0];
            dy = pptr[10] - pptr[ 1];
            dz = pptr[11] - pptr[ 2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 3: p2 - p1
            dx = pptr[ 6] - pptr[ 3];
            dy = pptr[ 7] - pptr[ 4];
            dz = pptr[ 8] - pptr[ 5];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 4: p3 - p1
            dx = pptr[ 9] - pptr[ 3];
            dy = pptr[10] - pptr[ 4];
            dz = pptr[11] - pptr[ 5];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 5: p3 - p2
            dx = pptr[ 9] - pptr[ 6];
            dy = pptr[10] - pptr[ 7];
            dz = pptr[11] - pptr[ 8];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            *(outptr++) = sqrt(ls);
        }
    }
    else if (presult == CTPT_P_AND_T)
    {
        npy_intp tidx;
        npy_double * p0; 
        npy_double * p1;
        npy_double * p2;
        npy_double * p3;
        while (tptr < tptr_max)
        {
            TPTR_TO_4PPTR
            
            // Edge 0: p1 - p0
            dx = p1[0] - p0[0];
            dy = p1[1] - p0[1];
            dz = p1[2] - p0[2];
            ls = (dx * dx) + (dy * dy) + (dz * dz);
            
            // Edge 1: p2 - p0
            dx = p2[0] - p0[0];
            dy = p2[1] - p0[1];
            dz = p2[2] - p0[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 2: p3 - p0
            dx = p3[0] - p0[0];
            dy = p3[1] - p0[1];
            dz = p3[2] - p0[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 3: p2 - p1
            dx = p2[0] - p1[0];
            dy = p2[1] - p1[1];
            dz = p2[2] - p1[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 4: p3 - p1
            dx = p3[0] - p1[0];
            dy = p3[1] - p1[1];
            dz = p3[2] - p1[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 5: p3 - p2
            dx = p3[0] - p2[0];
            dy = p3[1] - p2[1];
            dz = p3[2] - p2[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            *(outptr++) = sqrt(ls);
        }
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError, "Bug: code shouldn't reach here (0)");
        goto ret_err;
    }

    // Return output.
    Py_XDECREF(p);
    Py_XDECREF(t);
    return out;
    
ret_err:
    // Signal error.
    Py_XDECREF(p);
    Py_XDECREF(t);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PyDoc_STRVAR(tet_longestEdge2__doc__,
"Compute the length of the longest square edge of one or more \n\
tetrahedrons.\n\
\n\
PARAMETERS:\n\
    p\n\
        An N_pnt * 3 array of points, or a sequence object\n\
        that can be transformed into such an array.\n\
    t\n\
        An optional argument.\n\
        If specified, it is interpreted as an N_tet * 4 array of\n\
        integer indices into p (or any nested sequence object that\n\
        can be converted to such an array).\n\
        If not specified, N_pnt must be a multiple of four, meaning\n\
        that each consecutive block of four rows in p defines a\n\
        tetrahedron.\n\
\n\
RETURNS:\n\
    A 1D array of N_tet longest square edge lengths.\n\
\n\
SEE ALSO:\n\
    steps.math.tetrahedron.shortestEdge()\n\
    steps.math.tetrahedron.shortestEdge2()\n\
    steps.math.tetrahedron.longestEdge()\n");

static PyObject * tet_longestEdge2(PyObject * self, PyObject * args)
{
    // OVERVIEW:
    // Parse input arguments
    // Prepare output array (1D ntets, double)
    // Find longest edge
    //   If case1flag:
    //     Loop over all elements of p in blocks of 4*3=12 doubles
    //      Find square longest (r01, r02, r03, r12, r13, r23) 
    //   If case2flag:
    //     Loop over each row of t
    //       Check index of each column in each strip
    //       Find square longest 
    // Return output

    PyObject * out = 0;
    npy_intp out_dims[1] = { 0 };
    npy_double * outptr;
    double dx, dy, dz;
    double ls, ln;
    
    PyObject * p = 0;
    npy_intp npnts = 0;
    npy_double * pptr;
    npy_double * pptr_max;
    PyObject * t = 0;
    npy_intp ntets = 0;
    npy_intp * tptr;
    npy_intp * tptr_max;
    CallTypePT presult = 
        parse_input_CallTypePT(args, &p, &t, npnts, ntets, 
        &pptr, &pptr_max, &tptr, &tptr_max);
    if (presult == CTPT_ERR) goto ret_err;

    CREATE_OUTPUT_NDARRAY_DOUBLE(1, ntets)
    
    if (presult == CTPT_P_NOT_T)
    {
        for (;pptr < pptr_max; pptr += 12)
        {
            // Edge 0: p1 - p0
            dx = pptr[ 3] - pptr[ 0];
            dy = pptr[ 4] - pptr[ 1];
            dz = pptr[ 5] - pptr[ 2];
            ls = (dx * dx) + (dy * dy) + (dz * dz);
            
            // Edge 1: p2 - p0
            dx = pptr[ 6] - pptr[ 0];
            dy = pptr[ 7] - pptr[ 1];
            dz = pptr[ 8] - pptr[ 2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 2: p3 - p0
            dx = pptr[ 9] - pptr[ 0];
            dy = pptr[10] - pptr[ 1];
            dz = pptr[11] - pptr[ 2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 3: p2 - p1
            dx = pptr[ 6] - pptr[ 3];
            dy = pptr[ 7] - pptr[ 4];
            dz = pptr[ 8] - pptr[ 5];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 4: p3 - p1
            dx = pptr[ 9] - pptr[ 3];
            dy = pptr[10] - pptr[ 4];
            dz = pptr[11] - pptr[ 5];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 5: p3 - p2
            dx = pptr[ 9] - pptr[ 6];
            dy = pptr[10] - pptr[ 7];
            dz = pptr[11] - pptr[ 8];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            *(outptr++) = ls;
        }
    }
    else if (presult == CTPT_P_AND_T)
    {
        npy_intp tidx;
        npy_double * p0; 
        npy_double * p1;
        npy_double * p2;
        npy_double * p3;
        while (tptr < tptr_max)
        {
            TPTR_TO_4PPTR
            
            // Edge 0: p1 - p0
            dx = p1[0] - p0[0];
            dy = p1[1] - p0[1];
            dz = p1[2] - p0[2];
            ls = (dx * dx) + (dy * dy) + (dz * dz);
            
            // Edge 1: p2 - p0
            dx = p2[0] - p0[0];
            dy = p2[1] - p0[1];
            dz = p2[2] - p0[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 2: p3 - p0
            dx = p3[0] - p0[0];
            dy = p3[1] - p0[1];
            dz = p3[2] - p0[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 3: p2 - p1
            dx = p2[0] - p1[0];
            dy = p2[1] - p1[1];
            dz = p2[2] - p1[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 4: p3 - p1
            dx = p3[0] - p1[0];
            dy = p3[1] - p1[1];
            dz = p3[2] - p1[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            // Edge 5: p3 - p2
            dx = p3[0] - p2[0];
            dy = p3[1] - p2[1];
            dz = p3[2] - p2[2];
            ln = (dx * dx) + (dy * dy) + (dz * dz);
            if (ln > ls) ls = ln;
            
            *(outptr++) = ls;
        }
    }
    else
    {
        PyErr_Format(PyExc_RuntimeError, "Bug: code shouldn't reach here (0)");
        goto ret_err;
    }

    // Return output.
    Py_XDECREF(p);
    Py_XDECREF(t);
    return out;
    
ret_err:
    // Signal error.
    Py_XDECREF(p);
    Py_XDECREF(t);
    Py_XDECREF(out);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

static PyMethodDef tet_methods[] = 
{
    { "vol", (PyCFunction)tet_vol, METH_VARARGS, tet_vol__doc__ },
    { "barycenter", (PyCFunction)tet_barycenter, METH_VARARGS, tet_barycenter__doc__ },
    { "toBarycentric", (PyCFunction)tet_toBarycentric, METH_VARARGS, tet_toBarycentric__doc__ },
    { "inside", (PyCFunction)tet_inside, METH_VARARGS, tet_inside__doc__ },
    { "ranpnt", (PyCFunction)tet_ranpnt, METH_VARARGS, tet_ranpnt__doc__ },
    { "circumsphere", (PyCFunction)tet_circumsphere, METH_VARARGS, tet_circumsphere__doc__ },
    { "circumradius", (PyCFunction)tet_circumradius, METH_VARARGS, tet_circumradius__doc__ },
    { "circumradius2", (PyCFunction)tet_circumradius2, METH_VARARGS, tet_circumradius2__doc__ },
    { "shortestEdge", (PyCFunction)tet_shortestEdge, METH_VARARGS, tet_shortestEdge__doc__ },
    { "shortestEdge2", (PyCFunction)tet_shortestEdge2, METH_VARARGS, tet_shortestEdge2__doc__ },
    { "longestEdge", (PyCFunction)tet_longestEdge, METH_VARARGS, tet_longestEdge__doc__ },
    { "longestEdge2", (PyCFunction)tet_longestEdge2, METH_VARARGS, tet_longestEdge2__doc__ },
    { NULL, NULL, 0, NULL }
};

PyDoc_STRVAR(tet__doc__,
"A variety of auxiliary functions for dealing with tetrahedrons.\n\
\n\
The following sets of functionality are offered by this module:\n\
    - Computing tetrahedron volumes\n\
    - Transformation to/from barycentric coordinates\n\
    - Computing circumspheres\n\
    - Conputing shortest/longest edge lengths\n\
\n\
All of this functionality has been implemented in C++ for speed (an\n\
older implementation in Python proved to be a bottleneck while working\n\
with meshes).\n\
\n\
SEE ALSO:\n\
    steps.math.tetrahedron_test\n\
    steps.math.triangle\n");

extern "C"
{
    PyMODINIT_FUNC inittetrahedron(void);
}

PyMODINIT_FUNC inittetrahedron(void)
{
    Py_InitModule3("tetrahedron", tet_methods, tet__doc__);
    import_array();
}

////////////////////////////////////////////////////////////////////////////////

// END
