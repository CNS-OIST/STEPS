////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cmath>
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/linalg.hpp>

////////////////////////////////////////////////////////////////////////////////

int steps::math::linsolve(uint n, uint rhs_num, double a[])
{
    // Precompute n+ rhs_num
    uint n_plus_rhs_num = n + rhs_num;

    // Loop over all rows.
    for (uint j = 0; j < n; ++j)
    {
        //  Choose a pivot row: first we select j.
        uint ipivot = j;
        double apivot = a[j + j * n];
        // But we really want the largest.
        for (uint i = j + 1; i < n; ++i)
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
        for (uint i = 0; i < n_plus_rhs_num; ++i)
        {
            double temp          = a[ipivot + i * n];
            a[ipivot + i * n]    = a[j + i * n];
            a[j + i * n]         = temp;
        }

        // a[j,j] becomes 1.
        //a[j + j * n] = 1.0;
        for (uint k = j; k < n_plus_rhs_num; ++k)
        {
            a[j + k * n] = a[j + k * n] / apivot;
        }

        // a[i,j] becomes 0.
        for (uint i = 0; i < n; ++i)
        {
            if (i != j)
            {
                double factor = a[i + j * n];
                //a[i + j * n] = 0.0;
                for (uint k = j; k < n_plus_rhs_num; ++k)
                {
                    a[i + k * n] = a[i + k * n] - factor * a[j + k * n];
                }
            }
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
