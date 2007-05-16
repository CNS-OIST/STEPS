//
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2006 Stefan Wils.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//

////////////////////////////////////////////////////////////////////////////////

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
        //std::cout << a[j + j * n] << std::endl;

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
                //std::cout << a[i+j*n] << std::endl;
            }
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
