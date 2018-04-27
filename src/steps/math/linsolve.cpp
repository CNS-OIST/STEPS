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

// STL headers
#include <cmath>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/linsolve.hpp"

////////////////////////////////////////////////////////////////////////////////

int steps::math::linsolve(int n, int rhs_num, double a[])
{
    // SOURCE: I forgot. Probably Numerical Recipes.

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

// END
