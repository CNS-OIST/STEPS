/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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


#include <algorithm>
#include <cmath>
#include <utility>

#include "steps/common.h"
#include "steps/solver/efield/bdsystem.hpp"

#include "easylogging++.h"

#include <iostream>

namespace steps {
namespace solver {
namespace efield {

// LU and back substitution implementations are from Numerical Recipes, 2nd ed., section 2.4.

inline void swap_row(double *u,double *v,int n) {
    for (int j=0;j<n;++j) std::swap(u[j],v[j]);
}

void BDSystem::solve()
{
    constexpr double TINY = 1.0e-20;

    auto n = pN;
    auto h = pHalfBW;
    auto w = 2 * h + 1;
    double *a = pA.data(); // will hold U from LU decomposition
    double *l = h > 0? pL.data(): nullptr;

    // w is equivalent to pBW in the prop. 'a' dimensions are nrow ('n') x 'w'

    // 1. LU decomposition.

    // Shift initial h rows left in compact representation, overwriting the
    // unused elements and padding with zeroes on the right.
    double *ai = a; // ai corresponds to ith row of a

    for (auto i = 0u; i < h; ++i)
    {
        auto p = h - i;

        for (auto j = p; j < w; ++j) {
            ai[j-p] = ai[j];
}

        for (auto j = w - p; j < w; ++j) {
            ai[j] = 0.0;
}

        ai += w;
    }

    double *ak = a;
    double *lk = l;
    for (auto k = 0u; k < n; ++k)
    {
        // find pivot in the following h rows
        auto p = std::min(k + h + 1, n);
 
        double dum = std::abs(ak[0]);
        auto ipiv = k;
        double *aj = ak;
        for (auto j = k + 1; j < p; ++j)
        {
            aj += w;
            double abs_aj = std::abs(aj[0]);
            if (abs_aj > dum)
            {
                dum = abs_aj;
                ipiv = j;
            }
        }
        if (dum == 0.) { ak[0] = TINY;
}

        pp[k] = ipiv;
        if (ipiv != k) {
            swap_row(ak, a + ipiv * w, w);
        }

        // perform eliminiation
        double *aij = ak+w; // aij == a + (k+1)*w + 0
        double *lki = lk;

        for (auto i = k+1; i < p; ++i)
        {
            dum = aij[0]/ak[0];

            // lki == l + k*h + (i-k-1)
            *lki++ = dum;

            for (auto j = 1u; j < w; ++j)
            {
                ++aij; // aij == a + i*w + j
                aij[-1] =  aij[0] - ak[j]*dum;
            }
            aij[0] = 0.0; // aij == a + i*w + (w-1)
            ++aij; // aij == a + (i+1)*w 
        }
        ak += w;
        lk += h;
    }

    // 2. Forward substitution, b into x.
    std::copy(pb.begin(),pb.end(),px.begin());
    double *x = px.data();
    double *xk = x;
    lk = l;
    ak = a;
    for (auto k = 0u; k < n; ++k)
    {
        auto i = pp[k];
        if (i != static_cast<int>(k)) {
            std::swap(xk[0], x[i]);
        }
        
        auto p = std::min(h + 1, n - k);
        for (auto j = 1u; j < p; ++j) {
            xk[j] -= lk[j-1]*xk[0];
        }

        lk += h;
        ++xk;
    }

    // 3. Backward substitution on x
    ak = a+n*w;
    xk = x+n;
  {
    for (auto k = static_cast<int>(n) - 1; k >= 0; --k)
    {
      ak-=w;
      --xk;

      auto p = std::min(w,n-k);

      double d = xk[0];
      for (auto i = 1u; i < p; ++i) {
        d -= ak[i]*xk[i];
      }

      xk[0] = d/ak[0];
    }
  }
}

}  // namespace efield
}  // namespace solver
}  // namespace steps
