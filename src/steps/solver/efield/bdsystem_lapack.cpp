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


#include <algorithm>
#include <utility>

#include "steps/common.h"
#include "steps/solver/efield/bdsystem_lapack.hpp"

#include "easylogging++.h"

namespace steps {
namespace solver {
namespace efield {

extern "C" {
extern void dgbsv_(int *n,int *kl,int *ku,int *nrhs,double *ab,int *ldab,int *ipiv,double *b,int *ldb,int *info);
}

void BDSystemLapack::solve()
{
    int n=pN;
    int h=pHalfBW;
    int nrhs=1;
    int ldab=3*h+1;
    int info=0;

    std::copy(pb.begin(),pb.end(),px.begin());
    dgbsv_(&n,&h,&h,&nrhs,pA.data(),&ldab,&pwork[0],&px[0],&n,&info);
}

}  // namespace efield
}  // namespace solver
}  // namespace steps
