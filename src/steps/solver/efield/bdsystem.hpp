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


#ifndef STEPS_SOLVER_EFIELD_BDSYSTEM_HPP
#define STEPS_SOLVER_EFIELD_BDSYSTEM_HPP 1

#include <cstddef>
#include <vector>

#include "steps/common.h"
#include "steps/solver/efield/linsystem.hpp"

namespace steps {
namespace solver {
namespace efield {


class BDMatrix: public AMatrix {
public:
    BDMatrix(size_t n,size_t halfbw):
        pN(n), pData((2*halfbw+1)*n,0.0), p00(&pData[halfbw]), pStride(2*halfbw)
    {}

    size_t nRow() const override final { return pN; }
    size_t nCol() const override final { return pN; }

    double get(size_t row,size_t col) const override final { return p00[row*pStride+col]; }
    void set(size_t row,size_t col,double value) override final { p00[row*pStride+col]=value; }

    void zero() override final {
        pData.assign(pData.size(),0.0);
    }

    // more direct access via operator[]:

    const double *operator[](size_t row) const { return p00+row*pStride; }
    double *operator[](size_t row) { return p00+row*pStride; }

    // direct access to compact representation

    inline const double *data() const noexcept { return pData.data(); }
    inline double *data() noexcept { return pData.data(); }
    inline size_t data_size() const noexcept { return pData.size(); }

private:
    size_t pN;
    std::vector<double> pData;
    double *p00;
    size_t pStride;
};

class BDSystem
{
public:
    typedef BDMatrix matrix_type;
    typedef VVector vector_type;

    BDSystem(size_t n, size_t halfbw):
        pN(n),
        pHalfBW(halfbw),
        pA(n,halfbw),
        pb(n,0.0),
        px(n,0.0),
        pL(halfbw*n,0.0),
        pp(n,0),
        pb_view(n, pb.data()),
        px_view(n, px.data())
    {}

    const matrix_type &A() const { return pA; }
    matrix_type &A() { return pA; }

    const vector_type &b() const { return pb_view; }
    vector_type &b() { return pb_view; }

    const vector_type &x() const { return px_view; }

    void solve(); // destructive: overwrites pA

private:
    size_t pN,pHalfBW;

    BDMatrix pA; // will contain U after LU-decomposition
    std::vector<double> pb;
    std::vector<double> px;
    std::vector<double> pL;
    std::vector<int> pp; // permutation vector post LU-decomposition

    vector_type pb_view;
    vector_type px_view;
};


}}} // namespace steps::solver::efield

#endif // ndef STEPS_SOLVER_EFIELD_BDSYSTEM_HPP

