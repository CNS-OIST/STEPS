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


#ifndef STEPS_SOLVER_EFIELD_SLUSYSTEM_HPP
#define STEPS_SOLVER_EFIELD_SLUSYSTEM_HPP 1

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include <memory>

#include <iterator>
#include <algorithm>

#include <mpi.h>

#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/efield/linsystem.hpp"


// logging
#include "third_party/easyloggingpp/src/easylogging++.h"

namespace steps {
namespace solver {
namespace efield {

struct sparsity_template {
    explicit sparsity_template(size_t n): pN(n), cols(n), count(0) {}

    void insert(std::pair<int,int> ij) {
        cols.at(ij.second).push_back(ij.first);
        ++count;
    }

    template <typename Iter>
    void insert(Iter b,Iter e) {
        while (b != e) insert(*b++);
    }

    size_t size() const { return count; }
    size_t dim() const { return pN; }

    void clear() {
        std::vector<int> empty_col;
        cols.assign(pN,empty_col);
        count = 0;
    }

    size_t pN;
    std::vector<std::vector<int>> cols;
    size_t count;
};

struct supermatrix_nc_view;

struct SLU_NCMatrix: public AMatrix {
    explicit SLU_NCMatrix(const sparsity_template &S): pN(S.dim()), pNnz(S.size()), pValues(S.size()), pRidx(S.size()), pCoff(S.dim()+1)
    {
        int offset = 0;
        int *ridx = &pRidx.front();
        for (int col=0; col<pN; ++col) {
            const auto &scol = S.cols[col];

            pCoff[col] = offset;
            std::copy(scol.begin(), scol.end(), ridx+offset);

            int next_offset = offset + scol.size();
            std::sort(ridx+offset, ridx+next_offset);
            
            if (ridx[offset]<0 || ridx[next_offset-1]>=pN)
                ProgErrLog("out of range element in sparsity matrix");

            // handle duplicates gracefully
            const int *adj = std::unique(ridx+offset, ridx+next_offset);
            next_offset = adj-ridx;

            offset = next_offset;
        }
        pCoff[pN] = offset;
        if (offset < pNnz) {
            // there were some duplicates
            pNnz = offset;

            pValues.resize(pNnz);
            pValues.shrink_to_fit();

            pRidx.resize(pNnz);
            pValues.shrink_to_fit();
        }
    }

    size_t nRow() const override final { return pN; }
    size_t nCol() const override final { return pN; }

    void zero() override final {
        pValues.assign(pNnz, 0.0);
    }

    double get(size_t row, size_t col) const override final {
        int i=get_offset((int)row, (int)col);
        return i>=0?pValues[i]:0;
    }

    void set(size_t row, size_t col, double value) override final {
        int i=get_offset((int)row, (int)col);
        if (i<0) ArgErrLog("index not in sparse template");

        pValues[i]=value;
    }

    size_t nNz() const { return pNnz; }

    friend struct supermatrix_nc_view;

private:
    size_t pN, pNnz;

    std::vector<double> pValues;
    std::vector<int> pRidx;
    std::vector<int> pCoff;

    int get_offset(int i, int j) const {
        const int *r0 = &pRidx.front();
        const int *rb = r0+pCoff[j];
        const int *re = r0+pCoff[j+1];

        const int *r = std::lower_bound(rb, re, i);
        if (r==re || *r!=i)
            return -1;
        else return r-r0;
    }

};

struct SLUData;

// fill this out ...
struct SLUSolverStats {

};

class SLUSystem {
public:
    typedef SLU_NCMatrix matrix_type;
    typedef VVector vector_type;

    // ... plus add SuperLU options ...
    explicit SLUSystem(const sparsity_template &S, MPI_Comm mpi_comm);
    ~SLUSystem();

    const matrix_type &A() const { return pA; }
    matrix_type &A() { return pA; }

    const vector_type &b() const { return pb_view; }
    vector_type &b() { return pb_view; }

    const vector_type &x() const { return px_view; }

    void solve();
    
    // query solver stats, error
    double berr() const { return pBerr; }

    // (stub)
    SLUSolverStats solver_stats() const { return SLUSolverStats(); }

    // (stub)
    SLUSolverStats solver_stats_global() const { return SLUSolverStats(); }

private:
    int pN;
    matrix_type pA;
    std::vector<double> pb, px;
    vector_type pb_view, px_view;

    std::unique_ptr<SLUData> slu;
    double pBerr;
};


}}} // namespace steps::solver::efield

#endif // ndef STEPS_SOLVER_EFIELD_SLUSYSTEM_HPP
