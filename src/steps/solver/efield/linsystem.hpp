/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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


#ifndef STEPS_SOLVER_EFIELD_LINSYSTEM_HPP
#define STEPS_SOLVER_EFIELD_LINSYSTEM_HPP

#include <algorithm>
#include <cstddef>

/// Abstract base class for solving linear system; interface to
/// existing banded matrix solver and parallel implementations.

namespace steps {
namespace solver {
namespace efield {

// Abstract matrix and vector classes exist only to provide generic
// access through LinSystem; use concrete implementations directly
// for efficiency.

class AMatrix
{
public:
    virtual size_t nRow() const =0;
    virtual size_t nCol() const =0;

    virtual double get(size_t row,size_t col) const =0;
    virtual void set(size_t row,size_t col,double value) =0;
    virtual void zero() =0;

    virtual ~AMatrix() {}
};

class AVector
{
public:
    virtual size_t size() const =0;

    virtual double get(size_t i) const =0;
    virtual void set(size_t i,double value) =0;
    virtual void zero() =0;

    virtual ~AVector() {}
};

// Vector as contiguous view on memory

class VVector: public AVector {
public:
    VVector(size_t n,double *data): pN(n), pData(data) {}

    size_t size() const override final { return pN; }

    double get(size_t i) const override final { return pData[i]; }
    void set(size_t i,double value) override final { pData[i]=value; }
    void zero() override final { std::fill(pData,pData+pN,0.0); }

    // more direct access via operator[]:

    double operator[](size_t i) const { return *(pData+i); }
    double &operator[](size_t i) { return *(pData+i); }

private:
    size_t pN;
    double *pData;
};

class LinSystem
{
public:
    LinSystem(size_t nrow,size_t ncol): pNRow(nrow), pNCol(ncol) {}

    virtual ~LinSystem() = default;

    virtual const AMatrix &A() const =0;
    virtual AMatrix &A() =0;

    virtual const AVector &b() const =0;
    virtual AVector &b() =0;

    virtual const AVector &x() const =0;
    
    virtual void solve() =0;

    inline size_t nrow() const noexcept { return pNRow; }
    inline size_t ncol() const noexcept { return pNCol; }

private:
    size_t pNRow,pNCol;
};


}}}  // namespace steps::solver::efield


#endif // ndef  STEPS_SOLVER_EFIELD_LINSYSTEM_HPP
