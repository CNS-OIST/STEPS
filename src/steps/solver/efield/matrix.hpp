/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_SOLVER_EFIELD_MATRIX_HPP
#define STEPS_SOLVER_EFIELD_MATRIX_HPP 1

#include <fstream>
#include <iostream>

// STEPS headers.
#include "util/common.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace solver {
namespace efield {

////////////////////////////////////////////////////////////////////////////////

/// \todo Clean up (especially get rid of the pointer-to-pointer storage);
/// could be a useful addition for Boost STEPS.
///
/// \author Robert Cannon
///
class Matrix
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor that creates an uninitialized n0 * n0 matrix.
    ///
    Matrix(uint n0);

    /// Constructor that creates an nn * nn matrix and initializes it
    /// by copying the contents of da.
    ///
    Matrix(uint nn, double ** da);

    /// Destructor.
    ///
    ~Matrix();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // MATRIX OPERATIONS
    ////////////////////////////////////////////////////////////////////////

    /// Makes a deep copy of the matrix.
    ///
    Matrix * copy();

    /// Computes left-hand vector product.
    ///
    /// The resulting array needs to deallocated by the caller!
    ///
    double * lvprod(double * v);

    /// Returns the transpose of this matrix.
    ///
    Matrix * transpose();

    /// Computes the determinant of this matrix.
    ///
    double det();

    /// Returns the inverse of this matrix.
    ///
    Matrix * inverse();

    /// Compute the LU decomposition.
    ///
    void LU();

    ///
    double * lubksb(double*);

    /// NOT IMPLEMENTED????
    double * rvprod(double*);

    ////////////////////////////////////////////////////////////////////////

private:

    double       ** pA;
    double        * pWS;
    uint            pN;
    int           * pPerm;
    int             pSign;

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SOLVER_EFIELD_MATRIX_HPP

// END
