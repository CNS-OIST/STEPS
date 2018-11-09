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


#ifndef STEPS_SOLVER_EFIELD_TETCOUPLER_HPP
#define STEPS_SOLVER_EFIELD_TETCOUPLER_HPP 1

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/efield/tetmesh.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps{
 namespace solver {
 namespace efield {

////////////////////////////////////////////////////////////////////////////////

/// It is temporarily created in the constructor of class EField, after a
/// TetMesh has been (partially) constructed.
///
/// \author Robert Cannon
///
class TetCoupler
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor. Just copies the mesh pointer.
    ///
    TetCoupler(TetMesh * mesh);

    /// Destructor.
    ///
    ~TetCoupler();

    ////////////////////////////////////////////////////////////////////////

    /// The major method in this class... it couples a mesh!
    ///
    /// The coupling constants are stored in the VertexConnection
    /// objects stored in the mesh.
    ///
    void coupleMesh();

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////
    // AUXILIARY FUNCTIONS FOR COUPLEMESH()
    ////////////////////////////////////////////////////////////////////////

    /// Checks whether two doubles differ.
    ///
    bool dblsDiffer(double, double);

    /// Compute the corss product between two vectors
    void cross_product(double * a, double * b, double * c);

    /// Computes the actual flux coefficients.
    ///
    ///
    void fluxCoeficients(VertexElement*, VertexElement**, double * ret);

    ////////////////////////////////////////////////////////////////////////
    // DATA FIELDS
    ////////////////////////////////////////////////////////////////////////

    TetMesh *                   pMesh;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_SIM_EFIELD_TETCOUPLER_HPP

// END
