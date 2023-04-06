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


#ifndef STEPS_SOLVER_EFIELD_VERTEXCONNECTION_HPP
#define STEPS_SOLVER_EFIELD_VERTEXCONNECTION_HPP 1

// STL headers.
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

// STEPS headers.
#include "util/common.h"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace solver {
namespace efield {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class VertexElement;
class VertexConnection;
class Mesh;

// Auxiliary declarations.
typedef VertexConnection *                     VertexConnectionP;
typedef std::vector<VertexConnectionP>         VertexConnectionPVec;
typedef VertexConnectionPVec::iterator         VertexConnectionPVecI;
typedef VertexConnectionPVec::const_iterator   VertexConnectionPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Class VertexConnection only contains pointers to vertices, so when
/// vertex indices are changed, there is no need to update anything in
/// these objects.
///
class VertexConnection
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor. Copies the VertexElement pointers and notifies them
    /// that they are a part of this connection.
    ///
    VertexConnection(VertexElement * v1, VertexElement * v2);

    /// Destructor.
    ///
    ~VertexConnection();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////

    //bool isEdge();
    //bool hasInternalEnd();

    /// When called with a VertexElement that is part of this connection,
    /// this method returns the VertexElement on the other side of the
    /// connection. If the parameter is not a part of the connection,
    /// an assertion is raised.
    ///
    VertexElement * getOther(VertexElement *);

    inline void setGeomCouplingConstant(double d) noexcept
    { pGeomCC = d; }

    inline VertexElement * getA() const noexcept
    { return pVert1; }

    inline VertexElement * getB() const noexcept
    { return pVert2; }

    inline double getGeomCouplingConstant() const noexcept
    { return pGeomCC; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT DATA
    ////////////////////////////////////////////////////////////////////////

    ///@{
    /// Point to the vertices on this connection.
    VertexElement *             pVert1;
    VertexElement *             pVert2;
    ///@}

    /// Geometric coupling constant.
    double                      pGeomCC;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_EFIELD_VERTEXCONNECTION_HPP

// END
