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


// STL headers.
#include <iostream>
#include <sstream>
#include <cassert>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/efield/vertexconnection.hpp"
#include "steps/solver/efield/vertexelement.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace sefield = steps::solver::efield;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

sefield::VertexConnection::VertexConnection(sefield::VertexElement * v1, sefield::VertexElement * v2)
: pVert1(v1)
, pVert2(v2)
, pGeomCC(0.0)
{
    AssertLog(v1 != 0);
    AssertLog(v2 != 0);
    pVert1->addConnection(this);
    pVert2->addConnection(this);
}

////////////////////////////////////////////////////////////////////////////////

sefield::VertexConnection::~VertexConnection(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void sefield::VertexConnection::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&pGeomCC, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void sefield::VertexConnection::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&pGeomCC, sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

sefield::VertexElement * sefield::VertexConnection::getOther(sefield::VertexElement * element)
{
    VertexElement * ret;
    if (pVert1 == element)
    {
        ret = pVert2;
    }
    else if (pVert2 == element)
    {
        ret = pVert1;
    }
    else
    {
        AssertLog(0);
        ret = 0;
    }
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

//bool VertexConnection::hasInternalEnd(void)
//{
//    return (vea->isInternal() || veb->isInternal());
//}

////////////////////////////////////////////////////////////////////////////////

//bool VertexConnection::isEdge(void)
//{
//    return (vea->isEdge() && veb->isEdge());
//}

////////////////////////////////////////////////////////////////////////////////

// END
