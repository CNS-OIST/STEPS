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

#include "vertexconnection.hpp"

#include "util/checkpointing.hpp"
#include "util/error.hpp"
#include "vertexelement.hpp"

namespace steps::solver::efield {

VertexConnection::VertexConnection(VertexElement* v1, VertexElement* v2)
    : pVert1(v1)
    , pVert2(v2)
    , pGeomCC(0.0) {
    AssertLog(v1 != nullptr);
    AssertLog(v2 != nullptr);
    pVert1->addConnection(this);
    pVert2->addConnection(this);
}

////////////////////////////////////////////////////////////////////////////////

VertexConnection::~VertexConnection() = default;

////////////////////////////////////////////////////////////////////////////////

void VertexConnection::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pGeomCC);
}

////////////////////////////////////////////////////////////////////////////////

void VertexConnection::restore(std::fstream& cp_file) {
    util::compare(cp_file, pGeomCC, "Mismatched pGeomCC restore value.");
}

////////////////////////////////////////////////////////////////////////////////

VertexElement* VertexConnection::getOther(VertexElement* element) {
    VertexElement* ret;
    if (pVert1 == element) {
        ret = pVert2;
    } else if (pVert2 == element) {
        ret = pVert1;
    } else {
        AssertLog(false);
    }
    return ret;
}

}  // namespace steps::solver::efield
