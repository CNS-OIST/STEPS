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

#include "api.hpp"

#include "geom/tetmesh.hpp"
#include "util/error.hpp"

namespace steps::solver {

double API::getVertV(vertex_id_t vidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        if (vidx >= mesh->countVertices()) {
            std::ostringstream os;
            os << "Vertex index out of range.";
            ArgErrLog(os.str());
        }
        return _getVertV(vidx);
    }

    else {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setVertV(vertex_id_t vidx, double v) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        if (vidx >= mesh->countVertices()) {
            std::ostringstream os;
            os << "Vertex index out of range.";
            ArgErrLog(os.str());
        }
        _setVertV(vidx, v);
    }

    else {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getVertVClamped(vertex_id_t vidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        if (vidx >= mesh->countVertices()) {
            std::ostringstream os;
            os << "Vertex index out of range.";
            ArgErrLog(os.str());
        }
        return _getVertVClamped(vidx);
    }

    else {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setVertVClamped(vertex_id_t vidx, bool cl) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        if (vidx >= mesh->countVertices()) {
            std::ostringstream os;
            os << "Vertex index out of range.";
            ArgErrLog(os.str());
        }
        _setVertVClamped(vidx, cl);
    }

    else {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getVertIClamp(vertex_id_t vidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        if (vidx >= mesh->countVertices()) {
            std::ostringstream os;
            os << "Vertex index out of range.";
            ArgErrLog(os.str());
        }
        return _getVertIClamp(vidx);
    } else {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setVertIClamp(vertex_id_t vidx, double i) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        if (vidx >= mesh->countVertices()) {
            std::ostringstream os;
            os << "Vertex index out of range.";
            ArgErrLog(os.str());
        }
        _setVertIClamp(vidx, i);
    }

    else {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::_getVertV(vertex_id_t /*vidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setVertV(vertex_id_t /*vidx*/, double /*v*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getVertVClamped(vertex_id_t /*vidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setVertVClamped(vertex_id_t /*vidx*/, bool /*cl*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getVertIClamp(vertex_id_t /*vidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setVertIClamp(vertex_id_t /*vidx*/, double /*i*/) {
    NotImplErrLog("");
}

}  // namespace steps::solver
