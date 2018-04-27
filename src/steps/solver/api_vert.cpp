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
#include <string>
#include <sstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/api.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/geom/tetmesh.hpp"
// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

USING(std, string);
using namespace steps::solver;

////////////////////////////////////////////////////////////////////////////////

double API::getVertV(uint vidx) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (vidx >= mesh->countVertices())
        {
            std::ostringstream os;
            os << "Vertex index out of range.";
            ArgErrLog(os.str());
        }
        return _getVertV(vidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setVertV(uint vidx, double v)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (vidx >= mesh->countVertices())
        {
            std::ostringstream os;
            os << "Vertex index out of range.";
            ArgErrLog(os.str());
        }
        _setVertV(vidx, v);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getVertVClamped(uint vidx) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (vidx >= mesh->countVertices())
        {
            std::ostringstream os;
            os << "Vertex index out of range.";
            ArgErrLog(os.str());
        }
        _getVertVClamped(vidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setVertVClamped(uint vidx, bool cl)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (vidx >= mesh->countVertices())
        {
            std::ostringstream os;
            os << "Vertex index out of range.";
            ArgErrLog(os.str());
        }
        _setVertVClamped(vidx, cl);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setVertIClamp(uint vidx, double i)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (vidx >= mesh->countVertices())
        {
            std::ostringstream os;
            os << "Vertex index out of range.";
            ArgErrLog(os.str());
        }
        _setVertIClamp(vidx, i);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::_getVertV(uint vidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setVertV(uint vidx, double v)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getVertVClamped(uint vidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setVertVClamped(uint vidx, bool cl)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setVertIClamp(uint vidx, double i)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

// END
