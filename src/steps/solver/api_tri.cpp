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

double API::getTriArea(uint tidx) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        return _getTriArea(tidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriArea(uint tidx, double area)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // NOTE: the following method may never be implemented
        _setTriArea(tidx, area);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriCount(uint tidx, string const & s) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may raise exceptions if strings are unused
        uint sidx = pStatedef->getSpecIdx(s);

        return _getTriCount(tidx, sidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriSpecDefined(uint tidx, string const & s) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        return _getTriSpecDefined(tidx, sidx);
    }
    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriCount(uint tidx, string const & s, double n)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        if (n < 0.0)
        {
            std::ostringstream os;
            os << "Number of molecules cannot be negative.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is uknown
        uint sidx = pStatedef->getSpecIdx(s);

        _setTriCount(tidx, sidx, n);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriAmount(uint tidx, string const & s) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        return _getTriAmount(tidx, sidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriAmount(uint tidx, string const & s, double m)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        if (m < 0.0)
        {
            std::ostringstream os;
            os << "Amount of mols cannot be negative.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        _setTriAmount(tidx, sidx, m);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriClamped(uint tidx, string const & s) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        return _getTriClamped(tidx, sidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriClamped(uint tidx, string const & s, bool buf)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        _setTriClamped(tidx, sidx, buf);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacK(uint tidx, string const & r) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is unknown
        uint sridx = pStatedef->getSReacIdx(r);

        return _getTriSReacK(tidx, sridx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSReacK(uint tidx, string const & r, double kf)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        if (kf < 0.0)
        {
            std::ostringstream os;
            os << "Reaction constant cannot be negative.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is unknown
        uint sridx = pStatedef->getSReacIdx(r);

        _setTriSReacK(tidx, sridx, kf);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriSReacActive(uint tidx, string const & r) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is unknown
        uint sridx = pStatedef->getSReacIdx(r);

        return _getTriSReacActive(tidx, sridx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSReacActive(uint tidx, string const & r, bool act)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        // the following may raise exception if string is unknown
        uint sridx = pStatedef->getSReacIdx(r);

        _setTriSReacActive(tidx, sridx, act);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacH(uint tidx, string const & r) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is unknown
        uint sridx = pStatedef->getSReacIdx(r);

        return _getTriSReacH(tidx, sridx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacC(uint tidx, string const & r) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is unknown
        uint sridx = pStatedef->getSReacIdx(r);

        return _getTriSReacC(tidx, sridx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacA(uint tidx, string const & r) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is unknown
        uint sridx = pStatedef->getSReacIdx(r);

        return _getTriSReacA(tidx, sridx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriDiffD(uint tidx, string const & d, uint direction_tri) const
{
    return getTriSDiffD(tidx, d, direction_tri);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSDiffD(uint tidx, string const & d, uint direction_tri) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint didx = pStatedef->getSurfDiffIdx(d);
        
        return _getTriSDiffD(tidx, didx, direction_tri);
    }
    
    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}


////////////////////////////////////////////////////////////////////////////////

void API::setTriDiffD(uint tidx, string const & d, double dk, uint direction_tri) {
    setTriSDiffD(tidx, d, dk, direction_tri);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSDiffD(uint tidx, string const & d, double dk, uint direction_tri)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        if (direction_tri != std::numeric_limits<uint>::max() &&direction_tri >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Direction tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        if (dk < 0.0)
        {
            std::ostringstream os;
            os << "Diffusion constant cannot be negative.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint didx = pStatedef->getSurfDiffIdx(d);
        
        _setTriSDiffD(tidx, didx, dk, direction_tri);
    }
    
    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}
////////////////////////////////////////////////////////////////////////////////

double API::getTriV(uint tidx) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        return _getTriV(tidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriV(uint tidx, double v)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        _setTriV(tidx, v);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriVClamped(uint tidx) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        return _getTriVClamped(tidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriVClamped(uint tidx, bool cl)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        _setTriVClamped(tidx, cl);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriOhmicI(uint tidx)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        return _getTriOhmicI(tidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriOhmicI(uint tidx, string const & oc)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is unknown
        uint ocidx = pStatedef->getOhmicCurrIdx(oc);

        return _getTriOhmicI(tidx, ocidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriGHKI(uint tidx)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        return _getTriGHKI(tidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriGHKI(uint tidx, string const & ghk)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is unknown
        uint ghkidx = pStatedef->getGHKcurrIdx(ghk);

        return _getTriGHKI(tidx, ghkidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}
////////////////////////////////////////////////////////////////////////////////

double API::getTriI(uint tidx) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        return _getTriI(tidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriIClamp(uint tidx, double i)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        _setTriIClamp(tidx, i);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriCapac(uint tidx, double cm)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        _setTriCapac(tidx, cm);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriVDepSReacActive(uint tidx, string const & vsr) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }
        // the following may raise exception if string is unknown
        uint vsridx = pStatedef->getVDepSReacIdx(vsr);

        return _getTriVDepSReacActive(tidx, vsridx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriVDepSReacActive(uint tidx, string const & vsr, bool act)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        // the following may raise exception if string is unknown
        uint vsridx = pStatedef->getVDepSReacIdx(vsr);

        _setTriVDepSReacActive(tidx, vsridx, act);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriArea(uint tidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriArea(uint tidx, double area)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriSpecDefined(uint tidx, uint sidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriCount(uint tidx, uint sidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriCount(uint tidx, uint sidx, double n)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriAmount(uint tidx, uint sidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriAmount(uint tidx, uint sidx, double m)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriClamped(uint tidx, uint sidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriClamped(uint tidx, uint sidx, bool buf)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacK(uint tidx, uint ridx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSReacK(uint tidx, uint ridx, double kf)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriSReacActive(uint tidx, uint ridx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSReacActive(uint tidx, uint ridx, bool act)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSDiffD(uint tidx, uint didx, uint direction_tri) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSDiffD(uint tidx, uint didx, double dk, uint direction_tri)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacH(uint tidx, uint ridx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacC(uint tidx, uint ridx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacA(uint tidx, uint ridx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriV(uint tidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriV(uint tidx, double v)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriVClamped(uint tidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriVClamped(uint tidx, bool cl)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriOhmicI(uint tidx)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriOhmicI(uint tidx, uint ocidx)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriGHKI(uint tidx)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriGHKI(uint tidx, uint ghkidx)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriI(uint tidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriIClamp(uint tidx, double i)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriCapac(uint tidx, double cm)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriVDepSReacActive(uint tidx, uint vsridx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriVDepSReacActive(uint tidx, uint vsridx, bool act)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

// END
