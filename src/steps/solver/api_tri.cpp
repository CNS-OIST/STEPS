/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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
#include <sstream>
#include <string>


// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/geom/tetmesh.hpp"
#include "steps/solver/api.hpp"
#include "steps/solver/statedef.hpp"
// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace solver {

////////////////////////////////////////////////////////////////////////////////

double API::getTriArea(triangle_id_t tidx) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

void API::setTriArea(triangle_id_t tidx, double area)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriCount(triangle_id_t tidx, const std::string&  s) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

bool API::getTriSpecDefined(triangle_id_t tidx, const std::string&  s) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

void API::setTriCount(triangle_id_t tidx, const std::string&  s, double n)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriAmount(triangle_id_t tidx, const std::string&  s) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

void API::setTriAmount(triangle_id_t tidx, const std::string&  s, double m)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

bool API::getTriClamped(triangle_id_t tidx, const std::string&  s) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

void API::setTriClamped(triangle_id_t tidx, const std::string&  s, bool buf)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriSReacK(triangle_id_t tidx, const std::string&  r)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

void API::setTriSReacK(triangle_id_t tidx, const std::string&  r, double kf)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

bool API::getTriSReacActive(triangle_id_t tidx, const std::string&  r)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

void API::setTriSReacActive(triangle_id_t tidx, const std::string&  r, bool act)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriSReacH(triangle_id_t tidx, const std::string&  r)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriSReacC(triangle_id_t tidx, const std::string&  r)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriSReacA(triangle_id_t tidx, const std::string&  r)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriDiffD(triangle_id_t tidx, const std::string&  d, uint direction_tri)
{
    return getTriSDiffD(tidx, d, direction_tri);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSDiffD(triangle_id_t tidx, const std::string &d,
                         triangle_id_t direction_tri)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

void API::setTriDiffD(triangle_id_t tidx, const std::string&  d, double dk, triangle_id_t direction_tri) {
    setTriSDiffD(tidx, d, dk, direction_tri);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSDiffD(triangle_id_t tidx, const std::string&  d, double dk, triangle_id_t direction_tri)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        if (direction_tri != UNKNOWN_TRI &&direction_tri >= mesh->countTris())
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

double API::getTriV(triangle_id_t tidx) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

void API::setTriV(triangle_id_t tidx, double v)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

bool API::getTriVClamped(triangle_id_t tidx) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

void API::setTriVClamped(triangle_id_t tidx, bool cl)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriOhmicI(triangle_id_t tidx) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriOhmicI(triangle_id_t tidx, const std::string&  oc) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriGHKI(triangle_id_t tidx) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriGHKI(triangle_id_t tidx, const std::string&  ghk) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriI(triangle_id_t tidx) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::getTriIClamp(triangle_id_t tidx) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTris())
        {
            std::ostringstream os;
            os << "Triangle index out of range.";
            ArgErrLog(os.str());
        }

        return _getTriIClamp(tidx);
    }
    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriIClamp(triangle_id_t tidx, double i)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

void API::setTriCapac(triangle_id_t tidx, double cm)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

bool API::getTriVDepSReacActive(triangle_id_t tidx, const std::string&  vsr)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

void API::setTriVDepSReacActive(triangle_id_t tidx, const std::string&  vsr, bool act)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
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

double API::_getTriArea(triangle_id_t /*tidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriArea(triangle_id_t /*tidx*/, double /*area*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriSpecDefined(triangle_id_t /*tidx*/, uint /*sidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriCount(triangle_id_t /*tidx*/, uint /*sidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriCount(triangle_id_t /*tidx*/, uint /*sidx*/, double /*n*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriAmount(triangle_id_t /*tidx*/, uint /*sidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriAmount(triangle_id_t /*tidx*/, uint /*sidx*/, double /*m*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriClamped(triangle_id_t /*tidx*/, uint /*sidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriClamped(triangle_id_t /*tidx*/, uint /*sidx*/, bool /*buf*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacK(triangle_id_t /*tidx*/, uint /*ridx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSReacK(triangle_id_t /*tidx*/, uint /*ridx*/, double /*kf*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriSReacActive(triangle_id_t /*tidx*/, uint /*ridx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSReacActive(triangle_id_t /*tidx*/, uint /*ridx*/, bool /*act*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSDiffD(triangle_id_t /*tidx*/, uint /*didx*/, triangle_id_t /*direction_tri*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSDiffD(triangle_id_t /*tidx*/, uint /*didx*/, double /*dk*/,
                        triangle_id_t /*direction_tri*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacH(triangle_id_t /*tidx*/, uint /*ridx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacC(triangle_id_t /*tidx*/, uint /*ridx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacA(triangle_id_t /*tidx*/, uint /*ridx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriV(triangle_id_t /*tidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriV(triangle_id_t /*tidx*/, double /*v*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriVClamped(triangle_id_t /*tidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriVClamped(triangle_id_t /*tidx*/, bool /*cl*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriOhmicI(triangle_id_t /*tidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriOhmicI(triangle_id_t /*tidx*/, uint /*ocidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriGHKI(triangle_id_t /*tidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriGHKI(triangle_id_t /*tidx*/, uint /*ghkidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriI(triangle_id_t /*tidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriIClamp(triangle_id_t /*tidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriIClamp(triangle_id_t /*tidx*/, double /*i*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriCapac(triangle_id_t /*tidx*/, double /*cm*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriVDepSReacActive(triangle_id_t /*tidx*/, uint /*vsridx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriVDepSReacActive(triangle_id_t /*tidx*/, uint /*vsridx*/, bool /*act*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

// END

} // namespace solver
} // namespace steps
