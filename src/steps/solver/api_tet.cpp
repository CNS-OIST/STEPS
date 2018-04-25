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

double API::getTetVol(uint tidx) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }

        return _getTetVol(tidx);
    }
    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetVol(uint tidx, double vol)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // NOTE: the following method may never be implemented
        _setTetVol(tidx, vol);
    }
    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetSpecDefined(uint tidx, string const & s) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        return _getTetSpecDefined(tidx, sidx);
    }
    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////


double API::getTetCount(uint tidx, string const & s) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        return _getTetCount(tidx, sidx);
    }
    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetCount(uint tidx, string const & s, double n)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        if (n < 0.0)
        {
            std::ostringstream os;
            os << "Number of molecules cannot be negative.";
            ArgErrLog(os.str());
        }

    // the following may throw exception if string is unknown
    uint sidx = pStatedef->getSpecIdx(s);

    _setTetCount(tidx, sidx, n);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetAmount(uint tidx, string const & s) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        return _getTetAmount(tidx, sidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetAmount(uint tidx, string const & s, double m)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
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

        _setTetAmount(tidx, sidx, m);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetConc(uint tidx, string const & s) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        return _getTetConc(tidx, sidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetConc(uint tidx, string const & s, double c)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        if (c < 0.0)
        {
            std::ostringstream os;
            os << "Concentration cannot be negative.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        _setTetConc(tidx, sidx, c);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetClamped(uint tidx, string const & s) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        return _getTetClamped(tidx, sidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetClamped(uint tidx, string const & s, bool buf)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint sidx = pStatedef->getSpecIdx(s);

        _setTetClamped(tidx, sidx, buf);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacK(uint tidx, string const & r) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint ridx = pStatedef->getReacIdx(r);

        return _getTetReacK(tidx, ridx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetReacK(uint tidx, string const & r, double kf)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        if (kf < 0.0)
        {
            std::ostringstream os;
            os << "Reaction constant cannot be negative.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint ridx = pStatedef->getReacIdx(r);

        _setTetReacK(tidx, ridx, kf);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetReacActive(uint tidx, string const & r) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint ridx = pStatedef->getReacIdx(r);

        return _getTetReacActive(tidx, ridx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetReacActive(uint tidx, string const & r, bool act)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint ridx = pStatedef->getReacIdx(r);

        _setTetReacActive(tidx, ridx, act);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetDiffD(uint tidx, string const & d, uint direction_tet) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint didx = pStatedef->getDiffIdx(d);

        return _getTetDiffD(tidx, didx, direction_tet);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetDiffD(uint tidx, string const & d, double dk, uint direction_tet)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        if (direction_tet != std::numeric_limits<uint>::max() &&direction_tet >= mesh->countTets())
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
        uint didx = pStatedef->getDiffIdx(d);

        _setTetDiffD(tidx, didx, dk, direction_tet);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetDiffActive(uint tidx, string const & d) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint didx = pStatedef->getDiffIdx(d);

        return _getTetDiffActive(tidx, didx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetDiffActive(uint tidx, string const & d, bool act)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {

        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint didx = pStatedef->getDiffIdx(d);

        _setTetDiffActive(tidx, didx, act);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacH(uint tidx, string const & r) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {

        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint ridx = pStatedef->getReacIdx(r);

        return _getTetReacH(tidx, ridx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacC(uint tidx, string const & r) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint ridx = pStatedef->getReacIdx(r);

        return _getTetReacC(tidx, ridx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacA(uint tidx, string const & r) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {

        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint ridx = pStatedef->getReacIdx(r);

        return _getTetReacA(tidx, ridx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetDiffA(uint tidx, string const & d) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        // the following may throw exception if string is unknown
        uint didx = pStatedef->getDiffIdx(d);

        return _getTetDiffA(tidx, didx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetV(uint tidx) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }

        return _getTetV(tidx);
    }
    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetV(uint tidx, double v)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }

        _setTetV(tidx, v);
    }
    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetVClamped(uint tidx) const
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }

        return _getTetVClamped(tidx);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetVClamped(uint tidx, bool cl)
{
    if (steps::tetmesh::Tetmesh * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= mesh->countTets())
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }

        _setTetVClamped(tidx, cl);
    }

    else
    {
        std::ostringstream os;
        os << "Method not available for this solver.";
        NotImplErrLog("");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetVol(uint tidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetVol(uint tidx, double vol)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetSpecDefined(uint tidx, uint sidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetCount(uint tidx, uint sidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetCount(uint tidx, uint sidx, double n)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetAmount(uint tidx, uint sidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetAmount(uint tidx, uint sidx, double m)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetConc(uint tidx, uint sidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetConc(uint tidx, uint sidx, double c)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetClamped(uint tidx, uint sidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetClamped(uint tidx, uint sidx, bool buf)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacK(uint tidx, uint ridx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetReacK(uint tidx, uint ridx, double kf)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetReacActive(uint tidx, uint ridx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetReacActive(uint tidx, uint ridx, bool act)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetDiffD(uint tidx, uint didx, uint direction_tet) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetDiffD(uint tidx, uint didx, double dk, uint direction_tet)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetDiffActive(uint tidx, uint didx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetDiffActive(uint tidx, uint didx, bool act)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacH(uint tidx, uint ridx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacC(uint tidx, uint ridx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacA(uint tidx, uint ridx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetDiffA(uint tidx, uint didx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetV(uint tidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetV(uint tidx, double v)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetVClamped(uint tidx) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetVClamped(uint tidx, bool cl)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

// END

