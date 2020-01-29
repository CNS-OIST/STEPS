/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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

USING(std, string);

namespace steps {
namespace solver {

////////////////////////////////////////////////////////////////////////////////

double API::getTetVol(tetrahedron_id_t tidx) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

void API::setTetVol(tetrahedron_id_t tidx, double vol)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

bool API::getTetSpecDefined(tetrahedron_id_t tidx, string const & s) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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


double API::getTetCount(tetrahedron_id_t tidx, string const & s) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

void API::setTetCount(tetrahedron_id_t tidx, string const & s, double n)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

double API::getTetAmount(tetrahedron_id_t tidx, string const & s) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

void API::setTetAmount(tetrahedron_id_t tidx, string const & s, double m)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

double API::getTetConc(tetrahedron_id_t tidx, string const & s) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

void API::setTetConc(tetrahedron_id_t tidx, string const & s, double c)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

bool API::getTetClamped(tetrahedron_id_t tidx, string const & s) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

void API::setTetClamped(tetrahedron_id_t tidx, string const & s, bool buf)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

double API::getTetReacK(tetrahedron_id_t tidx, string const & r) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

void API::setTetReacK(tetrahedron_id_t tidx, string const & r, double kf)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

bool API::getTetReacActive(tetrahedron_id_t tidx, string const & r) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

void API::setTetReacActive(tetrahedron_id_t tidx, string const & r, bool act)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

double API::getTetDiffD(tetrahedron_id_t tidx, const string &d,
                        tetrahedron_id_t direction_tet) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

void API::setTetDiffD(tetrahedron_id_t tidx, const string &d, double dk,
                      tetrahedron_id_t direction_tet)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
        {
            std::ostringstream os;
            os << "Tetrahedron index out of range.";
            ArgErrLog(os.str());
        }
        if (direction_tet != UNKNOWN_TET && direction_tet >= static_cast<index_t>(mesh->countTets()))
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

bool API::getTetDiffActive(tetrahedron_id_t tidx, string const & d) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

void API::setTetDiffActive(tetrahedron_id_t tidx, string const & d, bool act)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {

        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

double API::getTetReacH(tetrahedron_id_t tidx, string const & r) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {

        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

double API::getTetReacC(tetrahedron_id_t tidx, string const & r) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

double API::getTetReacA(tetrahedron_id_t tidx, string const & r) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {

        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

double API::getTetDiffA(tetrahedron_id_t tidx, string const & d) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

double API::getTetV(tetrahedron_id_t tidx) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

void API::setTetV(tetrahedron_id_t tidx, double v)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

bool API::getTetVClamped(tetrahedron_id_t tidx) const
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

void API::setTetVClamped(tetrahedron_id_t tidx, bool cl)
{
    if (auto * mesh = dynamic_cast<steps::tetmesh::Tetmesh*>(geom()))
    {
        if (tidx >= static_cast<index_t>(mesh->countTets()))
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

double API::_getTetVol(tetrahedron_id_t /*tidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetVol(tetrahedron_id_t /*tidx*/, double /*vol*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetSpecDefined(tetrahedron_id_t /*tidx*/, uint /*sidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetCount(tetrahedron_id_t /*tidx*/, uint /*sidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetCount(tetrahedron_id_t /*tidx*/, uint /*sidx*/, double /*n*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetAmount(tetrahedron_id_t /*tidx*/, uint /*sidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetAmount(tetrahedron_id_t /*tidx*/, uint /*sidx*/, double /*m*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetConc(tetrahedron_id_t /*tidx*/, uint /*sidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetConc(tetrahedron_id_t /*tidx*/, uint /*sidx*/, double /*c*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetClamped(tetrahedron_id_t /*tidx*/, uint /*sidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetClamped(tetrahedron_id_t /*tidx*/, uint /*sidx*/, bool /*buf*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacK(tetrahedron_id_t /*tidx*/, uint /*ridx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetReacK(tetrahedron_id_t /*tidx*/, uint /*ridx*/, double /*kf*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetReacActive(tetrahedron_id_t /*tidx*/, uint /*ridx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetReacActive(tetrahedron_id_t /*tidx*/, uint /*ridx*/, bool /*act*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetDiffD(tetrahedron_id_t /*tidx*/, uint /*didx*/, tetrahedron_id_t /*direction_tet*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetDiffD(tetrahedron_id_t /*tidx*/, uint /*didx*/, double /*dk*/, tetrahedron_id_t /*direction_tet*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetDiffActive(tetrahedron_id_t /*tidx*/, uint /*didx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetDiffActive(tetrahedron_id_t /*tidx*/, uint /*didx*/, bool /*act*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacH(tetrahedron_id_t /*tidx*/, uint /*ridx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacC(tetrahedron_id_t /*tidx*/, uint /*ridx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacA(tetrahedron_id_t /*tidx*/, uint /*ridx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetDiffA(tetrahedron_id_t /*tidx*/, uint /*didx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetV(tetrahedron_id_t /*tidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetV(tetrahedron_id_t /*tidx*/, double /*v*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetVClamped(tetrahedron_id_t /*tidx*/) const
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetVClamped(tetrahedron_id_t /*tidx*/, bool /*cl*/)
{
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

} // namespace solver
} // namespace steps
