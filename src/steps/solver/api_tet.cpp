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
#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

////////////////////////////////////////////////////////////////////////////////

double API::getTetVol(tetrahedron_global_id tidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        return _getTetVol(tidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetVol(tetrahedron_global_id tidx, double vol) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // NOTE: the following method may never be implemented
        _setTetVol(tidx, vol);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetSpecDefined(tetrahedron_global_id tidx, std::string const& s) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        return _getTetSpecDefined(tidx, sidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetSpecCount(tetrahedron_global_id tidx, std::string const& s) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        return _getTetSpecCount(tidx, sidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetSpecCount(tetrahedron_global_id tidx, std::string const& s, double n) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");
        ArgErrLogIf(n < 0.0, "Number of molecules cannot be negative.");

        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        _setTetSpecCount(tidx, sidx, n);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetSpecAmount(tetrahedron_global_id tidx, std::string const& s) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");
        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        return _getTetSpecAmount(tidx, sidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetSpecAmount(tetrahedron_global_id tidx, std::string const& s, double m) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");
        ArgErrLogIf(m < 0.0, "Amount of mols cannot be negative.");

        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        _setTetSpecAmount(tidx, sidx, m);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetSpecConc(tetrahedron_global_id tidx, std::string const& s) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        return _getTetSpecConc(tidx, sidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetSpecConc(tetrahedron_global_id tidx, std::string const& s, double c) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");
        ArgErrLogIf(c < 0.0, "Concentration cannot be negative.");

        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        _setTetSpecConc(tidx, sidx, c);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetSpecClamped(tetrahedron_global_id tidx, std::string const& s) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        return _getTetSpecClamped(tidx, sidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetSpecClamped(tetrahedron_global_id tidx, std::string const& s, bool buf) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        _setTetSpecClamped(tidx, sidx, buf);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacK(tetrahedron_global_id tidx, std::string const& r) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        reac_global_id ridx = pStatedef->getReacIdx(r);

        return _getTetReacK(tidx, ridx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetReacK(tetrahedron_global_id tidx, std::string const& r, double kf) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");
        ArgErrLogIf(kf < 0.0, "Reaction constant cannot be negative.");

        // the following may throw exception if string is unknown
        reac_global_id ridx = pStatedef->getReacIdx(r);

        _setTetReacK(tidx, ridx, kf);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetReacActive(tetrahedron_global_id tidx, std::string const& r) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        reac_global_id ridx = pStatedef->getReacIdx(r);

        return _getTetReacActive(tidx, ridx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetReacActive(tetrahedron_global_id tidx, std::string const& r, bool act) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        reac_global_id ridx = pStatedef->getReacIdx(r);

        _setTetReacActive(tidx, ridx, act);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetDiffD(tetrahedron_global_id tidx,
                        const std::string& d,
                        tetrahedron_global_id direction_tet) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        diff_global_id didx = pStatedef->getDiffIdx(d);

        return _getTetDiffD(tidx, didx, direction_tet);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetDiffD(tetrahedron_global_id tidx,
                      const std::string& d,
                      double dk,
                      tetrahedron_global_id direction_tet) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");
        ArgErrLogIf(direction_tet.valid() &&
                        direction_tet >= static_cast<index_t>(mesh->countTets()),
                    "Direction tetrahedron index out of range.");
        ArgErrLogIf(dk < 0.0, "Diffusion constant cannot be negative.");

        // the following may throw exception if string is unknown
        diff_global_id didx = pStatedef->getDiffIdx(d);

        _setTetDiffD(tidx, didx, dk, direction_tet);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetDiffActive(tetrahedron_global_id tidx, std::string const& d) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        diff_global_id didx = pStatedef->getDiffIdx(d);

        return _getTetDiffActive(tidx, didx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetDiffActive(tetrahedron_global_id tidx, std::string const& d, bool act) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        diff_global_id didx = pStatedef->getDiffIdx(d);

        _setTetDiffActive(tidx, didx, act);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacH(tetrahedron_global_id tidx, std::string const& r) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        reac_global_id ridx = pStatedef->getReacIdx(r);

        return _getTetReacH(tidx, ridx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacC(tetrahedron_global_id tidx, std::string const& r) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        reac_global_id ridx = pStatedef->getReacIdx(r);

        return _getTetReacC(tidx, ridx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetReacA(tetrahedron_global_id tidx, std::string const& r) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        reac_global_id ridx = pStatedef->getReacIdx(r);

        return _getTetReacA(tidx, ridx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetDiffA(tetrahedron_global_id tidx, std::string const& d) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        // the following may throw exception if string is unknown
        diff_global_id didx = pStatedef->getDiffIdx(d);

        return _getTetDiffA(tidx, didx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetV(tetrahedron_global_id tidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        return _getTetV(tidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetV(tetrahedron_global_id tidx, double v) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        _setTetV(tidx, v);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetVClamped(tetrahedron_global_id tidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        return _getTetVClamped(tidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetVClamped(tetrahedron_global_id tidx, bool cl) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(&geom())) {
        ArgErrLogIf(tidx >= static_cast<index_t>(mesh->countTets()),
                    "Tetrahedron index out of range.");

        _setTetVClamped(tidx, cl);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetVol(tetrahedron_global_id /*tidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetVol(tetrahedron_global_id /*tidx*/, double /*vol*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetSpecDefined(tetrahedron_global_id /*tidx*/, spec_global_id /*sidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetSpecCount(tetrahedron_global_id /*tidx*/, spec_global_id /*sidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetSpecCount(tetrahedron_global_id /*tidx*/, spec_global_id /*sidx*/, double /*n*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetSpecAmount(tetrahedron_global_id /*tidx*/, spec_global_id /*sidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetSpecAmount(tetrahedron_global_id /*tidx*/, spec_global_id /*sidx*/, double /*m*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetSpecConc(tetrahedron_global_id /*tidx*/, spec_global_id /*sidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetSpecConc(tetrahedron_global_id /*tidx*/, spec_global_id /*sidx*/, double /*c*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetSpecClamped(tetrahedron_global_id /*tidx*/, spec_global_id /*sidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetSpecClamped(tetrahedron_global_id /*tidx*/,
                             spec_global_id /*sidx*/,
                             bool /*buf*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacK(tetrahedron_global_id /*tidx*/, reac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetReacK(tetrahedron_global_id /*tidx*/, reac_global_id /*ridx*/, double /*kf*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetReacActive(tetrahedron_global_id /*tidx*/, reac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetReacActive(tetrahedron_global_id /*tidx*/, reac_global_id /*ridx*/, bool /*act*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetDiffD(tetrahedron_global_id /*tidx*/,
                         diff_global_id /*didx*/,
                         tetrahedron_global_id /*direction_tet*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetDiffD(tetrahedron_global_id /*tidx*/,
                       diff_global_id /*didx*/,
                       double /*dk*/,
                       tetrahedron_global_id /*direction_tet*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetDiffActive(tetrahedron_global_id /*tidx*/, diff_global_id /*didx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetDiffActive(tetrahedron_global_id /*tidx*/, diff_global_id /*didx*/, bool /*act*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacH(tetrahedron_global_id /*tidx*/, reac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacC(tetrahedron_global_id /*tidx*/, reac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetReacA(tetrahedron_global_id /*tidx*/, reac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetDiffA(tetrahedron_global_id /*tidx*/, diff_global_id /*didx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTetV(tetrahedron_global_id /*tidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetV(tetrahedron_global_id /*tidx*/, double /*v*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTetVClamped(tetrahedron_global_id /*tidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTetVClamped(tetrahedron_global_id /*tidx*/, bool /*cl*/) {
    NotImplErrLog("");
}

}  // namespace steps::solver
