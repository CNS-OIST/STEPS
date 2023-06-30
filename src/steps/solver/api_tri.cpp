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

// STL headers.
#include <sstream>
#include <string>

// STEPS headers.
#include "api.hpp"
#include "geom/tetmesh.hpp"
#include "statedef.hpp"
#include "util/common.hpp"
#include "util/error.hpp"
// logging
#include <easylogging++.h>

namespace steps::solver {

double API::getTriArea(triangle_global_id tidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        return _getTriArea(tidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriArea(triangle_global_id tidx, double area) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // NOTE: the following method may never be implemented
        _setTriArea(tidx, area);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSpecCount(triangle_global_id tidx, const std::string& s) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exceptions if strings are unused
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        return _getTriSpecCount(tidx, sidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriSpecDefined(triangle_global_id tidx, const std::string& s) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        return _getTriSpecDefined(tidx, sidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSpecCount(triangle_global_id tidx, const std::string& s, double n) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");
        ArgErrLogIf(n < 0.0, "Number of molecules cannot be negative.");

        // the following may raise exception if string is uknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        _setTriSpecCount(tidx, sidx, n);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSpecAmount(triangle_global_id tidx, const std::string& s) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        return _getTriSpecAmount(tidx, sidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSpecAmount(triangle_global_id tidx, const std::string& s, double m) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");
        ArgErrLogIf(m < 0.0, "Amount of mols cannot be negative.");

        // the following may throw exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        _setTriSpecAmount(tidx, sidx, m);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriSpecClamped(triangle_global_id tidx, const std::string& s) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        return _getTriSpecClamped(tidx, sidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSpecClamped(triangle_global_id tidx, const std::string& s, bool buf) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        spec_global_id sidx = pStatedef->getSpecIdx(s);

        _setTriSpecClamped(tidx, sidx, buf);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacK(triangle_global_id tidx, const std::string& r) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        sreac_global_id sridx = pStatedef->getSReacIdx(r);

        return _getTriSReacK(tidx, sridx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSReacK(triangle_global_id tidx, const std::string& r, double kf) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");
        ArgErrLogIf(kf < 0.0, "Reaction constant cannot be negative.");

        // the following may raise exception if string is unknown
        sreac_global_id sridx = pStatedef->getSReacIdx(r);

        _setTriSReacK(tidx, sridx, kf);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriSReacActive(triangle_global_id tidx, const std::string& r) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        sreac_global_id sridx = pStatedef->getSReacIdx(r);

        return _getTriSReacActive(tidx, sridx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSReacActive(triangle_global_id tidx, const std::string& r, bool act) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        sreac_global_id sridx = pStatedef->getSReacIdx(r);

        _setTriSReacActive(tidx, sridx, act);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacH(triangle_global_id tidx, const std::string& r) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        sreac_global_id sridx = pStatedef->getSReacIdx(r);

        return _getTriSReacH(tidx, sridx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacC(triangle_global_id tidx, const std::string& r) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        sreac_global_id sridx = pStatedef->getSReacIdx(r);

        return _getTriSReacC(tidx, sridx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSReacA(triangle_global_id tidx, const std::string& r) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        sreac_global_id sridx = pStatedef->getSReacIdx(r);

        return _getTriSReacA(tidx, sridx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriDiffD(triangle_global_id tidx, const std::string& d, uint direction_tri) {
    return getTriSDiffD(tidx, d, triangle_global_id(direction_tri));
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriSDiffD(triangle_global_id tidx,
                         const std::string& d,
                         triangle_global_id direction_tri) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may throw exception if string is unknown
        surfdiff_global_id didx = pStatedef->getSurfDiffIdx(d);

        return _getTriSDiffD(tidx, didx, direction_tri);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriDiffD(triangle_global_id tidx,
                      const std::string& d,
                      double dk,
                      triangle_global_id direction_tri) {
    setTriSDiffD(tidx, d, dk, direction_tri);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriSDiffD(triangle_global_id tidx,
                       const std::string& d,
                       double dk,
                       triangle_global_id direction_tri) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Tetrahedron index out of range.");
        ArgErrLogIf(direction_tri.valid() && direction_tri >= mesh->countTris(),
                    "Direction tetrahedron index out of range.");
        ArgErrLogIf(dk < 0.0, "Diffusion constant cannot be negative.");

        // the following may throw exception if string is unknown
        surfdiff_global_id didx = pStatedef->getSurfDiffIdx(d);

        _setTriSDiffD(tidx, didx, dk, direction_tri);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}
////////////////////////////////////////////////////////////////////////////////

double API::getTriV(triangle_global_id tidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        return _getTriV(tidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriV(triangle_global_id tidx, double v) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        _setTriV(tidx, v);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriVClamped(triangle_global_id tidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        return _getTriVClamped(tidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriVClamped(triangle_global_id tidx, bool cl) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        _setTriVClamped(tidx, cl);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriOhmicErev(triangle_global_id tidx, const std::string& oc) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        ohmiccurr_global_id ocidx = pStatedef->getOhmicCurrIdx(oc);

        return _getTriOhmicErev(tidx, ocidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriOhmicErev(triangle_global_id tidx, const std::string& oc, double erev) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        ohmiccurr_global_id ocidx = pStatedef->getOhmicCurrIdx(oc);

        _setTriOhmicErev(tidx, ocidx, erev);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriOhmicI(triangle_global_id tidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        return _getTriOhmicI(tidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriOhmicI(triangle_global_id tidx, const std::string& oc) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        ohmiccurr_global_id ocidx = pStatedef->getOhmicCurrIdx(oc);

        return _getTriOhmicI(tidx, ocidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriGHKI(triangle_global_id tidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        return _getTriGHKI(tidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriGHKI(triangle_global_id tidx, const std::string& ghk) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        ghkcurr_global_id ghkidx = pStatedef->getGHKcurrIdx(ghk);

        return _getTriGHKI(tidx, ghkidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}
////////////////////////////////////////////////////////////////////////////////

double API::getTriI(triangle_global_id tidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        return _getTriI(tidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriIClamp(triangle_global_id tidx) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        return _getTriIClamp(tidx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriIClamp(triangle_global_id tidx, double i) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        _setTriIClamp(tidx, i);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriCapac(triangle_global_id tidx, double cm) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        _setTriCapac(tidx, cm);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriVDepSReacActive(triangle_global_id tidx, const std::string& vsr) const {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        vdepsreac_global_id vsridx = pStatedef->getVDepSReacIdx(vsr);

        return _getTriVDepSReacActive(tidx, vsridx);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriVDepSReacActive(triangle_global_id tidx, const std::string& vsr, bool act) {
    if (auto* mesh = dynamic_cast<tetmesh::Tetmesh*>(geom())) {
        ArgErrLogIf(tidx >= mesh->countTris(), "Triangle index out of range.");

        // the following may raise exception if string is unknown
        vdepsreac_global_id vsridx = pStatedef->getVDepSReacIdx(vsr);

        _setTriVDepSReacActive(tidx, vsridx, act);
    } else {
        NotImplErrLog("Method not available for this solver.");
    }
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriArea(triangle_global_id /*tidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriArea(triangle_global_id /*tidx*/, double /*area*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriSpecDefined(triangle_global_id /*tidx*/, spec_global_id /*sidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSpecCount(triangle_global_id /*tidx*/, spec_global_id /*sidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSpecCount(triangle_global_id /*tidx*/, spec_global_id /*sidx*/, double /*n*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSpecAmount(triangle_global_id /*tidx*/, spec_global_id /*sidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSpecAmount(triangle_global_id /*tidx*/, spec_global_id /*sidx*/, double /*m*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriSpecClamped(triangle_global_id /*tidx*/, spec_global_id /*sidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSpecClamped(triangle_global_id /*tidx*/, spec_global_id /*sidx*/, bool /*buf*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacK(triangle_global_id /*tidx*/, sreac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSReacK(triangle_global_id /*tidx*/, sreac_global_id /*ridx*/, double /*kf*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriSReacActive(triangle_global_id /*tidx*/, sreac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSReacActive(triangle_global_id /*tidx*/, sreac_global_id /*ridx*/, bool /*act*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSDiffD(triangle_global_id /*tidx*/,
                          surfdiff_global_id /*didx*/,
                          triangle_global_id /*direction_tri*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriSDiffD(triangle_global_id /*tidx*/,
                        surfdiff_global_id /*didx*/,
                        double /*dk*/,
                        triangle_global_id /*direction_tri*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacH(triangle_global_id /*tidx*/, sreac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacC(triangle_global_id /*tidx*/, sreac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriSReacA(triangle_global_id /*tidx*/, sreac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriV(triangle_global_id /*tidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriV(triangle_global_id /*tidx*/, double /*v*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriVClamped(triangle_global_id /*tidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriVClamped(triangle_global_id /*tidx*/, bool /*cl*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriOhmicErev(triangle_global_id /*tidx*/, ohmiccurr_global_id /*ocidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriOhmicErev(triangle_global_id /*tidx*/,
                           ohmiccurr_global_id /*ocidx*/,
                           double /*erev*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriOhmicI(triangle_global_id /*tidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriOhmicI(triangle_global_id /*tidx*/, ohmiccurr_global_id /*ocidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriGHKI(triangle_global_id /*tidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriGHKI(triangle_global_id /*tidx*/, ghkcurr_global_id /*ghkidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriI(triangle_global_id /*tidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getTriIClamp(triangle_global_id /*tidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriIClamp(triangle_global_id /*tidx*/, double /*i*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriCapac(triangle_global_id /*tidx*/, double /*cm*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getTriVDepSReacActive(triangle_global_id /*tidx*/,
                                 vdepsreac_global_id /*vsridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setTriVDepSReacActive(triangle_global_id /*tidx*/,
                                 vdepsreac_global_id /*vsridx*/,
                                 bool /*act*/) {
    NotImplErrLog("");
}

}  // namespace steps::solver
