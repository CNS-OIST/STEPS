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
#include "statedef.hpp"
#include "util/error.hpp"

// logging
#include <easylogging++.h>

namespace steps::solver {

////////////////////////////////////////////////////////////////////////////////

double API::getCompVol(std::string const& c) const {
    // the following may throw an exception if std::string is unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);

    return _getCompVol(cidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompVol(std::string const& c, double vol) {
    if (vol <= 0.0) {
        std::ostringstream os;
        os << "Volume cannot be negative or zero.";
        ArgErrLog(os.str());
    }
    // the following may throw an exception if string is unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);

    _setCompVol(cidx, vol);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompSpecCount(std::string const& c, std::string const& s) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getCompSpecCount(cidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompSpecCount(std::string const& c, std::string const& s, double n) {
    ArgErrLogIf(n < 0.0, "Number of molecules cannot be negative.");

    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    _setCompSpecCount(cidx, sidx, n);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompSpecAmount(std::string const& c, std::string const& s) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getCompSpecAmount(cidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompSpecAmount(std::string const& c, std::string const& s, double a) {
    ArgErrLogIf(a < 0.0, "Amount of mols cannot be negative.");

    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    _setCompSpecAmount(cidx, sidx, a);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompSpecConc(std::string const& c, std::string const& s) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getCompSpecConc(cidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompSpecConc(std::string const& c, std::string const& s, double conc) {
    ArgErrLogIf(conc < 0.0, "Concentration cannot be negative.");

    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    _setCompSpecConc(cidx, sidx, conc);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getCompSpecClamped(std::string const& c, std::string const& s) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getCompSpecClamped(cidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompSpecClamped(std::string const& c, std::string const& s, bool b) {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    _setCompSpecClamped(cidx, sidx, b);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompReacK(std::string const& c, std::string const& r) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    reac_global_id ridx = pStatedef->getReacIdx(r);

    return _getCompReacK(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompReacK(std::string const& c, std::string const& r, double kf) {
    ArgErrLogIf(kf < 0.0, "Reaction constant cannot be negative.");

    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    reac_global_id ridx = pStatedef->getReacIdx(r);

    _setCompReacK(cidx, ridx, kf);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getCompReacActive(std::string const& c, std::string const& r) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    reac_global_id ridx = pStatedef->getReacIdx(r);

    return _getCompReacActive(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompReacActive(std::string const& c, std::string const& r, bool a) {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    reac_global_id ridx = pStatedef->getReacIdx(r);

    _setCompReacActive(cidx, ridx, a);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompDiffD(std::string const& c, std::string const& d) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    diff_global_id didx = pStatedef->getDiffIdx(d);

    return _getCompDiffD(cidx, didx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompDiffD(std::string const& c, std::string const& d, double dcst) {
    ArgErrLogIf(dcst < 0.0, "Diffusion constant cannot be negative.");

    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    diff_global_id didx = pStatedef->getDiffIdx(d);

    _setCompDiffD(cidx, didx, dcst);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getCompDiffActive(std::string const& c, std::string const& d) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    diff_global_id didx = pStatedef->getDiffIdx(d);

    return _getCompDiffActive(cidx, didx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompDiffActive(std::string const& c, std::string const& d, bool act) {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    diff_global_id didx = pStatedef->getDiffIdx(d);

    _setCompDiffActive(cidx, didx, act);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompReacH(std::string const& c, std::string const& r) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    reac_global_id ridx = pStatedef->getReacIdx(r);

    return _getCompReacH(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompReacC(std::string const& c, std::string const& r) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    reac_global_id ridx = pStatedef->getReacIdx(r);

    return _getCompReacC(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompReacA(std::string const& c, std::string const& r) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    reac_global_id ridx = pStatedef->getReacIdx(r);

    return _getCompReacA(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long API::getCompReacExtent(std::string const& c, std::string const& r) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    reac_global_id ridx = pStatedef->getReacIdx(r);

    return _getCompReacExtent(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::resetCompReacExtent(std::string const& c, std::string const& r) {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    reac_global_id ridx = pStatedef->getReacIdx(r);

    _resetCompReacExtent(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompVol(comp_global_id /* cidx */, double /* vol */) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getCompDiffD(comp_global_id /*cidx*/, diff_global_id /*didx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompDiffD(comp_global_id /*cidx*/, diff_global_id /*didx*/, double /*dcst*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getCompDiffActive(comp_global_id /*cidx*/, diff_global_id /*didx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompDiffActive(comp_global_id /*cidx*/, diff_global_id /*didx*/, bool /*act*/) {
    NotImplErrLog("");
}
double API::_getCompReacH(comp_global_id /*cidx*/, reac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getCompReacC(comp_global_id /*cidx*/, reac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

long double API::_getCompReacA(comp_global_id /*cidx*/, reac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long API::_getCompReacExtent(comp_global_id /*cidx*/, reac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_resetCompReacExtent(comp_global_id /*cidx*/, reac_global_id /*ridx*/) {
    NotImplErrLog("");
}

}  // namespace steps::solver
