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
#include "model/complexevents.hpp"
#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

double API::getPatchArea(std::string const& p) const {
    // the following may raise an exception if string is unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);

    return _getPatchArea(pidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchArea(std::string const& p, double area) {
    ArgErrLogIf(area <= 0.0, "Area cannot be negative or zero.");

    // the following may raise an exception if string is unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);

    _setPatchArea(pidx, area);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSpecCount(std::string const& p, std::string const& s) const {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getPatchSpecCount(pidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchSpecCount(std::string const& p, std::string const& s, double n) {
    ArgErrLogIf(n < 0.0, "Number of molecules cannot be negative.");

    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    _setPatchSpecCount(pidx, sidx, n);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSpecAmount(std::string const& p, std::string const& s) const {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getPatchSpecAmount(pidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchSpecAmount(std::string const& p, std::string const& s, double a) {
    ArgErrLogIf(a < 0.0, "Amount of mols cannot be negative.");

    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    _setPatchSpecAmount(pidx, sidx, a);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getPatchSpecClamped(std::string const& p, std::string const& s) const {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getPatchSpecClamped(pidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchSpecClamped(std::string const& p, std::string const& s, bool buf) {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    _setPatchSpecClamped(pidx, sidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchComplexCount(
    std::string const& c,
    std::string const& s,
    const std::vector<std::vector<model::SubunitStateFilter>>& f) const {
    // the following may throw exceptions if strings are unknown
    patch_global_id cidx = pStatedef->getPatchIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    return _getPatchComplexCount(cidx, sidx, _convertComplexFilters(f));
}


////////////////////////////////////////////////////////////////////////////////

void API::setPatchComplexCount(std::string const& c,
                               std::string const& s,
                               const std::vector<std::vector<model::SubunitStateFilter>>& i,
                               double n) {
    if (n < 0.0) {
        std::ostringstream os;
        os << "Number of molecules cannot be negative.";
        ArgErrLog(os.str());
    }
    // the following may throw exceptions if strings are unknown
    patch_global_id cidx = pStatedef->getPatchIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    _setPatchComplexCount(cidx, sidx, _convertComplexState(i), n);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchComplexAmount(
    std::string const& c,
    std::string const& s,
    const std::vector<std::vector<model::SubunitStateFilter>>& f) const {
    // the following may throw exceptions if strings are unknown
    patch_global_id cidx = pStatedef->getPatchIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    return _getPatchComplexAmount(cidx, sidx, _convertComplexFilters(f));
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchComplexAmount(std::string const& c,
                                std::string const& s,
                                const std::vector<std::vector<model::SubunitStateFilter>>& i,
                                double a) {
    if (a < 0.0) {
        std::ostringstream os;
        os << "Amount of mols cannot be negative.";
        ArgErrLog(os.str());
    }
    // the following may throw exceptions if strings are unknown
    patch_global_id cidx = pStatedef->getPatchIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    _setPatchComplexAmount(cidx, sidx, _convertComplexState(i), a);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchComplexSUSCount(std::string const& c,
                                    std::string const& s,
                                    const std::vector<std::vector<model::SubunitStateFilter>>& f,
                                    uint m) const {
    // the following may throw exceptions if strings are unknown
    patch_global_id cidx = pStatedef->getPatchIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    return _getPatchComplexSUSCount(cidx, sidx, _convertComplexFilters(f), complex_substate_id(m));
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacK(std::string const& p, std::string const& sr) const {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    sreac_global_id sridx = pStatedef->getSReacIdx(sr);

    return _getPatchSReacK(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchSReacK(std::string const& p, std::string const& sr, double kf) {
    ArgErrLogIf(kf < 0.0, "Reaction constant cannot be negative.");

    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    sreac_global_id sridx = pStatedef->getSReacIdx(sr);

    _setPatchSReacK(pidx, sridx, kf);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getPatchSReacActive(std::string const& p, std::string const& sr) const {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    sreac_global_id sridx = pStatedef->getSReacIdx(sr);

    return _getPatchSReacActive(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchSReacActive(std::string const& p, std::string const& sr, bool a) {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    sreac_global_id sridx = pStatedef->getSReacIdx(sr);

    _setPatchSReacActive(pidx, sridx, a);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getPatchVDepSReacActive(std::string const& p, std::string const& vsr) const {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    vdepsreac_global_id vsridx = pStatedef->getVDepSReacIdx(vsr);

    return _getPatchVDepSReacActive(pidx, vsridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchVDepSReacActive(std::string const& p, std::string const& vsr, bool a) {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    vdepsreac_global_id vsridx = pStatedef->getVDepSReacIdx(vsr);

    _setPatchVDepSReacActive(pidx, vsridx, a);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacH(std::string const& p, std::string const& sr) const {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    sreac_global_id sridx = pStatedef->getSReacIdx(sr);

    return _getPatchSReacH(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacC(std::string const& p, std::string const& sr) const {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    sreac_global_id sridx = pStatedef->getSReacIdx(sr);

    return _getPatchSReacC(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchSReacA(std::string const& p, std::string const& sr) const {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    sreac_global_id sridx = pStatedef->getSReacIdx(sr);

    return _getPatchSReacA(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long API::getPatchSReacExtent(std::string const& p, std::string const& sr) const {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    sreac_global_id sridx = pStatedef->getSReacIdx(sr);

    return _getPatchSReacExtent(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::resetPatchSReacExtent(std::string const& p, std::string const& sr) {
    // the following may raise exceptions if strings are unused
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    sreac_global_id sridx = pStatedef->getSReacIdx(sr);

    _resetPatchSReacExtent(pidx, sridx);
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long API::getPatchComplexSReacExtent(std::string const& p,
                                                   std::string const& r) const {
    // the following may throw exceptions if strings are unknown
    patch_global_id pidx = pStatedef->getPatchIdx(p);
    complexsreac_global_id ridx = pStatedef->getComplexSReacIdx(r);

    return _getPatchComplexSReacExtent(pidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::_setPatchArea(patch_global_id /*pidx*/, double /*area*/) {
    NotImplErrLog("");
}


////////////////////////////////////////////////////////////////////////

double API::_getPatchComplexCount(
    patch_global_id /*pidx*/,
    complex_global_id /*cmplIdx*/,
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& /*f*/)
    const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////

void API::_setPatchComplexCount(patch_global_id /*pidx*/,
                                complex_global_id /*cmplIdx*/,
                                const util::strongid_vector<complex_substate_id, uint>& /*i*/,
                                double /*n*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////

double API::_getPatchComplexAmount(
    patch_global_id /*pidx*/,
    complex_global_id /*sidx*/,
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& /*f*/)
    const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////

void API::_setPatchComplexAmount(patch_global_id /*pidx*/,
                                 complex_global_id /*sidx*/,
                                 const util::strongid_vector<complex_substate_id, uint>& /*i*/,
                                 double /*a*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////

double API::_getPatchComplexSUSCount(
    patch_global_id /*pidx*/,
    complex_global_id /*sidx*/,
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& /*f*/,
    complex_substate_id /*m*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getPatchSReacH(patch_global_id /*pidx*/, sreac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getPatchSReacC(patch_global_id /*pidx*/, sreac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getPatchSReacA(patch_global_id /*pidx*/, sreac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

unsigned long long API::_getPatchSReacExtent(patch_global_id /*pidx*/,
                                             sreac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_resetPatchSReacExtent(patch_global_id /*pidx*/, sreac_global_id /*ridx*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////

unsigned long long API::_getPatchComplexSReacExtent(patch_global_id /*pidx*/,
                                                    complexsreac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getPatchVDepSReacActive(patch_global_id /*pidx*/
                                   ,
                                   vdepsreac_global_id /*vsridx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setPatchVDepSReacActive(patch_global_id /*pidx*/,
                                   vdepsreac_global_id /*vsridx*/,
                                   bool /*a*/) {
    NotImplErrLog("");
}

}  // namespace steps::solver
