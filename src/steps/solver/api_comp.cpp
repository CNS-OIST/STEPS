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

////////////////////////////////////////////////////////////////////////////////

std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>
_convertComplexFilters(const std::vector<std::vector<model::SubunitStateFilter>>& f) {
    std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>> filts;
    filts.reserve(f.size());
    for (auto& filt: f) {
        filts.emplace_back(filt);
    }
    return filts;
}

////////////////////////////////////////////////////////////////////////////////

util::strongid_vector<complex_substate_id, uint> _convertComplexState(
    const std::vector<std::vector<model::SubunitStateFilter>>& f) {
    util::strongid_vector<complex_substate_id, uint> state;
    AssertLog(f.size() == 1);
    state.container().reserve(f[0].size());
    for (const auto& susfilt: f[0]) {
        AssertLog(susfilt.min == susfilt.max);
        state.container().push_back(susfilt.min);
    }
    return state;
}

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

double API::getCompComplexCount(
    std::string const& c,
    std::string const& s,
    const std::vector<std::vector<model::SubunitStateFilter>>& f) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    return _getCompComplexCount(cidx, sidx, _convertComplexFilters(f));
}


////////////////////////////////////////////////////////////////////////////////

void API::setCompComplexCount(std::string const& c,
                              std::string const& s,
                              const std::vector<std::vector<model::SubunitStateFilter>>& i,
                              double n) {
    if (n < 0.0) {
        std::ostringstream os;
        os << "Number of molecules cannot be negative.";
        ArgErrLog(os.str());
    }
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    _setCompComplexCount(cidx, sidx, _convertComplexState(i), n);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompComplexAmount(
    std::string const& c,
    std::string const& s,
    const std::vector<std::vector<model::SubunitStateFilter>>& f) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    return _getCompComplexAmount(cidx, sidx, _convertComplexFilters(f));
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompComplexAmount(std::string const& c,
                               std::string const& s,
                               const std::vector<std::vector<model::SubunitStateFilter>>& i,
                               double a) {
    if (a < 0.0) {
        std::ostringstream os;
        os << "Amount of mols cannot be negative.";
        ArgErrLog(os.str());
    }
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    _setCompComplexAmount(cidx, sidx, _convertComplexState(i), a);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompComplexConc(std::string const& c,
                               std::string const& s,
                               const std::vector<std::vector<model::SubunitStateFilter>>& f) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    return _getCompComplexConc(cidx, sidx, _convertComplexFilters(f));
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompComplexConc(std::string const& c,
                             std::string const& s,
                             const std::vector<std::vector<model::SubunitStateFilter>>& i,
                             double conc) {
    if (conc < 0.0) {
        std::ostringstream os;
        os << "Concentration cannot be negative.";
        ArgErrLog(os.str());
    }
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    _setCompComplexConc(cidx, sidx, _convertComplexState(i), conc);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompComplexSUSCount(std::string const& c,
                                   std::string const& s,
                                   const std::vector<std::vector<model::SubunitStateFilter>>& f,
                                   uint m) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    return _getCompComplexSUSCount(cidx, sidx, _convertComplexFilters(f), complex_substate_id(m));
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompComplexSUSConc(std::string const& c,
                                  std::string const& s,
                                  const std::vector<std::vector<model::SubunitStateFilter>>& f,
                                  uint m) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    return _getCompComplexSUSConc(cidx, sidx, _convertComplexFilters(f), complex_substate_id(m));
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompComplexSUSAmount(std::string const& c,
                                    std::string const& s,
                                    const std::vector<std::vector<model::SubunitStateFilter>>& f,
                                    uint m) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    complex_global_id sidx = pStatedef->getComplexIdx(s);

    return _getCompComplexSUSAmount(cidx, sidx, _convertComplexFilters(f), complex_substate_id(m));
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

unsigned long long API::getCompComplexReacExtent(std::string const& c, std::string const& r) const {
    // the following may throw exceptions if strings are unknown
    comp_global_id cidx = pStatedef->getCompIdx(c);
    complexreac_global_id ridx = pStatedef->getComplexReacIdx(r);

    return _getCompComplexReacExtent(cidx, ridx);
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompVol(comp_global_id /* cidx */, double /* vol */) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getCompComplexCount(
    comp_global_id /*cidx*/,
    complex_global_id /*sidx*/,
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& /*f*/)
    const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompComplexCount(comp_global_id /*cidx*/,
                               complex_global_id /*sidx*/,
                               const util::strongid_vector<complex_substate_id, uint>& /*f*/,
                               double /*n*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getCompComplexAmount(
    comp_global_id /*cidx*/,
    complex_global_id /*sidx*/,
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& /*f*/)
    const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompComplexAmount(comp_global_id /*cidx*/,
                                complex_global_id /*sidx*/,
                                const util::strongid_vector<complex_substate_id, uint>& /*f*/,
                                double /*a*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////


double API::_getCompComplexConc(
    comp_global_id /*cidx*/,
    complex_global_id /*sidx*/,
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& /*f*/)
    const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setCompComplexConc(comp_global_id /*cidx*/,
                              complex_global_id /*sidx*/,
                              const util::strongid_vector<complex_substate_id, uint>& /*f*/,
                              double /*c*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getCompComplexSUSCount(
    comp_global_id /*cidx*/,
    complex_global_id /*sidx*/,
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& /*f*/,
    complex_substate_id /*m*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getCompComplexSUSConc(
    comp_global_id /*cidx*/,
    complex_global_id /*sidx*/,
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& /*f*/,
    complex_substate_id /*m*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

double API::_getCompComplexSUSAmount(
    comp_global_id /*cidx*/,
    complex_global_id /*sidx*/,
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& /*f*/,
    complex_substate_id /*m*/) const {
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

////////////////////////////////////////////////////////////////////////////////

unsigned long long API::_getCompComplexReacExtent(comp_global_id /*cidx*/,
                                                  complexreac_global_id /*ridx*/) const {
    NotImplErrLog("");
}

}  // namespace steps::solver
