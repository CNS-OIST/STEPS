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

#include "solver/api.hpp"

#include "solver/compdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/specdef.hpp"
#include "solver/statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

std::vector<double> API::getBatchTetCounts(const std::vector<index_t>& tets,
                                           std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getBatchTetCounts' method is deprecated and will be removed in STEPS 6, use "
           "'getBatchTetSpecCounts' instead."
        << std::endl;
    return getBatchTetSpecCounts(tets, s);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getBatchTetConcs(const std::vector<index_t>& tets,
                                          std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getBatchTetConcs' method is deprecated and will be removed in STEPS 6, use "
           "'getBatchTetSpecConcs' instead."
        << std::endl;
    return getBatchTetSpecConcs(tets, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setBatchTetConcs(const std::vector<index_t>& tets,
                           std::string const& s,
                           const std::vector<double>& concs) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setBatchTetConcs' method is deprecated and will be removed in STEPS 6, use "
           "'setBatchTetSpecConcs' instead."
        << std::endl;
    setBatchTetSpecConcs(tets, s, concs);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getBatchTriCounts(const std::vector<index_t>& tris,
                                           std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getBatchTriCounts' method is deprecated and will be removed in STEPS 6, use "
           "'getBatchTriSpecCounts' instead."
        << std::endl;
    return getBatchTriSpecCounts(tris, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::getBatchTetCountsNP(const index_t* indices,
                              size_t input_size,
                              std::string const& s,
                              double* counts,
                              size_t output_size) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getBatchTetCountsNP' method is deprecated and will be removed in STEPS 6, use "
           "'getBatchTetSpecCountsNP' instead."
        << std::endl;
    getBatchTetSpecCountsNP(indices, input_size, s, counts, output_size);
}

////////////////////////////////////////////////////////////////////////////////

void API::getBatchTetConcsNP(const index_t* indices,
                             size_t input_size,
                             std::string const& s,
                             double* counts,
                             size_t output_size) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getBatchTetConcsNP' method is deprecated and will be removed in STEPS 6, use "
           "'getBatchTetSpecConcsNP' instead."
        << std::endl;
    getBatchTetSpecConcsNP(indices, input_size, s, counts, output_size);
}

////////////////////////////////////////////////////////////////////////////////

void API::setBatchTetConcsNP(const index_t* indices,
                             size_t ntets,
                             std::string const& s,
                             const double* concs,
                             size_t output_size) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setBatchTetConcsNP' method is deprecated and will be removed in STEPS 6, use "
           "'setBatchTetSpecConcsNP' instead."
        << std::endl;
    setBatchTetSpecConcsNP(indices, ntets, s, concs, output_size);
}

////////////////////////////////////////////////////////////////////////////////

void API::getBatchTriCountsNP(const index_t* indices,
                              size_t input_size,
                              std::string const& s,
                              double* counts,
                              size_t output_size) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getBatchTriCountsNP' method is deprecated and will be removed in STEPS 6, use "
           "'getBatchTriSpecCountsNP' instead."
        << std::endl;
    return getBatchTriSpecCountsNP(indices, input_size, s, counts, output_size);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getROITetCounts(const std::string& ROI_id, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getROITetCounts' method is deprecated and will be removed in STEPS 6, use "
           "'getROITetSpecCounts' instead."
        << std::endl;
    return getROITetSpecCounts(ROI_id, s);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> API::getROITriCounts(const std::string& ROI_id, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getROITriCounts' method is deprecated and will be removed in STEPS 6, use "
           "'getROITriSpecCounts' instead."
        << std::endl;
    return getROITriSpecCounts(ROI_id, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::getROITetCountsNP(const std::string& ROI_id,
                            std::string const& s,
                            double* counts,
                            size_t output_size) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The getROITetCountsNP'' method is deprecated and will be removed in STEPS 6, use "
           "'getROITetSpecCountsNP' instead."
        << std::endl;
    getROITetSpecCountsNP(ROI_id, s, counts, output_size);
}

////////////////////////////////////////////////////////////////////////////////

void API::getROITriCountsNP(const std::string& ROI_id,
                            std::string const& s,
                            double* counts,
                            size_t output_size) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getROITriCountsNP' method is deprecated and will be removed in STEPS 6, use "
           "'getROITriSpecCountsNP' instead."
        << std::endl;
    getROITriSpecCountsNP(ROI_id, s, counts, output_size);
}

////////////////////////////////////////////////////////////////////////////////

double API::getROICount(const std::string& ROI_id, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'getROICount' method is deprecated and will be "
                                               "removed in STEPS 6, use 'getROISpecCount' instead."
                                            << std::endl;
    return getROISpecCount(ROI_id, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setROICount(const std::string& ROI_id, std::string const& s, double count) {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'setROICount' method is deprecated and will be "
                                               "removed in STEPS 6, use 'setROISpecCount' instead."
                                            << std::endl;
    setROISpecCount(ROI_id, s, count);
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIAmount(const std::string& ROI_id, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getROIAmount' method is deprecated and will be removed in STEPS 6, use "
           "'getROISpecAmount' instead."
        << std::endl;
    return getROISpecAmount(ROI_id, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIAmount(const std::string& ROI_id, std::string const& s, double amount) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setROIAmount' method is deprecated and will be removed in STEPS 6, use "
           "'setROISpecAmount' instead."
        << std::endl;
    setROISpecAmount(ROI_id, s, amount);
}

////////////////////////////////////////////////////////////////////////////////

double API::getROIConc(const std::string& ROI_id, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'getROIConc' method is deprecated and will be "
                                               "removed in STEPS 6, use 'getROISpecConc' instead."
                                            << std::endl;
    return getROISpecConc(ROI_id, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIConc(const std::string& ROI_id, std::string const& s, double conc) {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'setROIConc' method is deprecated and will be "
                                               "removed in STEPS 6, use 'setROISpecConc' instead."
                                            << std::endl;
    setROISpecConc(ROI_id, s, conc);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompCount(std::string const& c, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getCompCount' method is deprecated and will be removed in STEPS 6, use "
           "'getCompSpecCount' instead."
        << std::endl;
    return getCompSpecCount(c, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompCount(std::string const& c, std::string const& s, double n) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setCompCount' method is deprecated and will be removed in STEPS 6, use "
           "'setCompSpecCount' instead."
        << std::endl;
    setCompSpecCount(c, s, n);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompAmount(std::string const& c, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getCompAmount' method is deprecated and will be removed in STEPS 6, use "
           "'getCompSpecAmount' instead."
        << std::endl;
    return getCompSpecAmount(c, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompAmount(std::string const& c, std::string const& s, double a) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setCompAmount' method is deprecated and will be removed in STEPS 6, use "
           "'setCompSpecAmount' instead."
        << std::endl;
    setCompSpecAmount(c, s, a);
}

////////////////////////////////////////////////////////////////////////////////

double API::getCompConc(std::string const& c, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'getCompConc' method is deprecated and will be "
                                               "removed in STEPS 6, use 'getCompSpecConc' instead."
                                            << std::endl;
    return getCompSpecConc(c, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompConc(std::string const& c, std::string const& s, double conc) {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'setCompConc' method is deprecated and will be "
                                               "removed in STEPS 6, use 'setCompSpecConc' instead."
                                            << std::endl;
    setCompSpecConc(c, s, conc);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetCount(tetrahedron_global_id tidx, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'getTetCount' method is deprecated and will be "
                                               "removed in STEPS 6, use 'getTetSpecCount' instead."
                                            << std::endl;
    return getTetSpecCount(tidx, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetCount(tetrahedron_global_id tidx, std::string const& s, double n) {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'setTetCount' method is deprecated and will be "
                                               "removed in STEPS 6, use 'setTetSpecCount' instead."
                                            << std::endl;
    setTetSpecCount(tidx, s, n);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetAmount(tetrahedron_global_id tidx, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getTetAmount' method is deprecated and will be removed in STEPS 6, use "
           "'getTetSpecAmount' instead."
        << std::endl;
    return getTetSpecAmount(tidx, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetAmount(tetrahedron_global_id tidx, std::string const& s, double m) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setTetAmount' method is deprecated and will be removed in STEPS 6, use "
           "'setTetSpecAmount' instead."
        << std::endl;
    setTetSpecAmount(tidx, s, m);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTetConc(tetrahedron_global_id tidx, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'getTetConc' method is deprecated and will be "
                                               "removed in STEPS 6, use 'getTetSpecConc' instead."
                                            << std::endl;
    return getTetSpecConc(tidx, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetConc(tetrahedron_global_id tidx, std::string const& s, double c) {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'setTetConc' method is deprecated and will be "
                                               "removed in STEPS 6, use 'setTetSpecConc' instead."
                                            << std::endl;
    setTetSpecConc(tidx, s, c);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchCount(std::string const& p, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getPatchCount' method is deprecated and will be removed in STEPS 6, use "
           "'getPatchSpecCount' instead."
        << std::endl;
    return getPatchSpecCount(p, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchCount(std::string const& p, std::string const& s, double n) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setPatchCount' method is deprecated and will be removed in STEPS 6, use "
           "'setPatchSpecCount' instead."
        << std::endl;
    setPatchSpecCount(p, s, n);
}

////////////////////////////////////////////////////////////////////////////////

double API::getPatchAmount(std::string const& p, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getPatchAmount' method is deprecated and will be removed in STEPS 6, use "
           "'getPatchSpecAmount' instead."
        << std::endl;
    return getPatchSpecAmount(p, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchAmount(std::string const& p, std::string const& s, double a) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setPatchAmount' method is deprecated and will be removed in STEPS 6, use "
           "'setPatchSpecAmount' instead."
        << std::endl;
    setPatchSpecAmount(p, s, a);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriCount(triangle_global_id tidx, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'getTriCount' method is deprecated and will be "
                                               "removed in STEPS 6, use 'getTriSpecCount' instead."
                                            << std::endl;
    return getTriSpecCount(tidx, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriCount(triangle_global_id tidx, std::string const& s, double n) {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'setTriCount' method is deprecated and will be "
                                               "removed in STEPS 6, use 'setTriSpecCount' instead."
                                            << std::endl;
    setTriSpecCount(tidx, s, n);
}

////////////////////////////////////////////////////////////////////////////////

double API::getTriAmount(triangle_global_id tidx, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getTriAmount' method is deprecated and will be removed in STEPS 6, use "
           "'getTriSpecAmount' instead."
        << std::endl;
    return getTriSpecAmount(tidx, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriAmount(triangle_global_id tidx, std::string const& s, double m) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setTriAmount' method is deprecated and will be removed in STEPS 6, use "
           "'setTriSpecAmount' instead."
        << std::endl;
    setTriSpecAmount(tidx, s, m);
}

////////////////////////////////////////////////////////////////////////////////

void API::setROIClamped(const std::string& ROI_id, std::string const& s, bool b) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setROIClamped' method is deprecated and will be removed in STEPS 6, use "
           "'setROISpecClamped' instead."
        << std::endl;
    setROISpecClamped(ROI_id, s, b);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getCompClamped(std::string const& c, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getCompClamped' method is deprecated and will be removed in STEPS 6, use "
           "'getCompSpecClamped' instead."
        << std::endl;
    return getCompSpecClamped(c, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setCompClamped(std::string const& c, std::string const& s, bool b) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setCompClamped' method is deprecated and will be removed in STEPS 6, use "
           "'setCompSpecClamped' instead."
        << std::endl;
    setCompSpecClamped(c, s, b);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTetClamped(tetrahedron_global_id tidx, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getTetClamped' method is deprecated and will be removed in STEPS 6, use "
           "'getTetSpecClamped' instead."
        << std::endl;
    return getTetSpecClamped(tidx, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTetClamped(tetrahedron_global_id tidx, std::string const& s, bool buf) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setTetClamped' method is deprecated and will be removed in STEPS 6, use "
           "'setTetSpecClamped' instead."
        << std::endl;
    setTetSpecClamped(tidx, s, buf);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getPatchClamped(std::string const& p, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getPatchClamped' method is deprecated and will be removed in STEPS 6, use "
           "'getPatchSpecClamped' instead."
        << std::endl;
    return getPatchSpecClamped(p, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setPatchClamped(std::string const& p, std::string const& s, bool buf) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setPatchClamped' method is deprecated and will be removed in STEPS 6, use "
           "'setPatchSpecClamped' instead."
        << std::endl;
    setPatchSpecClamped(p, s, buf);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getTriClamped(triangle_global_id tidx, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'getTriClamped' method is deprecated and will be removed in STEPS 6, use "
           "'getTriSpecClamped' instead."
        << std::endl;
    return getTriSpecClamped(tidx, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setTriClamped(triangle_global_id tidx, std::string const& s, bool buf) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setTriClamped' method is deprecated and will be removed in STEPS 6, use "
           "'setTriSpecClamped' instead."
        << std::endl;
    setTriSpecClamped(tidx, s, buf);
}

////////////////////////////////////////////////////////////////////////////////

void API::setDiffBoundaryDiffusionActive(std::string const& db, std::string const& s, bool act) {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'setDiffBoundaryDiffusionActive' method is "
                                               "deprecated and will be removed in STEPS 6, "
                                               "use 'setDiffBoundarySpecDiffusionActive' instead."
                                            << std::endl;
    setDiffBoundarySpecDiffusionActive(db, s, act);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getDiffBoundaryDiffusionActive(std::string const& db, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'getDiffBoundaryDiffusionActive' method is "
                                               "deprecated and will be removed in STEPS 6, "
                                               "use 'getDiffBoundarySpecDiffusionActive' instead."
                                            << std::endl;
    return getDiffBoundarySpecDiffusionActive(db, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setDiffBoundaryDcst(std::string const& db,
                              std::string const& s,
                              double dcst,
                              std::string const& direction_comp) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setDiffBoundaryDcst' method is deprecated and will be removed in STEPS 6, use "
           "'setDiffBoundarySpecDcst' instead."
        << std::endl;
    setDiffBoundarySpecDcst(db, s, dcst, direction_comp);
}

////////////////////////////////////////////////////////////////////////////////

void API::setSDiffBoundaryDiffusionActive(std::string const& sdb, std::string const& s, bool act) {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'setSDiffBoundaryDiffusionActive' method is "
                                               "deprecated and will be removed in STEPS 6, "
                                               "use 'setSDiffBoundarySpecDiffusionActive' instead."
                                            << std::endl;
    setSDiffBoundarySpecDiffusionActive(sdb, s, act);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getSDiffBoundaryDiffusionActive(std::string const& sdb, std::string const& s) const {
    CLOG_N_TIMES(1, WARNING, "general_log") << "The 'getSDiffBoundaryDiffusionActive' method is "
                                               "deprecated and will be removed in STEPS 6, "
                                               "use 'getSDiffBoundarySpecDiffusionActive' instead."
                                            << std::endl;
    return getSDiffBoundarySpecDiffusionActive(sdb, s);
}

////////////////////////////////////////////////////////////////////////////////

void API::setSDiffBoundaryDcst(std::string const& sdb,
                               std::string const& s,
                               double dcst,
                               std::string const& direction_patch) {
    CLOG_N_TIMES(1, WARNING, "general_log")
        << "The 'setSDiffBoundaryDcst' method is deprecated and will be removed in STEPS 6, use "
           "'setSDiffBoundarySpecDcst' instead."
        << std::endl;
    setSDiffBoundarySpecDcst(sdb, s, dcst, direction_patch);
}

}  // namespace steps::solver
