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

#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

////////////////////////////////////////////////////////////////////////////////

void API::setSDiffBoundarySpecDiffusionActive(std::string const& sdb,
                                              std::string const& s,
                                              bool act) {
    sdiffboundary_global_id sdbidx = pStatedef->getSDiffBoundaryIdx(sdb);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _setSDiffBoundarySpecDiffusionActive(sdbidx, sidx, act);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getSDiffBoundarySpecDiffusionActive(std::string const& sdb, std::string const& s) const {
    sdiffboundary_global_id sdbidx = pStatedef->getSDiffBoundaryIdx(sdb);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getSDiffBoundarySpecDiffusionActive(sdbidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setSDiffBoundarySpecDcst(std::string const& sdb,
                                   std::string const& s,
                                   double dcst,
                                   std::string const& direction_patch) {
    sdiffboundary_global_id sdbidx = pStatedef->getSDiffBoundaryIdx(sdb);
    spec_global_id sidx = pStatedef->getSpecIdx(s);
    if (direction_patch.empty()) {
        _setSDiffBoundarySpecDcst(sdbidx, sidx, dcst);
    } else {
        patch_global_id pidx = pStatedef->getPatchIdx(direction_patch);
        _setSDiffBoundarySpecDcst(sdbidx, sidx, dcst, pidx);
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::_setSDiffBoundarySpecDiffusionActive(sdiffboundary_global_id /*sdbidx*/,
                                               spec_global_id /*sidx*/,
                                               bool /*act*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getSDiffBoundarySpecDiffusionActive(sdiffboundary_global_id /*sdbidx*/,
                                               spec_global_id /*sidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setSDiffBoundarySpecDcst(sdiffboundary_global_id /*sdbidx*/,
                                    spec_global_id /*sidx*/,
                                    double /*dcst*/,
                                    patch_global_id /*direction_patch*/) {
    NotImplErrLog("");
}

}  // namespace steps::solver
