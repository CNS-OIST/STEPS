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
// util
#include "util/error.hpp"
// logging
#include <easylogging++.h>

namespace steps::solver {

////////////////////////////////////////////////////////////////////////////////

void API::setDiffBoundarySpecDiffusionActive(std::string const& db,
                                             std::string const& s,
                                             bool act) {
    diffboundary_global_id dbidx = pStatedef->getDiffBoundaryIdx(db);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _setDiffBoundarySpecDiffusionActive(dbidx, sidx, act);
}

////////////////////////////////////////////////////////////////////////////////

bool API::getDiffBoundarySpecDiffusionActive(std::string const& db, std::string const& s) const {
    diffboundary_global_id dbidx = pStatedef->getDiffBoundaryIdx(db);
    spec_global_id sidx = pStatedef->getSpecIdx(s);

    return _getDiffBoundarySpecDiffusionActive(dbidx, sidx);
}

////////////////////////////////////////////////////////////////////////////////

void API::setDiffBoundarySpecDcst(std::string const& db,
                                  std::string const& s,
                                  double dcst,
                                  std::string const& direction_comp) {
    diffboundary_global_id dbidx = pStatedef->getDiffBoundaryIdx(db);
    spec_global_id sidx = pStatedef->getSpecIdx(s);
    if (direction_comp.empty()) {
        _setDiffBoundarySpecDcst(dbidx, sidx, dcst);
    } else {
        comp_global_id cidx = pStatedef->getCompIdx(direction_comp);
        _setDiffBoundarySpecDcst(dbidx, sidx, dcst, cidx);
    }
}

////////////////////////////////////////////////////////////////////////////////

void API::_setDiffBoundarySpecDiffusionActive(diffboundary_global_id /*dbidx*/,
                                              spec_global_id /*sidx*/,
                                              bool /*act*/) {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

bool API::_getDiffBoundarySpecDiffusionActive(diffboundary_global_id /*dbidx*/,
                                              spec_global_id /*sidx*/) const {
    NotImplErrLog("");
}

////////////////////////////////////////////////////////////////////////////////

void API::_setDiffBoundarySpecDcst(diffboundary_global_id /*dbidx*/,
                                   spec_global_id /*sidx*/,
                                   double /*dcst*/,
                                   comp_global_id /*direction_comp*/) {
    NotImplErrLog("");
}

}  // namespace steps::solver
