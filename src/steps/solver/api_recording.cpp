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
#include "compdef.hpp"
#include "patchdef.hpp"
#include "specdef.hpp"
#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

uint API::getNComps() const {
    return pStatedef->countComps();
}

////////////////////////////////////////////////////////////////////////////////

uint API::getNPatches() const {
    return pStatedef->countPatches();
}

////////////////////////////////////////////////////////////////////////////////

std::string API::getCompName(comp_global_id c_idx) const {
    return pStatedef->compdef(c_idx)->name();
}

////////////////////////////////////////////////////////////////////////////////

std::string API::getPatchName(patch_global_id p_idx) const {
    return pStatedef->patchdef(p_idx)->name();
}

////////////////////////////////////////////////////////////////////////////////

uint API::getNCompSpecs(comp_global_id c_idx) const {
    return pStatedef->compdef(c_idx)->countSpecs();
}

////////////////////////////////////////////////////////////////////////////////

uint API::getNPatchSpecs(patch_global_id p_idx) const {
    return pStatedef->patchdef(p_idx)->countSpecs();
}

////////////////////////////////////////////////////////////////////////////////

std::string API::getCompSpecName(comp_global_id c_idx, spec_local_id s_idx) const {
    return pStatedef->specdef(pStatedef->compdef(c_idx)->specL2G(s_idx))->name();
}

////////////////////////////////////////////////////////////////////////////////

std::string API::getPatchSpecName(patch_global_id p_idx, spec_local_id s_idx) const {
    return pStatedef->specdef(pStatedef->patchdef(p_idx)->specL2G(s_idx))->name();
}

}  // namespace steps::solver
