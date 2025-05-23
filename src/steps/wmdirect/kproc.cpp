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

// STEPS headers.
#include "kproc.hpp"
// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::wmdirect {

KProc::KProc() = default;

////////////////////////////////////////////////////////////////////////////////

KProc::~KProc() = default;

////////////////////////////////////////////////////////////////////////////////

void KProc::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, rExtent);
}

////////////////////////////////////////////////////////////////////////////////

void KProc::restore(std::fstream& cp_file) {
    util::restore(cp_file, rExtent);
}

////////////////////////////////////////////////////////////////////////////////

solver::Reacdef* KProc::defr() const {
    // Should only be called on derived object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

solver::ComplexReacdef& KProc::defcr() const {
    // Should only be called on derived object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

solver::SReacdef* KProc::defsr() const {
    // Should only be called on derived object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

solver::ComplexSReacdef& KProc::defcsr() const {
    // Should only be called on derived object
    AssertLog(false);
}

}  // namespace steps::wmdirect
