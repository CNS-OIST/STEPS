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

/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

#include <random>
#include <vector>

#include "kproc.hpp"
// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::mpi::tetopsplit {

////////////////////////////////////////////////////////////////////////////////

KProc::KProc()
    : rExtent(0)
    , pFlags(0)
    , pSchedIDX(0u) {}

////////////////////////////////////////////////////////////////////////////////

KProc::~KProc() = default;

////////////////////////////////////////////////////////////////////////////////

void KProc::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, rExtent);
    util::checkpoint(cp_file, pFlags);
    util::checkpoint(cp_file, crData);
}

////////////////////////////////////////////////////////////////////////////////

void KProc::restore(std::fstream& cp_file) {
    util::restore(cp_file, rExtent);
    util::restore(cp_file, pFlags);
    util::restore(cp_file, crData);
}

////////////////////////////////////////////////////////////////////////////////

void KProc::setActive(bool active) {
    if (active == true) {
        pFlags &= ~INACTIVATED;
    } else {
        pFlags |= INACTIVATED;
    }
}

////////////////////////////////////////////////////////////////////////////////

unsigned long long KProc::getExtent() const {
    return rExtent;
}

////////////////////////////////////////////////////////////////////////////////

void KProc::resetExtent() {
    rExtent = 0;
}
////////////////////////////////////////////////////////////////////////////////

void KProc::resetCcst() {
    // This should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

double KProc::c() const {
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

double KProc::h() {
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

int KProc::apply(const rng::RNGptr& /*rng*/) {
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

int KProc::apply(const rng::RNGptr& /*rng*/, uint /*nmolcs*/) {
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

void KProc::apply(const rng::RNGptr& /*rng*/,
                  double /*dt*/,
                  double /*simtime*/,
                  double /*period*/) {
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

void KProc::resetOccupancies() {
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& KProc::getLocalUpdVec(int /*direction*/) const {
    // Should never get called on base object
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::kproc_global_id> const& KProc::getRemoteUpdVec(int /*direction*/) const {
    // Should never get called on base object
    AssertLog(false);
}

}  // namespace steps::mpi::tetopsplit
