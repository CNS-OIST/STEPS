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

#pragma once

#include <cmath>
#include <iostream>

#include "kproc.hpp"

// logging
#include "util/error.hpp"
#include <easylogging++.h>

namespace steps::tetexact {

class KProc;

struct CRGroup {
    CRGroup(int power, uint init_size = 1024) {
        max = pow(2, power);
        sum = 0.0;
        capacity = init_size;
        size = 0;
        indices = static_cast<KProc**>(malloc(sizeof(KProc*) * init_size));
        if (indices == NULL)
            SysErrLog("DirectCR: unable to allocate memory for SSA group.");

#ifdef SSA_DEBUG
        CLOG(INFO, "general_log") << "SSA: CRGroup Created\n";
        CLOG(INFO, "general_log") << "power: " << power << "\n";
        CLOG(INFO, "general_log") << "max: " << max << "\n";
        CLOG(INFO, "general_log") << "capacity: " << capacity << "\n";
        CLOG(INFO, "general_log") << "--------------------------------------------------------\n";
#endif
    }

    void free_indices() {
        free(indices);
        indices = 0;
    }

    unsigned capacity;
    unsigned size;
    double max;
    double sum;
    KProc** indices;
};

struct CRKProcData {
    CRKProcData() {
        recorded = false;
        pow = 0;
        pos = 0;
        rate = 0.0;
    }

    bool recorded;
    int pow;
    unsigned pos;
    double rate;
};

}  // namespace steps::tetexact
