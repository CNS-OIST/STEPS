/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

#ifndef STEPS_MPI_TETOPSPLIT_CRSTRUCT_HPP
#define STEPS_MPI_TETOPSPLIT_CRSTRUCT_HPP 1

#include <iostream>
#include <cmath>

#include "steps/error.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
// logging
#include <easylogging++.h>
////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace mpi {
 namespace tetopsplit {

class KProc;

struct CRGroup {
    explicit CRGroup(int power, uint init_size = 1024)
    : capacity(init_size),
      max(std::pow(2, power))
    {
        indices = static_cast<KProc**> (malloc(sizeof (KProc *) * init_size));
        if (indices == nullptr)
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
        indices = nullptr;
    }

    unsigned                                capacity;
    unsigned                                size{0};
    double                                  max;
    double                                  sum{0};
    KProc**                                 indices;
};

struct CRKProcData {
    bool                                    recorded{false};
    int                                     pow{0};
    unsigned                                pos{0};
    double                                  rate{0.0};
};

////////////////////////////////////////////////////////////////////////////////

}
}
}

#endif

// STEPS_MPI_TETOPSPLIT_CRSTRUCT_HPP

// END
