/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace mpi {
 namespace tetopsplit {

class KProc;

struct CRGroup {
    CRGroup(int power, uint init_size = 1024) {
        max = pow(2, power);
        sum = 0.0;
        capacity = init_size;
        size = 0;
        indices = (KProc**)malloc(sizeof(KProc*) * init_size);
        if (indices == NULL)
            throw steps::SysErr("DirectCR: unable to allocate memory for SSA group.");

        #ifdef SSA_DEBUG
        std::cout << "SSA: CRGroup Created\n";
        std::cout << "power: " << power << "\n";
        std::cout << "max: " << max << "\n";
        std::cout << "capacity: " << capacity << "\n";
        std::cout << "--------------------------------------------------------\n";
        #endif
    }

    void free_indices() {
        free(indices);
        indices = 0;
    }

    unsigned                                capacity;
    unsigned                                size;
    double                                  max;
    double                                  sum;
    KProc**                                 indices;
};

struct CRKProcData {
    CRKProcData() {
        recorded = false;
        pow = 0;
        pos = 0;
        rate = 0.0;
    }

    bool                                    recorded;
    int                                     pow;
    unsigned                                pos;
    double                                  rate;
};

////////////////////////////////////////////////////////////////////////////////

}
}
}

#endif

// STEPS_MPI_TETOPSPLIT_CRSTRUCT_HPP

// END
