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

#include <string>

#include <mpi.h>

#include "steps/mpi/mpi_init.hpp"
#include "third_party/easylogging++.h"



void steps::mpi::mpiInit(void) {
    /* Initialize MPI */
    MPI_Init(NULL, NULL);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // MPI DEBUG Logger
    el::Configurations mpi_debug_conf;
    mpi_debug_conf.set(el::Level::Debug, el::ConfigurationType::Format, "[%datetime][%func][%loc]: %msg");
    mpi_debug_conf.set(el::Level::Debug,
                   el::ConfigurationType::ToFile, "true");
    mpi_debug_conf.set(el::Level::Debug,
                   el::ConfigurationType::ToStandardOutput, "false");
    std::string file = ".logs/mpi_debug_log_";
    file += std::to_string(rank);
    mpi_debug_conf.set(el::Level::Debug,
                   el::ConfigurationType::Filename, file);
    
    el::Loggers::getLogger("mpi_debug");
    el::Loggers::reconfigureLogger("mpi_debug", mpi_debug_conf);

    
    #ifdef MPI_DEBUG
    CLOG(DEBUG, "mpi_debug") << "######## SIMULATION START ############\n";
    #endif
    MPI_Barrier(MPI_COMM_WORLD);
}

int steps::mpi::getRank(void) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int steps::mpi::getNHosts(void) {
    int nhosts;
    MPI_Comm_size(MPI_COMM_WORLD, &nhosts);
    return nhosts;
}
