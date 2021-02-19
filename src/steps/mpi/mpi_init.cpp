/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

#include "easylogging++.h"
#include "steps/mpi/mpi_common.hpp"
#include "steps/mpi/mpi_init.hpp"

namespace steps {
namespace mpi {
    bool internally_initialized = false;
    void mpiInit() {
        /* Initialize MPI */
        {
            int flag;
            MPI_Initialized(&flag);
            if (!flag) {
                internally_initialized = true;
                MPI_Init(NULL, NULL);
            }
        }

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // This will replace the general log setup in serial init()
        // parallel Logger
        el::Configurations parallel_conf;

        // Global conf for the logger
        parallel_conf.set(el::Level::Global, el::ConfigurationType::Format, "[%datetime][%level][%loc][%func]: %msg");
        parallel_conf.set(el::Level::Global,
                el::ConfigurationType::ToStandardOutput, "false");
        parallel_conf.set(el::Level::Global,
                el::ConfigurationType::ToFile, "true");
        std::string file = ".logs/general_log_";
        file += std::to_string(rank);
        file += ".txt";
        parallel_conf.set(el::Level::Global,
                el::ConfigurationType::Filename, file);
        parallel_conf.set(el::Level::Global,
                el::ConfigurationType::MaxLogFileSize, "2097152");

        parallel_conf.set(el::Level::Fatal, el::ConfigurationType::ToStandardOutput, "true");
        parallel_conf.set(el::Level::Error, el::ConfigurationType::ToStandardOutput, "true");
        parallel_conf.set(el::Level::Warning, el::ConfigurationType::ToStandardOutput, "true");

        el::Loggers::getLogger("general_log");
        el::Loggers::reconfigureLogger("general_log", parallel_conf);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    int getRank() {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank;
    }

    int getNHosts() {
        int nhosts;
        MPI_Comm_size(MPI_COMM_WORLD, &nhosts);
        return nhosts;
    }
}
}
