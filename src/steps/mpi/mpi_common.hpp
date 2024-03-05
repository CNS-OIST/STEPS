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

// Rationale:
// Store common definitions for MPI, such as communication tags

#pragma once

#include <limits.h>
#include <stdint.h>

namespace steps::mpi {

extern bool internally_initialized;

enum MsgTag {
    MPI_CONDITIONAL_BCAST = 1000,
    OPSPLIT_MOLECULE_CHANGE = 10000,
    OPSPLIT_MOLECULE_CHANGE_COMPLETE = 10001,
    OPSPLIT_COUNT_SYNC_INFO = 10100,
    OPSPLIT_COUNT_SYNC_DATA = 10101,
    OPSPLIT_SYNC_COMPLETE = 10102,
    OPSPLIT_KPROC_UPD = 10200,
    OPSPLIT_UPD_COMPLETE = 10201,

    RDEF_VESRAFT_POOL_SYNC = 20000
};

#ifdef STEPS_USE_64BITS_INDICES
#define MPI_STEPS_INDEX MPI_UNSIGNED_LONG
#else
#define MPI_STEPS_INDEX MPI_UNSIGNED
#endif

#if SIZE_MAX == UCHAR_MAX
#define MPI_STD_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define MPI_STD_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define MPI_STD_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define MPI_STD_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define MPI_STD_SIZE_T MPI_UNSIGNED_LONG_LONG
#endif

}  // namespace steps::mpi
