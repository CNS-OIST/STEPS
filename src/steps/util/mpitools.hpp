#pragma once

#include <mpi.h>

#ifdef USE_PETSC
#include <petscsys.h>
#endif  // USE_PETSC

namespace steps::util {


// this is duplicated code, see steps/mpi/mpi_init.hpp
inline int mpi_comm_rank(MPI_Comm comm) {
    int rank{};
    MPI_Comm_rank(comm, &rank);
    return rank;
}

inline int mpi_comm_size(MPI_Comm comm) {
    int size{};
    MPI_Comm_size(comm, &size);
    return size;
}

#ifdef USE_PETSC
struct PetscFixture {
    PetscFixture(int* argc, char*** argv, const char file[], const char help[]) {
        PetscErrorCode ierr = PetscInitialize(argc, argv, file, help);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
    }
    ~PetscFixture() {
        PetscErrorCode ierr = PetscFinalize();
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
    }
};
#endif  // USE_PETSC

}  // namespace steps::util
