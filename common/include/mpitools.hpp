#pragma once

#include <mpi.h>
#include <petscsys.h>

namespace zee {

inline int mpi_comm_rank(MPI_Comm comm = MPI_COMM_WORLD) {
    int rank{};
    MPI_Comm_rank(comm, &rank);
    return rank;
}

inline int mpi_comm_size(MPI_Comm comm = MPI_COMM_WORLD) {
    int size{};
    MPI_Comm_size(comm, &size);
    return size;
}

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

}  // namespace zee
