#pragma once

#include <mpi.h>

#include "common.hpp"
#include "type_traits.hpp"

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

template <typename T>
constexpr MPI_Datatype mpi_get_type() noexcept {
    if constexpr (std::is_same_v<T, double>) {
        return MPI_DOUBLE;
    } else if constexpr (std::is_same_v<T, osh::I64>) {
        return MPI_INT64_T;
    } else if constexpr (std::is_same_v<T, bool>) {
        return MPI_CHAR;
    } else {
        static_assert(util::always_false_v<T>, "unmanaged entity type");
    }
    return MPI_DATATYPE_NULL;
}

}  // namespace steps::util
