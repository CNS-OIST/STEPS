#include "debug.hpp"

#include <climits>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <unistd.h>

#ifndef HOST_NAME_MAX
#if defined(_POSIX_HOST_NAME_MAX)
#define HOST_NAME_MAX _POSIX_HOST_NAME_MAX
#elif defined(MAXHOSTNAMELEN)
#define HOST_NAME_MAX MAXHOSTNAMELEN
#endif
#endif /* HOST_NAME_MAX */



namespace steps {
namespace util {

void wait_for_gdb() {
  int i = 0;
  char hostname[HOST_NAME_MAX];
  if (gethostname(hostname, sizeof(hostname)) != 0) {
    std::cerr << std::strerror(errno) << '\n';
    std::exit(EXIT_FAILURE);
  }
  std::cout << "PID " << getpid() << " on " << hostname
            << " ready for GDB to attach.\n"
               "When attached, pause execution and update variable i to exit "
               "the infinite loop."
            << std::endl;
  while (0 == i) {
    sleep(5);
  }
}

} // namespace util
} // namespace steps

#if USE_PETSC

std::ostream& print_comm_stats(std::ostream& os, const PetscObject& obj) {
    MPI_Comm comm;
    PetscObjectGetComm(obj, &comm);

    int rank{};
    MPI_Comm_rank(comm, &rank);
    int size{};
    MPI_Comm_size(comm, &size);
    os << "rank/size: " << rank << '/' << size;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Vec& v) {
    VecAssemblyBegin(v);
    VecAssemblyEnd(v);

    PetscInt N;
    VecGetSize(v, &N);
    os << '(' << N << "): \n";
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_DENSE);
    VecView(v, PETSC_VIEWER_STDOUT_WORLD);
    return os << '\n';
}

std::ostream& operator<<(std::ostream& os, const Mat& m) {
    MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);

    PetscInt M, N;
    MatGetSize(m, &M, &N);
    os << '(' << M << 'x' << N << "): \n";
    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_DENSE);
    MatView(m, PETSC_VIEWER_STDOUT_WORLD);
    return os << '\n';
}

#endif  // USE_PETSC
