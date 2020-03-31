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

namespace zee {
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

}  // namespace zee
