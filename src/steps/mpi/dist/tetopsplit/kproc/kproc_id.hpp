#pragma once

#include "mpi/dist/tetopsplit/kproc/fwd.hpp"
#include "util/type_id.hpp"

namespace steps {
namespace dist {

using namespace util;

namespace kproc {

/**
 * \brief Class encapsulates the Kinetic process state of a STEPS simulation.
 *
 * This class encapsulates all Kinetic Processes (KProc) such as Reac
 * in the simulation.
 *
 * This class has no corresponding class in STEPS as KProc state is stored in
 * TetOpSplit/Tetexact solver class without encapsulation.
 */
constexpr unsigned num_bits(unsigned x) {
  return x < 2 ? x : 1 + num_bits(x >> 1);
}
static_assert(num_kproc_types() >= 1,
              "There must be at least one type of kin. process.");
using KProcID =
    CompactTypeId<KProcType, num_bits(num_kproc_types() - 1u), unsigned>;

} // namespace kproc
} // namespace dist
} // namespace steps
