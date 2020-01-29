#ifndef STEPS_GEOM_FWD_HPP
#define STEPS_GEOM_FWD_HPP

#include <limits>
#include <iosfwd>
#include <type_traits>

#include <steps/util/strong_id.hpp>

namespace steps {

#ifdef STEPS_USE_64BITS_INDICES
using index_t = unsigned long;
#else
using index_t = unsigned int;
#endif

using tetrahedron_id_t = strong_id<index_t, struct tetrahedron_id_trait>;
using triangle_id_t = strong_id<index_t, struct triangle_id_trait>;
using vertex_id_t = strong_id<index_t, struct vertex_id_trait>;
using bar_id_t = strong_id<index_t, struct bar_id_trait>;
/// TODO TCL add host_id_t primitive type

static const tetrahedron_id_t UNKNOWN_TET(std::numeric_limits<index_t>::max());
static const triangle_id_t UNKNOWN_TRI(std::numeric_limits<index_t>::max());
static const vertex_id_t UNKNOWN_VER(std::numeric_limits<index_t>::max());
static const bar_id_t UNKNOWN_BAR(std::numeric_limits<index_t>::max());

namespace tetmesh {

// Forward declarations.
class Tetmesh;
class Memb;
class TmPatch;

} // namespace tetmesh

} // namespace steps

#endif //!STEPS_GEOM_FWD_HPP
