#pragma once

#include <limits>
#include <type_traits>

#include "math/point.hpp"
#include "util/common.hpp"
#include "util/strong_id.hpp"
#include "util/strong_ra.hpp"

namespace steps {

inline constexpr index_t UNKNOWN_INDEX(std::numeric_limits<index_t>::max());

namespace solver {

template <typename StrongId, typename ValueType>
using strongid_span = util::strong_random_access<StrongId, gsl::span<ValueType>>;

using spec_global_id = util::strong_id<index_t, struct spec_global_id_trait>;
using spec_local_id = util::strong_id<index_t, struct spec_local_id_trait>;
using chan_global_id = util::strong_id<index_t, struct chan_global_id_trait>;
using chan_local_id = util::strong_id<index_t, struct chan_local_id_trait>;

using complex_global_id = util::strong_id<index_t, struct complex_global_id_trait>;
using complex_local_id = util::strong_id<index_t, struct complex_local_id_trait>;
using complex_substate_id = util::strong_id<index_t, struct complex_substate_id_trait>;
struct complex_individual_id_trait {};
using complex_individual_id = util::strong_id<index_t, complex_individual_id_trait>;
using complex_filter_id = util::strong_id<index_t, struct complex_filter_id_trait>;

using kproc_global_id = util::strong_id<index_t, struct kproc_global_id_trait>;

using reac_global_id = util::strong_id<index_t, struct reac_global_id_trait>;
using reac_local_id = util::strong_id<index_t, struct reac_local_id_trait>;
using complexreac_global_id = util::strong_id<index_t, struct complexreac_global_id_trait>;
using complexreac_local_id = util::strong_id<index_t, struct complexreac_local_id_trait>;
using diff_global_id = util::strong_id<index_t, struct diff_global_id_trait>;
using diff_local_id = util::strong_id<index_t, struct diff_local_id_trait>;
using surfdiff_global_id = util::strong_id<index_t, struct surfdiff_global_id_trait>;
using surfdiff_local_id = util::strong_id<index_t, struct surfdiff_local_id_trait>;
using sreac_global_id = util::strong_id<index_t, struct sreac_global_id_trait>;
using sreac_local_id = util::strong_id<index_t, struct sreac_local_id_trait>;
using complexsreac_global_id = util::strong_id<index_t, struct complexsreac_global_id_trait>;
using complexsreac_local_id = util::strong_id<index_t, struct complexsreac_local_id_trait>;
using vdepsreac_global_id = util::strong_id<index_t, struct vdepsreac_global_id_trait>;
using vdepsreac_local_id = util::strong_id<index_t, struct vdepsreac_local_id_trait>;
using ohmiccurr_global_id = util::strong_id<index_t, struct ohmiccurr_global_id_trait>;
using ohmiccurr_local_id = util::strong_id<index_t, struct ohmiccurr_local_id_trait>;
using ghkcurr_global_id = util::strong_id<index_t, struct ghkcurr_global_id_trait>;
using ghkcurr_local_id = util::strong_id<index_t, struct ghkcurr_local_id_trait>;

using vesicle_global_id = util::strong_id<index_t, struct vesicle_global_id_trait>;
// using vesicle_local_id = util::strong_id<index_t, struct vesicle_local_id_trait>;
using vesproxy_global_id = util::strong_id<index_t, struct vesproxy_global_id_trait>;
using vesicle_individual_id = util::strong_id<index_t, struct vesicle_individual_id_trait>;
using linkspec_global_id = util::strong_id<index_t, struct linkspec_global_id_trait>;
using linkspec_local_id = util::strong_id<index_t, struct linkspec_local_id_trait>;
using linkspec_individual_id = util::strong_id<index_t, struct linkspec_individual_id_trait>;
using pointspec_individual_id = util::strong_id<index_t, struct pointspec_individual_id_trait>;

using vesbind_global_id = util::strong_id<index_t, struct vesbind_global_id_trait>;
using vesbind_local_id = util::strong_id<index_t, struct vesbind_local_id_trait>;
using vesunbind_global_id = util::strong_id<index_t, struct vesunbind_global_id_trait>;
using vesunbind_local_id = util::strong_id<index_t, struct vesunbind_local_id_trait>;
using endocytosis_global_id = util::strong_id<index_t, struct endocytosis_global_id_trait>;
using endocytosis_local_id = util::strong_id<index_t, struct endocytosis_local_id_trait>;
using exocytosis_global_id = util::strong_id<index_t, struct exocytosis_global_id_trait>;
using exocytosis_local_id = util::strong_id<index_t, struct exocytosis_local_id_trait>;
using vessreac_global_id = util::strong_id<index_t, struct vessreac_global_id_trait>;
using vessreac_local_id = util::strong_id<index_t, struct vessreac_local_id_trait>;
using vessdiff_global_id = util::strong_id<index_t, struct vessdiff_global_id_trait>;
using vessdiff_local_id = util::strong_id<index_t, struct vessdiff_local_id_trait>;

using raft_global_id = util::strong_id<index_t, struct raft_global_id_trait>;
using raft_individual_id = util::strong_id<index_t, struct raft_individual_id_trait>;
using raftsreac_global_id = util::strong_id<index_t, struct raftsreac_global_id_trait>;
using raftsreac_local_id = util::strong_id<index_t, struct raftsreac_local_id_trait>;
using raftendocytosis_global_id = util::strong_id<index_t, struct raftendocytosis_global_id_trait>;
using raftendocytosis_local_id = util::strong_id<index_t, struct raftendocytosis_local_id_trait>;
using raftdis_global_id = util::strong_id<index_t, struct raftdis_global_id_trait>;
using raftdis_local_id = util::strong_id<index_t, struct raftdis_local_id_trait>;
using raftgen_global_id = util::strong_id<index_t, struct raftgen_global_id_trait>;
using raftgen_local_id = util::strong_id<index_t, struct raftgen_local_id_trait>;

using comp_global_id = util::strong_id<index_t, struct comp_global_id_trait>;
using patch_global_id = util::strong_id<index_t, struct patch_global_id_trait>;
using membrane_global_id = util::strong_id<index_t, struct membrane_global_id_trait>;

using diffboundary_global_id = util::strong_id<index_t, struct diffboundary_global_id_trait>;
using sdiffboundary_global_id = util::strong_id<index_t, struct sdiffboundary_global_id_trait>;


typedef std::vector<spec_global_id> spec_global_id_vec;
typedef spec_global_id_vec::const_iterator spec_global_id_vecCI;

typedef std::vector<linkspec_global_id> linkspec_global_id_vec;

typedef int depT;

static const depT DEP_NONE = 0;
static const depT DEP_STOICH = 1;
static const depT DEP_RATE = 2;

// Used for near-zero volumes in TetVesicle solver
static const double TINY_VOLUME = 1e-100;


// Forwards declarations
class Statedef;

class Compdef;
class ComplexFilter;
class ComplexState;
class Patchdef;
class EndocyticZonedef;
class Specdef;
class Complexdef;
class Reacdef;
class ComplexReacdef;
class SReacdef;
class ComplexSReacdef;
class Chandef;
class VDepSReacdef;
class OhmicCurrdef;
class GHKcurrdef;
class DiffBoundarydef;
class SDiffBoundarydef;

template <typename GlobalId>
class MetaDiffdef;
using Diffdef = MetaDiffdef<diff_global_id>;
using SurfDiffdef = MetaDiffdef<surfdiff_global_id>;

class Vesicledef;
class LinkSpecdef;
class Endocytosisdef;
class EndocytosisEvent;
class Exocytosisdef;
class ExocytosisEvent;
class VesBinddef;
class VesUnbinddef;
class VesSDiffdef;
class VesSReacdef;
class Raftdef;
class RaftSReacdef;
class RaftEndocytosisdef;
class RaftEndocytosisEvent;
class RaftGendef;
class RaftDisdef;

}  // namespace solver

}  // namespace steps
