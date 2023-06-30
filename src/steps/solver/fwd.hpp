#pragma once

#include <iosfwd>
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

using kproc_global_id = util::strong_id<index_t, struct kproc_global_id_trait>;

using reac_global_id = util::strong_id<index_t, struct reac_global_id_trait>;
using reac_local_id = util::strong_id<index_t, struct reac_local_id_trait>;
using diff_global_id = util::strong_id<index_t, struct diff_global_id_trait>;
using diff_local_id = util::strong_id<index_t, struct diff_local_id_trait>;
using surfdiff_global_id = util::strong_id<index_t, struct surfdiff_global_id_trait>;
using surfdiff_local_id = util::strong_id<index_t, struct surfdiff_local_id_trait>;
using sreac_global_id = util::strong_id<index_t, struct sreac_global_id_trait>;
using sreac_local_id = util::strong_id<index_t, struct sreac_local_id_trait>;
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
typedef spec_global_id_vec::iterator spec_global_id_vecI;
typedef spec_global_id_vec::const_iterator spec_global_id_vecCI;

typedef std::vector<linkspec_global_id> linkspec_global_id_vec;
typedef linkspec_global_id_vec::iterator linkspec_global_id_vecI;
typedef linkspec_global_id_vec::const_iterator linkspec_global_id_vecCI;

typedef int depT;
typedef std::vector<depT> depTVec;
typedef depTVec::iterator depTVecI;
typedef depTVec::const_iterator depTVecCI;

static const depT DEP_NONE = 0;
static const depT DEP_STOICH = 1;
static const depT DEP_RATE = 2;

// Used for near-zero volumes in TetVesicle solver
static const double TINY_VOLUME = 1e-100;


// Forwards declarations
class Statedef;

class Compdef;
class Patchdef;

class Specdef;
class Reacdef;
class SReacdef;
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

// Auxiliary declarations.

typedef Specdef* SpecdefP;
typedef std::vector<SpecdefP> SpecdefPVec;
typedef SpecdefPVec::iterator SpecdefPVecI;
typedef SpecdefPVec::const_iterator SpecdefPVecCI;

typedef Compdef* CompdefP;
typedef std::vector<CompdefP> CompdefPVec;
typedef CompdefPVec::iterator CompdefPVecI;
typedef CompdefPVec::const_iterator CompdefPVecCI;

typedef Patchdef* PatchdefP;
typedef std::vector<PatchdefP> PatchdefPVec;
typedef PatchdefPVec::iterator PatchdefPVecI;
typedef PatchdefPVec::const_iterator PatchdefPVecCI;

typedef Reacdef* ReacdefP;
typedef std::vector<ReacdefP> ReacdefPVec;
typedef ReacdefPVec::iterator ReacdefPVecI;
typedef ReacdefPVec::const_iterator ReacdefPVecCI;

typedef SReacdef* SReacdefP;
typedef std::vector<SReacdefP> SReacdefPVec;
typedef SReacdefPVec::iterator SReacdefPVecI;
typedef SReacdefPVec::const_iterator SReacdefPVecCI;

typedef Chandef* ChandefP;
typedef std::vector<ChandefP> ChandefPVec;
typedef ChandefPVec::iterator ChandefPVecI;
typedef ChandefPVec::const_iterator ChandefPVecCI;

typedef VDepSReacdef* VDepSReacdefP;
typedef std::vector<VDepSReacdefP> VDepSReacdefPVec;
typedef VDepSReacdefPVec::iterator VDepSReacdefPVecI;
typedef VDepSReacdefPVec::const_iterator VDepSReacdefPVecCI;

typedef OhmicCurrdef* OhmicCurrdefP;
typedef std::vector<OhmicCurrdefP> OhmicCurrdefPVec;
typedef OhmicCurrdefPVec::iterator OhmicCurrdefPVecI;
typedef OhmicCurrdefPVec::const_iterator OhmicCurrdefPVecCI;

typedef GHKcurrdef* GHKcurrdefP;
typedef std::vector<GHKcurrdefP> GHKcurrdefPVec;
typedef GHKcurrdefPVec::iterator GHKcurrdefPVecI;
typedef GHKcurrdefPVec::const_iterator GHKcurrdefPVecCI;

typedef DiffBoundarydef* DiffBoundarydefP;
typedef std::vector<DiffBoundarydefP> DiffBoundarydefPVec;
typedef DiffBoundarydefPVec::iterator DiffBoundarydefPVecI;
typedef DiffBoundarydefPVec::const_iterator DiffBoundarydefPVecCI;

typedef SDiffBoundarydef* SDiffBoundarydefP;
typedef std::vector<SDiffBoundarydefP> SDiffBoundarydefPVec;
typedef SDiffBoundarydefPVec::iterator SDiffBoundarydefPVecI;
typedef SDiffBoundarydefPVec::const_iterator SDiffBoundarydefPVecCI;

typedef LinkSpecdef* LinkSpecdefP;
typedef std::vector<LinkSpecdefP> LinkSpecdefPVec;
typedef LinkSpecdefPVec::iterator LinkSpecdefPVecI;
typedef LinkSpecdefPVec::const_iterator LinkSpecdefPVecCI;

typedef Vesicledef* VesicledefP;
typedef std::vector<VesicledefP> VesicledefPVec;
typedef VesicledefPVec::iterator VesicledefPVecI;
typedef VesicledefPVec::const_iterator VesicledefPVecCI;

typedef Endocytosisdef* EndocytosisdefP;
typedef std::vector<EndocytosisdefP> EndocytosisdefPVec;
typedef EndocytosisdefPVec::iterator EndocytosisdefPVecI;
typedef EndocytosisdefPVec::const_iterator EndocytosisdefPVecCI;

typedef Exocytosisdef* ExocytosisdefP;
typedef std::vector<ExocytosisdefP> ExocytosisdefPVec;
typedef ExocytosisdefPVec::iterator ExocytosisdefPVecI;
typedef ExocytosisdefPVec::const_iterator ExocytosisdefPVecCI;

typedef VesBinddef* VesBinddefP;
typedef std::vector<VesBinddefP> VesBinddefPVec;
typedef VesBinddefPVec::iterator VesBinddefPVecI;
typedef VesBinddefPVec::const_iterator VesBinddefPVecCI;

typedef VesUnbinddef* VesUnbinddefP;
typedef std::vector<VesUnbinddefP> VesUnbinddefPVec;
typedef VesUnbinddefPVec::iterator VesUnbinddefPVecI;
typedef VesUnbinddefPVec::const_iterator VesUnbinddefPVecCI;

typedef VesSDiffdef* VesSDiffdefP;
typedef std::vector<VesSDiffdefP> VesSDiffdefPVec;
typedef VesSDiffdefPVec::iterator VesSDiffdefPVecI;
typedef VesSDiffdefPVec::const_iterator VesSDiffdefPVecCI;

typedef VesSReacdef* VesSReacdefP;
typedef std::vector<VesSReacdefP> VesSReacdefPVec;
typedef VesSReacdefPVec::iterator VesSReacdefPVecI;
typedef VesSReacdefPVec::const_iterator VesSReacdefPVecCI;

typedef Raftdef* RaftdefP;
typedef std::vector<RaftdefP> RaftdefPVec;
typedef RaftdefPVec::iterator RaftdefPVecI;
typedef RaftdefPVec::const_iterator RaftdefPVecCI;

typedef RaftSReacdef* RaftSReacdefP;
typedef std::vector<RaftSReacdefP> RaftSReacdefPVec;
typedef RaftSReacdefPVec::iterator RaftSReacdefPVecI;
typedef RaftSReacdefPVec::const_iterator RaftSReacdefPVecCI;

typedef RaftEndocytosisdef* RaftEndocytosisdefP;
typedef std::vector<RaftEndocytosisdefP> RaftEndocytosisdefPVec;
typedef RaftEndocytosisdefPVec::iterator RaftEndocytosisdefPVecI;
typedef RaftEndocytosisdefPVec::const_iterator RaftEndocytosisdefPVecCI;

typedef RaftGendef* RaftGendefP;
typedef std::vector<RaftGendefP> RaftGendefPVec;
typedef RaftGendefPVec::iterator RaftGendefPVecI;
typedef RaftGendefPVec::const_iterator RaftGendefPVecCI;

typedef RaftDisdef* RaftDisdefP;
typedef std::vector<RaftDisdefP> RaftDisdefPVec;
typedef RaftDisdefPVec::iterator RaftDisdefPVecI;
typedef RaftDisdefPVec::const_iterator RaftDisdefPVecCI;

typedef DiffBoundarydef* DiffBoundaryDefP;
typedef std::vector<DiffBoundaryDefP> DiffBoundaryDefPVec;
typedef DiffBoundaryDefPVec::iterator DiffBoundaryDefPVecI;
typedef DiffBoundaryDefPVec::const_iterator DiffBoundaryDefPVecCI;

typedef SDiffBoundarydef* SDiffBoundaryDefP;
typedef std::vector<SDiffBoundaryDefP> SDiffBoundaryDefPVec;
typedef SDiffBoundaryDefPVec::iterator SDiffBoundaryDefPVecI;
typedef SDiffBoundaryDefPVec::const_iterator SDiffBoundaryDefPVecCI;

}  // namespace solver

}  // namespace steps
