#include "geom/dist/distcomp.hpp"
#include "geom/dist/distmesh.hpp"
#include "model/model.hpp"
#include "model/spec.hpp"
#include "model/volsys.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "mpi/dist/tetopsplit/kproc/fwd.hpp"
#include "mpi/dist/tetopsplit/kproc/kproc_id.hpp"
#include "mpi/dist/tetopsplit/kproc/kproc_state.hpp"

#include "mpi/dist/tetopsplit/mol_state.hpp"
#include "test_common.hpp"
#include "util/vocabulary.hpp"

#include <catch2/catch_test_macros.hpp>

namespace osh = Omega_h;

auto lib = Omega_h::Library();

TEST_CASE("KProc_state_independent_group") {
    steps::model::Model mdl;
    steps::model::Spec SA{"SA", mdl};
    steps::model::Spec SB{"SB", mdl};
    steps::model::Spec SC{"SC", mdl};
    steps::model::Spec SD{"SD", mdl};

    steps::model::Volsys vsys{"vsys", mdl};

    // SA + SB -> SC
    std::vector<steps::model::Spec*> lhs_1, rhs_1;
    lhs_1.emplace_back(&SA);
    lhs_1.emplace_back(&SB);
    rhs_1.emplace_back(&SC);
    steps::model::Reac r1{"R1", vsys, lhs_1, rhs_1};

    // SC -> SA + SB
    std::vector<steps::model::Spec*> lhs_2, rhs_2;
    rhs_2.emplace_back(&SA);
    rhs_2.emplace_back(&SB);
    lhs_2.emplace_back(&SC);
    steps::model::Reac r2{"R2", vsys, lhs_2, rhs_2};

    // SD -> None
    std::vector<steps::model::Spec*> lhs_3, rhs_3;
    lhs_3.emplace_back(&SD);
    steps::model::Reac r3{"R3", vsys, lhs_3, rhs_3};

    const auto mesh_file = Omega_h::filesystem::path(STEPS_SOURCE_DIR) / "test" / "mesh" /
                           "3_tets.msh";
    steps::dist::DistMesh mesh(lib, mesh_file.string());
    steps::dist::DistComp comp{steps::dist::mesh::compartment_name("__MESH__"), mesh};
    comp.addVolsys("vsys");

    steps::dist::Statedef sd{mdl, mesh};

    osh::LOs spec_per_elem{4, 4, 4};
    osh::LOs cmplx_per_elem{0, 0, 0};
    steps::dist::RankInfo info{0, 1};
    steps::dist::MolState mol_state{mesh, sd, spec_per_elem, cmplx_per_elem, info};

    // Without independent kprocs
    {
        steps::dist::kproc::KProcState kps_1{sd, mesh, mol_state, false};

        REQUIRE(kps_1.groups().size() == 1);
    }

    // With independent kprocs
    {
        steps::dist::kproc::KProcState kps_2{sd, mesh, mol_state, true};

        REQUIRE(kps_2.groups().size() == 6);  // 3 tets, 2 groups per tet because R3 is independent
                                              // from R1 and R2

        // Check that each Kproc only appears in one group
        std::set<osh::LO> kproc_ids;
        for (const auto& group: kps_2.groups()) {
            for (const auto kid: group) {
                REQUIRE(kproc_ids.insert(kid).second);
            }
        }

        // Check that adding an outdated KProcID only adds it for one propensity group
        steps::dist::kproc::KProcID kpid_0{steps::dist::kproc::KProcType::Reac, 0};
        std::vector<steps::dist::kproc::KProcID> _outdated;
        _outdated.emplace_back(kpid_0);
        kps_2.add_outdated_kprocs(_outdated);
        for (uint i = 0; i < 6; ++i) {
            const auto& outdated = kps_2.get_outdated_kprocs(i);
            if (i == 0) {
                REQUIRE(outdated.size() == 1);
                REQUIRE(outdated[0].data() == kpid_0.data());
            } else {
                REQUIRE(outdated.empty());
            }
        }
    }
}
