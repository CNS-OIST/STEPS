#include <map>
#include <vector>

#include <hadoken/format/format.hpp>

#include "opsplit/compdef.hpp"
#include "opsplit/diffdef.hpp"
#include "opsplit/mesh.hpp"
#include "opsplit/patchdef.hpp"
#include "opsplit/reacdef.hpp"
#include "opsplit/sreacdef.hpp"
#include "opsplit/test/simdef.hpp"
#include "opsplit/test/simulation.hpp"
#include "sreac_validation.hpp"

namespace zee {

using hadoken::scat;

/// Sampling time-step
static const PetscScalar DT = 0.01;
/// Simulation  end time
static const PetscScalar INT = 1.1;

static inline bool isTolerable(PetscScalar val_first,
                               PetscScalar val_second,
                               PetscScalar tolerance) {
    return std::abs(2 * (val_first - val_second) / (val_first + val_second)) <= tolerance;
}

SReacValidationTest::SReacValidationTest(const ScenarioInput& t_input)
    : Scenario("sreac validation", "parallel surface reaction validations", t_input)
    , tpnts(arange<PetscScalar>(0.0, INT, DT))
    , res_m_so2d(static_cast<size_t>(SurfaceReactionsValidationSimdef::NITER_so2d()),
                 Vec2<PetscScalar>(tpnts.size(), Vec1<PetscScalar>(3, 0.0)))

{}

std::unique_ptr<Statedef> SReacValidationTest::createStatedef() const {
    SurfaceReactionsValidationSimdef simdef;
    return std::move(simdef.getStatedef());
}

void SReacValidationTest::register_compartments(DistMesh& mesh) const {
    mesh.addComp("comp1", model::compartment_label(1));
}

void SReacValidationTest::fill_compartments(simulation_t& /* simulation */) const {}

void SReacValidationTest::run_simulation_impl(simulation_t& simulation) {
    const SurfaceReactionsValidationSimdef simdef;

    const auto num_tpnts = tpnts.size();

    for (auto iter = 0u; iter < SurfaceReactionsValidationSimdef::NITER_so2d(); iter++) {
        simulation.reset();
        simulation.setPatchCount(simdef.getso2d_PatchCount());

        for (auto t = 0u; t < num_tpnts; t++) {
            simulation.run(tpnts[t]);
            res_m_so2d[iter][t][0] = simulation.getPatchCount("patch1", "A_so2d");
            res_m_so2d[iter][t][1] = simulation.getPatchCount("patch1", "B_so2d");
        }
    }
}

int SReacValidationTest::check_and_log_results_impl(const simulation_t& simulation) const {
    // generate mean result
    Vec2<PetscScalar> mean_res_so2d;
    iterMeanVec3(res_m_so2d, mean_res_so2d);

    bool test_passed = true;
    int status{};
    PetscScalar C_so2d = SurfaceReactionsValidationSimdef::COUNTA_so2d() -
                         SurfaceReactionsValidationSimdef::COUNTB_so2d();
    /// TODO: Remove the const_cast and warm up DistMesh
    PetscScalar patch_area =
        const_cast<DistMesh*>(&simulation.getMesh())->getTotalPatchArea("patch1");
    PetscScalar CCST_so2d = SurfaceReactionsValidationSimdef::KCST_so2d() /
                            (6.02214179e23 * patch_area);
    for (size_t i = 0; i < tpnts.size(); i++) {
        PetscScalar A_so2d = mean_res_so2d[i][0];
        PetscScalar B_so2d = mean_res_so2d[i][1];
        PetscScalar lnBA_so2d = std::log(B_so2d / A_so2d);
        PetscScalar lineAB_so2d = std::log(SurfaceReactionsValidationSimdef::COUNTB_so2d() /
                                           SurfaceReactionsValidationSimdef::COUNTA_so2d()) -
                                  C_so2d * CCST_so2d * tpnts[i];

        if (!isTolerable(lnBA_so2d,
                         lineAB_so2d,
                         SurfaceReactionsValidationSimdef::tolerance_so2d())) {
            simulation.log_once("Validation Failed: Reaction - Second-order, irreversible, 2D");
            test_passed = false;
            status += 1;
            break;
        }
    }
    if (test_passed) {
        simulation.log_once("Validation passed: Reaction - Second-order, irreversible, 2D");
    }
    return status;
}

}  // namespace zee
