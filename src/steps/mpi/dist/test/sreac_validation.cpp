#include "sreac_validation.hpp"

#include <map>
#include <vector>

#include "geom/dist/distmesh.hpp"
#include "model/model.hpp"
#include "mpi/dist/test/simulation.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "util/vec3.hpp"

namespace steps {
namespace dist {

static inline bool isTolerable(osh::Real val_first, osh::Real val_second,
                               osh::Real tolerance) {
  return std::abs(2 * (val_first - val_second) / (val_first + val_second)) <=
         tolerance;
}

SReacValidationTest::SReacValidationTest(const ScenarioInput &t_input)
    : Scenario("sreac validation", "parallel surface reaction validations",
               t_input),
      tpnts(util::arange<osh::Real>(sampling_dt(), sim_end_time(),
                                    sampling_dt())),
      res_m_so2d(
          static_cast<size_t>(SurfaceReactionsValidationSimdef::NITER_so2d()),
          util::Vec2<osh::Real>(tpnts.size(), util::Vec1<osh::Real>(3, 0.0)))

{}

std::unique_ptr<Statedef>
SReacValidationTest::createStatedef(const simulation_t &simulation) const {
  steps::model::Model model;
  SurfaceReactionsValidationSimdef simdef(model, simulation.getMesh());
  return std::move(simdef.getStatedef());
}

void SReacValidationTest::register_compartments(DistMesh &mesh) const {
  mesh.addComp("comp1", model::compartment_label(1));
}

void SReacValidationTest::fill_compartments(
    simulation_t & /* simulation */) const {}

void SReacValidationTest::run_simulation_impl(simulation_t &simulation) {
  steps::model::Model model;
  const SurfaceReactionsValidationSimdef simdef(model, simulation.getMesh());

  const auto num_tpnts = tpnts.size();

  for (auto iter = 0u; iter < SurfaceReactionsValidationSimdef::NITER_so2d();
       iter++) {
    simulation.reset();
    simulation.setPatchCount(simdef.getso2d_PatchCount());

    for (auto t = 0u; t < num_tpnts; t++) {
      simulation.run(tpnts[t]);

      res_m_so2d[iter][t][0] = simulation.getPatchCount("patch1", "A_so2d");
      res_m_so2d[iter][t][1] = simulation.getPatchCount("patch1", "B_so2d");
    }
    simulation.log_progress(iter,
                            SurfaceReactionsValidationSimdef::NITER_so2d());
  }
}

int SReacValidationTest::check_and_log_results_impl(
    simulation_t &simulation) const {
  // generate mean result
  util::Vec2<osh::Real> mean_res_so2d;
  util::iterMeanVec3(res_m_so2d, mean_res_so2d);

  bool test_passed = true;
  int status{};
  osh::Real C_so2d = SurfaceReactionsValidationSimdef::COUNTA_so2d() -
                     SurfaceReactionsValidationSimdef::COUNTB_so2d();
  osh::Real patch_area = simulation.getMesh().total_measure(model::patch_id("patch1"));
  osh::Real CCST_so2d = SurfaceReactionsValidationSimdef::KCST_so2d() /
                        (6.02214179e23 * patch_area);
  for (size_t i = 0; i < tpnts.size(); i++) {
    osh::Real A_so2d = mean_res_so2d[i][0];
    osh::Real B_so2d = mean_res_so2d[i][1];
    osh::Real lnBA_so2d = std::log(B_so2d / A_so2d);
    osh::Real lineAB_so2d =
        std::log(SurfaceReactionsValidationSimdef::COUNTB_so2d() /
                 SurfaceReactionsValidationSimdef::COUNTA_so2d()) -
        C_so2d * CCST_so2d * tpnts[i];

    if (!isTolerable(lnBA_so2d, lineAB_so2d,
                     SurfaceReactionsValidationSimdef::tolerance_so2d())) {
      simulation.log_once(
          "Validation Failed: Reaction - Second-order, irreversible, 2D");
      test_passed = false;
      status += 1;
      break;
    }
  }
  if (test_passed) {
    simulation.log_once(
        "Validation passed: Reaction - Second-order, irreversible, 2D");
  }
  return status;
}

} // namespace dist
} // namespace steps
