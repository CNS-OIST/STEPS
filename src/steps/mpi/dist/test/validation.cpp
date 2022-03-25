#include "validation.hpp"

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

/// Sampling time-step
static const osh::Real DT = 0.1;
/// Simulation  end time
static const osh::Real INT = 1.1;
static const size_t NITER_max = 100000;

static inline bool isTolerable(osh::Real val_first, osh::Real val_second,
                               osh::Real tolerance) {
  return std::abs(2 * (val_first - val_second) / (val_first + val_second)) <=
         tolerance;
}

ValidationTest::ValidationTest(const ScenarioInput &t_input)
    : Scenario("validation", "parallel validations", t_input),
      tpnts(util::arange<osh::Real>(0.0, INT, DT)),
      res_m_foi(
          static_cast<size_t>(ValidationSimdef::NITER_foi()),
          util::Vec2<osh::Real>(tpnts.size(), util::Vec1<osh::Real>(1, 0.0))),
      res_m_for(
          static_cast<size_t>(ValidationSimdef::NITER_for()),
          util::Vec2<osh::Real>(tpnts.size(), util::Vec1<osh::Real>(2, 0.0)))

      ,
      res_m_soA2(
          static_cast<size_t>(ValidationSimdef::NITER_soA2()),
          util::Vec2<osh::Real>(tpnts.size(), util::Vec1<osh::Real>(2, 0.0)))

      ,
      res_m_soAA(
          static_cast<size_t>(ValidationSimdef::NITER_soAA()),
          util::Vec2<osh::Real>(tpnts.size(), util::Vec1<osh::Real>(3, 0.0)))

      ,
      res_m_soAB(
          static_cast<size_t>(ValidationSimdef::NITER_soAB()),
          util::Vec2<osh::Real>(tpnts.size(), util::Vec1<osh::Real>(3, 0.0)))

      ,
      res_m_toA3(
          static_cast<size_t>(ValidationSimdef::NITER_toA3()),
          util::Vec2<osh::Real>(tpnts.size(), util::Vec1<osh::Real>(2, 0.0))),
      res_m_toA2B(
          static_cast<size_t>(ValidationSimdef::NITER_toA2B()),
          util::Vec2<osh::Real>(tpnts.size(), util::Vec1<osh::Real>(3, 0.0)))
      // SReac validations
      ,
      res_m_so2d(
          static_cast<size_t>(ValidationSimdef::NITER_so2d()),
          util::Vec2<osh::Real>(tpnts.size(), util::Vec1<osh::Real>(3, 0.0)))

{}

std::unique_ptr<Statedef>
ValidationTest::createStatedef(const simulation_t &simulation) const {
  steps::model::Model model;
  ValidationSimdef simdef(model, simulation.getMesh());
  return std::move(simdef.getStatedef());
}

void ValidationTest::register_compartments(DistMesh &mesh) const {
  mesh.addComp("comp1", model::compartment_label(1));
}

void ValidationTest::fill_compartments(simulation_t & /* simulation */) const {}

void ValidationTest::run_simulation_impl(simulation_t &simulation) {
  steps::model::Model model;
  const ValidationSimdef simdef(model, simulation.getMesh());

  const auto num_tpnts = tpnts.size();

  for (auto iter = 0u; iter < NITER_max; iter++) {
    simulation.reset();

    simulation.setCompCount(
        simdef.getfoi_CompartementCounts(!(iter < ValidationSimdef::NITER_foi())));
    simulation.setCompCount(
        simdef.getfor_CompartementCounts(!(iter < ValidationSimdef::NITER_for())));
    simulation.setCompConc(
        simdef.getsoA2_CompartementConc(!(iter < ValidationSimdef::NITER_soA2())));
    simulation.setCompConc(
        simdef.getsoAA_CompartementConc(!(iter < ValidationSimdef::NITER_soAA())));
    simulation.setCompConc(
        simdef.getsoAB_CompartementConc(!(iter < ValidationSimdef::NITER_soAB())));
    simulation.setCompConc(
        simdef.gettoA3_CompartementConc(!(iter < ValidationSimdef::NITER_toA3())));
    simulation.setCompConc(
        simdef.gettoA2B_CompartementConc(!(iter < ValidationSimdef::NITER_toA2B())));

    for (auto t = 0u; t < num_tpnts; t++) {
      simulation.run(tpnts[t]);

      if (iter < ValidationSimdef::NITER_foi()) {
        res_m_foi[iter][t][0] = simulation.getCompCount("comp1", "A_foi");
      }

      if (iter < ValidationSimdef::NITER_for()) {
        res_m_for[iter][t][0] = simulation.getCompConc("comp1", "A_for") * 1e6;
        res_m_for[iter][t][1] = simulation.getCompConc("comp1", "B_for") * 1e6;
      }
      if (iter < ValidationSimdef::NITER_soA2()) {
        res_m_soA2[iter][t][0] = simulation.getCompConc("comp1", "A_soA2");
      }
      if (iter < ValidationSimdef::NITER_soAA()) {
        res_m_soAA[iter][t][0] = simulation.getCompConc("comp1", "A_soAA");
        res_m_soAA[iter][t][1] = simulation.getCompConc("comp1", "B_soAA");
      }
      if (iter < ValidationSimdef::NITER_soAB()) {
        res_m_soAB[iter][t][0] = simulation.getCompConc("comp1", "A_soAB");
        res_m_soAB[iter][t][1] = simulation.getCompConc("comp1", "B_soAB");
      }
      if (iter < ValidationSimdef::NITER_toA3()) {
        res_m_toA3[iter][t][0] = simulation.getCompConc("comp1", "A_toA3");
      }
      if (iter < ValidationSimdef::NITER_toA2B()) {
        res_m_toA2B[iter][t][0] = simulation.getCompConc("comp1", "A_toA2B");
        res_m_toA2B[iter][t][1] = simulation.getCompConc("comp1", "B_toA2B");
        res_m_toA2B[iter][t][2] = simulation.getCompConc("comp1", "C_toA2B");
      }
    }
    simulation.log_progress(iter, NITER_max);
  }
}

int ValidationTest::check_and_log_results_impl(simulation_t &simulation) const {
  // generate mean result
  util::Vec2<osh::Real> mean_res_foi;
  util::iterMeanVec3(res_m_foi, mean_res_foi);
  util::Vec2<osh::Real> std_res_foi;
  util::iterSTDVec3(res_m_foi, std_res_foi);

  util::Vec2<osh::Real> mean_res_for;
  util::iterMeanVec3(res_m_for, mean_res_for);

  util::Vec2<osh::Real> mean_res_soA2;
  util::iterMeanVec3(res_m_soA2, mean_res_soA2);

  util::Vec2<osh::Real> mean_res_soAA;
  util::iterMeanVec3(res_m_soAA, mean_res_soAA);

  util::Vec2<osh::Real> mean_res_soAB;
  util::iterMeanVec3(res_m_soAB, mean_res_soAB);

  util::Vec2<osh::Real> mean_res_toA3;
  util::iterMeanVec3(res_m_toA3, mean_res_toA3);

  util::Vec2<osh::Real> mean_res_toA2B;
  util::iterMeanVec3(res_m_toA2B, mean_res_toA2B);

  util::Vec2<osh::Real> mean_res_so2d;
  util::iterMeanVec3(res_m_so2d, mean_res_so2d);

  bool test_passed = true;
  int status{};
  for (size_t i = 1; i < tpnts.size(); i++) {
    osh::Real analy = ValidationSimdef::N_foi() *
                      std::exp(-ValidationSimdef::KCST_foi() * tpnts[i]);
    osh::Real std =
        std::pow((ValidationSimdef::N_foi() *
                  (std::exp(-ValidationSimdef::KCST_foi() * tpnts[i])) *
                  (1 - (std::exp(-ValidationSimdef::KCST_foi() * tpnts[i])))),
                 0.5);
    if (!(isTolerable(analy, mean_res_foi[i][0],
                      ValidationSimdef::tolerance_foi()) &&
          isTolerable(std, std_res_foi[i][0],
                      ValidationSimdef::tolerance_foi()))) {
      simulation.log_once(
          "Validation Failed: Reaction - First order, irreversible");
      test_passed = false;
      status += 1;
      break;
    }
  }
  if (test_passed) {
    simulation.log_once(
        "Validation passed: Reaction - First order, irreversible");
  }

  const osh::Real VOL =
      simulation.getMesh().total_measure(model::compartment_id("comp1"));

  test_passed = true;
  osh::Real Aeq =
      ValidationSimdef::COUNT_for() *
      (ValidationSimdef::KCST_b_for() / ValidationSimdef::KCST_f_for()) /
      (1 + (ValidationSimdef::KCST_b_for() / ValidationSimdef::KCST_f_for())) /
      (VOL * 6.0221415e26) * 1e6;

  osh::Real Beq =
      (ValidationSimdef::COUNT_for() / (VOL * 6.0221415e26)) * 1e6 - Aeq;

  for (size_t i = 7; i < tpnts.size(); i++) {
    if (!(isTolerable(mean_res_for[i][0], Aeq,
                      ValidationSimdef::tolerance_for()) &&
          isTolerable(mean_res_for[i][1], Beq,
                      ValidationSimdef::tolerance_for()))) {
      simulation.log_once(
          "Validation Failed: Reaction - First order, reversible");
      test_passed = false;
      status += 1;
      break;
    }
  }
  if (test_passed) {
    simulation.log_once(
        "Validation passed: Reaction - First order, reversible");
  }

  test_passed = true;
  for (size_t i = 0; i < tpnts.size(); i++) {
    osh::Real invA = (1.0 / mean_res_soA2[i][0]);
    osh::Real lineA = (1.0 / ValidationSimdef::CONCA_soA2() +
                       ((tpnts[i] * 2 * ValidationSimdef::KCST_soA2())));
    if (!isTolerable(invA, lineA, ValidationSimdef::tolerance_soA2())) {
      simulation.log_once(
          "Validation Failed: Reaction - Second order, irreversible, 2A->C");
      test_passed = false;
      status += 1;
      break;
    }
  }
  if (test_passed) {
    simulation.log_once(
        "Validation passed: Reaction - Second order, irreversible, 2A->C");
  }

  test_passed = true;
  for (size_t i = 0; i < tpnts.size(); i++) {
    osh::Real inv2A = (1.0 / std::pow(mean_res_toA3[i][0], 2));
    osh::Real line2A = (1.0 / std::pow(ValidationSimdef::CONCA_toA3(), 2) +
                        ((tpnts[i] * 6 * ValidationSimdef::KCST_toA3())));
    if (!isTolerable(inv2A, line2A, ValidationSimdef::tolerance_toA3())) {
      simulation.log_once(
          "Validation Failed: Reaction - Third order, irreversible, 3A->C");
      test_passed = false;
      status += 1;
      break;
    }
  }
  if (test_passed) {
    simulation.log_once(
        "Validation passed: Reaction - Third order, irreversible, 3A->C");
  }

  test_passed = true;
  osh::Real A0 = ValidationSimdef::CONCA_toA2B();
  osh::Real B0 = ValidationSimdef::CONCB_toA2B();
  osh::Real delta_AB = A0 - 2 * B0;
  osh::Real delta_BA = 2 * B0 - A0;

  for (size_t i = 1; i < tpnts.size(); i++) {
    osh::Real A = mean_res_toA2B[i][0];
    osh::Real B = mean_res_toA2B[i][1];
    osh::Real lineA =
        (-1.0 / delta_AB) * ((-1.0 / delta_BA) * std::log((B / A) / (B0 / A0)) +
                             1.0 / A - 1.0 / A0);
    osh::Real kt = (tpnts[i] * ValidationSimdef::KCST_toA2B());
    if (!isTolerable(kt, lineA, ValidationSimdef::tolerance_toA2B())) {
      simulation.log_once(
          "Validation Failed: Reaction - Third order, irreversible, 2A+B->C");
      test_passed = false;
      status += 1;
      break;
    }
  }
  if (test_passed) {
    simulation.log_once(
        "Validation passed: Reaction - Third order, irreversible, 2A+B->C");
  }

  test_passed = true;
  for (size_t i = 0; i < tpnts.size(); i++) {
    osh::Real invA = (1.0 / mean_res_soAA[i][0]);
    osh::Real invB = (1.0 / mean_res_soAA[i][1]);
    osh::Real lineA = (1.0 / ValidationSimdef::CONCA_soAA() +
                       ((tpnts[i] * ValidationSimdef::KCST_soAA())));
    osh::Real lineB = (1.0 / ValidationSimdef::CONCB_soAA() +
                       ((tpnts[i] * ValidationSimdef::KCST_soAA())));

    if (!(isTolerable(invA, lineA, ValidationSimdef::tolerance_soAA()) &&
          isTolerable(invB, lineB, ValidationSimdef::tolerance_soAA()))) {
      simulation.log_once(
          "Validation Failed: Reaction - Second order, irreversible, A0=B0");
      test_passed = false;
      status += 1;
      break;
    }
  }
  if (test_passed) {
    simulation.log_once(
        "Validation passed: Reaction - Second order, irreversible, A0=B0");
  }

  test_passed = true;
  osh::Real C_soAB =
      ValidationSimdef::CONCA_soAB() - ValidationSimdef::CONCB_soAB();

  for (size_t i = 0; i < tpnts.size(); i++) {
    osh::Real A_soAB = mean_res_soAB[i][0];
    osh::Real B_soAB = mean_res_soAB[i][1];
    osh::Real lnBA_soAB = std::log(B_soAB / A_soAB);
    osh::Real lineAB_soAB = std::log(ValidationSimdef::CONCB_soAB() /
                                     ValidationSimdef::CONCA_soAB()) -
                            C_soAB * ValidationSimdef::KCST_soAB() * tpnts[i];

    if (!isTolerable(lnBA_soAB, lineAB_soAB,
                     ValidationSimdef::tolerance_soAB())) {
      simulation.log_once(
          "Validation Failed: Reaction - Second order, irreversible, A0!=B0");
      test_passed = false;
      status += 1;
      break;
    }
  }
  if (test_passed) {
    simulation.log_once(
        "Validation passed: Reaction - Second order, irreversible, A0!=B0");
  }

  return status;
}

} // namespace dist
} // namespace steps
