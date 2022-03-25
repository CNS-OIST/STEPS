#include "sreac_unittest.hpp"

#include "geom/dist/distmesh.hpp"
#include "model/model.hpp"
#include "mpi/dist/test/simulation.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"

namespace steps {
namespace dist {

SReacUnitTest::SReacUnitTest(const ScenarioInput &t_input)
    : Scenario("SReacUnitTest", "Unit test for surface reactions", t_input) {}

std::unique_ptr<Statedef>
SReacUnitTest::createStatedef(const simulation_t &simulation) const {
  steps::model::Model model;
  SReacUnitTestSimdef simdef(model, simulation.getMesh());
  return std::move(simdef.getStatedef());
}

void SReacUnitTest::register_compartments(DistMesh &mesh) const {
  mesh.addComp("comp_i", model::compartment_label(1));
  mesh.addComp("comp_o", model::compartment_label(2));
}

void SReacUnitTest::fill_compartments(simulation_t &simulation) const {
  steps::model::Model model;
  SReacUnitTestSimdef simdef(model, simulation.getMesh());

  simulation.setCompCount(simdef.getCompartementCounts(input.num_mols_factor));

  for (const auto &initializer : simdef.getPatchCounts()) {
      simulation.setPatchCount(initializer.patch,
                               initializer.species,
                               initializer.num_mols * input.num_mols_factor);
  }
}

void SReacUnitTest::run_simulation_impl(simulation_t &simulation) {
  simulation.run(input.end_time);
}

int SReacUnitTest::check_and_log_results_impl(simulation_t &simulation) const {
  osh::Real A = simulation.getCompCount("comp_i", "A");
  osh::Real B = simulation.getCompCount("comp_i", "B");
  osh::Real C = simulation.getCompCount("comp_o", "C");
  osh::Real D = simulation.getPatchCount("patch", "D");
  simulation.log_once("A: " + std::to_string(A) + ", B: " + std::to_string(B) +
                      ", C: " + std::to_string(C) +
                      ", D: " + std::to_string(D));
  bool status = EXIT_SUCCESS;
  if (static_cast<osh::I64>(A) != 0) {
    simulation.log_once("Error: expected concentration of A to be 0 but got " +
                        std::to_string(A));
    status = EXIT_FAILURE;
  }
  if (static_cast<osh::I64>(B) != 0) {
    simulation.log_once("Error: expected concentration of B to be 0 but got " +
                        std::to_string(B));
    status = EXIT_FAILURE;
  }
  if (static_cast<osh::I64>(C) != 0) {
    simulation.log_once("Error: expected concentration of C to be 0 but got " +
                        std::to_string(C));
    status = EXIT_FAILURE;
  }
  if (static_cast<osh::I64>(D) != 3000) {
    simulation.log_once("Error: expected 3000 mols D but got " +
                        std::to_string(D));
    status = EXIT_FAILURE;
  }
  return status;
}

} // namespace dist
} // namespace steps
