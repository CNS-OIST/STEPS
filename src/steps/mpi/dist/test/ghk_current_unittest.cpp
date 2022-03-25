#include "ghk_current_unittest.hpp"

#include <iostream>

#include "geom/dist/distmesh.hpp"
#include "math/constants.hpp"
#include "model/model.hpp"
#include "mpi/dist/test/simulation.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "mpi/dist/tetopsplit/kproc/surface_reactions.hpp"

namespace steps {
namespace dist {

GHKCurrentUnitTest::GHKCurrentUnitTest(const ScenarioInput &t_input)
    : Scenario("GHKCurrentUnitTest", "Unit test for GHK current", t_input) {}

std::unique_ptr<Statedef>
GHKCurrentUnitTest::createStatedef(const simulation_t &simulation) const {
  steps::model::Model model;
  GHKCurrentUnitTestSimdef simdef(model, simulation.getMesh());
  return std::move(simdef.getStatedef());
}

void GHKCurrentUnitTest::register_compartments(DistMesh &mesh) const {
  mesh.addComp("comp_i", model::compartment_label(1));
  mesh.addComp("comp_o", model::compartment_label(2));
}

void GHKCurrentUnitTest::fill_compartments(simulation_t &simulation) const {
  steps::model::Model model;
  GHKCurrentUnitTestSimdef simdef(model, simulation.getMesh());
  simulation.setCompCount(simdef.getCompartementCounts(input.num_mols_factor));

  for (const auto &initializer : simdef.getPatchCounts()) {
      simulation.setPatchCount(initializer.patch,
                               initializer.species,
                               initializer.num_mols * input.num_mols_factor);
  }
}

void GHKCurrentUnitTest::run_simulation_impl(simulation_t &simulation) {
  const osh::Real step{1.0e-8};
  const size_t ns{1000};
  res.resize(ns);
  conc_i = simulation.getCompConc("comp_i", "Na") * 1.0e3; //[mol/m^3]
  conc_o = simulation.getCompConc("comp_o", "Na") * 1.0e3; //[mol/m^3]
  for (size_t k = 0; k < ns; ++k) {
    simulation.reset();
    simulation.setPotential(equilibrium_potential);
    fill_compartments(simulation);
    simulation.run(step);
    res[k] = simulation.getGHKCurrents()[0] + simulation.getGHKCurrents()[1];
    simulation.log_once("test n: " + std::to_string(k));
  }
  const size_t n_steps = 10000;
  simulation.log_once("Last cumulative test start. N steps: " +
                      std::to_string(n_steps));

  simulation.run(step * n_steps);
  simulation.log_once("Last cumulative test end.");
}

int GHKCurrentUnitTest::check_and_log_results_impl(
    simulation_t &simulation) const {
  auto AI = simulation.getCompConc("comp_i", "Na") * 1.0e3; //[mol/m^3]
  auto AO = simulation.getCompConc("comp_o", "Na") * 1.0e3; //[mol/m^3]
  int status = 0;
  auto statedef = createStatedef(simulation);
  auto ghk_info = statedef->getPatchdef(model::patch_id("patch")).ghkSReacdefs().front()->getInfo();
  double concentration_ratio = std::exp(-equilibrium_potential * math::FARADAY *
                                        static_cast<double>(ghk_info.valence) / math::GAS_CONSTANT /
                                        statedef->getTemp());
  std::clog << "Checking the concentration ratio inside/outside the membrane "
            << AI / AO
            << " with respect to Nernst equation : " << concentration_ratio
            << "\n";
  if (std::abs(AI / AO - concentration_ratio) > 0.01 * concentration_ratio) {
    // Error out if simulated concentration ratio is more than 1% off analytic
    // value
    status = 1;
  }
  double c = std::accumulate(res.begin(), res.end(), 0.0) / res.size();
  // init_current is the value of the current as computed from the GHK formula
  double init_current = ghk_info.permeability * ghk_info.valence * math::FARADAY *
                          (-std::log(concentration_ratio)) /
                          (1.0 - concentration_ratio) *
                          (conc_i - 
				concentration_ratio * conc_o);
  std::clog << "GHK theoretical current  " << init_current << " vs average simulated current: " << c;  
    // Check that simulated current does not differ from  init_current by more
    // than 1%.
    if (std::abs(init_current - c) >
        0.01 * std::abs(init_current)) {
      std::clog << " : FAIL\n";
      status = 1;
    } else {
      std::clog << " : PASS\n";
    }
  return status;
}

} // namespace dist
} // namespace steps
