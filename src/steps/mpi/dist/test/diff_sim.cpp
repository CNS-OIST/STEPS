#include "diff_sim.hpp"

#include <map>

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

DiffusionOnly::DiffusionOnly(const ScenarioInput &t_input)
    : Scenario("diffusionOnly", "Only diffuse one species", t_input) {}

std::unique_ptr<Statedef>
DiffusionOnly::createStatedef(const simulation_t &simulation) const {
  steps::model::Model model;
  auto statedef = std::make_unique<Statedef>(model, simulation.getMesh());
  statedef->addComp("comp");
  std::vector<model::species_name> species = {"A"};
  statedef->addCompSpecs("comp", species);

  statedef->addCompDiff("comp", "A", dcst);

  return statedef;
}

void DiffusionOnly::register_compartments(DistMesh &mesh) const {
  mesh.addComp("comp", comp_tag);
}

void DiffusionOnly::fill_compartments(simulation_t &simulation) const {
  if (input.do_inject_in_elt == -1) {
      simulation.setCompCount("comp", "A", NMols);
  } else {
    simulation.setOwnedElementCount("comp", input.do_inject_in_elt, "A", NMols);
  }
}

void DiffusionOnly::run_simulation_impl(simulation_t &simulation) {
  int iter = 0;

  for (osh::Real t = 0.0; t < input.end_time;
       t += input.do_interval) { // NOLINT
    simulation.log_once("time_step=" + std::to_string(t));
    simulation.exportMolStateToVTK("count_" + std::to_string(iter) + ".vtk");
    simulation.run(t);
    iter++;
  }
}

} // namespace dist
} // namespace steps
