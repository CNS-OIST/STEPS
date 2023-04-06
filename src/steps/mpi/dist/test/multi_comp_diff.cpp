#include "multi_comp_diff.hpp"

#include <complex>
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

MultipleCompartmentDiff::MultipleCompartmentDiff(const ScenarioInput &t_input)
    : Scenario("MultipleCompartmentDiff", "Simulation occurs in 2 compartments",
               t_input) {}

std::unique_ptr<Statedef>
MultipleCompartmentDiff::createStatedef(const simulation_t &simulation) const {
  steps::model::Model model;
  MultiCompartmentDiffSimdef simdef(model, simulation.getMesh());
  return std::move(simdef.getStatedef());
}

void MultipleCompartmentDiff::register_compartments(DistMesh &mesh) const {
  mesh.addComp("comp_i", model::compartment_label(1));
  mesh.addComp("comp_o", model::compartment_label(2));
}

void MultipleCompartmentDiff::fill_compartments(
    simulation_t &simulation) const {
  steps::model::Model model;
  simulation.setCompCount(
      MultiCompartmentDiffSimdef(model, simulation.getMesh()).getCompartementCounts());
}

void MultipleCompartmentDiff::run_simulation_impl(simulation_t &simulation) {
  osh::Real interval = 0.00001;
  simulation.getMesh().addDiffusionBoundary("0", "comp_i", "comp_o",
                                            std::nullopt);
  C_i_init = simulation.getCompCount("comp_i", "C");
  C_o_init = simulation.getCompCount("comp_o", "C");
  simulation.run(interval);
  C_i_nodiff = simulation.getCompCount("comp_i", "C");
  C_o_nodiff = simulation.getCompCount("comp_o", "C");
  simulation.setDiffusionBoundaryActive("0", "C", true);
  simulation.run(100.0 * interval);
  C_i_bdiff = simulation.getCompCount("comp_i", "C");
  C_o_bdiff = simulation.getCompCount("comp_o", "C");
}

int MultipleCompartmentDiff::check_and_log_results_impl(
    simulation_t &simulation) const {
  int status = 0;
  simulation.log_once("Init: " + std::to_string(C_i_init) + " " +
                      std::to_string(C_o_init));
  simulation.log_once("No diff: " + std::to_string(C_i_nodiff) + " " +
                      std::to_string(C_o_nodiff));
  simulation.log_once("B. diff: " + std::to_string(C_i_bdiff) + " " +
                      std::to_string(C_o_bdiff));
  // Concentrations should be equal if no diffusion boundary is defined
  if (C_i_init != C_i_nodiff || C_o_init != C_o_init) {
    status = 1;
  }
  // Concentration inside should be half that outside
  if (std::abs(C_i_bdiff - (C_i_init + C_o_init) / 3.0) > 0.1 * C_i_bdiff ||
      std::abs(C_o_bdiff - (C_i_init + C_o_init) * 2.0 / 3.0) >
          0.1 * C_o_bdiff) {
    status = 1;
  }
  return status;
}

} // namespace dist
} // namespace steps
