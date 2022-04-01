#include "multi_comp.hpp"

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

MultipleCompartment::MultipleCompartment(const ScenarioInput &t_input)
    : Scenario("MultipleCompartment", "Simulation occurs in 2 compartments",
               t_input) {}

std::unique_ptr<Statedef>
MultipleCompartment::createStatedef(const simulation_t &simulation) const {
  steps::model::Model model;
  MultiCompartmentSimdef simdef(model, simulation.getMesh());
  return std::move(simdef.getStatedef());
}

void MultipleCompartment::register_compartments(DistMesh &mesh) const {
  mesh.addComp("Left", model::compartment_label(1));
  mesh.addComp("Right", model::compartment_label(2));
}

void MultipleCompartment::fill_compartments(simulation_t &simulation) const {
  steps::model::Model model;
  MultiCompartmentSimdef simdef(model, simulation.getMesh());
  simulation.setCompCount(simdef.getCompartementCounts());
}

void MultipleCompartment::log_num_molecules(
    const simulation_t &simulation) const {
  osh::Real A_Left = simulation.getCompCount("Left", "A");
  osh::Real B_Left = simulation.getCompCount("Left", "B");
  osh::Real C_Left = simulation.getCompCount("Left", "C");
  osh::Real C_Right = simulation.getCompCount("Right", "C");
  osh::Real D_Right = simulation.getCompCount("Right", "D");
  osh::Real E_Right = simulation.getCompCount("Right", "E");

  simulation.log_once("Comp Left has " + std::to_string(std::real(A_Left)) +
                      " of spec A");
  simulation.log_once("Comp Left has " + std::to_string(std::real(B_Left)) +
                      "of spec B");
  simulation.log_once("Comp Left has " + std::to_string(std::real(C_Left)) +
                      "of spec C");
  simulation.log_once("Comp Right has " + std::to_string(std::real(C_Right)) +
                      "of spec C");
  simulation.log_once("Comp Right has " + std::to_string(std::real(D_Right)) +
                      "of spec D");
  simulation.log_once("Comp Right has " + std::to_string(std::real(E_Right)) +
                      "of spec E");
}

void MultipleCompartment::run_simulation_impl(simulation_t &simulation) {
  osh::Real interval = 0.00001;
  simulation.log_once("Binomial Threshold: " + std::to_string(input.threshold));
  log_num_molecules(simulation);
  for (auto i = 1; i < 100; i++) {
    const std::string output_name = "conc_" + std::to_string(i) + ".vtk";
    simulation.exportMolStateToVTK(output_name);
    simulation.run(interval * i);
    log_num_molecules(simulation);
    simulation.log_once("Number of iterations: " +
                        std::to_string(simulation.getNIterations()));
  }
}

} // namespace dist
} // namespace steps
