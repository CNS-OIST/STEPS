#include "simple_model.hpp"

#include <map>

#include "geom/dist/distcomp.hpp"
#include "geom/dist/distmesh.hpp"
#include "model/diff.hpp"
#include "model/model.hpp"
#include "model/reac.hpp"
#include "model/spec.hpp"
#include "model/volsys.hpp"
#include "mpi/dist/test/simulation.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "util/profile/profiler_interface.h"

namespace steps {
namespace dist {

SimpleModel2::SimpleModel2(const ScenarioInput &t_input)
    : Scenario("SimpleModel", "Simple container in Chen_FNeuroinf_2017",
               t_input) {}

std::unique_ptr<Statedef>
SimpleModel2::createStatedef(const simulation_t &simulation) const {
  Instrumentor::phase p("SimpleModel2::createStatedef");
  steps::model::Model model;
  SimpleSimdef simdef(model, simulation.getMesh());
  return std::move(simdef.getStatedef());
}

void SimpleModel2::register_compartments(DistMesh &mesh) const {
  Instrumentor::phase p("SimpleModel2::register_compartments");
  mesh.addComp("comp1", input.comp_label);
}

void SimpleModel2::fill_compartments(simulation_t &simulation) const {
  Instrumentor::phase p("SimpleModel2::fill_compartments");
  steps::model::Model model;
  SimpleSimdef simdef(model, simulation.getMesh());
  simulation.setCompCount(simdef.getCompartementCounts(input.num_mols_factor));
}

void SimpleModel2::run_simulation_impl(simulation_t &simulation) {
  Instrumentor::phase p("SimpleModel2::run_simulation_impl");
  simulation.run(input.end_time);
}

int SimpleModel2::check_and_log_results_impl(simulation_t &simulation) const {
  Instrumentor::phase p("SimpleModel2::check_and_log_results_impl");
  osh::Real A = simulation.getCompCount("comp1", "A");
  osh::Real B = simulation.getCompCount("comp1", "B");
  osh::Real C = simulation.getCompCount("comp1", "C");
  osh::Real D = simulation.getCompCount("comp1", "D");
  osh::Real E = simulation.getCompCount("comp1", "E");
  osh::Real F = simulation.getCompCount("comp1", "F");
  osh::Real G = simulation.getCompCount("comp1", "G");
  osh::Real H = simulation.getCompCount("comp1", "H");
  osh::Real I = simulation.getCompCount("comp1", "I");
  osh::Real J = simulation.getCompCount("comp1", "J");
  simulation.log_once(
      "A: " + std::to_string(A) + ", B: " + std::to_string(B) +
      ", C: " + std::to_string(C) + ", D: " + std::to_string(D) +
      ", E: " + std::to_string(E) + ", F: " + std::to_string(F) +
      ", G: " + std::to_string(G) + ", H: " + std::to_string(H) +
      ", I: " + std::to_string(I) + ", J: " + std::to_string(J));
  return EXIT_SUCCESS;
}


SimpleModel3::SimpleModel3(const ScenarioInput &t_input)
    : Scenario("SimpleModel", "Simple container in Chen_FNeuroinf_2017",
               t_input) {}

std::unique_ptr<Statedef>
SimpleModel3::createStatedef(const simulation_t &simulation) const {
  steps::model::Model model;
  
  // Build STEPS model object
  auto *A = new steps::model::Spec("A", &model);
  auto *B = new steps::model::Spec("B", &model);
  auto *C = new steps::model::Spec("C", &model);
  auto *D = new steps::model::Spec("D", &model);
  auto *E = new steps::model::Spec("E", &model);
  auto *F = new steps::model::Spec("F", &model);
  auto *G = new steps::model::Spec("G", &model);
  auto *H = new steps::model::Spec("H", &model);
  auto *I = new steps::model::Spec("I", &model);
  auto *J = new steps::model::Spec("J", &model);

  auto *vsys = new steps::model::Volsys("vsys", &model);

  new steps::model::Reac("R1", vsys, {A, B}, {C}, 1000.0e6);
  new steps::model::Reac("R2", vsys, {C}, {A, B}, 100);
  new steps::model::Reac("R3", vsys, {C, D}, {E}, 100e6);
  new steps::model::Reac("R4", vsys, {E}, {C, D}, 10);
  new steps::model::Reac("R5", vsys, {F, G}, {H}, 10e6);
  new steps::model::Reac("R6", vsys, {H}, {F, G}, 1);
  new steps::model::Reac("R7", vsys, {H, I}, {J}, 1e6);
  new steps::model::Reac("R8", vsys, {J}, {H, I}, 0.1 * 10);

  new steps::model::Diff("D1", vsys, A, 100e-12);
  new steps::model::Diff("D2", vsys, B, 90e-12);
  new steps::model::Diff("D3", vsys, C, 80e-12);
  new steps::model::Diff("D4", vsys, D, 70e-12);
  new steps::model::Diff("D5", vsys, E, 60e-12);
  new steps::model::Diff("D6", vsys, F, 50e-12);
  new steps::model::Diff("D7", vsys, G, 40e-12);
  new steps::model::Diff("D8", vsys, H, 30e-12);
  new steps::model::Diff("D9", vsys, I, 20e-12);
  new steps::model::Diff("D10", vsys, J, 10e-12);

  // Use a SimDef, not a SimpleSimDef
  Simdef simdef(model, simulation.getMesh());
  return std::move(simdef.getStatedef());
}

void SimpleModel3::register_compartments(DistMesh &mesh) const {
  auto *comp1 = new steps::dist::DistComp("comp1", mesh);
  comp1->addVolsys("vsys");

  mesh.addComp("comp1", input.comp_label);
}

void SimpleModel3::fill_compartments(simulation_t &simulation) const {
  const osh::Real N0A = 1000.0 * input.num_mols_factor;
  const osh::Real N0B = 2000.0 * input.num_mols_factor;
  const osh::Real N0C = 3000.0 * input.num_mols_factor;
  const osh::Real N0D = 4000.0 * input.num_mols_factor;
  const osh::Real N0E = 5000.0 * input.num_mols_factor;
  const osh::Real N0F = 6000.0 * input.num_mols_factor;
  const osh::Real N0G = 7000.0 * input.num_mols_factor;
  const osh::Real N0H = 8000.0 * input.num_mols_factor;
  const osh::Real N0I = 9000.0 * input.num_mols_factor;
  const osh::Real N0J = 10000.0 * input.num_mols_factor;

  Simdef::compartment_counts_t compartmentCounts{{"comp1",
                        {
                            {"A", N0A},
                            {"B", N0B},
                            {"C", N0C},
                            {"D", N0D},
                            {"E", N0E},
                            {"F", N0F},
                            {"G", N0G},
                            {"H", N0H},
                            {"I", N0I},
                            {"J", N0J},
                        }}};
  simulation.setCompCount(compartmentCounts);
}

void SimpleModel3::run_simulation_impl(simulation_t &simulation) {
  simulation.run(input.end_time);
}

int SimpleModel3::check_and_log_results_impl(simulation_t &simulation) const {
  osh::Real A = simulation.getCompCount("comp1", "A");
  osh::Real B = simulation.getCompCount("comp1", "B");
  osh::Real C = simulation.getCompCount("comp1", "C");
  osh::Real D = simulation.getCompCount("comp1", "D");
  osh::Real E = simulation.getCompCount("comp1", "E");
  osh::Real F = simulation.getCompCount("comp1", "F");
  osh::Real G = simulation.getCompCount("comp1", "G");
  osh::Real H = simulation.getCompCount("comp1", "H");
  osh::Real I = simulation.getCompCount("comp1", "I");
  osh::Real J = simulation.getCompCount("comp1", "J");
  simulation.log_once(
      "A: " + std::to_string(A) + ", B: " + std::to_string(B) +
      ", C: " + std::to_string(C) + ", D: " + std::to_string(D) +
      ", E: " + std::to_string(E) + ", F: " + std::to_string(F) +
      ", G: " + std::to_string(G) + ", H: " + std::to_string(H) +
      ", I: " + std::to_string(I) + ", J: " + std::to_string(J));
  return EXIT_SUCCESS;
}

} // namespace dist
} // namespace steps
