#include "scenario.hpp"

#include <algorithm>

#include "geom/dist/distmesh.hpp"
#include "mpi/dist/test/ca_burst_background.hpp"
#include "mpi/dist/test/ca_burst_full.hpp"
#include "mpi/dist/test/ca_burst_integration_test.hpp"
#include "mpi/dist/test/diff_sim.hpp"
#include "mpi/dist/test/ghk_current_unittest.hpp"
#include "mpi/dist/test/multi_comp.hpp"
#include "mpi/dist/test/multi_comp_diff.hpp"
#include "mpi/dist/test/rallpack3.hpp"
#include "mpi/dist/test/simple_model.hpp"
#include "mpi/dist/test/simulation.hpp"
#include "mpi/dist/test/single_comp_dist.hpp"
#include "mpi/dist/test/sreac_unittest.hpp"
#include "mpi/dist/test/sreac_validation.hpp"
#include "mpi/dist/test/validation.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "util/tracker/peak_rss.hpp"

namespace steps {
namespace dist {

ScenarioResult ScenarioResult::invalid(int status) {
  return {status, 0, boost::none, 0, 0, 0, {}, 0, {}, 0, {}, {}};
}

template <typename RNG>
Scenario<RNG>::Scenario(const std::string &t_name,
                        const std::string &t_description,
                        const ScenarioInput &t_input)
    : name(t_name), description(t_description), input(t_input) {}

template <typename RNG> Scenario<RNG>::~Scenario() noexcept = default;

template <typename RNG>
ScenarioResult Scenario<RNG>::execute(simulation_t &simulation) {

  simulation.log_once("Initializing simulation");
  timer_initialization.start();
  initialize(simulation);
  timer_initialization.stop();
  const auto &result = run_simulation(simulation);
  return result;
}

template <typename RNG>
void Scenario<RNG>::initialize(simulation_t &simulation) {
  if (!initialized) {
    load_mesh(simulation);
    apply_statedef(simulation);
    initialized = true;
  }
}

template <typename RNG>
void Scenario<RNG>::load_mesh(Simulation<RNG> &simulation) {
  simulation.log_once("Loading mesh file \"" + input.mesh_file + "\"");
  this->register_compartments(simulation.getMesh());
  simulation.log_once("Mesh file loaded");
}

template <typename RNG>
void Scenario<RNG>::log_owned_elements(
    const Simulation<RNG> &simulation) const {
    simulation.log_all("Number of owned elements: " +
                       std::to_string(simulation.getMesh().owned_elems().size()));
    auto err = MPI_Barrier(simulation.comm());
    if (err != MPI_SUCCESS) {
        MPI_Abort(simulation.comm(), err);
  }
}

template <typename RNG>
void Scenario<RNG>::apply_statedef(Simulation<RNG> &simulation) {
  auto statedef = this->createStatedef(simulation);
  if (input.log_statedef_report) {
    simulation.log_once(statedef->createReport());
  }
  simulation.init(std::move(statedef));
  this->fill_compartments(simulation);
  this->fill_patches(simulation);
  simulation.setDiffOpBinomialThreshold(input.threshold);
  if (input.log_state_report) {
    simulation.log_once(simulation.createStateReport());
  }
}

template <typename RNG>
void ScenarioResult::log_summary(const Simulation<RNG> &simulation,
                                 const ScenarioInput &input) const {
  // simulation.log_once("hpcbench_meta_zee_revision=", Zee_COMMIT));
  simulation.log_once("hpcbench_meta_mesh_backend=" +
                      simulation.getMesh().getBackend());
  simulation.log_once("hpcbench_meta_mesh_file=" + input.mesh_file);
  simulation.log_once("hpcbench_meta_num_processes=" +
                      std::to_string(simulation.comm_size));
  simulation.log_once("hpcbench_meta_num_elements=" +
                      std::to_string(num_elements));
  simulation.log_once("hpcbench_meta_num_elements_per_process=" +
                      std::to_string(num_elements / simulation.comm_size));
  simulation.log_once("hpcbench_meta_threshold=" +
                      std::to_string(input.threshold));
  simulation.log_once("hpcbench_meta_scale=" + std::to_string(input.scale));
  simulation.log_once("hpcbench_meta_num_mols_factor=" +
                      std::to_string(input.num_mols_factor));
  simulation.log_once("hpcbench_meta_end_time_sec=" +
                      std::to_string(input.end_time));
  simulation.log_once("hpcbench_metric_num_iterations=" +
                      std::to_string(num_iterations));
  simulation.log_once("hpcbench_metric_time_step_sec=" +
                      std::to_string(time_step));
  if (efield_dt)
    simulation.log_once("hpcbench_metric_efield_dt_sec=" +
                        std::to_string(*efield_dt));
  simulation.log_once(
      "hpcbench_metric_elapsed_initialization_sec=" +
      std::to_string(elapsed_initialization_sec));
  simulation.log_once("hpcbench_metric_reaction_method=" +
                      std::to_string(simulation.ssaMethod()));
  simulation.log_once("hpcbench_metric_num_reactions=" +
                      std::to_string(num_reactions));
  simulation.log_once(
      "hpcbench_metric_elapsed_reactions_sec=" +
      std::to_string(elapsed_reactions));
  if (elapsed_efield) {
      simulation.log_once(
          "hpcbench_metric_elapsed_efield_sec=" +
          std::to_string(*elapsed_efield));
  }
  simulation.log_once("hpcbench_metric_num_diffusions=" +
                      std::to_string(num_diffusions));
  simulation.log_once(
      "hpcbench_metric_elapsed_diffusions_sec=" +
      std::to_string(elapsed_diffusions));
  simulation.log_once("hpcbench_metric_elapsed_sim_sec=" +
                      std::to_string(elapsed_sim_sec));
  simulation.log_once("hpcbench_metric_status=" + std::to_string(status));

}

template <typename RNG>
ScenarioResult Scenario<RNG>::run_simulation(Simulation<RNG> &simulation) {
  simulation.log_once("Starting simulation...");
  util::TimeTracker timer_sim;
  timer_sim.start();
  this->run_simulation_impl(simulation);
  timer_sim.stop();
  const auto status = this->check_and_log_results(simulation);

  // scenario will be removed, no need to spend extra time to pass this
  // temporary measurement around
  simulation.log_once("hpcbench_metric_sim_peak_rss_MB=" +
                      std::to_string(util::peak_rss()*1.e-6));

  boost::optional<osh::Real> efield_dt(false, 0.);
#if USE_PETSC
  efield_dt = simulation.getEfieldDt();
#endif // USE_PETSC

  return {status,
          simulation.getIterationTimeStep(),
          efield_dt,
          simulation.getNIterations(),
          timer_sim.diff(),
          simulation.getMesh().total_num_elems(),
          timer_initialization.diff(),
          simulation.getSSAOpExtent(),
          simulation.getElapsedSSA(),
          simulation.getDiffOpExtent(),
          simulation.getElapsedDiff(),
          simulation.getElapsedEField()};
}

template <typename RNG>
ScenarioResult test_splitting_operator(Simulation<RNG> &simulation,
                                       int scenario,
                                       const ScenarioInput &input) {
    switch (scenario) {
    case 0:
        return ValidationTest(input).execute(simulation);
    case 1:
        return SimpleModel2(input).execute(simulation);
    case 2:
        return MultipleCompartment(input).execute(simulation);
    case 3:
        return DiffusionOnly(input).execute(simulation);
    case 4:
        return SReacValidationTest(input).execute(simulation);
    case 5:
        return SReacUnitTest(input).execute(simulation);
    case 6:
        return CaBurstBackground(input).execute(simulation);
    case 7:
        return Rallpack3(input).execute(simulation);
    case 8:
        return GHKCurrentUnitTest(input).execute(simulation);
    case 9:
        return MultipleCompartmentDiff(input).execute(simulation);
    case 10:
        return CaBurstFull(input).execute(simulation);
    case 11:
        return SingleCompDist(input).execute(simulation);
    case 12:
        return SimpleModel3(input).execute(simulation);
    case 13:
        return CaBurstIntegrationTest(input).execute(simulation);
    default:
        simulation.log_once("Unknown test id: " + std::to_string(scenario));
        return ScenarioResult::invalid(-1);
    }
}

template <typename RNG>
int Scenario<RNG>::check_and_log_results(Simulation<RNG> &simulation) const {
  auto status = check_and_log_results_impl(simulation);
  const auto num_diffusions = simulation.getDiffOpExtent();
  const auto num_reactions = simulation.getSSAOpExtent();
  if (simulation.comm_rank == 0) {
    if (input.expected_num_diffusions >= 0) {
      const auto min_diffusions = static_cast<osh::I64>(
          static_cast<osh::Real>(input.expected_num_diffusions) *
          (1. - input.expected_num_diffusions_tolerance));
      const auto max_diffusions = static_cast<osh::I64>(
          static_cast<osh::Real>(input.expected_num_diffusions) *
          (1. + input.expected_num_diffusions_tolerance));
      if (num_diffusions < min_diffusions) {
        simulation.log_once(
            "Error: too few number of diffusions, expected at most " +
            std::to_string(min_diffusions));
        status = EXIT_FAILURE;
      } else if (num_diffusions > max_diffusions) {
        simulation.log_once(
            "Error: too many number of diffusions, expected at most " +
            std::to_string(max_diffusions));
        status = EXIT_FAILURE;
      }
    }
    if (input.expected_num_reactions >= 0) {
      const auto min_reactions = static_cast<osh::I64>(
          static_cast<osh::Real>(input.expected_num_reactions) *
          (1. - input.expected_num_reactions_tolerance));
      const auto max_reactions = static_cast<osh::I64>(
          static_cast<osh::Real>(input.expected_num_reactions) *
          (1. + input.expected_num_reactions_tolerance));
      if (num_reactions < min_reactions) {
        simulation.log_once(
            "Error: too few number of reactions, expected at least " +
            std::to_string(min_reactions));
        status = EXIT_FAILURE;
      } else if (num_reactions > max_reactions) {
        simulation.log_once(
            "Error: too many number of reactions, expected at most " +
            std::to_string(max_reactions));
        status = EXIT_FAILURE;
      }
    }
  }
  return status;
}

// explicit template instantiation definitions
template class Scenario<std::mt19937>;
template void
ScenarioResult::log_summary(const Simulation<std::mt19937> &simulation,
                            const ScenarioInput &input) const;
template ScenarioResult
test_splitting_operator(Simulation<std::mt19937> &simulation, int scenario,
                        const ScenarioInput &input);

} // namespace dist
} // namespace steps
