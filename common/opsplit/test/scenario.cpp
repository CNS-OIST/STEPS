#include <algorithm>
#include <hadoken/format/format.hpp>

#include "Ca_burst_background.hpp"
#include "diff_sim.hpp"
#include "multi_comp.hpp"
#include "simple_model.hpp"
#include "sreac_unittest.hpp"
#include "sreac_validation.hpp"
#include "validation.hpp"
#include "version.hpp"
#include <mpitools.hpp>
#include <opsplit/compdef.hpp>
#include <opsplit/diffdef.hpp>
#include <opsplit/mesh.hpp>
#include <opsplit/patchdef.hpp>
#include <opsplit/reacdef.hpp>
#include <opsplit/sreacdef.hpp>
#include <opsplit/test/scenario.hpp>
#include <opsplit/test/simulation.hpp>

namespace zee {

using hadoken::scat;

ScenarioResult ScenarioResult::invalid(int status) {
    return {status, 0, 0, 0, 0, 0, {}, 0, {}};
}

template <typename RNG>
Scenario<RNG>::Scenario(const std::string& t_name,
                        const std::string& t_description,
                        const ScenarioInput& t_input)
    : name(t_name)
    , description(t_description)
    , rank(mpi_comm_rank())
    , num_ranks(mpi_comm_size())
    , input(t_input) {}

template <typename RNG>
Scenario<RNG>::~Scenario() noexcept = default;

template <typename RNG>
ScenarioResult Scenario<RNG>::execute(simulation_t& simulation) {
    timemory_fixture ti_exec_timer("Scenario::execute");
    simulation.log_once("Initializing simulation");
    load_mesh(simulation);
    apply_statedef(simulation);
    const auto& result = run_simulation(simulation);
    result.ti_execution_timer = ti_exec_timer.stop();
    return result;
}

template <typename RNG>
PetscLogDouble Scenario<RNG>::max_runtime(PetscLogDouble runtime) {
    decltype(runtime) max_runtime{};
    {
        auto err =
            MPI_Reduce(&runtime, &max_runtime, 1, MPIU_PETSCLOGDOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (err != 0) {
            MPI_Abort(MPI_COMM_WORLD, err);
        }
    }
    return max_runtime;
}

template <typename RNG>
void Scenario<RNG>::load_mesh(Simulation<RNG>& simulation) {
    simulation.log_once(scat("Loading mesh file \"", input.mesh_file, '\"'));
    this->register_compartments(simulation.getMesh());
    if (input.log_mesh_report) {
        simulation.log_once(simulation.getMesh().createReport());
        this->log_owned_elements(simulation);
    }
}

template <typename RNG>
void Scenario<RNG>::log_owned_elements(const Simulation<RNG>& simulation) const {
    simulation.log_all(
        scat("Number of owned elements: ", simulation.getMesh().getOwnedNumElements()));
    MPI_Barrier(MPI_COMM_WORLD);
    simulation.log_once("");
}

template <typename RNG>
void Scenario<RNG>::apply_statedef(Simulation<RNG>& simulation) {
    auto statedef = this->createStatedef();
    if (input.log_statedef_report) {
        simulation.log_once(statedef->createReport());
    }
    simulation.init(std::move(statedef));
    this->fill_compartments(simulation);
    simulation.setDiffOpBinomialThreshold(input.threshold);
    if (input.log_state_report) {
        simulation.log_once(simulation.createStateReport());
    }
}

template <typename RNG>
void ScenarioResult::log_summary(const Simulation<RNG>& simulation,
                                 const ScenarioInput& input) const {
    simulation.log_once(scat("hpcbench_meta_zee_revision=", Zee_COMMIT));
    simulation.log_once(scat("hpcbench_meta_mesh_backend=", simulation.getMesh().getBackend()));
    simulation.log_once(scat("hpcbench_meta_mesh_file=", input.mesh_file));
    simulation.log_once(scat("hpcbench_meta_num_processes=", simulation.comm_size));
    simulation.log_once(scat("hpcbench_meta_num_elements=", num_elements));
    simulation.log_once(scat("hpcbench_meta_threshold=", input.threshold));
    simulation.log_once(scat("hpcbench_meta_scale=", input.scale));
    simulation.log_once(scat("hpcbench_meta_num_mols_factor=", input.num_mols_factor));
    simulation.log_once(scat("hpcbench_meta_end_time_sec=", input.end_time));
    simulation.log_once(scat("hpcbench_metric_num_iterations=", num_iterations));
    simulation.log_once(scat("hpcbench_metric_time_step_sec=", time_step));
    simulation.log_once(scat("hpcbench_metric_num_reactions=", num_reactions));
    simulation.log_once(scat("hpcbench_metric_elapsed_reactions_sec=",
                             static_cast<double>(elapsed_reactions.wall) / 1e9));
    simulation.log_once(scat("hpcbench_metric_num_diffusions=", num_diffusions));
    simulation.log_once(scat("hpcbench_metric_elapsed_diffusions_sec=",
                             static_cast<double>(elapsed_diffusions.wall) / 1e9));
    simulation.log_once(scat("hpcbench_metric_max_runtime=", max_runtime));
    simulation.log_once(scat("hpcbench_metric_status=", status));
#ifdef Zee_USE_TIMEMORY
    simulation.log_once(scat("hpcbench_metric_timemory_execute_peak_rss_MB=",
                             ti_execution_timer.get<timc::peak_rss>().get()));
    simulation.log_once(scat("hpcbench_metric_timemory_execute_data_rss_MB=",
                             ti_execution_timer.get<timc::data_rss>().get()));
    simulation.log_once(scat("hpcbench_metric_timemory_execute_real_time_sec=",
                             ti_execution_timer.get<timc::real_clock>().get()));
    simulation.log_once(scat("hpcbench_metric_timemory_init_exec_peak_rss_MB=",
                             ti_init_exec_timer.get<timc::peak_rss>().get()));
    simulation.log_once(scat("hpcbench_metric_timemory_init_exec_data_rss_MB=",
                             ti_init_exec_timer.get<timc::data_rss>().get()));
    simulation.log_once(scat("hpcbench_metric_timemory_init_exec_real_time_sec=",
                             ti_init_exec_timer.get<timc::real_clock>().get()));
    simulation.log_once(scat("hpcbench_metric_timemory_full_peak_rss_MB=",
                             ti_full_timer.get<timc::peak_rss>().get()));
    simulation.log_once(scat("hpcbench_metric_timemory_full_data_rss_MB=",
                             ti_full_timer.get<timc::data_rss>().get()));
    simulation.log_once(scat("hpcbench_metric_timemory_full_real_time_sec=",
                             ti_full_timer.get<timc::real_clock>().get()));
#endif  // Zee_USE_TIMEMORY
}

template <typename RNG>
ScenarioResult Scenario<RNG>::run_simulation(Simulation<RNG>& simulation) {
    PetscLogDouble time_beginning, time_end;
    simulation.log_once("Starting simulation...");
    PetscTime(&time_beginning);
    this->run_simulation_impl(simulation);
    PetscTime(&time_end);
    const auto status = this->check_and_log_results(simulation);
    return {status,
            simulation.getIterationTimeStep(),
            simulation.getNIterations(),
            this->max_runtime(time_end - time_beginning),
            simulation.getMesh().getTotalNumElements(),
            simulation.getSSAOpExtent(),
            simulation.getElapsedSSA(),
            simulation.getDiffOpExtent(),
            simulation.getElapsedDiff()};
}

template <typename RNG>
ScenarioResult test_splitting_operator(Simulation<RNG>& simulation,
                                       int scenario,
                                       const ScenarioInput& input) {
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
    default:
        simulation.log_once(scat("Unknown test id: ", scenario));
        return ScenarioResult::invalid(1);
    }
}

template <typename RNG>
int Scenario<RNG>::check_and_log_results(const Simulation<RNG>& simulation) const {
    auto status = check_and_log_results_impl(simulation);
    const auto num_diffusions = simulation.getDiffOpExtent();
    const auto num_reactions = simulation.getSSAOpExtent();
    if (rank == 0) {
        if (input.expected_num_diffusions >= 0) {
            const auto min_diffusions = static_cast<PetscInt>(
                static_cast<PetscScalar>(input.expected_num_diffusions) *
                (1. - input.expected_num_diffusions_tolerance));
            const auto max_diffusions = static_cast<PetscInt>(
                static_cast<PetscScalar>(input.expected_num_diffusions) *
                (1. + input.expected_num_diffusions_tolerance));
            if (num_diffusions < min_diffusions) {
                simulation.log_once(
                    scat("Error: too few number of diffusions, expected at most ", min_diffusions));
                status = EXIT_FAILURE;
            } else if (num_diffusions > max_diffusions) {
                simulation.log_once(scat("Error: too many number of diffusions, expected at most ",
                                         max_diffusions));
                status = EXIT_FAILURE;
            }
        }
        if (input.expected_num_reactions >= 0) {
            const auto min_reactions = static_cast<PetscInt>(
                static_cast<PetscScalar>(input.expected_num_reactions) *
                (1. - input.expected_num_reactions_tolerance));
            const auto max_reactions = static_cast<PetscInt>(
                static_cast<PetscScalar>(input.expected_num_reactions) *
                (1. + input.expected_num_reactions_tolerance));
            if (num_reactions < min_reactions) {
                simulation.log_once(
                    scat("Error: too few number of reactions, expected at least ", min_reactions));
                status = EXIT_FAILURE;
            } else if (num_reactions > max_reactions) {
                simulation.log_once(
                    scat("Error: too many number of reactions, expected at most ", max_reactions));
                status = EXIT_FAILURE;
            }
        }
    }
    return status;
}

// explicit template instantiation definitions
template class Scenario<std::mt19937>;
template void ScenarioResult::log_summary(const Simulation<std::mt19937>& simulation,
                                          const ScenarioInput& input) const;
template ScenarioResult test_splitting_operator(Simulation<std::mt19937>& simulation,
                                                int scenario,
                                                const ScenarioInput& input);


}  // namespace zee
