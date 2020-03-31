#pragma once

#include <memory>
#include <random>

#include <boost/timer/timer.hpp>
#include <memutils.hpp>
#include <opsplit/fwd.hpp>
#include <opsplit/vocabulary.hpp>

namespace zee {

using boost::timer::cpu_times;

struct ScenarioInput {
    ScenarioInput(std::string t_mesh_file,
                  PetscScalar t_threshold,
                  PetscScalar t_num_mols_factor,
                  PetscScalar t_scale,
                  model::compartment_label t_comp_label,
                  PetscScalar t_do_interval,
                  PetscScalar t_end_time,
                  PetscInt t_expected_num_reactions,
                  PetscScalar t_expected_num_reactions_tolerance,
                  PetscInt t_expected_num_diffusions,
                  PetscScalar t_expected_num_diffusions_tolerance,
                  int t_inject_in_elt,
                  bool t_log_mesh_report,
                  bool t_log_statedef_report,
                  bool t_log_state_report,
                  bool t_molecules_pools_force_dist_for_variable_sized)
        : mesh_file(std::move(t_mesh_file))
        , threshold(t_threshold)
        , num_mols_factor(t_num_mols_factor)
        , scale(t_scale)
        , comp_label(t_comp_label)
        , do_interval(t_do_interval)
        , end_time(t_end_time)
        , expected_num_reactions(t_expected_num_reactions)
        , expected_num_reactions_tolerance(t_expected_num_reactions_tolerance)
        , expected_num_diffusions(t_expected_num_diffusions)
        , expected_num_diffusions_tolerance(t_expected_num_diffusions_tolerance)
        , do_inject_in_elt(t_inject_in_elt)
        , log_mesh_report(t_log_mesh_report)
        , log_statedef_report(t_log_statedef_report)
        , log_state_report(t_log_state_report)
        , molecules_pools_force_dist_for_variable_sized(
              t_molecules_pools_force_dist_for_variable_sized) {}

    /// the mesh file to use during the simulation
    const std::string mesh_file;
    /// diffusion operator molecules threshold
    const PetscScalar threshold;
    /// Multiplication factor of the initial number of molecules
    const PetscScalar num_mols_factor;
    /// The import scale from the Gmsh unit to STEPS unit.
    /// e.g. if 1 unit in Gmsh represents 1 um, then scale=1e-6
    const PetscScalar scale;
    /// For DMPlex solutions. Provide the mesh label for
    /// tetrahedrons in a compartment.
    const model::compartment_label comp_label;
    /// step interval to use for the "diffusionOnly" scenario
    const PetscScalar do_interval;
    /// simulation end time
    const PetscScalar end_time;

    /// Expected number of reactions during the simulation
    /// the program will fail if number of reactions is not in
    /// the following interval
    /// [expected_num_reactions * (1 - expected_num_reactions_tolerance),
    ///  expected_num_reactions * (1 + expected_num_reactions_tolerance)]
    const PetscInt expected_num_reactions;
    /// tolerance between 0 and 1
    const PetscScalar expected_num_reactions_tolerance;
    /// Expected number of reactions during the simulation
    /// the program will fail if number of reactions is not in
    /// the following interval
    /// [expected_num_reactions * (1 - expected_num_reactions_tolerance),
    ///  expected_num_reactions * (1 + expected_num_reactions_tolerance)]
    const PetscInt expected_num_diffusions;
    /// tolerance between 0 and 1
    const PetscScalar expected_num_diffusions_tolerance;

    /// For the "diffusionOnly" scenario. Tells in which element to
    /// inject molecules. Default is -1 meaning that molecules are
    /// spread in all elements.
    const int do_inject_in_elt;
    /// Log the mesh report before starting the simulation
    bool log_mesh_report;
    /// log the statedef report before starting the simulation
    bool log_statedef_report;
    /// log the state report before starting the simulation
    bool log_state_report;

    /// if true, force usage of this Omega_h synchronization method
    /// even if all elements have the same number of species
    bool molecules_pools_force_dist_for_variable_sized;
};

struct ScenarioResult {
    ScenarioResult(int t_status,
                   PetscScalar t_time_step,
                   PetscInt64 t_num_iterations,
                   PetscScalar t_max_runtime,
                   PetscInt t_num_elements,
                   PetscInt64 t_num_reactions,
                   cpu_times t_elapsed_reactions,
                   PetscInt64 t_num_diffusions,
                   cpu_times t_elapsed_diffusions)
        : status(t_status)
        , time_step(t_time_step)
        , num_iterations(t_num_iterations)
        , max_runtime(t_max_runtime)
        , num_elements(t_num_elements)
        , num_reactions(t_num_reactions)
        , elapsed_reactions(t_elapsed_reactions)
        , num_diffusions(t_num_diffusions)
        , elapsed_diffusions(t_elapsed_diffusions) {}

    static ScenarioResult invalid(int status);

    /// status code the process may return
    const int status;
    /// simulation time step
    PetscScalar time_step;
    /// number of iterations that occured during the simulation
    PetscInt64 num_iterations;
    /// maximum simulation runtime among the ranks
    const PetscScalar max_runtime;
    /// total number of elements in the mesh
    PetscInt num_elements;
    /// total number of reactions that occurred during the simulation
    PetscInt64 num_reactions;
    /// total time spent performing reactions
    cpu_times elapsed_reactions;
    /// total number of diffusions that occurred during the simulation
    PetscInt64 num_diffusions;
    /// total time spent performing diffusions
    cpu_times elapsed_diffusions;

    template <typename RNG>
    void log_summary(const Simulation<RNG>& simulation, const ScenarioInput& input) const;

    /// timemory timer to measure Scenario::execute
    mutable measurement_t ti_execution_timer;
    /// timemory timer to measure test_splitting_operator
    mutable measurement_t ti_init_exec_timer;
    /// timemory timer to measure client main function
    mutable measurement_t ti_full_timer;
};

template <typename RNG>
class Scenario {
  public:
    using simulation_t = Simulation<RNG>;
    Scenario(const std::string& name,
             const std::string& t_description,
             const ScenarioInput& t_input);
    virtual ~Scenario() noexcept;
    /// Execute the given simulation
    ScenarioResult execute(simulation_t& simulation);

    /// brief name of the scenario
    const std::string& name;
    /// longer description of the simulation
    const std::string& description;

  protected:
    /** \brief The following member functions are called in this particular order
     * by \a execute member function
     * */

    /**
     * \return the container to use for this simulation
     */
    virtual std::unique_ptr<Statedef> createStatedef() const = 0;

    /**
     * add the proper compartments to the Mesh object
     */
    virtual void register_compartments(DistMesh& mesh) const = 0;

    /**
     * Add molecules in the compartments
     */
    virtual void fill_compartments(simulation_t& simulation) const = 0;

    /**
     * Main loop
     */
    virtual void run_simulation_impl(simulation_t& simulation) = 0;

    /**
     * Check that simulation ran properly, and log results if any.
     * \return the status code the process may return
     */
    virtual int check_and_log_results_impl(const simulation_t& /* simulation */) const {
        return EXIT_SUCCESS;
    }

    /// MPI rank of the current process
    const int rank;
    /// Number of MPI ranks involved in the simulation
    const int num_ranks;
    /// Simulation input
    const ScenarioInput& input;

  private:
    void load_mesh(simulation_t& simulation);
    void log_owned_elements(const simulation_t& simulation) const;
    void apply_statedef(simulation_t& simulation);
    PetscLogDouble max_runtime(PetscLogDouble runtime);
    ScenarioResult run_simulation(simulation_t& simulation);
    int check_and_log_results(const simulation_t& simulation) const;
};

template <typename RNG>
ScenarioResult test_splitting_operator(Simulation<RNG>& simulation,
                                       int scenario,
                                       const ScenarioInput& input);

// explicit template instantiation declarations
extern template class Scenario<std::mt19937>;
extern template void ScenarioResult::log_summary(const Simulation<std::mt19937>& simulation,
                                                 const ScenarioInput& input) const;
extern template ScenarioResult test_splitting_operator(Simulation<std::mt19937>& simulation,
                                                       int scenario,
                                                       const ScenarioInput& input);

}  // namespace zee
