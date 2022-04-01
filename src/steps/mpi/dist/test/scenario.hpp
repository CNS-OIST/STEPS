#pragma once

#include <memory>
#include <random>

#include "geom/dist/fwd.hpp"
#include "mpi/dist/tetopsplit/definition/fwd.hpp"
#include "util/vocabulary.hpp"
#include "util/tracker/time_tracker.hpp"

namespace steps {
namespace dist {

struct ScenarioInput {
    ScenarioInput(std::string t_mesh_file,
                  osh::Real t_threshold,
                  osh::Real t_num_mols_factor,
                  osh::Real t_efield_dt,
                  osh::Real t_scale,
                  model::compartment_label t_comp_label,
                  osh::Real t_do_interval,
                  osh::Real t_end_time,
                  osh::I64 t_expected_num_reactions,
                  osh::Real t_expected_num_reactions_tolerance,
                  osh::I64 t_expected_num_diffusions,
                  osh::Real t_expected_num_diffusions_tolerance,
                  int t_inject_in_elt,
                  bool t_log_statedef_report,
                  bool t_log_state_report,
                  const std::string& t_logfile,
                  const std::string& t_depgraphfile,
                  bool t_molecules_pools_force_dist_for_variable_sized,
                  bool t_sim_indep_k_procs)
        : mesh_file(std::move(t_mesh_file))
        , threshold(t_threshold)
        , num_mols_factor(t_num_mols_factor)
        , efield_dt(t_efield_dt)
        , scale(t_scale)
        , comp_label(t_comp_label)
        , do_interval(t_do_interval)
        , end_time(t_end_time)
        , expected_num_reactions(t_expected_num_reactions)
        , expected_num_reactions_tolerance(t_expected_num_reactions_tolerance)
        , expected_num_diffusions(t_expected_num_diffusions)
        , expected_num_diffusions_tolerance(t_expected_num_diffusions_tolerance)
        , do_inject_in_elt(t_inject_in_elt)
        , log_statedef_report(t_log_statedef_report)
        , log_state_report(t_log_state_report)
        , logfile(t_logfile)
        , depgraphfile(t_depgraphfile)
        , molecules_pools_force_dist_for_variable_sized(
              t_molecules_pools_force_dist_for_variable_sized)
        , sim_indep_k_procs(t_sim_indep_k_procs) {}

    ScenarioInput(std::string t_mesh_file, osh::Real t_scale = 1.0e-6)
        : mesh_file(std::move(t_mesh_file))
        , scale(t_scale) {}

    /// the mesh file to use during the simulation
    const std::string mesh_file;
    /// diffusion operator molecules threshold
    const osh::Real threshold{10.0};
    /// Multiplication factor of the initial number of molecules
    const osh::Real num_mols_factor{1.0};
    /// E-field time step. Assigned by CLI
    const osh::Real efield_dt{1.0e-5};
    /// The import scale from the Gmsh unit to STEPS unit.
    /// e.g. if 1 unit in Gmsh represents 1 um, then scale=1e-6
    const osh::Real scale{1.0e-6};
    /// For DMPlex solutions. Provide the mesh label for
    /// tetrahedrons in a compartment.
    const model::compartment_label comp_label{100};
    /// step interval to use for the "diffusionOnly" scenario
    const osh::Real do_interval{1e-7};
    /// simulation end time
    const osh::Real end_time{20.0};

    /// Expected number of reactions during the simulation
    /// the program will fail if number of reactions is not in
    /// the following interval
    /// [expected_num_reactions * (1 - expected_num_reactions_tolerance),
    ///  expected_num_reactions * (1 + expected_num_reactions_tolerance)]
    const osh::I64 expected_num_reactions{-1};
    /// tolerance between 0 and 1
    const osh::Real expected_num_reactions_tolerance{0.01};
    /// Expected number of reactions during the simulation
    /// the program will fail if number of reactions is not in
    /// the following interval
    /// [expected_num_reactions * (1 - expected_num_reactions_tolerance),
    ///  expected_num_reactions * (1 + expected_num_reactions_tolerance)]
    const osh::I64 expected_num_diffusions{-1};
    /// tolerance between 0 and 1
    const osh::Real expected_num_diffusions_tolerance{0.02};

    /// For the "diffusionOnly" scenario. Tells in which element to
    /// inject molecules. Default is -1 meaning that molecules are
    /// spread in all elements.
    const int do_inject_in_elt{-1};
    /// log the statedef report before starting the simulation
    bool log_statedef_report{false};
    /// log the state report before starting the simulation
    bool log_state_report{false};
    /// file where the log is printed. If no file is provided, it is printed to clog
    std::string logfile;
    /// file where the dependency graph is printed. It should end in .dot
    std::string depgraphfile;

    /// if true, force usage of this Omega_h synchronization method
    /// even if all elements have the same number of species

    bool molecules_pools_force_dist_for_variable_sized{false};
    /// simulate independently groups of independent kinetic processes
    bool sim_indep_k_procs;
};

struct ScenarioResult {
    ScenarioResult(int t_status,
                   osh::Real t_time_step,
                   boost::optional<osh::Real> t_efield_dt,
                   osh::I64 t_num_iterations,
                   double t_elapsed_sim_sec,
                   osh::I64 t_num_elements,
                   double t_elapsed_initialization_sec,
                   osh::I64 t_num_reactions,
                   double t_elapsed_reactions,
                   osh::I64 t_num_diffusions,
                   double t_elapsed_diffusions,
                   boost::optional<double> t_elapsed_efield)
        : status(t_status)
        , time_step(t_time_step)
        , efield_dt(t_efield_dt)
        , num_iterations(t_num_iterations)
        , elapsed_sim_sec(t_elapsed_sim_sec)
        , num_elements(t_num_elements)
        , elapsed_initialization_sec(t_elapsed_initialization_sec)
        , num_reactions(t_num_reactions)
        , elapsed_reactions(t_elapsed_reactions)
        , num_diffusions(t_num_diffusions)
        , elapsed_diffusions(t_elapsed_diffusions)
        , elapsed_efield(t_elapsed_efield) {}

    static ScenarioResult invalid(int status);

    /// status code the process may return
    const int status;
    /// simulation time step
    osh::Real time_step;
    /// max E-field dt
    const boost::optional<osh::Real> efield_dt;
    /// number of iterations that occured during the simulation
    osh::I64 num_iterations;
    /// total time spent by the simulations
    const osh::Real elapsed_sim_sec;
    /// total number of elements in the mesh
    osh::I64 num_elements;
    /// time spent during the initialization of the simulation
    double elapsed_initialization_sec;
    /// total number of reactions that occurred during the simulation
    osh::I64 num_reactions;
    /// total time in seconds spent performing reactions
    double elapsed_reactions;
    /// total number of diffusions that occurred during the simulation
    osh::I64 num_diffusions;
    /// total time in seconds spent performing diffusions
    double elapsed_diffusions;
    /// total time in seconds spent performing efield
    boost::optional<double> elapsed_efield;

    template <typename RNG>
    void log_summary(const Simulation<RNG>& simulation, const ScenarioInput& input) const;
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

    /// Initialize and fill the mesh with molecules
    void initialize(simulation_t& simulation);

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
    virtual std::unique_ptr<Statedef> createStatedef(const simulation_t& simulation) const = 0;

    /**
     * add the proper compartments to the Mesh object
     */
    virtual void register_compartments(DistMesh& mesh) const = 0;

    /**
     * Add molecules in the compartments
     */
    virtual void fill_compartments(simulation_t& simulation) const = 0;

    /**
     * Add molecules in patches
     */
    virtual void fill_patches(simulation_t& /* simulation */) const {}

    /**
     * Main loop
     */
    virtual void run_simulation_impl(simulation_t& simulation) = 0;

    /**
     * Check that simulation ran properly, and log results if any.
     * \return the status code the process may return
     */
    virtual int check_and_log_results_impl(simulation_t& /* simulation */) const {
        return EXIT_SUCCESS;
    }

    /// Simulation input
    const ScenarioInput& input;
    /// whether method \a initialize has been called or not
    bool initialized{false};

    util::TimeTracker timer_initialization;

  private:
    void load_mesh(simulation_t& simulation);
    void log_owned_elements(const simulation_t& simulation) const;
    void apply_statedef(simulation_t& simulation);
    osh::Real max_runtime(const simulation_t&, osh::Real runtime);
    ScenarioResult run_simulation(simulation_t& simulation);
    int check_and_log_results(simulation_t& simulation) const;
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

}  // namespace dist
} // namespace steps
