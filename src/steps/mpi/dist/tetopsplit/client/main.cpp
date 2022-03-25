#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_shape.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graphviz.hpp>

#ifdef OMEGA_H_USE_GMSH
#include <gmsh.h>
#endif // OMEGA_H_USE_GMSH

// clang-format off
#include "mpi/dist/test/scenario.hpp"
// clang-format on

#include "geom/dist/distmesh.hpp"
#include "mpi/dist/tetopsplit/simulation.hpp"
#include "util/debug.hpp"
#include "util/mesh.hpp"
#include "util/mpitools.hpp"
#include "util/profile/profiler_interface.h"
//#include "test_rssa.hpp"

namespace steps {
namespace dist {

// see README.md for description

// Execution example on 8 cores. 4 physical CPU, 1 socket with OpenMP 3.0
// export OMP_NUM_THREADS=2
// export OMP_PROC_BIND=true
// mpirun --bind-to core -np 4 ./OpSplitOmega_h

#undef MPI_Allreduce

/**
 * A wrapper that provides a convenient RAII-style mechanism to manipulate GMSH
 * library.
 */
struct GmshFixture {
  GmshFixture() {
#ifdef OMEGA_H_USE_GMSH
    ::gmsh::initialize();
#endif // OMEGA_H_USE_GMSH
  }

  ~GmshFixture() {
#ifdef OMEGA_H_USE_GMSH
    ::gmsh::finalize();
#endif // OMEGA_H_USE_GMSH
  }
};

template <osh::Int Dim> struct MeshTraits {};

template <> struct MeshTraits<2> {
  static const int boundary_type = osh::EDGE;
};

template <> struct MeshTraits<3> {
  static const int boundary_type = osh::FACE;
};

/**
 * Provide coordinates of the point where molecules are inserted at the
 * beginning of the simulation \tparam Dim mesh dimension
 */
template <osh::Int Dim> osh::Vector<Dim> seed_point() {
  return {util::initializer_list<Dim>(0.5)};
}


static int splitting_operator(osh::Mesh &mesh, Omega_h::CmdLine &cmdline) {
  {
    int debugRank{-1};
    if (cmdline.parsed("--debug-rank")) {
      debugRank = cmdline.get<int>("--debug-rank", "value");
      if (debugRank >= 0) {
        const auto rank = util::mpi_comm_rank(MPI_COMM_WORLD);
        if (rank == debugRank) {
          util::wait_for_gdb();
        }
      }
    }
  }

  const auto mesh_file = cmdline.get<std::string>("square.msh");
  const bool use_big_int_for_num_mols =
      cmdline.parsed("--use-big-int-for-num-mols");
  int scenario = -1;
  if (cmdline.parsed("--test")) {
    scenario = cmdline.get<int>("--test", "number");
  }

  osh::Real threshold = 10.0;
  if (cmdline.parsed("--threshold")) {
    threshold = cmdline.get<int>("--threshold", "value");
  }

  osh::Real scale{1.0e-6};
  if (cmdline.parsed("--scale")) {
    scale = cmdline.get<osh::Real>("--scale", "value");
  }

  int comp_label{100};
  if (cmdline.parsed("--comp-label")) {
    comp_label = cmdline.get<int>("--comp-label", "value");
  }

  osh::Real num_mols_factor{1};
  if (cmdline.parsed("--num-mols-factor")) {
    num_mols_factor = cmdline.get<osh::Real>("--num-mols-factor", "value");
  }

  osh::Real efield_dt{1.0e-5};
  if (cmdline.parsed("--efield-dt")) {
    efield_dt = cmdline.get<osh::Real>("--efield-dt", "value");
  }

  bool log_statedef_report{false};
  if (cmdline.parsed("--log-statedef-report")) {
    log_statedef_report = true;
  }

  bool log_state_report{false};
  if (cmdline.parsed("--log-state-report")) {
    log_state_report = true;
  }

  if (cmdline.parsed("--log-all-reports")) {
    log_statedef_report = true;
    log_state_report = true;
  }

  std::string logfile = "";
  if (cmdline.parsed("--logfile")) {
    logfile = cmdline.get<std::string>("--logfile", "string");
  }

  std::string depgraphfile = "";
  if (cmdline.parsed("--depgraphfile")) {
      depgraphfile = cmdline.get<std::string>("--depgraphfile", "string");
  }

  bool molecules_pools_force_dist_for_variable_sized = false;
  if (cmdline.parsed("--molecules-pools-force-dist-for-variable-sized")) {
    molecules_pools_force_dist_for_variable_sized = true;
  }

  unsigned rng_seed = std::random_device()();
  if (cmdline.parsed("--rng-seed")) {
    rng_seed = static_cast<unsigned>(cmdline.get<int>("--rng-seed", "value"));
  }

  if (cmdline.parsed("--rng-seed-add-process-rank")) {
    rng_seed += mesh.comm()->rank();
  }

  osh::Real do_interval{1e-7};
  if (cmdline.parsed("--do-interval")) {
    do_interval = cmdline.get<osh::Real>("--do-interval", "value");
  }

  osh::Real end_time{0};
  if (cmdline.parsed("--end-time")) {
    end_time = cmdline.get<osh::Real>("--end-time", "value");
  }

  int do_inject_in_elt = -1;
  if (cmdline.parsed("--do-inject-in-elt")) {
    do_inject_in_elt = cmdline.get<int>("--do-inject-in-elt", "value");
  }

  int expected_num_reactions{-1};
  if (cmdline.parsed("--expected-reactions")) {
    expected_num_reactions = cmdline.get<int>("--expected-reactions", "value");
  }

  osh::Real expected_num_reactions_tolerance{0.01};
  if (cmdline.parsed("--expected-reactions-tolerance")) {
    expected_num_reactions_tolerance =
        cmdline.get<osh::Real>("--expected-reactions-tolerance", "value");
  }

  int expected_num_diffusions{-1};
  if (cmdline.parsed("--expected-diffusions")) {
    expected_num_diffusions =
        cmdline.get<int>("--expected-diffusions", "value");
  }

  osh::Real expected_num_diffusions_tolerance{0.02};
  if (cmdline.parsed("--expected-diffusions-tolerance")) {
    expected_num_diffusions_tolerance =
        cmdline.get<osh::Real>("--expected-diffusions-tolerance", "value");
  }

  DistMesh mesh_wrapper(mesh, mesh_file, scale);

  const ScenarioInput input(mesh_file,
                            threshold,
                            num_mols_factor,
                            efield_dt,
                            scale,
                            model::compartment_label(comp_label),
                            do_interval,
                            end_time,
                            expected_num_reactions,
                            expected_num_reactions_tolerance,
                            expected_num_diffusions,
                            expected_num_diffusions_tolerance,
                            do_inject_in_elt,
                            log_statedef_report,
                            log_state_report,
                            logfile,
                            depgraphfile,
                            molecules_pools_force_dist_for_variable_sized,
                            cmdline.parsed("--independent-kprocs"));

  const auto run_sim = [&](auto &sim) {
    const auto &result = test_splitting_operator(sim, scenario, input);
    result.log_summary(sim, input);
    return result.status;
  };

  std::mt19937 rng(rng_seed);
  std::ofstream outfile;
  if (!input.logfile.empty()) {
    outfile.open(input.logfile, std::ofstream::out);
  }
  std::ostream &outstream = outfile.is_open() ? outfile : std::clog;

  if (cmdline.parsed("--use-rssa")) {
    if (cmdline.parsed("--use-gibson-bruck")) {
      throw std::logic_error("Can't use RSSA with Gibson-Bruck");
    } else {
      if (use_big_int_for_num_mols) {
        OmegaHSimulation<SSAMethod::RSSA, decltype(rng), osh::I64,
                         NextEventSearchMethod::Direct>
            simulation(input, mesh_wrapper, rng, outstream);
        return run_sim(simulation);
      } else {
        OmegaHSimulation<SSAMethod::RSSA, decltype(rng), osh::I32,
                         NextEventSearchMethod::Direct>
            simulation(input, mesh_wrapper, rng, outstream);
        return run_sim(simulation);
      }
    }

  } else {
    if (cmdline.parsed("--use-gibson-bruck")) {
      if (use_big_int_for_num_mols) {
        OmegaHSimulation<SSAMethod::SSA, decltype(rng), osh::I64,
                         NextEventSearchMethod ::GibsonBruck>
            simulation(input, mesh_wrapper, rng, outstream);
        return run_sim(simulation);
      } else {
        OmegaHSimulation<SSAMethod::SSA, decltype(rng), osh::I32,
                         NextEventSearchMethod ::GibsonBruck>
            simulation(input, mesh_wrapper, rng, outstream);
        return run_sim(simulation);
      }
    } else {
      if (use_big_int_for_num_mols) {
        OmegaHSimulation<SSAMethod::SSA, decltype(rng), osh::I64,
                         NextEventSearchMethod ::Direct>
            simulation(input, mesh_wrapper, rng, outstream);
        return run_sim(simulation);
      } else {
        OmegaHSimulation<SSAMethod::SSA, decltype(rng), osh::I32,
                         NextEventSearchMethod ::Direct>
            simulation(input, mesh_wrapper, rng, outstream);
        return run_sim(simulation);
      }
    }
  }
}

} // namespace dist
} // namespace steps

static inline bool is_power_of_2(int value) {
  return !(value == 0) && !(value & (value - 1));
}

/**
 *  Create Omega_h mesh object from a GMSH mesh
 *  \param filename use parallel import if the prefix (without leading '_')
 *  of a multi-part mesh is given, let Omega_h do the partitioning otherwise.
 *  \param world Omega_h MPI wrapper
 *  \return an Omega_h mesh object
 */
static Omega_h::Mesh load_mesh(const std::string &filename,
                               const Omega_h::CommPtr &world) {
  const auto rank0 = world->rank() == 0;
#ifdef OMEGA_H_USE_GMSH
  {
    using namespace Omega_h::filesystem;
    std::ostringstream part_name;
    part_name << filename << '_' << world->rank() + 1 << ".msh";
    if (exists(path(part_name.str()))) {
      if (rank0) {
        std::clog << "Creating Omega_h mesh from GMSH mesh partition "
                  << part_name.str() << '\n';
      }
      return Omega_h::gmsh::read_parallel(filename, world);
    }
  }
#endif // OMEGA_H_USE_GMSH
  if (rank0) {
    std::clog << "Creating Omega_h mesh from GMSH mesh " << filename << '\n';
  }
  return Omega_h::gmsh::read(filename, world);
}

int main(int argc, char **argv) {
  steps::Instrumentor::init_profile();

  auto lib = steps::osh::Library(&argc, &argv);

  if (lib.world()->size() > 1 && !is_power_of_2(lib.world()->size())) {
    if (lib.world()->rank() == 0) {
      std::cerr << "Warning: Omega_h requires the number of processes to be a power of 2 to split the mesh. It must "
                   "be already split.\n";
    }
  }
  //    auto mesh = osh::gmsh::read("square.msh", lib.world());
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("square.msh");

  cmdline.add_flag("--use-big-int-for-num-mols",
                   "Use I64 instead of I32 for MolNumbers");

  auto &test_flag = cmdline.add_flag("--test", "scenario");
  test_flag.add_arg<int>("number");

  auto &threshold_flag = cmdline.add_flag("--threshold", "threshold");
  threshold_flag.add_arg<int>("value");

  auto &scale_flag = cmdline.add_flag("--scale", "scale");
  scale_flag.add_arg<steps::osh::Real>("value");

  auto &comp_label_flag = cmdline.add_flag("--comp-label", "compLabel");
  comp_label_flag.add_arg<int>("value");

  auto &num_mols_factor_flag =
      cmdline.add_flag("--num-mols-factor", "numMolsFactor");
  num_mols_factor_flag.add_arg<steps::osh::Real>("value");

  auto &efield_dt_flag = cmdline.add_flag("--efield-dt", "EfieldDtMax");
  efield_dt_flag.add_arg<steps::osh::Real>("value");

  cmdline.add_flag("--log-mesh-report", "logging");
  cmdline.add_flag("--log-statedef-report", "logging");
  cmdline.add_flag("--log-state-report", "logging");
  cmdline.add_flag("--log-all-reports", "logging");

  auto &logfile_flag = cmdline.add_flag("--logfile", "logging");
  logfile_flag.add_arg<std::string>("string");

  auto& depgraphfile_flag = cmdline.add_flag("--depgraphfile",
                                             "Print the dependency graph in the provided path (it"
                                             " should end in .dot)");
  depgraphfile_flag.add_arg<std::string>("string");

  cmdline.add_flag("--molecules-pools-force-dist-for-variable-sized", "mpi");

  auto &num_iterations = cmdline.add_flag("--num-iterations", "numIterations");
  num_iterations.add_arg<int>("value");

  auto &rng_seed_flag = cmdline.add_flag("--rng-seed", "RngSeed");
  rng_seed_flag.add_arg<int>("value");
  cmdline.add_flag("--rng-seed-add-process-rank",
                   "add process rank to rng seed");

  auto &expected_num_diffusions_flag =
      cmdline.add_flag("--expected-diffusions", "numDiffusions");
  expected_num_diffusions_flag.add_arg<int>("value");

  auto &expected_num_diffusions_tolerance_flag =
      cmdline.add_flag("--expected-diffusions-tolerance", "tolerance");
  expected_num_diffusions_tolerance_flag.add_arg<steps::osh::Real>("value");

  auto &expected_num_reactions_flag =
      cmdline.add_flag("--expected-reactions", "numReactions");
  expected_num_reactions_flag.add_arg<int>("value");

  auto &expected_num_reactions_tolerance_flag =
      cmdline.add_flag("--expected-reactions-tolerance", "tolerance");
  expected_num_reactions_tolerance_flag.add_arg<steps::osh::Real>("value");

  // diffusionOnly parameters
  auto &do_interval_flag =
      cmdline.add_flag("--do-interval", "diffusionOnlyInterval");
  do_interval_flag.add_arg<steps::osh::Real>("value");

  auto &endtime_flag = cmdline.add_flag("--end-time", "simulationEndTime");
  endtime_flag.add_arg<steps::osh::Real>("value");

  auto &do_inject_in_elt_flag =
      cmdline.add_flag("--do-inject-in-elt", "InjectInElement");
  do_inject_in_elt_flag.add_arg<int>("value");

  cmdline.add_flag("--use-rssa", "Use RSSA");
  cmdline.add_flag("--use-gibson-bruck",
                   "use Gibson-Bruck method to determine next kinetic event");
  cmdline.add_flag("--independent-kprocs",
                   "Simulate unrelated subsets of kprocs independently");

  auto &debug_rank_flag = cmdline.add_flag("--debug-rank", "mpiRank");
  debug_rank_flag.add_arg<int>("value");

  // end diffusionOnly parameters

  steps::dist::GmshFixture gmsh_context;
#if USE_PETSC
  // linear solver options
  cmdline.add_flag("-ksp_view",
                   "PETSC option (efield operator), print krylov solver info");
  cmdline.add_flag("-ksp_monitor",
                   "PETSC option (efield operator), print residual");
  cmdline.add_flag("-ksp_monitor_short",
                   "PETSC option (efield operator), print residual");
  auto &ksp_type = cmdline.add_flag(
      "-ksp_type", "PETSC option (efield operator), linear solver type");
  ksp_type.add_arg<std::string>("value");
  auto &ksp_rtol = cmdline.add_flag(
      "-ksp_rtol", "PETSC option (efield operator), set relative tol");
  ksp_rtol.add_arg<double>("value");
  auto &ksp_atol = cmdline.add_flag(
      "-ksp_atol", "PETSC option (efield operator), set absolute tol");
  ksp_atol.add_arg<double>("value");
  //-ksp_divtol <dtol>, and-ksp_max_it
  auto &ksp_pc_side = cmdline.add_flag(
      "-ksp_pc_side",
      "PETSC option (efield operator), preconditioner side [left|right]");
  ksp_pc_side.add_arg<std::string>("value");
  auto &ksp_norm_type = cmdline.add_flag(
      "-ksp_norm_type",
      "PETSC option (efield operator), norm [unpreconditioned]");
  ksp_norm_type.add_arg<std::string>("value");
  auto &ksp_initial_guess_nonzero =
      cmdline.add_flag("-ksp_initial_guess_nonzero",
                       "PETSC option (efield operator), [true|false]");
  ksp_initial_guess_nonzero.add_arg<std::string>("value");
  // GMRES
  auto &ksp_gmres_restart =
      cmdline.add_flag("-ksp_gmres_restart",
                       "PETSC option (efield operator), set gmres restart");
  ksp_gmres_restart.add_arg<int>("value");
  // preconditioner options
  auto &pc_type = cmdline.add_flag(
      "-pc_type", "PETSC option (efield operator), preconditioner type");
  pc_type.add_arg<std::string>("value");
  //
  steps::util::PetscFixture petsc_fixture(&argc, &argv, nullptr,
                                          "OpSplit solver");
#endif // USE_PETSC

  if (!cmdline.parse_final(lib.world(), &argc, argv)) {
    return EXIT_FAILURE;
  }

  const auto mesh_file = cmdline.get<std::string>("square.msh");
  auto mesh = load_mesh(mesh_file, lib.world());

  auto ret = steps::dist::splitting_operator(mesh, cmdline);
  
  steps::Instrumentor::finalize_profile();
  
  return ret;
}
