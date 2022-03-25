#include "steps/mpi/dist/tetopsplit/kproc/kproc_id.hpp"

#define CATCH_CONFIG_RUNNER
#include <Omega_h_file.hpp>
#include <catch2/catch.hpp>
#include <mpi.h>

#include "test_common.hpp"

#include "steps/geom/dist/distmesh.hpp"
#include "steps/mpi/dist/test/ca_burst_background.hpp"
#include "steps/mpi/dist/test/multi_comp.hpp"
#include "steps/mpi/dist/tetopsplit/simulation.hpp"
#include "steps/util/finish.hpp"
#include "steps/util/init.hpp"

enum class My2Enum { VAL1 = 0, VAL2 = 1 };
enum class My3Enum { VAL1 = 0, VAL2 = 1, VAL3 = 3 };

/**
 * Singleton class used by test cases to access global variables
 */
class Context {
public:
  using simulation_t =
      steps::dist::OmegaHSimulation<steps::dist::SSAMethod::SSA, std::mt19937>;

  Context(int *argc, char ***argv) : library_(argc, argv, MPI_COMM_WORLD) {
    steps::init();
  }

  ~Context() { steps::finish(); }

  /// internal class used to reset the simulation pointer when not required
  /// anymore
  struct Fixture {
    explicit Fixture(Context &context) : context_(&context) {}
    Fixture(Fixture &&other) noexcept { std::swap(context_, other.context_); }
    Fixture(const Fixture &other) = delete;
    ~Fixture() {
      if (context_ != nullptr) {
        context_->simulation_ = nullptr;
      }
    }

  private:
    Context *context_{};
  };

  /**
   * globally set the simulation object pointer. The pointer is reset when the
   * Fixture object is destructed.
   */
  Fixture simulation(simulation_t &simulation) {
    simulation_ = &simulation;
    return Fixture(*this);
  }

  /// Get global simulation instance.
  simulation_t &simulation() const {
    assert(simulation_ != nullptr);
    return *simulation_;
  }

  /// Get reference to the project source directory
  Omega_h::filesystem::path source_dir() const noexcept {
    return STEPS_SOURCE_DIR;
  }

  /// Get reference to the project binary directory
  Omega_h::filesystem::path binary_dir() const noexcept {
    return STEPS_BINARY_DIR;
  }

  /// Get Omega_h root object
  Omega_h::Library &library() noexcept { return library_; }

private:
  /// simulation object used by test cases
  simulation_t *simulation_{};
  /// top Omega_h object
  Omega_h::Library library_;
};

/// global context used by test cases
static Context *context;

TEST_CASE("ca_burst_mol_count", "[ca_burst_background]") {
  auto &sim = context->simulation();
  {
    // check number of molecules in elements
    std::map<std::string, steps::osh::Real> mols_counts{
        {"Mg", 4.9854e+08},      {"iCBsf", 2.34094e+07},
        {"iCBCaf", 2.22839e+06}, {"iCBsCa", 1.27998e+06},
        {"iCBCaCa", 121846},     {"CBsf", 9.36409e+07},
        {"CBCaf", 8.91372e+06},  {"CBsCa", 5.12017e+06},
        {"CBCaCa", 487403},      {"PV", 2.70952e+06},
        {"PVCa", 1.37327e+07},   {"PVMg", 5.11561e+07}};
    for (const auto &molecule : mols_counts) {
      const auto mols_count = sim.getCompCount(
          "__MESH__", steps::dist::model::species_name(molecule.first));
      REQUIRE(mols_count == Approx(molecule.second));
    }
  }
  {
    // check number of molecules in boundaries
    REQUIRE(sim.getPatchCount("smooth.__BOUNDARY__", "Pump") == Approx(2067.0));
    REQUIRE(sim.getPatchCount("spiny.__BOUNDARY__", "Pump") == Approx(22131.0));
    REQUIRE(sim.getPatchCount("smooth.__BOUNDARY__", "CaPump") == Approx(0.0));
    REQUIRE(sim.getPatchCount("spiny.__BOUNDARY__", "CaPump") == Approx(0.0));
  }
}

TEST_CASE("multi_comp_mol_count", "[two_comp_cylinder]") {
  auto &sim = context->simulation();
  REQUIRE(sim.getCompCount("Left", "A") == Approx(10.));
  REQUIRE(sim.getCompCount("Left", "B") == Approx(10.));
  REQUIRE(sim.getCompCount("Left", "C") == Approx(0.));
  REQUIRE(sim.getCompCount("Right", "C") == Approx(0.));
  REQUIRE(sim.getCompCount("Right", "D") == Approx(0.));
  REQUIRE(sim.getCompCount("Right", "E") == Approx(10.));

  std::array<Omega_h::Real, 6> local_counts{
      sim.getOwnedCompCount("Left", "A"),  sim.getOwnedCompCount("Left", "B"),
      sim.getOwnedCompCount("Left", "C"),  sim.getOwnedCompCount("Right", "C"),
      sim.getOwnedCompCount("Right", "D"), sim.getOwnedCompCount("Right", "E"),
  };

  std::array<Omega_h::Real, 6> counts{};
  MPI_Reduce(local_counts.data(),
             counts.data(),
             local_counts.size(),
             MPI_REAL8,
             MPI_SUM,
             0,
             sim.getMesh().comm_impl());
  if (sim.comm_rank == 0) {
    REQUIRE(counts[0] == Approx(10));
    REQUIRE(counts[1] == Approx(10));
    REQUIRE(counts[2] == Approx(0));
    REQUIRE(counts[3] == Approx(0));
    REQUIRE(counts[4] == Approx(0));
    REQUIRE(counts[5] == Approx(10));
  }
}

TEST_CASE("multi_comp_measure", "[two_comp_cylinder]") {
  auto &sim = context->simulation();
  static const Omega_h::Real mesh_measure = 57.188949391;
  static const Omega_h::Real compartment_measure = 28.5945;

  auto &mesh = sim.mesh_impl();
  const auto &measure = mesh.getMeasure();

  REQUIRE(measure.mesh_measure() == Approx(mesh_measure));

  auto left = mesh.getCompID("Left");
  auto right = mesh.getCompID("Right");
  REQUIRE(measure.mesh_measure(left) == Approx(compartment_measure));
  REQUIRE(measure.mesh_measure(right) == Approx(compartment_measure));

  std::array<Omega_h::Real, 3> local_measures{
      measure.rank_measure(),
      measure.rank_measure(left),
      measure.rank_measure(right),
  };
  std::array<Omega_h::Real, 3> measures{};

  MPI_Reduce(local_measures.data(),
             measures.data(),
             local_measures.size(),
             MPI_REAL8,
             MPI_SUM,
             0,
             sim.getMesh().comm_impl());
  if (sim.comm_rank == 0) {
    REQUIRE(measures[0] == Approx(mesh_measure));
    REQUIRE(measures[1] == Approx(compartment_measure));
    REQUIRE(measures[2] == Approx(compartment_measure));
  }
}

TEST_CASE("SignedCompactTypeId", "[type_id]") {
  using TId = steps::util::CompactTypeId<My2Enum, 1, std::int32_t>;
  {
    REQUIRE(TId::fits(My2Enum::VAL1));
    REQUIRE(TId::fits(My2Enum::VAL2));
    REQUIRE(TId::fits(static_cast<My2Enum>(0)));
    REQUIRE(TId::fits(static_cast<My2Enum>(1)));
    REQUIRE_FALSE(TId::fits(static_cast<My2Enum>(2)));
    REQUIRE_FALSE(TId::fits(static_cast<My2Enum>(3)));
  }
  {
    REQUIRE(TId::fits(0));
    REQUIRE(TId::fits(std::numeric_limits<int>::max()));
    REQUIRE_FALSE(TId::fits(-1));
    REQUIRE_FALSE(TId::fits(-3));
    REQUIRE_FALSE(TId::fits(std::numeric_limits<int>::min()));
  }

  {
    TId id(My2Enum::VAL1, 0);
    REQUIRE(id.type() == My2Enum::VAL1);
    REQUIRE(id.id() == 0);
  }
  {
    auto max_supported_id = std::numeric_limits<int>::max();
    TId id(My2Enum::VAL2, max_supported_id);
    REQUIRE(id.type() == My2Enum::VAL2);
    REQUIRE(id.id() == max_supported_id);
  }
}

TEST_CASE("2BitsSignedCompactTypeId", "[type_id]") {
  using TId = steps::util::CompactTypeId<My3Enum, 2, std::int32_t>;
  {
    REQUIRE(TId::fits(My3Enum::VAL1));
    REQUIRE(TId::fits(My3Enum::VAL2));
    REQUIRE(TId::fits(My3Enum::VAL3));
    REQUIRE(TId::fits(static_cast<My3Enum>(0)));
    REQUIRE(TId::fits(static_cast<My3Enum>(1)));
    REQUIRE(TId::fits(static_cast<My3Enum>(2)));
    REQUIRE(TId::fits(static_cast<My3Enum>(3)));
    REQUIRE_FALSE(TId::fits(static_cast<My3Enum>(4)));
  }
  {
    REQUIRE(TId::fits(0));
    REQUIRE(TId::fits(std::numeric_limits<int>::max() >> 1));
    REQUIRE_FALSE(TId::fits(std::numeric_limits<int>::max()));
    REQUIRE_FALSE(TId::fits(-1));
    REQUIRE_FALSE(TId::fits(-3));
    REQUIRE_FALSE(TId::fits(std::numeric_limits<int>::min()));
  }

  {
    TId id(My3Enum::VAL1, 0);
    REQUIRE(id.type() == My3Enum::VAL1);
    REQUIRE(id.id() == 0);
  }
  {
    TId id(My3Enum::VAL2, 2);
    REQUIRE(id.type() == My3Enum::VAL2);
    REQUIRE(id.id() == 2);
  }
  {
    auto max_supported_id = std::numeric_limits<int>::max() >> 1;
    TId id(My3Enum::VAL3, max_supported_id);
    REQUIRE(id.type() == My3Enum::VAL3);
    REQUIRE(id.id() == max_supported_id);
  }
}

TEST_CASE("UnsignedCompactTypeId", "[type_id]") {
  using TId = steps::util::CompactTypeId<My2Enum, 1, unsigned>;
  {
    REQUIRE(TId::fits(My2Enum::VAL1));
    REQUIRE(TId::fits(My2Enum::VAL2));
    REQUIRE(TId::fits(static_cast<My2Enum>(0)));
    REQUIRE(TId::fits(static_cast<My2Enum>(1)));
    REQUIRE_FALSE(TId::fits(static_cast<My2Enum>(2)));
    REQUIRE_FALSE(TId::fits(static_cast<My2Enum>(3)));
  }
  {
    REQUIRE(TId::fits(0));
    REQUIRE(TId::fits(std::numeric_limits<int>::max()));
    REQUIRE_FALSE(TId::fits(static_cast<unsigned>(-1)));
    REQUIRE_FALSE(
        TId::fits(static_cast<unsigned>(std::numeric_limits<int>::max()) + 1));
    REQUIRE_FALSE(
        TId::fits(static_cast<unsigned>(std::numeric_limits<int>::max()) + 2));
  }

  {
    TId id(My2Enum::VAL1, 0);
    REQUIRE(id.type() == My2Enum::VAL1);
    REQUIRE(id.id() == 0);
  }
  {
    TId id(My2Enum::VAL2, 2);
    REQUIRE(id.type() == My2Enum::VAL2);
    REQUIRE(id.id() == 2);
  }
}

TEST_CASE("entities", "[mesh]") {
  const auto mesh_file = context->source_dir() / "test" / "mesh" / "box.msh";
  steps::dist::DistMesh mesh(context->library(), mesh_file.string());

  {
    // check the count of elements
    constexpr Omega_h::GO expected_num_elements = 4662;

    REQUIRE(mesh.total_num_elems() == expected_num_elements);
    const Omega_h::GO num_owned_elements = mesh.num_elems();
    Omega_h::GO total_num_owned_elements{};
    MPI_Allreduce(
        &num_owned_elements, &total_num_owned_elements, 1, MPI_INT64_T, MPI_SUM, mesh.comm_impl());
    REQUIRE(total_num_owned_elements == expected_num_elements);
  }

  {
      // check the count of boundaries
      constexpr Omega_h::GO expected_num_bounds = 10052;

      REQUIRE(mesh.total_num_bounds() == expected_num_bounds);
      const Omega_h::GO num_owned_bounds = mesh.num_bounds();
      Omega_h::GO total_num_owned_bounds{};
      MPI_Allreduce(
          &num_owned_bounds, &total_num_owned_bounds, 1, MPI_INT64_T, MPI_SUM, mesh.comm_impl());
      REQUIRE(total_num_owned_bounds == expected_num_bounds);
  }

  {
      // check the count of vertices
      constexpr Omega_h::GO expected_num_verts = 1149;

      REQUIRE(mesh.total_num_verts() == expected_num_verts);
      const Omega_h::GO num_owned_verts = mesh.num_verts();
      Omega_h::GO total_num_owned_verts{};
      MPI_Allreduce(
          &num_owned_verts, &total_num_owned_verts, 1, MPI_INT64_T, MPI_SUM, mesh.comm_impl());
      REQUIRE(total_num_owned_verts == expected_num_verts);
  }
}

TEST_CASE("KProcID", "[type_id]") {
  using steps::dist::kproc::KProcID;
  using steps::dist::kproc::KProcType;

  {
    KProcID id(KProcType::Reac, 42);
    REQUIRE(id.type() == KProcType::Reac);
    REQUIRE(id.id() == 42);
  }
  {
    KProcID id(KProcType::Diff, 42);
    REQUIRE(id.type() == KProcType::Diff);
    REQUIRE(id.id() == 42);
  }
  {
    KProcID id(KProcType::SReac, 42);
    REQUIRE(id.type() == KProcType::SReac);
    REQUIRE(id.id() == 42);
  }
  static_assert(sizeof(KProcID) == sizeof(unsigned),
                "Unexpected datastructure size");
}

TEST_CASE("optional_num", "[type_id]") {
  {
    using optint = steps::util::OptionalNum<int>;
    static_assert(std::is_same<int, optint::value_type>::value,
                  "Missing type definition 'value_type'");
    static_assert(optint::none_value() ==
                      std::numeric_limits<optint::value_type>::max(),
                  "Cannot call 'none_value' at compile time");

    {
      optint val{};
      REQUIRE(!val); // operator bool
      REQUIRE_THROWS(val.value());
      REQUIRE(val.value_or(42) == 42);
    }
    {
      optint val;
      REQUIRE(!val); // operator bool
      REQUIRE_THROWS(val.value());
      REQUIRE(val.value_or(42) == 42);
    }
    {
      optint val(42);
      REQUIRE(*val == 42);
      REQUIRE(val.has_value());
      REQUIRE(val); // operator bool
      REQUIRE_NOTHROW(val.value() == 42);

      val = steps::util::nothing;
      REQUIRE(!val); // operator bool
      REQUIRE(!val.has_value());
      REQUIRE_THROWS(val.value());
    }
    {
      optint val(steps::util::nothing);
      REQUIRE(!val.has_value());
      REQUIRE(!val); // operator bool
      REQUIRE_THROWS(val.value());
      REQUIRE(val.value_or(42) == 42);
    }
    {
      optint val;
      auto val42 = *(val = 42);
      REQUIRE(*val == 42);
      REQUIRE(val42 == 42);
      REQUIRE(val); // operator bool
      REQUIRE_NOTHROW(val.value() == 42);
    }
  }
  {
    using optint = steps::util::OptionalNum<int, 42>;
    optint val(steps::util::nothing);
    val = std::numeric_limits<optint::value_type>::max();
    REQUIRE(val);
    REQUIRE(val.has_value());
    REQUIRE_NOTHROW(val.value() ==
                    std::numeric_limits<optint::value_type>::max());
  }
}

static int run_catch_test_or_tag(Catch::Session &session,
                                 const std::string &name) {
  session.configData().testsOrTags.push_back(name);
  // force reset of configuration to forget previous runs
  session.useConfigData(session.configData());
  auto failed_tests = session.run();
  session.configData().testsOrTags.pop_back();
  return failed_tests;
}

template <typename Scenario>
int run_simulation_catch_test(Catch::Session &session,
                              const steps::dist::ScenarioInput &input,
                              const std::string &catch2_tag) {
  std::mt19937 rng;
  using simulation_type = steps::dist::OmegaHSimulation<>;

  // create mesh and opsplit simulation objects
  simulation_type::mesh_type mesh(context->library(), input.mesh_file,
                                  input.scale);

  std::ofstream outfile;
  if (!input.logfile.empty()) {
    outfile.open(input.logfile, std::ofstream::out);
  }
  std::ostream &outstream = outfile.is_open() ? outfile : std::clog;
  simulation_type simulation(input, mesh, rng, outstream);

  // load the simulation with 2 compartments
  Scenario scenario(input);
  scenario.initialize(simulation);

  // assign global simulation object as long as `fixture` instance exists
  const auto fixture = context->simulation(simulation);
  // execute only involved tests
  return run_catch_test_or_tag(session, catch2_tag);
}
static int run_two_comp_cylinder_catch_tests(Catch::Session &session) {
  const auto mesh_file =
      context->source_dir() / "test" / "mesh" / "two_comp_cyl.msh";
  return run_simulation_catch_test<steps::dist::MultipleCompartment>(
      session, {mesh_file.string(), 1.0}, "[two_comp_cylinder]");
}

static int run_caburst_catch_tests(Catch::Session &session) {
  const auto mesh_file = context->source_dir() / "test" / "mesh" / "CaBurst" /
                         "branch_labeledV4.msh";
  return run_simulation_catch_test<steps::dist::CaBurstBackground>(
      session, {mesh_file.string()}, "[ca_burst_background]");
}

int main(int argc, char *argv[]) {
  auto ctx = std::make_unique<Context>(&argc, &argv);
  // set global variable
  context = ctx.get();

  // top catch2 object
  Catch::Session session;

  // execute tests
  int failed_tests{};
  failed_tests += run_two_comp_cylinder_catch_tests(session);
  failed_tests += run_caburst_catch_tests(session);
  failed_tests += run_catch_test_or_tag(session, "[type_id]");
  failed_tests += run_catch_test_or_tag(session, "[mesh]");

  return failed_tests;
}
