#pragma once

#include "mpi/dist/test/scenario.hpp"
#include "util/vec3.hpp"

namespace steps {
namespace dist {

class SReacValidationTest : public Scenario<std::mt19937> {
public:
  explicit SReacValidationTest(const ScenarioInput &input);

private:
  std::unique_ptr<Statedef>
  createStatedef(const simulation_t &simulation) const override;
  void register_compartments(DistMesh &mesh) const override;
  void fill_compartments(simulation_t &simulation) const override;
  void run_simulation_impl(simulation_t &simulation) override;
  int check_and_log_results_impl(simulation_t &simulation) const override;

  const std::vector<osh::Real> tpnts;

  util::Vec3<osh::Real> res_m_so2d;

  /// Sampling time-step
  static constexpr osh::Real sampling_dt() { return 0.01; }
  /// Simulation  end time
  static constexpr osh::Real sim_end_time() { return 1.1; }
};

} // namespace dist
} // namespace steps
