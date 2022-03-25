#pragma once

#include "mpi/dist/test/scenario.hpp"

namespace steps {
namespace dist {

class GHKCurrentUnitTest : public Scenario<std::mt19937> {
public:
  explicit GHKCurrentUnitTest(const ScenarioInput &input);

private:
  std::unique_ptr<Statedef>
  createStatedef(const simulation_t &simulation) const override;
  void register_compartments(DistMesh &mesh) const override;
  void fill_compartments(simulation_t &simulation) const override;
  void run_simulation_impl(simulation_t &simulation) override;
  int check_and_log_results_impl(simulation_t &simulation) const override;
  
  double conc_i, conc_o;
  std::vector<double> res;
  static constexpr double equilibrium_potential = -65.0e-3;
};

} // namespace dist
} // namespace steps
