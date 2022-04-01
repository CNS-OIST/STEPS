#pragma once

#include "mpi/dist/test/scenario.hpp"

namespace steps {
namespace dist {

class SReacUnitTest : public Scenario<std::mt19937> {
public:
  explicit SReacUnitTest(const ScenarioInput &input);

private:
  std::unique_ptr<Statedef>
  createStatedef(const simulation_t &simulation) const override;
  void register_compartments(DistMesh &mesh) const override;
  void fill_compartments(simulation_t &simulation) const override;
  void run_simulation_impl(simulation_t &simulation) override;
  int check_and_log_results_impl(simulation_t &simulation) const override;
};

} // namespace dist
} // namespace steps
