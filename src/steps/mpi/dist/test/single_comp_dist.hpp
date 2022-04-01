#pragma once

#include "mpi/dist/test/scenario.hpp"

namespace steps {
namespace dist {

class SingleCompDist : public Scenario<std::mt19937> {

public:
  explicit SingleCompDist(const ScenarioInput &input);

private:
  std::unique_ptr<Statedef>
  createStatedef(const simulation_t &simulation) const override;
  void register_compartments(DistMesh &mesh) const override;
  void fill_compartments(simulation_t &simulation) const override;
  void run_simulation_impl(simulation_t &) override {}
  int check_and_log_results_impl(simulation_t &simulation) const override;
};

} // namespace dist
} // namespace steps
