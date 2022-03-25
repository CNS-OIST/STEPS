#pragma once

#include "mpi/dist/test/scenario.hpp"

namespace steps {
namespace dist {

class MultipleCompartmentDiff : public Scenario<std::mt19937> {
public:
  explicit MultipleCompartmentDiff(const ScenarioInput &input);

private:
  std::unique_ptr<Statedef>
  createStatedef(const simulation_t &simulation) const override;
  void register_compartments(DistMesh &mesh) const override;
  void fill_compartments(simulation_t &simulation) const override;
  void run_simulation_impl(simulation_t &simulation) override;
  int check_and_log_results_impl(simulation_t &simulation) const override;

  osh::Real C_i_init, C_o_init;
  osh::Real C_i_nodiff, C_o_nodiff;
  osh::Real C_i_bdiff, C_o_bdiff;
};

} // namespace dist
} // namespace steps
