#pragma once

#include "mpi/dist/test/scenario.hpp"

namespace steps {
namespace dist {

class MultipleCompartment : public Scenario<std::mt19937> {
public:
  explicit MultipleCompartment(const ScenarioInput &input);

private:
  std::unique_ptr<Statedef>
  createStatedef(const simulation_t &simulation) const override;
  void register_compartments(DistMesh &mesh) const override;
  void fill_compartments(simulation_t &simulation) const override;
  void run_simulation_impl(simulation_t &simulation) override;
  void log_num_molecules(const simulation_t &simulation) const;
};

} // namespace dist
} // namespace steps
