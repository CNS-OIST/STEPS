#pragma once

#include "mpi/dist/test/scenario.hpp"

namespace steps {
namespace dist {

class DiffusionOnly : public Scenario<std::mt19937> {
public:
  explicit DiffusionOnly(const ScenarioInput &input);

private:
  std::unique_ptr<Statedef>
  createStatedef(const simulation_t &simulation) const override;
  void register_compartments(DistMesh &mesh) const override;
  void fill_compartments(simulation_t &simulation) const override;
  void run_simulation_impl(simulation_t &simulation) override;

  model::compartment_label comp_tag{};
  osh::Real dcst{1e-8};
  osh::Real NMols{100000};
};

} // namespace dist
} // namespace steps
