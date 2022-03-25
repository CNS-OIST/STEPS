#pragma once

#include "mpi/dist/test/scenario.hpp"

namespace steps {
namespace dist {

class Rallpack3 : public Scenario<std::mt19937> {
public:
  explicit Rallpack3(const ScenarioInput &input);

private:
  std::unique_ptr<Statedef>
  createStatedef(const simulation_t &simulation) const override;
  void register_compartments(DistMesh &mesh) const override;
  void register_patches(DistMesh &mesh) const;
  void fill_compartments(simulation_t &simulation) const override;
  void fill_patches(simulation_t &simulation) const override;
  void run_simulation_impl(simulation_t &simulation) override;

  model::compartment_label cyto_tag{};
  model::patch_label memb_tag{};
};

} // namespace dist
} // namespace steps
