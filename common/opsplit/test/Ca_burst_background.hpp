#pragma once

#include <opsplit/test/scenario.hpp>

namespace zee {

class CaBurstBackground: public Scenario<std::mt19937> {
  public:
    explicit CaBurstBackground(const ScenarioInput& input);

  private:
    std::unique_ptr<Statedef> createStatedef() const override;
    void register_compartments(DistMesh& mesh) const override;
    void register_patches(DistMesh& mesh) const;
    void fill_compartments(simulation_t& simulation) const override;
    void fill_patches(simulation_t& simulation) const;
    void run_simulation_impl(simulation_t& simulation) override;

    model::compartment_label cyto_tag{};
    model::patch_label memb_tag{};
};

}  // namespace zee
