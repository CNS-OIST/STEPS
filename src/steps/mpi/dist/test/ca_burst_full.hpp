#pragma once

#include <vector>

#include "mpi/dist/test/scenario.hpp"

namespace steps {
namespace dist {

class CaBurstFull: public Scenario<std::mt19937> {
  public:
    explicit CaBurstFull(const ScenarioInput& input);

  private:
    std::unique_ptr<Statedef> createStatedef(const simulation_t& simulation) const override;
    void register_compartments(DistMesh& mesh) const override;
    void register_patches(DistMesh& mesh) const;
    void fill_compartments(simulation_t& simulation) const override;
    void fill_patches(simulation_t& simulation) const override;
    void run_simulation_impl(simulation_t& simulation) override;

    model::compartment_label cyto_tag{};
    model::patch_label smooth_tag{};
    model::patch_label spiny_tag{};
    uint n_timepoints{3001};
    double time_point_interval{2.0e-5};
};

}  // namespace dist
}  // namespace steps
