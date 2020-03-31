#pragma once

#include <opsplit/test/scenario.hpp>

namespace zee {

class SimpleModel2: public Scenario<std::mt19937> {
  public:
    explicit SimpleModel2(const ScenarioInput& input);

  private:
    std::unique_ptr<Statedef> createStatedef() const override;
    void register_compartments(DistMesh& mesh) const override;
    void fill_compartments(simulation_t& simulation) const override;
    void run_simulation_impl(simulation_t& simulation) override;
    int check_and_log_results_impl(const simulation_t& simulation) const override;
};

}  // namespace zee
