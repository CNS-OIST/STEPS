#pragma once

#include "mpi/dist/test/scenario.hpp"

namespace steps {
namespace dist {

class CaBurstIntegrationTest: public Scenario<std::mt19937> {
  public:
    explicit CaBurstIntegrationTest(const ScenarioInput& input);

  private:
    std::unique_ptr<Statedef> createStatedef(const simulation_t& simulation) const override;
    void register_compartments(DistMesh& mesh) const override;
    void fill_compartments(simulation_t& simulation) const override;
    void fill_patches(simulation_t& simulation) const override;
    void run_simulation_impl(simulation_t& simulation) override;
    int check_and_log_results_impl(simulation_t& simulation) const override;

    std::vector<double> A_count;
    std::vector<double> B_count;
    std::vector<double> C_count;
    std::vector<double> C_surf_count;
    std::vector<double> V_min;
    std::vector<double> V_max;
    std::vector<double> i_ghk;
    std::vector<double> t;
};

}  // namespace dist
}  // namespace steps
