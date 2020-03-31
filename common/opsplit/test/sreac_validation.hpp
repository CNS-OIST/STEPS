#pragma once

#include "opsplit/test/scenario.hpp"
#include "vec3.hpp"

namespace zee {

class SReacValidationTest: public Scenario<std::mt19937> {
  public:
    explicit SReacValidationTest(const ScenarioInput& input);

  private:
    std::unique_ptr<Statedef> createStatedef() const override;
    void register_compartments(DistMesh& mesh) const override;
    void fill_compartments(simulation_t& simulation) const override;
    void run_simulation_impl(simulation_t& simulation) override;
    int check_and_log_results_impl(const simulation_t& simulation) const override;

    const std::vector<PetscScalar> tpnts;

    Vec3<PetscScalar> res_m_so2d;
};

}  // namespace zee
