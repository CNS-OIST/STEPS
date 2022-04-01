#pragma once

#include "mpi/dist/test/scenario.hpp"
#include "util/vec3.hpp"

namespace steps {
namespace dist {

class ValidationTest : public Scenario<std::mt19937> {
public:
  explicit ValidationTest(const ScenarioInput &input);

private:
  std::unique_ptr<Statedef>
  createStatedef(const simulation_t &simulation) const override;
  void register_compartments(DistMesh &mesh) const override;
  void fill_compartments(simulation_t &simulation) const override;
  void run_simulation_impl(simulation_t &simulation) override;
  int check_and_log_results_impl(simulation_t &simulation) const override;

  const std::vector<osh::Real> tpnts;

  util::Vec3<osh::Real> res_m_foi;
  util::Vec3<osh::Real> res_m_for;
  util::Vec3<osh::Real> res_m_soA2;
  util::Vec3<osh::Real> res_m_soAA;
  util::Vec3<osh::Real> res_m_soAB;
  util::Vec3<osh::Real> res_m_toA3;
  util::Vec3<osh::Real> res_m_toA2B;
  util::Vec3<osh::Real> res_m_so2d;
};

} // namespace dist
} // namespace steps
