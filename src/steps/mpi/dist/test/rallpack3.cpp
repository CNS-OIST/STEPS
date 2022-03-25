#include "rallpack3.hpp"

#include <iostream>
#include <map>

#include "geom/dist/distmesh.hpp"
#include "model/model.hpp"
#include "mpi/dist/test/simulation.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"

namespace steps {
namespace dist {

Rallpack3::Rallpack3(const ScenarioInput &t_input)
    : Scenario("Rallpack ", "Background buffers only", t_input) {}

std::unique_ptr<Statedef>
Rallpack3::createStatedef(const simulation_t &simulation) const {
  steps::model::Model model;
  Rallpack3Simdef simdef(model, simulation.getMesh());
  return std::move(simdef.getStatedef());
}

void Rallpack3::register_compartments(DistMesh &mesh) const {
  mesh.addComp("__MESH__", cyto_tag);
}

void Rallpack3::register_patches(DistMesh & /*mesh*/) const {
  /*mesh.addPatch("__MESH_BOUNDARY__", memb_tag);*/
}

void Rallpack3::fill_compartments(simulation_t & /*simulation*/) const {}

void Rallpack3::fill_patches(simulation_t &simulation) const {
  osh::Real surfarea_cyl{1.0e-3 * 1.0e-6 * 3.141592653589};
  simulation.setPatchCount(
      "memb", "Na_m0h1", (Rallpack3Simdef::Na_ro() * surfarea_cyl * Rallpack3Simdef::NA_FACS()[0]));
  simulation.setPatchCount(
      "memb", "Na_m1h1", (Rallpack3Simdef::Na_ro() * surfarea_cyl * Rallpack3Simdef::NA_FACS()[1]));
  simulation.setPatchCount(
      "memb", "Na_m2h1", (Rallpack3Simdef::Na_ro() * surfarea_cyl * Rallpack3Simdef::NA_FACS()[2]));
  simulation.setPatchCount(
      "memb", "Na_m3h1", (Rallpack3Simdef::Na_ro() * surfarea_cyl * Rallpack3Simdef::NA_FACS()[3]));
  simulation.setPatchCount(
      "memb", "Na_m0h0", (Rallpack3Simdef::Na_ro() * surfarea_cyl * Rallpack3Simdef::NA_FACS()[4]));
  simulation.setPatchCount(
      "memb", "Na_m1h0", (Rallpack3Simdef::Na_ro() * surfarea_cyl * Rallpack3Simdef::NA_FACS()[5]));
  simulation.setPatchCount(
      "memb", "Na_m2h0", (Rallpack3Simdef::Na_ro() * surfarea_cyl * Rallpack3Simdef::NA_FACS()[6]));
  simulation.setPatchCount(
      "memb", "Na_m3h0", (Rallpack3Simdef::Na_ro() * surfarea_cyl * Rallpack3Simdef::NA_FACS()[7]));
  simulation.setPatchCount("memb",
                           "K_n0",
                           (Rallpack3Simdef::K_ro() * surfarea_cyl * Rallpack3Simdef::K_FACS()[0]));
  simulation.setPatchCount("memb",
                           "K_n1",
                           (Rallpack3Simdef::K_ro() * surfarea_cyl * Rallpack3Simdef::K_FACS()[1]));
  simulation.setPatchCount("memb",
                           "K_n2",
                           (Rallpack3Simdef::K_ro() * surfarea_cyl * Rallpack3Simdef::K_FACS()[2]));
  simulation.setPatchCount("memb",
                           "K_n3",
                           (Rallpack3Simdef::K_ro() * surfarea_cyl * Rallpack3Simdef::K_FACS()[3]));
  simulation.setPatchCount("memb",
                           "K_n4",
                           (Rallpack3Simdef::K_ro() * surfarea_cyl * Rallpack3Simdef::K_FACS()[4]));
}

void Rallpack3::run_simulation_impl(simulation_t &simulation) {
  auto delta = 5e-6;
  simulation.setPotential(-65.0e-3);

  std::stringstream s;
  s << "t V_max_on_verts_at_z_min V_max_on_verts_at_z_max "
       "V_max_on_tris_at_z_min V_max_on_tris_at_z_max\n";

  s << 0.0 << ' ' << simulation.getMaxPotentialOnVertices("z_min") << ' '
    << simulation.getMaxPotentialOnVertices("z_max") << '\n';
  simulation.log_once(s.str());
  for (int k = 1; k * delta <= input.end_time; ++k) {
    const auto t = k * delta;
    simulation.run(t);

    // results
    s.str("");
    s << k * delta << " " << simulation.getMaxPotentialOnVertices("z_min")
      << ' ' << simulation.getMaxPotentialOnVertices("z_max") << '\n';
    simulation.log_once(s.str());

    simulation.log_progress(t, input.end_time);
  }
}

} // namespace dist
} // namespace steps
