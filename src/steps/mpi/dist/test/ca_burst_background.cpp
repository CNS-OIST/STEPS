#include "mpi/dist/test/ca_burst_background.hpp"

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

static const osh::Real SIM_TIME = 30.0e-5;

CaBurstBackground::CaBurstBackground(const ScenarioInput &t_input)
    : Scenario("Ca Burst Background", "Background buffers only", t_input) {}

std::unique_ptr<Statedef>
CaBurstBackground::createStatedef(const simulation_t &simulation) const {
  steps::model::Model model;
  CaBurstSimdef simdef(model, simulation.getMesh());
  return std::move(simdef.getStatedef());
}

void CaBurstBackground::register_compartments(DistMesh &mesh) const {
  mesh.addComp("__MESH__", cyto_tag);
}

void CaBurstBackground::register_patches(DistMesh & /*mesh*/) const {
  // mesh.addPatch("__MESH_BOUNDARY__", memb_tag);
}

void CaBurstBackground::fill_compartments(simulation_t &simulation) const {
  CompartmentConc mg("Mg", CaBurstSimdef::Mg_conc());
  CompartmentConc iCBsf("iCBsf", CaBurstSimdef::iCBsf_conc());
  CompartmentConc iCBCaf("iCBCaf", CaBurstSimdef::iCBCaf_conc());
  CompartmentConc iCBsCa("iCBsCa", CaBurstSimdef::iCBsCa_conc());
  CompartmentConc iCBCaCa("iCBCaCa", CaBurstSimdef::iCBCaCa_conc());
  CompartmentConc CBsf("CBsf", CaBurstSimdef::CBsf_conc());
  CompartmentConc CBCaf("CBCaf", CaBurstSimdef::CBCaf_conc());
  CompartmentConc CBsCa("CBsCa", CaBurstSimdef::CBsCa_conc());
  CompartmentConc CBCaCa("CBCaCa", CaBurstSimdef::CBCaCa_conc());
  CompartmentConc PV("PV", CaBurstSimdef::PV_conc());
  CompartmentConc PVCa("PVCa", CaBurstSimdef::PVCa_conc());
  CompartmentConc PVMg("PVMg", CaBurstSimdef::PVMg_conc());

  simulation.setCompConc(
      "__MESH__", {mg, iCBsf, iCBCaf, iCBsCa, iCBCaCa, CBsf, CBCaf, CBsCa, CBCaCa, PV, PVCa, PVMg});
}

void CaBurstBackground::fill_patches(simulation_t &simulation) const {

  double spiny_area =
      simulation.getMesh().total_measure(model::patch_id("spiny.__BOUNDARY__"));
  assert(spiny_area > 0.0);
  double smooth_area = simulation.getMesh().total_measure(
      model::patch_id("smooth.__BOUNDARY__"));
  ;
  assert(smooth_area > 0.0);

  PatchCount pump_smooth(
      "smooth.__BOUNDARY__", "Pump",
      std::round(CaBurstSimdef::pumpnbs_per_area() * smooth_area));
  PatchCount ca_pump_smooth("smooth.__BOUNDARY__", "CaPump", 0.0);

  PatchCount pump_spiny(
      "spiny.__BOUNDARY__", "Pump",
      std::round(CaBurstSimdef::pumpnbs_per_area() * spiny_area));
  PatchCount ca_pump_spiny("spiny.__BOUNDARY__", "CaPump", 0.0);

  simulation.setPatchCount({pump_smooth, pump_spiny, ca_pump_smooth, ca_pump_spiny});
}

void CaBurstBackground::run_simulation_impl(simulation_t &simulation) {
  simulation.run(SIM_TIME);
}

} // namespace dist
} // namespace steps
