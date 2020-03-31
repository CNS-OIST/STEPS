#include <map>
#include <vector>

#include <hadoken/format/format.hpp>
#include <petsctime.h>

#include "Ca_burst_background.hpp"
#include "opsplit/compdef.hpp"
#include "opsplit/diffdef.hpp"
#include "opsplit/mesh.hpp"
#include "opsplit/patchdef.hpp"
#include "opsplit/reacdef.hpp"
#include "opsplit/sreacdef.hpp"
#include "opsplit/statedef.hpp"

#include "opsplit/test/simulation.hpp"

namespace zee {

using hadoken::scat;

static const PetscScalar SIM_TIME = 30.0e-5;

CaBurstBackground::CaBurstBackground(const ScenarioInput& t_input)
    : Scenario("Ca Burst Background", "Background buffers only", t_input) {}


std::unique_ptr<Statedef> CaBurstBackground::createStatedef() const {
    CaBurstSimdef simdef;
    return std::move(simdef.getStatedef());
}

void CaBurstBackground::register_compartments(DistMesh& mesh) const {
    mesh.addComp("__MESH__", cyto_tag);
}

void CaBurstBackground::register_patches(DistMesh& /*mesh*/) const {
    // mesh.addPatch("__MESH_BOUNDARY__", memb_tag);
}

void CaBurstBackground::fill_compartments(simulation_t& simulation) const {
    CompartmentConc mg("__MESH__", "Mg", CaBurstSimdef::Mg_conc());
    CompartmentConc iCBsf("__MESH__", "iCBsf", CaBurstSimdef::iCBsf_conc());
    CompartmentConc iCBCaf("__MESH__", "iCBCaf", CaBurstSimdef::iCBCaf_conc());
    CompartmentConc iCBsCa("__MESH__", "iCBsCa", CaBurstSimdef::iCBsCa_conc());
    CompartmentConc iCBCaCa("__MESH__", "iCBCaCa", CaBurstSimdef::iCBCaCa_conc());
    CompartmentConc CBsf("__MESH__", "CBsf", CaBurstSimdef::CBsf_conc());
    CompartmentConc CBCaf("__MESH__", "CBCaf", CaBurstSimdef::CBCaf_conc());
    CompartmentConc CBsCa("__MESH__", "CBsCa", CaBurstSimdef::CBsCa_conc());
    CompartmentConc CBCaCa("__MESH__", "CBCaCa", CaBurstSimdef::CBCaCa_conc());
    CompartmentConc PV("__MESH__", "PV", CaBurstSimdef::PV_conc());
    CompartmentConc PVCa("__MESH__", "PVCa", CaBurstSimdef::PVCa_conc());
    CompartmentConc PVMg("__MESH__", "PVMg", CaBurstSimdef::PVMg_conc());
    simulation.setCompConc(
        {mg, iCBsf, iCBCaf, iCBsCa, iCBCaCa, CBsf, CBCaf, CBsCa, CBCaCa, PV, PVCa, PVMg});
}

void CaBurstBackground::fill_patches(simulation_t& simulation) const {
    PetscScalar patch_area = simulation.getMesh().getTotalPatchArea("__MESH_BOUNDARY__");
    PetscScalar pumpnbs = 6.022141e12 * patch_area;
    PatchCount pump("__MESH_BOUNDARY__", "Pump", round(pumpnbs));
    PatchCount ca_pump("__MESH_BOUNDARY__", "CaPump", 0.0);
    simulation.setPatchCount({pump, ca_pump});
}

void CaBurstBackground::run_simulation_impl(simulation_t& simulation) {
    simulation.run(SIM_TIME);
}

}  // namespace zee
