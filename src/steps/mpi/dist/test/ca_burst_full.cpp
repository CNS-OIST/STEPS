#include "mpi/dist/test/ca_burst_full.hpp"

#include <iomanip>
#include <iostream>
#include <map>

#include "geom/dist/distmesh.hpp"
#include "model/model.hpp"
#include "mpi/dist/test/glut_pulse_mvr.hpp"
#include "mpi/dist/test/simulation.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"

namespace steps {
namespace dist {

CaBurstFull::CaBurstFull(const ScenarioInput& t_input)
    : Scenario("Ca Burst Full", "Full Ca burst model", t_input) {}

std::unique_ptr<Statedef> CaBurstFull::createStatedef(const simulation_t& simulation) const {
    steps::model::Model model;
    CaBurstFullSimdef simdef(model, simulation.getMesh());
    return std::move(simdef.getStatedef());
}

void CaBurstFull::register_compartments(DistMesh& mesh) const {
    mesh.addComp("__MESH__", cyto_tag);
}

void CaBurstFull::register_patches(DistMesh& /* mesh */) const {
    /*
    mesh.addPatch("smooth", smooth_tag);
    mesh.addPatch("spiny", spiny_tag);
    */
}

void CaBurstFull::fill_compartments(simulation_t& simulation) const {
    CompartmentConc ca("Ca", CaBurstFullSimdef::Ca_iconc());
    CompartmentConc mg("Mg", CaBurstFullSimdef::Mg_conc());
    CompartmentConc iCBsf("iCBsf", CaBurstFullSimdef::iCBsf_conc());
    CompartmentConc iCBCaf("iCBCaf", CaBurstFullSimdef::iCBCaf_conc());
    CompartmentConc iCBsCa("iCBsCa", CaBurstFullSimdef::iCBsCa_conc());
    CompartmentConc iCBCaCa("iCBCaCa", CaBurstFullSimdef::iCBCaCa_conc());
    CompartmentConc CBsf("CBsf", CaBurstFullSimdef::CBsf_conc());
    CompartmentConc CBCaf("CBCaf", CaBurstFullSimdef::CBCaf_conc());
    CompartmentConc CBsCa("CBsCa", CaBurstFullSimdef::CBsCa_conc());
    CompartmentConc CBCaCa("CBCaCa", CaBurstFullSimdef::CBCaCa_conc());
    CompartmentConc PV("PV", CaBurstFullSimdef::PV_conc());
    CompartmentConc PVCa("PVCa", CaBurstFullSimdef::PVCa_conc());
    CompartmentConc PVMg("PVMg", CaBurstFullSimdef::PVMg_conc());
    simulation.setCompConc(
        "__MESH__",
        {ca, mg, iCBsf, iCBCaf, iCBsCa, iCBCaCa, CBsf, CBCaf, CBsCa, CBCaCa, PV, PVCa, PVMg});
}

void CaBurstFull::fill_patches(simulation_t& simulation) const {
    /// TODO
    double spiny_area = simulation.getMesh().total_measure(model::patch_id("spiny.__BOUNDARY__"));
    double smooth_area = simulation.getMesh().total_measure(model::patch_id("smooth.__BOUNDARY__")); ;

    osh::Real pumpnbs_smooth =
        CaBurstFullSimdef::pumpnbs_per_area() * smooth_area;
    PatchCount pump_smooth("smooth.__BOUNDARY__", "Pump", std::round(pumpnbs_smooth));

    PatchCount cap_m0_smooth("smooth.__BOUNDARY__",
                             "CaP_m0",
                             std::round(CaBurstFullSimdef::CaP_ro() * smooth_area *
                                        CaBurstFullSimdef::CaP_m0_p()));
    PatchCount cap_m1_smooth("smooth.__BOUNDARY__",
                             "CaP_m1",
                             std::round(CaBurstFullSimdef::CaP_ro() * smooth_area *
                                        CaBurstFullSimdef::CaP_m1_p()));
    PatchCount cap_m2_smooth("smooth.__BOUNDARY__",
                             "CaP_m2",
                             std::round(CaBurstFullSimdef::CaP_ro() * smooth_area *
                                        CaBurstFullSimdef::CaP_m2_p()));
    PatchCount cap_m3_smooth("smooth.__BOUNDARY__",
                             "CaP_m3",
                             std::round(CaBurstFullSimdef::CaP_ro() * smooth_area *
                                        CaBurstFullSimdef::CaP_m3_p()));
    PatchCount bk_c0_smooth("smooth.__BOUNDARY__",
                            "BK_C0",
                            std::round(CaBurstFullSimdef::BK_ro() * smooth_area *
                                       CaBurstFullSimdef::BK_C0_p()));
    PatchCount bk_c1_smooth("smooth.__BOUNDARY__",
                            "BK_C1",
                            std::round(CaBurstFullSimdef::BK_ro() * smooth_area *
                                       CaBurstFullSimdef::BK_C1_p()));
    PatchCount bk_c2_smooth("smooth.__BOUNDARY__",
                            "BK_C2",
                            std::round(CaBurstFullSimdef::BK_ro() * smooth_area *
                                       CaBurstFullSimdef::BK_C2_p()));
    PatchCount bk_c3_smooth("smooth.__BOUNDARY__",
                            "BK_C3",
                            std::round(CaBurstFullSimdef::BK_ro() * smooth_area *
                                       CaBurstFullSimdef::BK_C3_p()));
    PatchCount bk_c4_smooth("smooth.__BOUNDARY__",
                            "BK_C4",
                            std::round(CaBurstFullSimdef::BK_ro() * smooth_area *
                                       CaBurstFullSimdef::BK_C4_p()));
    PatchCount bk_o0_smooth("smooth.__BOUNDARY__",
                            "BK_O0",
                            std::round(CaBurstFullSimdef::BK_ro() * smooth_area *
                                       CaBurstFullSimdef::BK_O0_p()));
    PatchCount bk_o1_smooth("smooth.__BOUNDARY__",
                            "BK_O1",
                            std::round(CaBurstFullSimdef::BK_ro() * smooth_area *
                                       CaBurstFullSimdef::BK_O1_p()));
    PatchCount bk_o2_smooth("smooth.__BOUNDARY__",
                            "BK_O2",
                            std::round(CaBurstFullSimdef::BK_ro() * smooth_area *
                                       CaBurstFullSimdef::BK_O2_p()));
    PatchCount bk_o3_smooth("smooth.__BOUNDARY__",
                            "BK_O3",
                            std::round(CaBurstFullSimdef::BK_ro() * smooth_area *
                                       CaBurstFullSimdef::BK_O3_p()));
    PatchCount bk_o4_smooth("smooth.__BOUNDARY__",
                            "BK_O4",
                            std::round(CaBurstFullSimdef::BK_ro() * smooth_area *
                                       CaBurstFullSimdef::BK_O4_p()));

    PatchCount sk_c1_smooth("smooth.__BOUNDARY__",
                            "SK_C1",
                            std::round(CaBurstFullSimdef::SK_ro() * smooth_area *
                                       CaBurstFullSimdef::SK_C1_p()));
    PatchCount sk_c2_smooth("smooth.__BOUNDARY__",
                            "SK_C2",
                            std::round(CaBurstFullSimdef::SK_ro() * smooth_area *
                                       CaBurstFullSimdef::SK_C2_p()));
    PatchCount sk_c3_smooth("smooth.__BOUNDARY__",
                            "SK_C3",
                            std::round(CaBurstFullSimdef::SK_ro() * smooth_area *
                                       CaBurstFullSimdef::SK_C3_p()));
    PatchCount sk_c4_smooth("smooth.__BOUNDARY__",
                            "SK_C4",
                            std::round(CaBurstFullSimdef::SK_ro() * smooth_area *
                                       CaBurstFullSimdef::SK_C4_p()));

    PatchCount sk_o1_smooth("smooth.__BOUNDARY__",
                            "SK_O1",
                            std::round(CaBurstFullSimdef::SK_ro() * smooth_area *
                                       CaBurstFullSimdef::SK_O1_p()));
    PatchCount sk_o2_smooth("smooth.__BOUNDARY__",
                            "SK_O2",
                            std::round(CaBurstFullSimdef::SK_ro() * smooth_area *
                                       CaBurstFullSimdef::SK_O2_p()));

    PatchCount ampa_c_smooth("smooth.__BOUNDARY__", "AMPA_C", std::round(CaBurstFullSimdef::AMPA_receptors()));
    PatchCount leak_smooth("smooth.__BOUNDARY__",
                           "Leak",
                           std::round(CaBurstFullSimdef::L_ro_proximal() * smooth_area));

    osh::Real pumpnbs_spiny =
        CaBurstFullSimdef::pumpnbs_per_area() * spiny_area;
    PatchCount pump_spiny("spiny.__BOUNDARY__", "Pump", std::round(pumpnbs_spiny));

    PatchCount cap_m0_spiny("spiny.__BOUNDARY__",
                            "CaP_m0",
                            std::round(CaBurstFullSimdef::CaP_ro() * spiny_area *
                                       CaBurstFullSimdef::CaP_m0_p()));
    PatchCount cap_m1_spiny("spiny.__BOUNDARY__",
                            "CaP_m1",
                            std::round(CaBurstFullSimdef::CaP_ro() * spiny_area *
                                       CaBurstFullSimdef::CaP_m1_p()));
    PatchCount cap_m2_spiny("spiny.__BOUNDARY__",
                            "CaP_m2",
                            std::round(CaBurstFullSimdef::CaP_ro() * spiny_area *
                                       CaBurstFullSimdef::CaP_m2_p()));
    PatchCount cap_m3_spiny("spiny.__BOUNDARY__",
                            "CaP_m3",
                            std::round(CaBurstFullSimdef::CaP_ro() * spiny_area *
                                       CaBurstFullSimdef::CaP_m3_p()));
    PatchCount bk_c0_spiny("spiny.__BOUNDARY__",
                           "BK_C0",
                           std::round(CaBurstFullSimdef::BK_ro() * spiny_area *
                                      CaBurstFullSimdef::BK_C0_p()));
    PatchCount bk_c1_spiny("spiny.__BOUNDARY__",
                           "BK_C1",
                           std::round(CaBurstFullSimdef::BK_ro() * spiny_area *
                                      CaBurstFullSimdef::BK_C1_p()));
    PatchCount bk_c2_spiny("spiny.__BOUNDARY__",
                           "BK_C2",
                           std::round(CaBurstFullSimdef::BK_ro() * spiny_area *
                                      CaBurstFullSimdef::BK_C2_p()));
    PatchCount bk_c3_spiny("spiny.__BOUNDARY__",
                           "BK_C3",
                           std::round(CaBurstFullSimdef::BK_ro() * spiny_area *
                                      CaBurstFullSimdef::BK_C3_p()));
    PatchCount bk_c4_spiny("spiny.__BOUNDARY__",
                           "BK_C4",
                           std::round(CaBurstFullSimdef::BK_ro() * spiny_area *
                                      CaBurstFullSimdef::BK_C4_p()));
    PatchCount bk_o0_spiny("spiny.__BOUNDARY__",
                           "BK_O0",
                           std::round(CaBurstFullSimdef::BK_ro() * spiny_area *
                                      CaBurstFullSimdef::BK_O0_p()));
    PatchCount bk_o1_spiny("spiny.__BOUNDARY__",
                           "BK_O1",
                           std::round(CaBurstFullSimdef::BK_ro() * spiny_area *
                                      CaBurstFullSimdef::BK_O1_p()));
    PatchCount bk_o2_spiny("spiny.__BOUNDARY__",
                           "BK_O2",
                           std::round(CaBurstFullSimdef::BK_ro() * spiny_area *
                                      CaBurstFullSimdef::BK_O2_p()));
    PatchCount bk_o3_spiny("spiny.__BOUNDARY__",
                           "BK_O3",
                           std::round(CaBurstFullSimdef::BK_ro() * spiny_area *
                                      CaBurstFullSimdef::BK_O3_p()));
    PatchCount bk_o4_spiny("spiny.__BOUNDARY__",
                           "BK_O4",
                           std::round(CaBurstFullSimdef::BK_ro() * spiny_area *
                                      CaBurstFullSimdef::BK_O4_p()));

    PatchCount sk_c1_spiny("spiny.__BOUNDARY__",
                           "SK_C1",
                           std::round(CaBurstFullSimdef::SK_ro() * spiny_area *
                                      CaBurstFullSimdef::SK_C1_p()));
    PatchCount sk_c2_spiny("spiny.__BOUNDARY__",
                           "SK_C2",
                           std::round(CaBurstFullSimdef::SK_ro() * spiny_area *
                                      CaBurstFullSimdef::SK_C2_p()));
    PatchCount sk_c3_spiny("spiny.__BOUNDARY__",
                           "SK_C3",
                           std::round(CaBurstFullSimdef::SK_ro() * spiny_area *
                                      CaBurstFullSimdef::SK_C3_p()));
    PatchCount sk_c4_spiny("spiny.__BOUNDARY__",
                           "SK_C4",
                           std::round(CaBurstFullSimdef::SK_ro() * spiny_area *
                                      CaBurstFullSimdef::SK_C4_p()));

    PatchCount sk_o1_spiny("spiny.__BOUNDARY__",
                           "SK_O1",
                           std::round(CaBurstFullSimdef::SK_ro() * spiny_area *
                                      CaBurstFullSimdef::SK_O1_p()));
    PatchCount sk_o2_spiny("spiny.__BOUNDARY__",
                           "SK_O2",
                           std::round(CaBurstFullSimdef::SK_ro() * spiny_area *
                                      CaBurstFullSimdef::SK_O2_p()));

    PatchCount leak_spiny("spiny.__BOUNDARY__",
                          "Leak",
                          std::round(CaBurstFullSimdef::L_ro_spiny() * spiny_area));

    simulation.setPatchCount(
        {pump_smooth,  cap_m0_smooth, cap_m1_smooth, cap_m2_smooth, cap_m3_smooth, bk_c0_smooth,
         bk_c1_smooth, bk_c2_smooth,  bk_c3_smooth,  bk_c4_smooth,  bk_o0_smooth,  bk_o1_smooth,
         bk_o2_smooth, bk_o3_smooth,  bk_o4_smooth,  sk_c1_smooth,  sk_c2_smooth,  sk_c3_smooth,
         sk_c4_smooth, sk_o1_smooth,  sk_o2_smooth,  ampa_c_smooth, leak_smooth,   pump_spiny,
         cap_m0_spiny, cap_m1_spiny,  cap_m2_spiny,  cap_m3_spiny,  bk_c0_spiny,   bk_c1_spiny,
         bk_c2_spiny,  bk_c3_spiny,   bk_c4_spiny,   bk_o0_spiny,   bk_o1_spiny,   bk_o2_spiny,
         bk_o3_spiny,  bk_o4_spiny,   sk_c1_spiny,   sk_c2_spiny,   sk_c3_spiny,   sk_c4_spiny,
         sk_o1_spiny,  sk_o2_spiny,   leak_spiny});
}

void CaBurstFull::run_simulation_impl(simulation_t& simulation) {
    auto sim_time = 0.0f;
    std::stringstream s;
    steps::model::Model model;
    CaBurstFullSimdef simdef(model, simulation.getMesh());

    simulation.setPotential(simdef.init_pot());

    const std::vector<std::pair<std::string, std::string>> channels = {
        {"smooth_memb", "BKchan"},
        {"spiny_memb", "BKchan"},
        {"smooth_memb", "SKchan"},
        {"spiny_memb", "SKchan"},
        {"smooth_memb", "AMPA"}};

    s << "it t "
      << "smooth_max_V_on_verts "
      << "smooth_min_V_on_verts "
      << "spiny_max_V_on_verts "
      << "spiny_min_V_on_verts "
      << "smooth_max_V_on_tris "
      << "smooth_min_V_on_tris "
      << "spiny_max_V_on_tris "
      << "spiny_min_V_on_tris "
      << "smooth_CaP_m0 "
      << "smooth_CaP_m1 "
      << "smooth_CaP_m2 "
      << "smooth_CaP_m3 "
      << "spiny_CaP_m0 "
      << "spiny_CaP_m1 "
      << "spiny_CaP_m2 "
      << "spiny_CaP_m3 "
      << "tot_GHK_curr";
    for (auto p : channels) {
      s << ' ' << p.first << '_' << p.second << "_ohm_curr";
    }
    s << '\n';

    simulation.log_once(s.str());

    n_timepoints = simulation.getScenario().end_time / time_point_interval;

    s << "running for " << n_timepoints << " iterations\n";

    simulation.log_once(s.str());

    for (auto i = 0u; i < n_timepoints; i++) {
      sim_time = time_point_interval * i;
      simulation.setKCstSReac("smooth.__BOUNDARY__", simdef.ampacc1,
                              1.0e-3 * CaBurstFullSimdef::rb() *
                                  glut_pulse_mvr[i + 2000]);
      simulation.setKCstSReac("smooth.__BOUNDARY__", simdef.ampac1c2,
                              1.0e-3 * CaBurstFullSimdef::rb() *
                                  glut_pulse_mvr[i + 2000]);
      simulation.run(osh::Real(sim_time));

      if (i % 10 == 0) {
        s.str("");
        s << i << ' ' << sim_time << std::setprecision(8) << ' '
          << simulation.getMaxPotentialOnVertices("smooth.__BOUNDARY__") << ' '
          << simulation.getMinPotentialOnVertices("smooth.__BOUNDARY__") << ' '
          << simulation.getMaxPotentialOnVertices("spiny.__BOUNDARY__") << ' '
          << simulation.getMinPotentialOnVertices("spiny.__BOUNDARY__")

          << ' ' << simulation.getMaxPotentialOnTriangles("smooth.__BOUNDARY__")
          << ' ' << simulation.getMinPotentialOnTriangles("smooth.__BOUNDARY__")
          << ' ' << simulation.getMaxPotentialOnTriangles("spiny.__BOUNDARY__")
          << ' ' << simulation.getMinPotentialOnTriangles("spiny.__BOUNDARY__");

        //          << ' ' << simulation.getPatchCount("smooth.__BOUNDARY__",
        //          "CaP_m0")
        //          << ' ' << simulation.getPatchCount("smooth.__BOUNDARY__",
        //          "CaP_m1")
        //          << ' ' << simulation.getPatchCount("smooth.__BOUNDARY__",
        //          "CaP_m2")
        //          << ' ' << simulation.getPatchCount("smooth.__BOUNDARY__",
        //          "CaP_m3")
        //
        //          << ' ' << simulation.getPatchCount("spiny.__BOUNDARY__",
        //          "CaP_m0")
        //          << ' ' << simulation.getPatchCount("spiny.__BOUNDARY__",
        //          "CaP_m1")
        //          << ' ' << simulation.getPatchCount("spiny.__BOUNDARY__",
        //          "CaP_m2")
        //          << ' ' << simulation.getPatchCount("spiny.__BOUNDARY__",
        //          "CaP_m3")
        //
        //          << ' ' << simulation.getTotalGHKCurrent();

        //        for (auto mem_chan : channels) {
        //          s << ' '
        //            << simulation.getTotalOhmicCurrent(mem_chan.first.data(),
        //                                               mem_chan.second.data());
        //        }
        s << '\n';

        simulation.log_once(s.str());
      }
      simulation.log_progress(i, n_timepoints);
    }
}

}  // namespace dist
}  // namespace steps
