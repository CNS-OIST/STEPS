#include "simdef.hpp"

#include <vector>

#include "geom/dist/distmesh.hpp"
#include "mpi/dist/test/simulation.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"

namespace steps {
namespace dist {

Simdef::Simdef(const steps::model::Model &model,
               const steps::dist::DistMesh &mesh)
    : statedef(new Statedef(model, mesh)) {}

SimpleSimdef::SimpleSimdef(const steps::model::Model &model,
                           const steps::dist::DistMesh &mesh)
    : Simdef(model, mesh) {
  statedef->addComp("comp1");

  const std::vector<model::species_name> species = {"A", "B", "C", "D", "E",
                                                    "F", "G", "H", "I", "J"};
  statedef->addCompSpecs("comp1", species);
  const std::vector<model::species_name> R1_lhs = {"A", "B"};
  const std::vector<model::species_name> R1_rhs = {"C"};
  statedef->addCompReac("comp1", R1_lhs, R1_rhs, 1000.0e6);

  const std::vector<model::species_name> R2_lhs = {"C"};
  const std::vector<model::species_name> R2_rhs = {"A", "B"};
  statedef->addCompReac("comp1", R2_lhs, R2_rhs, 100);

  const std::vector<model::species_name> R3_lhs = {"C", "D"};
  const std::vector<model::species_name> R3_rhs = {"E"};
  statedef->addCompReac("comp1", R3_lhs, R3_rhs, 100e6);

  const std::vector<model::species_name> R4_lhs = {"E"};
  const std::vector<model::species_name> R4_rhs = {"C", "D"};
  statedef->addCompReac("comp1", R4_lhs, R4_rhs, 10);

  const std::vector<model::species_name> R5_lhs = {"F", "G"};
  const std::vector<model::species_name> R5_rhs = {"H"};
  statedef->addCompReac("comp1", R5_lhs, R5_rhs, 10e6);

  const std::vector<model::species_name> R6_lhs = {"H"};
  const std::vector<model::species_name> R6_rhs = {"F", "G"};
  statedef->addCompReac("comp1", R6_lhs, R6_rhs, 1);

  const std::vector<model::species_name> R7_lhs = {"H", "I"};
  const std::vector<model::species_name> R7_rhs = {"J"};
  statedef->addCompReac("comp1", R7_lhs, R7_rhs, 1e6);

  const std::vector<model::species_name> R8_lhs = {"J"};
  const std::vector<model::species_name> R8_rhs = {"H", "I"};
  statedef->addCompReac("comp1", R8_lhs, R8_rhs, 0.1 * 10);

  statedef->addCompDiff("comp1", "A", 100e-12);
  statedef->addCompDiff("comp1", "B", 90e-12);
  statedef->addCompDiff("comp1", "C", 80e-12);
  statedef->addCompDiff("comp1", "D", 70e-12);
  statedef->addCompDiff("comp1", "E", 60e-12);
  statedef->addCompDiff("comp1", "F", 50e-12);
  statedef->addCompDiff("comp1", "G", 40e-12);
  statedef->addCompDiff("comp1", "H", 30e-12);
  statedef->addCompDiff("comp1", "I", 20e-12);
  statedef->addCompDiff("comp1", "J", 10e-12);

  const osh::Real N0A = 1000.0;
  const osh::Real N0B = 2000.0;
  const osh::Real N0C = 3000.0;
  const osh::Real N0D = 4000.0;
  const osh::Real N0E = 5000.0;
  const osh::Real N0F = 6000.0;
  const osh::Real N0G = 7000.0;
  const osh::Real N0H = 8000.0;
  const osh::Real N0I = 9000.0;
  const osh::Real N0J = 10000.0;

  compartmentCounts = {{"comp1",
                        {
                            {"A", N0A},
                            {"B", N0B},
                            {"C", N0C},
                            {"D", N0D},
                            {"E", N0E},
                            {"F", N0F},
                            {"G", N0G},
                            {"H", N0H},
                            {"I", N0I},
                            {"J", N0J},
                        }}};
}

SingleCompDistSimdef::SingleCompDistSimdef(const steps::model::Model &model,
                                           const steps::dist::DistMesh &mesh)
    : Simdef(model, mesh) {

  statedef->addComp("comp1");
  const std::vector<model::species_name> species = {"C"};
  statedef->addCompSpecs("comp1", species);

  compartmentCounts = {{"comp1", {{"C", 100000}}}};
}

MultiCompartmentDiffSimdef::MultiCompartmentDiffSimdef(
    const steps::model::Model &model, const steps::dist::DistMesh &mesh)
    : Simdef(model, mesh) {
  osh::Real dcst = 1.0e-8;

  // PetscOptionsGetScalar(nullptr, nullptr, "-d", &dcst, nullptr);
  statedef->addComp("comp_i");
  const std::vector<model::species_name> left_species = {"C"};
  statedef->addCompSpecs("comp_i", left_species);

  // add reaction2
  statedef->addCompDiff("comp_i", "C", dcst);

  // Right compartment
  statedef->addComp("comp_o");
  const std::vector<model::species_name> right_species = {"D", "E", "C"};
  statedef->addCompSpecs("comp_o", right_species);
  statedef->addCompDiff("comp_o", "C", dcst);
  compartmentCounts = {{"comp_i", {{"C", 1000000}}},
                       {"comp_o", {{"C", 1000000}}}};
}

MultiCompartmentSimdef::MultiCompartmentSimdef(
    const steps::model::Model &model, const steps::dist::DistMesh &mesh)
    : Simdef(model, mesh) {
  osh::Real const r1_kcst = 1000000.0e6;
  osh::Real const r2_kcst = 0.0005e6;
  osh::Real dcst = 1e-8;

  statedef->addComp("Left");
  const std::vector<model::species_name> left_species = {"A", "B", "C"};
  statedef->addCompSpecs("Left", left_species);

  // add reaction1
  const std::vector<model::species_name> left_r1_lhs = {"A", "B"};
  const std::vector<model::species_name> left_r1_rhs = {"C"};
  statedef->addCompReac("Left", left_r1_lhs, left_r1_rhs, r1_kcst);

  // add reaction2
  const std::vector<model::species_name> left_r2_lhs = {"C"};
  const std::vector<model::species_name> left_r2_rhs = {"A", "B"};
  statedef->addCompReac("Left", left_r2_lhs, left_r2_rhs, r2_kcst);

  statedef->addCompDiff("Left", "A", dcst);
  statedef->addCompDiff("Left", "B", dcst);
  statedef->addCompDiff("Left", "C", dcst);

  // Right compartment

  statedef->addComp("Right");
  const std::vector<model::species_name> right_species = {"C", "D", "E"};
  statedef->addCompSpecs("Right", right_species);

  // add reaction1
  const std::vector<model::species_name> right_r1_lhs = {"C", "D"};
  const std::vector<model::species_name> right_r1_rhs = {"E"};
  statedef->addCompReac("Right", right_r1_lhs, right_r1_rhs, r1_kcst);

  // add reaction2
  const std::vector<model::species_name> right_r2_lhs = {"E"};
  const std::vector<model::species_name> right_r2_rhs = {"C", "D"};
  statedef->addCompReac("Right", right_r2_lhs, right_r2_rhs, r2_kcst);
  statedef->addCompDiff("Right", "C", dcst);
  statedef->addCompDiff("Right", "D", dcst);
  statedef->addCompDiff("Right", "E", dcst);

  compartmentCounts = {{"Left", {{"A", 10}, {"B", 10}}},
                       {"Right", {{"E", 10}}}};
}

SReacUnitTestSimdef::SReacUnitTestSimdef(const steps::model::Model &model,
                                         const steps::dist::DistMesh &mesh)
    : Simdef(model, mesh) {
  statedef->addComp("comp_i");
  statedef->addComp("comp_o");
  statedef->addPatch("patch", "comp_i", model::compartment_id("comp_o"));

  statedef->addCompSpecs("comp_i", {"A", "B"});
  statedef->addCompSpecs("comp_o", {"C"});
  statedef->addPatchSpecs("patch", {"D"});

  const std::vector<model::species_name> R1_lhs = {"A", "B"};
  statedef->addCompReac("comp_i", R1_lhs, {}, 1e8);

  const std::vector<model::species_name> R2_lhs_p = {"D", "D", "D"};
  const std::vector<model::species_name> R2_lhs_o = {"C", "C"};
  const std::vector<model::species_name> R2_rhs_i = {"B", "B", "B", "B"};
  statedef->addSurfReac("patch", {}, R2_lhs_p, R2_lhs_o, R2_rhs_i, {}, {},
                        1e21);

  const osh::Real N0A = 4000.0;
  const osh::Real N0B = 0.0;
  const osh::Real N0C = 2000.0;
  const osh::Real N0D = 6000.0;

  compartmentCounts = {{"comp_i", {{"A", N0A}, {"B", N0B}}},
                       {"comp_o", {{"C", N0C}}}};
  patchCounts.emplace_back("patch", "D", N0D);
}

GHKCurrentUnitTestSimdef::GHKCurrentUnitTestSimdef(
    const steps::model::Model &model, const steps::dist::DistMesh &mesh)
    : Simdef(model, mesh) {
  statedef->disableEField();
  statedef->setTemp(273 + 30);
  statedef->addComp("comp_i");
  statedef->addComp("comp_o");
  statedef->addPatch("patch", "comp_i", model::compartment_id("comp_o"));
  statedef->addMembrane("membrane", "patch", 0.0);
  statedef->addChannel("membrane", "NaChan", {"NaChan_state_0"});
  statedef->addCompartmentConductivity("comp_i", 1.0);
  statedef->addCompartmentConductivity("comp_o", 1.0);
  statedef->addCompSpecs("comp_i", {"Na"});
  statedef->addCompSpecs("comp_o", {"Na"});
  statedef->addPatchSpecs("patch", {"NaChan_state_0"});
  statedef->addGHKCurrentSurfReac(
      "NaCurr", "membrane", "NaChan", "NaChan_state_0", "Na", 1.0e-14, 2);
  const osh::Real N0AI = 20e6;
  const osh::Real N0AO = 40e6;
  const osh::Real N0B = 1.0;

  compartmentCounts = {{"comp_i", {{"Na", N0AI}}}, {"comp_o", {{"Na", N0AO}}}};
  patchCounts.emplace_back("patch", "NaChan_state_0", N0B);
}

ValidationSimdef::ValidationSimdef(const steps::model::Model &model,
                                   const steps::dist::DistMesh &mesh)
    : Simdef(model, mesh) {
  statedef->addComp("comp1");
  const std::vector<model::species_name> comp_species = {
      "A_foi",  "A_for",   "B_for",   "A_soA2",  "C_soA2", "A_soAA",
      "B_soAA", "C_soAA",  "A_soAB",  "B_soAB",  "C_soAB", "A_toA3",
      "C_toA3", "A_toA2B", "B_toA2B", "C_toA2B",
  };
  statedef->addCompSpecs("comp1", comp_species);
  // First order irreversible

  const std::vector<model::species_name> R1_foi_lhs = {"A_foi"};
  const std::vector<model::species_name> R1_foi_rhs = {};
  statedef->addCompReac("comp1", R1_foi_lhs, R1_foi_rhs, KCST_foi());

  // First order reversible
  const std::vector<model::species_name> R1_for_lhs = {"A_for"};
  const std::vector<model::species_name> R1_for_rhs = {"B_for"};
  statedef->addCompReac("comp1", R1_for_lhs, R1_for_rhs, KCST_f_for());

  const std::vector<model::species_name> R2_for_lhs = {"B_for"};
  const std::vector<model::species_name> R2_for_rhs = {"A_for"};
  statedef->addCompReac("comp1", R2_for_lhs, R2_for_rhs, KCST_b_for());

  // Second order irreversible A2
  const std::vector<model::species_name> R1_soA2_lhs = {"A_soA2", "A_soA2"};
  const std::vector<model::species_name> R1_soA2_rhs = {"C_soA2"};
  statedef->addCompReac("comp1", R1_soA2_lhs, R1_soA2_rhs, KCST_soA2());

  // Second order irreversible AA
  const std::vector<model::species_name> R1_soAA_lhs = {"A_soAA", "B_soAA"};
  const std::vector<model::species_name> R1_soAA_rhs = {"C_soAA"};
  statedef->addCompReac("comp1", R1_soAA_lhs, R1_soAA_rhs, KCST_soAA());

  // Second order irreversible AB
  const std::vector<model::species_name> R1_soAB_lhs = {"A_soAB", "B_soAB"};
  const std::vector<model::species_name> R1_soAB_rhs = {"C_soAB"};
  statedef->addCompReac("comp1", R1_soAB_lhs, R1_soAB_rhs, KCST_soAB());

  // Third order irreversible A3
  const std::vector<model::species_name> R1_toA3_lhs = {"A_toA3", "A_toA3",
                                                        "A_toA3"};
  const std::vector<model::species_name> R1_toA3_rhs = {"C_toA3"};
  statedef->addCompReac("comp1", R1_toA3_lhs, R1_toA3_rhs, KCST_toA3());

  // Third order irreversible A2B
  const std::vector<model::species_name> R1_toA2B_lhs = {"A_toA2B", "A_toA2B",
                                                         "B_toA2B"};
  const std::vector<model::species_name> R1_toA2B_rhs = {"C_toA2B"};
  statedef->addCompReac("comp1", R1_toA2B_lhs, R1_toA2B_rhs, KCST_toA2B());

  statedef->addCompDiff("comp1", "A_foi", 0.01e-12);
  statedef->addCompDiff("comp1", "A_for", 0.01e-12);
  statedef->addCompDiff("comp1", "B_for", 0.01e-12);
  statedef->addCompDiff("comp1", "A_soA2", 1e-12);

  statedef->addCompDiff("comp1", "A_soAA", 0.2e-12);
  statedef->addCompDiff("comp1", "B_soAA", 0.2e-12);
  statedef->addCompDiff("comp1", "A_soAB", 0.1e-12);
  statedef->addCompDiff("comp1", "B_soAB", 0.1e-12);
  statedef->addCompDiff("comp1", "A_toA3", 0.2e-12);
  statedef->addCompDiff("comp1", "A_toA2B", 0.1e-12);
  statedef->addCompDiff("comp1", "B_toA2B", 0.1e-12);
}

ValidationSimdef::compartment_counts_t
ValidationSimdef::getfoi_CompartementCounts(bool zero) const {
  return {{"comp1", {{"A_foi", zero ? 0 : N_foi()}}}};
}

ValidationSimdef::compartment_counts_t
ValidationSimdef::getfor_CompartementCounts(bool zero) const {
  return {{"comp1", {{"A_for", zero ? 0 : COUNT_for()}, {"B_for", 0}}}};
}

ValidationSimdef::compartment_concs_t
ValidationSimdef::getsoA2_CompartementConc(bool zero) const {
  return {{"comp1", {{"A_soA2", zero ? 0 : CONCA_soA2()}}}};
}

ValidationSimdef::compartment_concs_t
ValidationSimdef::getsoAA_CompartementConc(bool zero) const {
  return {{"comp1",
           {{"A_soAA", zero ? 0 : CONCA_soAA()},
            {"B_soAA", zero ? 0 : CONCB_soAA()}}}};
}

ValidationSimdef::compartment_concs_t
ValidationSimdef::getsoAB_CompartementConc(bool zero) const {
  return {{"comp1",
           {{"A_soAB", zero ? 0 : CONCA_soAB()},
            {"B_soAB", zero ? 0 : CONCB_soAB()}}}};
}

ValidationSimdef::compartment_concs_t
ValidationSimdef::gettoA3_CompartementConc(bool zero) const {
  return {{"comp1", {{"A_toA3", zero ? 0 : CONCA_toA3()}}}};
}

ValidationSimdef::compartment_concs_t
ValidationSimdef::gettoA2B_CompartementConc(bool zero) const {
  return {{"comp1",
           {{"A_toA2B", zero ? 0 : CONCA_toA2B()},
            {"B_toA2B", zero ? 0 : CONCB_toA2B()}}}};
}

SurfaceReactionsValidationSimdef::SurfaceReactionsValidationSimdef(
    const steps::model::Model &model, const steps::dist::DistMesh &mesh)
    : Simdef(model, mesh) {
  statedef->disableEField();
  statedef->addComp("comp1");
  const std::vector<model::species_name> comp_species = {
      "A_foi",  "A_for",   "B_for",   "A_soA2",  "C_soA2", "A_soAA",
      "B_soAA", "C_soAA",  "A_soAB",  "B_soAB",  "C_soAB", "A_toA3",
      "C_toA3", "A_toA2B", "B_toA2B", "C_toA2B",
  };
  statedef->addCompSpecs("comp1", comp_species);
  statedef->addPatch("patch1", "comp1");
  const std::vector<model::species_name> patch_species = {
      "A_so2d",
      "B_so2d",
      "C_so2d",
  };
  statedef->addPatchSpecs("patch1", patch_species);
  // First order irreversible
  const std::vector<model::species_name> R1_foi_lhs = {"A_foi"};
  const std::vector<model::species_name> R1_foi_rhs = {};
  statedef->addCompReac("comp1", R1_foi_lhs, R1_foi_rhs, KCST_foi());

  // First order reversible
  const std::vector<model::species_name> R1_for_lhs = {"A_for"};
  const std::vector<model::species_name> R1_for_rhs = {"B_for"};
  statedef->addCompReac("comp1", R1_for_lhs, R1_for_rhs, KCST_f_for());

  const std::vector<model::species_name> R2_for_lhs = {"B_for"};
  const std::vector<model::species_name> R2_for_rhs = {"A_for"};
  statedef->addCompReac("comp1", R2_for_lhs, R2_for_rhs, KCST_b_for());

  // Second order irreversible A2
  const std::vector<model::species_name> R1_soA2_lhs = {"A_soA2", "A_soA2"};
  const std::vector<model::species_name> R1_soA2_rhs = {"C_soA2"};
  statedef->addCompReac("comp1", R1_soA2_lhs, R1_soA2_rhs, KCST_soA2());

  // Second order irreversible AA
  const std::vector<model::species_name> R1_soAA_lhs = {"A_soAA", "B_soAA"};
  const std::vector<model::species_name> R1_soAA_rhs = {"C_soAA"};
  statedef->addCompReac("comp1", R1_soAA_lhs, R1_soAA_rhs, KCST_soAA());

  // Second order irreversible AB
  const std::vector<model::species_name> R1_soAB_lhs = {"A_soAB", "B_soAB"};
  const std::vector<model::species_name> R1_soAB_rhs = {"C_soAB"};
  statedef->addCompReac("comp1", R1_soAB_lhs, R1_soAB_rhs, KCST_soAB());

  // Third order irreversible A3
  const std::vector<model::species_name> R1_toA3_lhs = {"A_toA3", "A_toA3",
                                                        "A_toA3"};
  const std::vector<model::species_name> R1_toA3_rhs = {"C_toA3"};
  statedef->addCompReac("comp1", R1_toA3_lhs, R1_toA3_rhs, KCST_toA3());

  // Third order irreversible A2B
  const std::vector<model::species_name> R1_toA2B_lhs = {"A_toA2B", "A_toA2B",
                                                         "B_toA2B"};
  const std::vector<model::species_name> R1_toA2B_rhs = {"C_toA2B"};
  statedef->addCompReac("comp1", R1_toA2B_lhs, R1_toA2B_rhs, KCST_toA2B());

  // Second order irreversible 2D
  const std::vector<model::species_name> R1_so2d_lhs = {"A_so2d", "B_so2d"};
  const std::vector<model::species_name> R1_so2d_rhs = {"C_so2d"};
  statedef->addSurfReac("patch1", {}, R1_so2d_lhs, {}, {}, R1_so2d_rhs, {},
                        KCST_so2d());

  statedef->addCompDiff("comp1", "A_foi", 0.01e-12);
  statedef->addCompDiff("comp1", "A_for", 0.01e-12);
  statedef->addCompDiff("comp1", "B_for", 0.01e-12);
  statedef->addCompDiff("comp1", "A_soA2", 1e-12);
  statedef->addCompDiff("comp1", "A_soAA", 0.2e-12);
  statedef->addCompDiff("comp1", "B_soAA", 0.2e-12);
  statedef->addCompDiff("comp1", "A_soAB", 0.1e-12);
  statedef->addCompDiff("comp1", "B_soAB", 0.1e-12);
  statedef->addCompDiff("comp1", "A_toA3", 0.2e-12);
  statedef->addCompDiff("comp1", "A_toA2B", 0.1e-12);
  statedef->addCompDiff("comp1", "B_toA2B", 0.1e-12);
}

SurfaceReactionsValidationSimdef::compartment_counts_t
SurfaceReactionsValidationSimdef::getfoi_CompartementCounts() const {
  return {{"comp1", {{"A_foi", N_foi()}}}};
}

SurfaceReactionsValidationSimdef::compartment_counts_t
SurfaceReactionsValidationSimdef::getfor_CompartementCounts() const {
  return {{"comp1", {{"A_for", COUNT_for()}, {"B_for", 0.0}}}};
}

SurfaceReactionsValidationSimdef::compartment_concs_t
SurfaceReactionsValidationSimdef::getsoA2_CompartementConc() const {
  return {{"comp1", {{"A_soA2", CONCA_soA2()}}}};
}

SurfaceReactionsValidationSimdef::compartment_concs_t
SurfaceReactionsValidationSimdef::getsoAA_CompartementConc() const {
  return {{"comp1", {{"A_soAA", CONCA_soAA()}, {"B_soAA", CONCB_soAA()}}}};
}

SurfaceReactionsValidationSimdef::compartment_concs_t
SurfaceReactionsValidationSimdef::getsoAB_CompartementConc() const {
  return {{"comp1", {{"A_soAB", CONCA_soAB()}, {"B_soAB", CONCB_soAB()}}}};
}

SurfaceReactionsValidationSimdef::compartment_concs_t
SurfaceReactionsValidationSimdef::gettoA3_CompartementConc() const {
  return {{"comp1", {{"A_toA3", CONCA_toA3()}}}};
}

SurfaceReactionsValidationSimdef::compartment_concs_t
SurfaceReactionsValidationSimdef::gettoA2B_CompartementConc() const {
  return {{"comp1", {{"A_toA2B", CONCA_toA2B()}, {"B_toA2B", CONCB_toA2B()}}}};
}

SurfaceReactionsValidationSimdef::patch_counts_t
SurfaceReactionsValidationSimdef::getso2d_PatchCount() const {
  patch_counts_t counts;

  counts.emplace_back("patch1", "A_so2d", COUNTA_so2d());
  counts.emplace_back("patch1", "B_so2d", COUNTB_so2d());

  return counts;
}

CaBurstSimdef::CaBurstSimdef(const steps::model::Model &model,
                             const steps::dist::DistMesh &mesh)
    : Simdef(model, mesh) {
    statedef->setTemp(273.15 + TEMPERATURE());
    statedef->addComp("__MESH__");
    const std::vector<model::species_name> comp_species = {"Ca",
                                                           "Pump",
                                                           "CaPump",
                                                           "iCBsf",
                                                           "iCBsCa",
                                                           "iCBCaf",
                                                           "iCBCaCa",
                                                           "CBsf",
                                                           "CBsCa",
                                                           "CBCaf",
                                                           "CBCaCa",
                                                           "PV",
                                                           "PVMg",
                                                           "PVCa",
                                                           "Mg"};
    statedef->addCompSpecs("__MESH__", comp_species);

    const std::vector<model::species_name> iCBsf1_f_lhs = {"Ca", "iCBsf"};
    const std::vector<model::species_name> iCBsf1_f_rhs = {"iCBsCa"};
    statedef->addCompReac("__MESH__", iCBsf1_f_lhs, iCBsf1_f_rhs, iCBsf1_f_kcst());
    statedef->addCompReac("__MESH__", iCBsf1_f_rhs, iCBsf1_f_lhs, iCBsf1_b_kcst());

    const std::vector<model::species_name> iCBsCa_f_lhs = {"Ca", "iCBsCa"};
    const std::vector<model::species_name> iCBsCa_f_rhs = {"iCBCaCa"};
    statedef->addCompReac("__MESH__", iCBsCa_f_lhs, iCBsCa_f_rhs, iCBsCa_f_kcst());
    statedef->addCompReac("__MESH__", iCBsCa_f_rhs, iCBsCa_f_lhs, iCBsCa_b_kcst());

    const std::vector<model::species_name> iCBsf2_f_lhs = {"Ca", "iCBsf"};
    const std::vector<model::species_name> iCBsf2_f_rhs = {"iCBCaf"};
    statedef->addCompReac("__MESH__", iCBsf2_f_lhs, iCBsf2_f_rhs, iCBsf2_f_kcst());
    statedef->addCompReac("__MESH__", iCBsf2_f_rhs, iCBsf2_f_lhs, iCBsf2_b_kcst());

    const std::vector<model::species_name> iCBCaf_f_lhs = {"Ca", "iCBCaf"};
    const std::vector<model::species_name> iCBCaf_f_rhs = {"iCBCaCa"};
    statedef->addCompReac("__MESH__", iCBCaf_f_lhs, iCBCaf_f_rhs, iCBCaf_f_kcst());
    statedef->addCompReac("__MESH__", iCBCaf_f_rhs, iCBCaf_f_lhs, iCBCaf_b_kcst());

    const std::vector<model::species_name> CBsf1_f_lhs = {"Ca", "CBsf"};
    const std::vector<model::species_name> CBsf1_f_rhs = {"CBsCa"};
    statedef->addCompReac("__MESH__", CBsf1_f_lhs, CBsf1_f_rhs, CBsf1_f_kcst());
    statedef->addCompReac("__MESH__", CBsf1_f_rhs, CBsf1_f_lhs, CBsf1_b_kcst());

    const std::vector<model::species_name> CBsCa_f_lhs = {"Ca", "CBsCa"};
    const std::vector<model::species_name> CBsCa_f_rhs = {"CBCaCa"};
    statedef->addCompReac("__MESH__", CBsCa_f_lhs, CBsCa_f_rhs, CBsCa_f_kcst());
    statedef->addCompReac("__MESH__", CBsCa_f_rhs, CBsCa_f_lhs, CBsCa_b_kcst());

    const std::vector<model::species_name> CBsf2_f_lhs = {"Ca", "CBsf"};
    const std::vector<model::species_name> CBsf2_f_rhs = {"CBCaf"};
    statedef->addCompReac("__MESH__", CBsf2_f_lhs, CBsf2_f_rhs, CBsf2_f_kcst());
    statedef->addCompReac("__MESH__", CBsf2_f_rhs, CBsf2_f_lhs, CBsf2_b_kcst());

    const std::vector<model::species_name> CBCaf_f_lhs = {"Ca", "CBCaf"};
    const std::vector<model::species_name> CBCaf_f_rhs = {"CBCaCa"};
    statedef->addCompReac("__MESH__", CBCaf_f_lhs, CBCaf_f_rhs, CBCaf_f_kcst());
    statedef->addCompReac("__MESH__", CBCaf_f_rhs, CBCaf_f_lhs, CBCaf_b_kcst());

    const std::vector<model::species_name> PVca_f_lhs = {"Ca", "PV"};
    const std::vector<model::species_name> PVca_f_rhs = {"PVCa"};
    statedef->addCompReac("__MESH__", PVca_f_lhs, PVca_f_rhs, PVca_f_kcst());
    statedef->addCompReac("__MESH__", PVca_f_rhs, PVca_f_lhs, PVca_b_kcst());

    const std::vector<model::species_name> PVmg_f_lhs = {"Mg", "PV"};
    const std::vector<model::species_name> PVmg_f_rhs = {"PVMg"};
    statedef->addCompReac("__MESH__", PVmg_f_lhs, PVmg_f_rhs, PVmg_f_kcst());
    statedef->addCompReac("__MESH__", PVmg_f_rhs, PVmg_f_lhs, PVmg_b_kcst());

    statedef->addCompDiff("__MESH__", "Ca", DCST());
    statedef->addCompDiff("__MESH__", "CBsf", DCB());
    statedef->addCompDiff("__MESH__", "CBsCa", DCB());
    statedef->addCompDiff("__MESH__", "CBCaf", DCB());
    statedef->addCompDiff("__MESH__", "PV", DPV());
    statedef->addCompDiff("__MESH__", "PVCa", DPV());
    statedef->addCompDiff("__MESH__", "PVMg", DPV());

    statedef->addPatch("smooth.__BOUNDARY__", "__MESH__");
    statedef->addPatch("spiny.__BOUNDARY__", "__MESH__");
    const std::vector<model::species_name> patch_species = {
        "Pump",
        "CaPump",
    };
    statedef->addPatchSpecs("smooth.__BOUNDARY__", patch_species);
    statedef->addPatchSpecs("spiny.__BOUNDARY__", patch_species);

    const std::vector<model::species_name> PumpD_f_ilhs = {"Ca"};
    const std::vector<model::species_name> PumpD_f_slhs = {"Pump"};
    const std::vector<model::species_name> PumpD_f_srhs = {"CaPump"};
    statedef->addSurfReac(
        "smooth.__BOUNDARY__", PumpD_f_ilhs, PumpD_f_slhs, {}, {}, PumpD_f_srhs, {}, P_f_kcst());
    statedef->addSurfReac(
        "spiny.__BOUNDARY__", PumpD_f_ilhs, PumpD_f_slhs, {}, {}, PumpD_f_srhs, {}, P_f_kcst());
    const std::vector<model::species_name> PumpD_b_slhs = {"CaPump"};
    const std::vector<model::species_name> PumpD_b_irhs = {"Ca"};
    const std::vector<model::species_name> PumpD_b_srhs = {"Pump"};
    statedef->addSurfReac(
        "smooth.__BOUNDARY__", {}, PumpD_b_slhs, {}, PumpD_b_irhs, PumpD_b_srhs, {}, P_b_kcst());
    statedef->addSurfReac(
        "spiny.__BOUNDARY__", {}, PumpD_b_slhs, {}, PumpD_b_irhs, PumpD_b_srhs, {}, P_b_kcst());
    const std::vector<model::species_name> PumpD_k_slhs = {"CaPump"};
    const std::vector<model::species_name> PumpD_k_srhs = {"Pump"};
    statedef->addSurfReac(
        "smooth.__BOUNDARY__", {}, PumpD_k_slhs, {}, {}, PumpD_k_srhs, {}, P_k_kcst());
    statedef->addSurfReac(
        "spiny.__BOUNDARY__", {}, PumpD_k_slhs, {}, {}, PumpD_k_srhs, {}, P_k_kcst());
}

CaBurstFullSimdef::CaBurstFullSimdef(const steps::model::Model &model,
                                     const steps::dist::DistMesh &mesh)
    : CaBurstSimdef(model, mesh) {
  statedef->addCompDiff("__MESH__", "CBCaCa", DCB());
  statedef->addCompartmentConductivity("__MESH__", 1.0 / Ra());

  statedef->addMembrane("smooth_memb", "smooth.__BOUNDARY__",
                        memb_capac_proximal());
  statedef->addPatchSpecs(
      "smooth.__BOUNDARY__",
      {"CaP_m0",  "CaP_m1",  "CaP_m2",  "CaP_m3",  "BK_C0",  "BK_C1", "BK_C2",
       "BK_C3",   "BK_C4",   "BK_O0",   "BK_O1",   "BK_O2",  "BK_O3", "BK_O4",
       "SK_C1",   "SK_C2",   "SK_C3",   "SK_C4",   "SK_O1",  "SK_O2", "AMPA_C",
       "AMPA_C1", "AMPA_C2", "AMPA_D1", "AMPA_D2", "AMPA_O", "Leak"});

  statedef->addChannel("smooth_memb", "CaPchan",
                       {"CaP_m0", "CaP_m1", "CaP_m2", "CaP_m3"});
  statedef->addVDepSurfReac(
      "smooth.__BOUNDARY__", {}, {"CaP_m0"}, {}, {}, {"CaP_m1"}, {},
      [&](auto V) { return 1.0e3 * 3. * alpha_cap(V * 1.0e3) * Qt(); });
  statedef->addVDepSurfReac(
      "smooth.__BOUNDARY__", {}, {"CaP_m1"}, {}, {}, {"CaP_m2"}, {},
      [&](auto V) { return 1.0e3 * 2. * alpha_cap(V * 1.0e3) * Qt(); });
  statedef->addVDepSurfReac(
      "smooth.__BOUNDARY__", {}, {"CaP_m2"}, {}, {}, {"CaP_m3"}, {},
      [&](auto V) { return 1.0e3 * 1. * alpha_cap(V * 1.0e3) * Qt(); });
  statedef->addVDepSurfReac(
      "smooth.__BOUNDARY__", {}, {"CaP_m3"}, {}, {}, {"CaP_m2"}, {},
      [&](auto V) { return 1.0e3 * 3. * beta_cap(V * 1.0e3) * Qt(); });
  statedef->addVDepSurfReac(
      "smooth.__BOUNDARY__", {}, {"CaP_m2"}, {}, {}, {"CaP_m1"}, {},
      [&](auto V) { return 1.0e3 * 2. * beta_cap(V * 1.0e3) * Qt(); });
  statedef->addVDepSurfReac(
      "smooth.__BOUNDARY__", {}, {"CaP_m1"}, {}, {}, {"CaP_m0"}, {},
      [&](auto V) { return 1.0e3 * 1. * beta_cap(V * 1.0e3) * Qt(); });
  statedef->addGHKCurrentSurfReac(
      "CaPCurrSmooth", "smooth_memb", "CaPchan", "CaP_m3", "Ca", CaP_P(), 2, Ca_oconc());

  statedef->addChannel("smooth_memb", "BKchan",
                       {"BK_C0", "BK_C1", "BK_C2", "BK_C3", "BK_C4", "BK_O0",
                        "BK_O1", "BK_O2", "BK_O3", "BK_O4"});
  statedef->addSurfReac("smooth.__BOUNDARY__", {"Ca"}, {"BK_C0"}, {}, {},
                        {"BK_C1"}, {}, c_01());
  statedef->addSurfReac("smooth.__BOUNDARY__", {"Ca"}, {"BK_C1"}, {}, {},
                        {"BK_C2"}, {}, c_12());
  statedef->addSurfReac("smooth.__BOUNDARY__", {"Ca"}, {"BK_C2"}, {}, {},
                        {"BK_C3"}, {}, c_23());
  statedef->addSurfReac("smooth.__BOUNDARY__", {"Ca"}, {"BK_C3"}, {}, {},
                        {"BK_C4"}, {}, c_34());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"BK_C1"}, {}, {"Ca"},
                        {"BK_C0"}, {}, c_10());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"BK_C2"}, {}, {"Ca"},
                        {"BK_C1"}, {}, c_21());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"BK_C3"}, {}, {"Ca"},
                        {"BK_C2"}, {}, c_32());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"BK_C4"}, {}, {"Ca"},
                        {"BK_C3"}, {}, c_43());

  statedef->addSurfReac("smooth.__BOUNDARY__", {"Ca"}, {"BK_O0"}, {}, {},
                        {"BK_O1"}, {}, o_01());
  statedef->addSurfReac("smooth.__BOUNDARY__", {"Ca"}, {"BK_O1"}, {}, {},
                        {"BK_O2"}, {}, o_12());
  statedef->addSurfReac("smooth.__BOUNDARY__", {"Ca"}, {"BK_O2"}, {}, {},
                        {"BK_O3"}, {}, o_23());
  statedef->addSurfReac("smooth.__BOUNDARY__", {"Ca"}, {"BK_O3"}, {}, {},
                        {"BK_O4"}, {}, o_34());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"BK_O1"}, {}, {"Ca"},
                        {"BK_O0"}, {}, o_10());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"BK_O2"}, {}, {"Ca"},
                        {"BK_O1"}, {}, o_21());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"BK_O3"}, {}, {"Ca"},
                        {"BK_O2"}, {}, o_32());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"BK_O4"}, {}, {"Ca"},
                        {"BK_O3"}, {}, o_43());

  statedef->addVDepSurfReac("smooth.__BOUNDARY__", {}, {"BK_C0"}, {}, {},
                            {"BK_O0"}, {}, [&](auto V) { return f_0(V); });
  statedef->addVDepSurfReac("smooth.__BOUNDARY__", {}, {"BK_C1"}, {}, {},
                            {"BK_O1"}, {}, [&](auto V) { return f_1(V); });
  statedef->addVDepSurfReac("smooth.__BOUNDARY__", {}, {"BK_C2"}, {}, {},
                            {"BK_O2"}, {}, [&](auto V) { return f_2(V); });
  statedef->addVDepSurfReac("smooth.__BOUNDARY__", {}, {"BK_C3"}, {}, {},
                            {"BK_O3"}, {}, [&](auto V) { return f_3(V); });
  statedef->addVDepSurfReac("smooth.__BOUNDARY__", {}, {"BK_C4"}, {}, {},
                            {"BK_O4"}, {}, [&](auto V) { return f_4(V); });

  statedef->addVDepSurfReac("smooth.__BOUNDARY__", {}, {"BK_O0"}, {}, {},
                            {"BK_C0"}, {}, [&](auto V) { return b_0(V); });
  statedef->addVDepSurfReac("smooth.__BOUNDARY__", {}, {"BK_O1"}, {}, {},
                            {"BK_C1"}, {}, [&](auto V) { return b_1(V); });
  statedef->addVDepSurfReac("smooth.__BOUNDARY__", {}, {"BK_O2"}, {}, {},
                            {"BK_C2"}, {}, [&](auto V) { return b_2(V); });
  statedef->addVDepSurfReac("smooth.__BOUNDARY__", {}, {"BK_O3"}, {}, {},
                            {"BK_C3"}, {}, [&](auto V) { return b_3(V); });
  statedef->addVDepSurfReac("smooth.__BOUNDARY__", {}, {"BK_O4"}, {}, {},
                            {"BK_C4"}, {}, [&](auto V) { return b_4(V); });

  // BK_G is single channel conductance
  // the same for all others currents in this model
  statedef->addOhmicCurrent("BK_O0SmoothCurr", "smooth_memb", "BKchan",
                            model::species_name("BK_O0"), BK_G(), BK_rev());
  statedef->addOhmicCurrent("BK_O1SmoothCurr", "smooth_memb", "BKchan",
                            model::species_name("BK_O1"), BK_G(), BK_rev());
  statedef->addOhmicCurrent("BK_O2SmoothCurr", "smooth_memb", "BKchan",
                            model::species_name("BK_O2"), BK_G(), BK_rev());
  statedef->addOhmicCurrent("BK_O3SmoothCurr", "smooth_memb", "BKchan",
                            model::species_name("BK_O3"), BK_G(), BK_rev());
  statedef->addOhmicCurrent("BK_O4SmoothCurr", "smooth_memb", "BKchan",
                            model::species_name("BK_O4"), BK_G(), BK_rev());

  // SK channel
  statedef->addChannel("smooth_memb", "SKchan",
                       {"SK_C1", "SK_C2", "SK_C3", "SK_C4", "SK_O1", "SK_O2"});
  statedef->addSurfReac("smooth.__BOUNDARY__", {"Ca"}, {"SK_C1"}, {}, {},
                        {"SK_C2"}, {}, dirc2_t());
  statedef->addSurfReac("smooth.__BOUNDARY__", {"Ca"}, {"SK_C2"}, {}, {},
                        {"SK_C3"}, {}, dirc3_t());
  statedef->addSurfReac("smooth.__BOUNDARY__", {"Ca"}, {"SK_C3"}, {}, {},
                        {"SK_C4"}, {}, dirc4_t());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"SK_C2"}, {}, {"Ca"},
                        {"SK_C1"}, {}, invc1_t());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"SK_C3"}, {}, {"Ca"},
                        {"SK_C2"}, {}, invc2_t());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"SK_C4"}, {}, {"Ca"},
                        {"SK_C3"}, {}, invc3_t());

  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"SK_C3"}, {}, {}, {"SK_O1"},
                        {}, diro1_t());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"SK_C4"}, {}, {}, {"SK_O2"},
                        {}, diro2_t());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"SK_O1"}, {}, {}, {"SK_C3"},
                        {}, invo1_t());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"SK_O2"}, {}, {}, {"SK_C4"},
                        {}, invo2_t());

  statedef->addOhmicCurrent("SK_O1SmoothCurr", "smooth_memb", "SKchan",
                            model::species_name("SK_O1"), SK_G(), SK_rev());
  statedef->addOhmicCurrent("SK_O2SmoothCurr", "smooth_memb", "SKchan",
                            model::species_name("SK_O2"), SK_G(), SK_rev());

  // AMPA channel
  statedef->addChannel(
      "smooth_memb", "AMPA",
      {"AMPA_C", "AMPA_C1", "AMPA_C2", "AMPA_D1", "AMPA_D2", "AMPA_O"});

  ampacc1 = statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"AMPA_C"}, {}, {},
                                  {"AMPA_C1"}, {}, 0.0);
  ampac1c2 = statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"AMPA_C1"}, {},
                                   {}, {"AMPA_C2"}, {}, 0.0);

  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"AMPA_C2"}, {}, {},
                        {"AMPA_O"}, {}, ro());

  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"AMPA_C1"}, {}, {},
                        {"AMPA_C"}, {}, ru1());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"AMPA_C2"}, {}, {},
                        {"AMPA_C1"}, {}, ru2());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"AMPA_O"}, {}, {},
                        {"AMPA_C2"}, {}, rc());

  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"AMPA_C1"}, {}, {},
                        {"AMPA_D1"}, {}, rd());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"AMPA_C2"}, {}, {},
                        {"AMPA_D2"}, {}, rd());

  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"AMPA_D1"}, {}, {},
                        {"AMPA_C1"}, {}, rr());
  statedef->addSurfReac("smooth.__BOUNDARY__", {}, {"AMPA_D2"}, {}, {},
                        {"AMPA_C2"}, {}, rr());

  statedef->addOhmicCurrent("AMPA_OSmoothCurr", "smooth_memb", "AMPA",
                            model::species_name("AMPA_O"), AMPA_G(),
                            AMPA_rev());

  // Leak channel
  statedef->addChannel("smooth_memb", "L", {"Leak"});
  statedef->addOhmicCurrent("LeakSmoothCurr", "smooth_memb", "L", model::species_name("Leak"),
                            L_G(true), L_rev());

  // setup spiny membrane
  statedef->addMembrane("spiny_memb", "spiny.__BOUNDARY__", memb_capac_spiny());
  statedef->addPatchSpecs(
      "spiny.__BOUNDARY__",
      {"CaP_m0",  "CaP_m1",  "CaP_m2",  "CaP_m3",  "BK_C0",  "BK_C1", "BK_C2",
       "BK_C3",   "BK_C4",   "BK_O0",   "BK_O1",   "BK_O2",  "BK_O3", "BK_O4",
       "SK_C1",   "SK_C2",   "SK_C3",   "SK_C4",   "SK_O1",  "SK_O2", "AMPA_C",
       "AMPA_C1", "AMPA_C2", "AMPA_D1", "AMPA_D2", "AMPA_O", "Leak"});

  statedef->addChannel("spiny_memb", "CaPchan",
                       {"CaP_m0", "CaP_m1", "CaP_m2", "CaP_m3"});
  statedef->addVDepSurfReac(
      "spiny.__BOUNDARY__", {}, {"CaP_m0"}, {}, {}, {"CaP_m1"}, {},
      [&](auto V) { return 1.0e3 * 3. * alpha_cap(V * 1.0e3) * Qt(); });
  statedef->addVDepSurfReac(
      "spiny.__BOUNDARY__", {}, {"CaP_m1"}, {}, {}, {"CaP_m2"}, {},
      [&](auto V) { return 1.0e3 * 2. * alpha_cap(V * 1.0e3) * Qt(); });
  statedef->addVDepSurfReac(
      "spiny.__BOUNDARY__", {}, {"CaP_m2"}, {}, {}, {"CaP_m3"}, {},
      [&](auto V) { return 1.0e3 * 1. * alpha_cap(V * 1.0e3) * Qt(); });
  statedef->addVDepSurfReac(
      "spiny.__BOUNDARY__", {}, {"CaP_m3"}, {}, {}, {"CaP_m2"}, {},
      [&](auto V) { return 1.0e3 * 3. * beta_cap(V * 1.0e3) * Qt(); });
  statedef->addVDepSurfReac(
      "spiny.__BOUNDARY__", {}, {"CaP_m2"}, {}, {}, {"CaP_m1"}, {},
      [&](auto V) { return 1.0e3 * 2. * beta_cap(V * 1.0e3) * Qt(); });
  statedef->addVDepSurfReac(
      "spiny.__BOUNDARY__", {}, {"CaP_m1"}, {}, {}, {"CaP_m0"}, {},
      [&](auto V) { return 1.0e3 * 1. * beta_cap(V * 1.0e3) * Qt(); });
  statedef->addGHKCurrentSurfReac(
      "CaPCurrSpiny", "spiny_memb", "CaPchan", "CaP_m3", "Ca", CaP_P(), 2, Ca_oconc());

  statedef->addChannel("spiny_memb", "BKchan",
                       {"BK_C0", "BK_C1", "BK_C2", "BK_C3", "BK_C4", "BK_O0",
                        "BK_O1", "BK_O2", "BK_O3", "BK_O4"});
  statedef->addSurfReac("spiny.__BOUNDARY__", {"Ca"}, {"BK_C0"}, {}, {},
                        {"BK_C1"}, {}, c_01());
  statedef->addSurfReac("spiny.__BOUNDARY__", {"Ca"}, {"BK_C1"}, {}, {},
                        {"BK_C2"}, {}, c_12());
  statedef->addSurfReac("spiny.__BOUNDARY__", {"Ca"}, {"BK_C2"}, {}, {},
                        {"BK_C3"}, {}, c_23());
  statedef->addSurfReac("spiny.__BOUNDARY__", {"Ca"}, {"BK_C3"}, {}, {},
                        {"BK_C4"}, {}, c_34());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"BK_C1"}, {}, {"Ca"},
                        {"BK_C0"}, {}, c_10());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"BK_C2"}, {}, {"Ca"},
                        {"BK_C1"}, {}, c_21());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"BK_C3"}, {}, {"Ca"},
                        {"BK_C2"}, {}, c_32());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"BK_C4"}, {}, {"Ca"},
                        {"BK_C3"}, {}, c_43());

  statedef->addSurfReac("spiny.__BOUNDARY__", {"Ca"}, {"BK_O0"}, {}, {},
                        {"BK_O1"}, {}, o_01());
  statedef->addSurfReac("spiny.__BOUNDARY__", {"Ca"}, {"BK_O1"}, {}, {},
                        {"BK_O2"}, {}, o_12());
  statedef->addSurfReac("spiny.__BOUNDARY__", {"Ca"}, {"BK_O2"}, {}, {},
                        {"BK_O3"}, {}, o_23());
  statedef->addSurfReac("spiny.__BOUNDARY__", {"Ca"}, {"BK_O3"}, {}, {},
                        {"BK_O4"}, {}, o_34());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"BK_O1"}, {}, {"Ca"},
                        {"BK_O0"}, {}, o_10());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"BK_O2"}, {}, {"Ca"},
                        {"BK_O1"}, {}, o_21());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"BK_O3"}, {}, {"Ca"},
                        {"BK_O2"}, {}, o_32());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"BK_O4"}, {}, {"Ca"},
                        {"BK_O3"}, {}, o_43());

  statedef->addVDepSurfReac("spiny.__BOUNDARY__", {}, {"BK_C0"}, {}, {},
                            {"BK_O0"}, {}, [&](auto V) { return f_0(V); });
  statedef->addVDepSurfReac("spiny.__BOUNDARY__", {}, {"BK_C1"}, {}, {},
                            {"BK_O1"}, {}, [&](auto V) { return f_1(V); });
  statedef->addVDepSurfReac("spiny.__BOUNDARY__", {}, {"BK_C2"}, {}, {},
                            {"BK_O2"}, {}, [&](auto V) { return f_2(V); });
  statedef->addVDepSurfReac("spiny.__BOUNDARY__", {}, {"BK_C3"}, {}, {},
                            {"BK_O3"}, {}, [&](auto V) { return f_3(V); });
  statedef->addVDepSurfReac("spiny.__BOUNDARY__", {}, {"BK_C4"}, {}, {},
                            {"BK_O4"}, {}, [&](auto V) { return f_4(V); });

  statedef->addVDepSurfReac("spiny.__BOUNDARY__", {}, {"BK_O0"}, {}, {},
                            {"BK_C0"}, {}, [&](auto V) { return b_0(V); });
  statedef->addVDepSurfReac("spiny.__BOUNDARY__", {}, {"BK_O1"}, {}, {},
                            {"BK_C1"}, {}, [&](auto V) { return b_1(V); });
  statedef->addVDepSurfReac("spiny.__BOUNDARY__", {}, {"BK_O2"}, {}, {},
                            {"BK_C2"}, {}, [&](auto V) { return b_2(V); });
  statedef->addVDepSurfReac("spiny.__BOUNDARY__", {}, {"BK_O3"}, {}, {},
                            {"BK_C3"}, {}, [&](auto V) { return b_3(V); });
  statedef->addVDepSurfReac("spiny.__BOUNDARY__", {}, {"BK_O4"}, {}, {},
                            {"BK_C4"}, {}, [&](auto V) { return b_4(V); });

  // BK_G is single channel conductance
  // the same for all others currents in this model
  statedef->addOhmicCurrent("BK_O0SpinyCurr", "spiny_memb", "BKchan",
                            model::species_name("BK_O0"), BK_G(), BK_rev());
  statedef->addOhmicCurrent("BK_O1SpinyCurr", "spiny_memb", "BKchan",
                            model::species_name("BK_O1"), BK_G(), BK_rev());
  statedef->addOhmicCurrent("BK_O2SpinyCurr", "spiny_memb", "BKchan",
                            model::species_name("BK_O2"), BK_G(), BK_rev());
  statedef->addOhmicCurrent("BK_O3SpinyCurr", "spiny_memb", "BKchan",
                            model::species_name("BK_O3"), BK_G(), BK_rev());
  statedef->addOhmicCurrent("BK_O4SpinyCurr", "spiny_memb", "BKchan",
                            model::species_name("BK_O4"), BK_G(), BK_rev());

  // SK channel
  statedef->addChannel("spiny_memb", "SKchan",
                       {"SK_C1", "SK_C2", "SK_C3", "SK_C4", "SK_O1", "SK_O2"});
  statedef->addSurfReac("spiny.__BOUNDARY__", {"Ca"}, {"SK_C1"}, {}, {},
                        {"SK_C2"}, {}, dirc2_t());
  statedef->addSurfReac("spiny.__BOUNDARY__", {"Ca"}, {"SK_C2"}, {}, {},
                        {"SK_C3"}, {}, dirc3_t());
  statedef->addSurfReac("spiny.__BOUNDARY__", {"Ca"}, {"SK_C3"}, {}, {},
                        {"SK_C4"}, {}, dirc4_t());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"SK_C2"}, {}, {"Ca"},
                        {"SK_C1"}, {}, invc1_t());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"SK_C3"}, {}, {"Ca"},
                        {"SK_C2"}, {}, invc2_t());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"SK_C4"}, {}, {"Ca"},
                        {"SK_C3"}, {}, invc3_t());

  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"SK_C3"}, {}, {}, {"SK_O1"},
                        {}, diro1_t());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"SK_C4"}, {}, {}, {"SK_O2"},
                        {}, diro2_t());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"SK_O1"}, {}, {}, {"SK_C3"},
                        {}, invo1_t());
  statedef->addSurfReac("spiny.__BOUNDARY__", {}, {"SK_O2"}, {}, {}, {"SK_C4"},
                        {}, invo2_t());

  statedef->addOhmicCurrent("SK_O1SpinyCurr", "spiny_memb", "SKchan",
                            model::species_name("SK_O1"), SK_G(), SK_rev());
  statedef->addOhmicCurrent("SK_O2SpinyCurr", "spiny_memb", "SKchan",
                            model::species_name("SK_O2"), SK_G(), SK_rev());

  // Leak channel
  statedef->addChannel("spiny_memb", "L", {"Leak"});
  statedef->addOhmicCurrent("LeakSpinyCurr", "spiny_memb", "L", model::species_name("Leak"),
                            L_G(true), L_rev());
}

CaBurstIntegrationTestSimdef::CaBurstIntegrationTestSimdef(const steps::model::Model& model,
                                                           const steps::dist::DistMesh& mesh)
    : Simdef(model, mesh) {
    statedef->setTemp(273.0 + 30.0);
    // add compartment
    statedef->addComp("comp1");
    statedef->addCompartmentConductivity("comp1", 1.0);
    statedef->addCompSpecs("comp1", {"A", "B", "C"});

    statedef->addComp("comp2");
    statedef->addCompartmentConductivity("comp2", 1.0);
    statedef->addCompSpecs("comp2", {"C"});

    // add patch
    statedef->addPatch("patch1", "comp1");
    statedef->addPatch("patchInBetween", "comp1", boost::optional<model::compartment_id>("comp2"));
    // add membrane
    statedef->addMembrane("memb1", "patch1", 1);

    statedef->addMembrane("membInBetween", "patchInBetween", 1);
    // add compartment reactions
    statedef->addCompReac("comp1", {"C"}, {"A"}, 1e3);
    statedef->addCompReac("comp1", {"A"}, {"B"}, 1e3);
    statedef->addCompReac("comp1", {"B"}, {"C"}, 1e3);

    statedef->addCompReac("comp2", {"C"}, {}, 1e3);

    // add compartment diffusion
    statedef->addCompDiff("comp1", "C", 1e-7);
    statedef->addCompDiff("comp2", "C", 1e-7);

    // add GHK surface reaction
    statedef->addChannel("memb1", "C_chan", {"C_chan_0"});
    statedef->addPatchSpecs("patch1", {"C_chan_0"});
    statedef->addGHKCurrentSurfReac("C_GHKcurr",
                                    "memb1",
                                    "C_chan",
                                    "C_chan_0",
                                    "C",
                                    1.0e-20,
                                    2,
                                    10 * CaBurstIntegrationTestSimdef::compConc());

    statedef->addChannel("membInBetween", "C_chan_inBetween", {"C_chan_0_inBetween"});
    statedef->addPatchSpecs("patchInBetween", {"C_chan_0_inBetween"});
    statedef->addGHKCurrentSurfReac("C_GHKcurr_inBetween",
                                    "membInBetween",
                                    "C_chan_inBetween",
                                    "C_chan_0_inBetween",
                                    "C",
                                    1.0e-20,
                                    2);

    // add open channels on the patch
    patchCounts.emplace_back("patch1", "C_chan_0", 10);
    patchCounts.emplace_back("patchInBetween", "C_chan_0_inBetween", 10);

    // Leak channel
    statedef->addChannel("memb1", "L_chan", {"L_chan_0"});
    statedef->addOhmicCurrent(
        "L_OHMcurr", "memb1", "L_chan", model::species_name("L_chan_0"), 1 / 1e9, potential() * 1);

    //   add open channels on the patch
    patchCounts.emplace_back("patch1", "L_chan_0", 1);

    // add surface reaction
    statedef->addPatchSpecs("patch1", {"C"});
    statedef->addSurfReac("patch1", {}, {"C"}, {}, {"C"}, {}, {}, 1e2);
    //    // add surface reaction relevant patch count
    patchCounts.emplace_back("patch1", "C", 0);

    // add voltage dependent surface reac
    statedef->addVDepSurfReac("patch1", {"C"}, {}, {}, {}, {"C"}, {}, [&](auto V) {
        double y1 = 1e1;
        double y2 = 1e2;
        double x1 = -0.065;
        double x2 = -0.0645;
        double a = (y1 - y2) / (x1 - x2);
        double q = y1 - a * x1;
        return a * V + q;
    });
}

Rallpack3Simdef::Rallpack3Simdef(const steps::model::Model &model,
                                 const steps::dist::DistMesh &mesh)
    : Simdef(model, mesh) {
  statedef->addComp("__MESH__");
  statedef->addPatch("z_min", "__MESH__");
  statedef->addPatch("memb", "__MESH__");
  auto corr_fac_area = 4 * 866e-9 * 1.0e-3 / (3.1415926535 * 1.0e-9);
  auto corr_fac_vol = 866e-9 * 866e-9 * 1.0e-3 / (3.1415926535 * 0.25 * 1e-15);
  statedef->addPatchSpecs("memb", {"K_n0", "K_n1", "K_n2", "K_n3", "K_n4",
                                   "Na_m0h1", "Na_m1h1", "Na_m2h1", "Na_m3h1",
                                   "Na_m0h0", "Na_m1h0", "Na_m2h0", "Na_m3h0"});
  statedef->addMembrane("membrane", "memb", 0.01 / corr_fac_area);
  statedef->addChannel("membrane", "K",
                       {"K_n0", "K_n1", "K_n2", "K_n3", "K_n4"});
  statedef->addChannel("membrane", "Na",
                       {"Na_m0h1", "Na_m1h1", "Na_m2h1", "Na_m3h1", "Na_m0h0",
                        "Na_m1h0", "Na_m2h0", "Na_m3h0"});
  statedef->addChannel("membrane", "Leak_Channel", {});
  statedef->addMembrane("current_injection", "z_min", 0.0);
  statedef->addCompartmentConductivity("__MESH__", 1.0 / Ra() / corr_fac_vol);
  statedef->addVDepSurfReac(
      "memb", {}, {"K_n0"}, {}, {}, {"K_n1"}, {},
      [&](auto V) { return 1.0e3 * 4.0 * _a_n(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"K_n1"}, {}, {}, {"K_n2"}, {},
      [&](auto V) { return 1.0e3 * 3.0 * _a_n(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"K_n2"}, {}, {}, {"K_n3"}, {},
      [&](auto V) { return 1.0e3 * 2.0 * _a_n(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"K_n3"}, {}, {}, {"K_n4"}, {},
      [&](auto V) { return 1.0e3 * 1.0 * _a_n(V * 1.0e3); });

  statedef->addVDepSurfReac(
      "memb", {}, {"K_n4"}, {}, {}, {"K_n3"}, {},
      [&](auto V) { return 1.0e3 * 4.0 * _b_n(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"K_n3"}, {}, {}, {"K_n2"}, {},
      [&](auto V) { return 1.0e3 * 3.0 * _b_n(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"K_n2"}, {}, {}, {"K_n1"}, {},
      [&](auto V) { return 1.0e3 * 2.0 * _b_n(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"K_n1"}, {}, {}, {"K_n0"}, {},
      [&](auto V) { return 1.0e3 * 1.0 * _b_n(V * 1.0e3); });

  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m0h1"}, {}, {}, {"Na_m1h1"}, {},
      [&](auto V) { return 1.0e3 * 3.0 * _a_m(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m1h1"}, {}, {}, {"Na_m2h1"}, {},
      [&](auto V) { return 1.0e3 * 2.0 * _a_m(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m2h1"}, {}, {}, {"Na_m3h1"}, {},
      [&](auto V) { return 1.0e3 * 1.0 * _a_m(V * 1.0e3); });

  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m3h1"}, {}, {}, {"Na_m2h1"}, {},
      [&](auto V) { return 1.0e3 * 3.0 * _b_m(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m2h1"}, {}, {}, {"Na_m1h1"}, {},
      [&](auto V) { return 1.0e3 * 2.0 * _b_m(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m1h1"}, {}, {}, {"Na_m0h1"}, {},
      [&](auto V) { return 1.0e3 * 1.0 * _b_m(V * 1.0e3); });

  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m0h0"}, {}, {}, {"Na_m1h0"}, {},
      [&](auto V) { return 1.0e3 * 3.0 * _a_m(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m1h0"}, {}, {}, {"Na_m2h0"}, {},
      [&](auto V) { return 1.0e3 * 2.0 * _a_m(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m2h0"}, {}, {}, {"Na_m3h0"}, {},
      [&](auto V) { return 1.0e3 * 1.0 * _a_m(V * 1.0e3); });

  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m3h0"}, {}, {}, {"Na_m2h0"}, {},
      [&](auto V) { return 1.0e3 * 3.0 * _b_m(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m2h0"}, {}, {}, {"Na_m1h0"}, {},
      [&](auto V) { return 1.0e3 * 2.0 * _b_m(V * 1.0e3); });
  statedef->addVDepSurfReac(
      "memb", {}, {"Na_m1h0"}, {}, {}, {"Na_m0h0"}, {},
      [&](auto V) { return 1.0e3 * 1.0 * _b_m(V * 1.0e3); });

  statedef->addVDepSurfReac("memb", {}, {"Na_m0h1"}, {}, {}, {"Na_m0h0"}, {},
                            [&](auto V) { return 1.0e3 * _a_h(V * 1.0e3); });
  statedef->addVDepSurfReac("memb", {}, {"Na_m1h1"}, {}, {}, {"Na_m1h0"}, {},
                            [&](auto V) { return 1.0e3 * _a_h(V * 1.0e3); });
  statedef->addVDepSurfReac("memb", {}, {"Na_m2h1"}, {}, {}, {"Na_m2h0"}, {},
                            [&](auto V) { return 1.0e3 * _a_h(V * 1.0e3); });
  statedef->addVDepSurfReac("memb", {}, {"Na_m3h1"}, {}, {}, {"Na_m3h0"}, {},
                            [&](auto V) { return 1.0e3 * _a_h(V * 1.0e3); });

  statedef->addVDepSurfReac("memb", {}, {"Na_m0h0"}, {}, {}, {"Na_m0h1"}, {},
                            [&](auto V) { return 1.0e3 * _b_h(V * 1.0e3); });
  statedef->addVDepSurfReac("memb", {}, {"Na_m1h0"}, {}, {}, {"Na_m1h1"}, {},
                            [&](auto V) { return 1.0e3 * _b_h(V * 1.0e3); });
  statedef->addVDepSurfReac("memb", {}, {"Na_m2h0"}, {}, {}, {"Na_m2h1"}, {},
                            [&](auto V) { return 1.0e3 * _b_h(V * 1.0e3); });
  statedef->addVDepSurfReac("memb", {}, {"Na_m3h0"}, {}, {}, {"Na_m3h1"}, {},
                            [&](auto V) { return 1.0e3 * _b_h(V * 1.0e3); });
  statedef->addOhmicCurrent("KCurr", "membrane", "K", model::species_name("K_n4"),
                            K_G() / K_ro(), K_rev());
  statedef->addOhmicCurrent("NaCurr", "membrane", "Na", model::species_name("Na_m3h0"),
                            Na_G() / Na_ro(), Na_rev());
  statedef->addOhmicCurrent("LeakCurr", "membrane", "Leak_Channel", boost::none,
                            L_G() / corr_fac_area, leak_rev());
  statedef->setStimulus("current_injection", Iinj());
}

} // namespace dist
} // namespace steps
