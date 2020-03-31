#include <string>
#include <vector>

#include "opsplit/compdef.hpp"
#include "opsplit/diffdef.hpp"
#include "opsplit/patchdef.hpp"
#include "opsplit/reacdef.hpp"
#include "opsplit/sreacdef.hpp"
#include "opsplit/test/simdef.hpp"

namespace zee {

Simdef::Simdef()
    : statedef(new Statedef()) {}

SimpleSimdef::SimpleSimdef()
    : Simdef() {
    statedef->addComp("comp");

    const std::vector<model::specie_name> species = {
        "A", "B", "C", "D", "E", "F", "G", "H", "I", "J"};
    statedef->addCompSpecs("comp", species);
    const std::vector<model::specie_name> R1_lhs = {"A", "B"};
    const std::vector<model::specie_name> R1_rhs = {"C"};
    statedef->addCompReac("comp", R1_lhs, R1_rhs, 1000.0e6);

    const std::vector<model::specie_name> R2_lhs = {"C"};
    const std::vector<model::specie_name> R2_rhs = {"A", "B"};
    statedef->addCompReac("comp", R2_lhs, R2_rhs, 100);

    const std::vector<model::specie_name> R3_lhs = {"C", "D"};
    const std::vector<model::specie_name> R3_rhs = {"E"};
    statedef->addCompReac("comp", R3_lhs, R3_rhs, 100e6);

    const std::vector<model::specie_name> R4_lhs = {"E"};
    const std::vector<model::specie_name> R4_rhs = {"C", "D"};
    statedef->addCompReac("comp", R4_lhs, R4_rhs, 10);

    const std::vector<model::specie_name> R5_lhs = {"F", "G"};
    const std::vector<model::specie_name> R5_rhs = {"H"};
    statedef->addCompReac("comp", R5_lhs, R5_rhs, 10e6);

    const std::vector<model::specie_name> R6_lhs = {"H"};
    const std::vector<model::specie_name> R6_rhs = {"F", "G"};
    statedef->addCompReac("comp", R6_lhs, R6_rhs, 1);

    const std::vector<model::specie_name> R7_lhs = {"H", "I"};
    const std::vector<model::specie_name> R7_rhs = {"J"};
    statedef->addCompReac("comp", R7_lhs, R7_rhs, 1e6);

    const std::vector<model::specie_name> R8_lhs = {"J"};
    const std::vector<model::specie_name> R8_rhs = {"H", "I"};
    statedef->addCompReac("comp", R8_lhs, R8_rhs, 0.1 * 10);

    statedef->addCompDiff("comp", "A", 100e-12);
    statedef->addCompDiff("comp", "B", 90e-12);
    statedef->addCompDiff("comp", "C", 80e-12);
    statedef->addCompDiff("comp", "D", 70e-12);
    statedef->addCompDiff("comp", "E", 60e-12);
    statedef->addCompDiff("comp", "F", 50e-12);
    statedef->addCompDiff("comp", "G", 40e-12);
    statedef->addCompDiff("comp", "H", 30e-12);
    statedef->addCompDiff("comp", "I", 20e-12);
    statedef->addCompDiff("comp", "J", 10e-12);

    const PetscScalar N0A = 1000.0;
    const PetscScalar N0B = 2000.0;
    const PetscScalar N0C = 3000.0;
    const PetscScalar N0D = 4000.0;
    const PetscScalar N0E = 5000.0;
    const PetscScalar N0F = 6000.0;
    const PetscScalar N0G = 7000.0;
    const PetscScalar N0H = 8000.0;
    const PetscScalar N0I = 9000.0;
    const PetscScalar N0J = 10000.0;

    compartmentCounts.emplace_back("comp", "A", N0A);
    compartmentCounts.emplace_back("comp", "B", N0B);
    compartmentCounts.emplace_back("comp", "C", N0C);
    compartmentCounts.emplace_back("comp", "D", N0D);
    compartmentCounts.emplace_back("comp", "E", N0E);
    compartmentCounts.emplace_back("comp", "F", N0F);
    compartmentCounts.emplace_back("comp", "G", N0G);
    compartmentCounts.emplace_back("comp", "H", N0H);
    compartmentCounts.emplace_back("comp", "I", N0I);
    compartmentCounts.emplace_back("comp", "J", N0J);
}

MultiCompartmentSimdef::MultiCompartmentSimdef() {
    PetscScalar const r1_kcst = 1000000.0e6;
    PetscScalar const r2_kcst = 0.0005e6;
    PetscScalar dcst = 1e-8;

    // PetscOptionsGetScalar(nullptr, nullptr, "-d", &dcst, nullptr);

    statedef->addComp("Left");
    const std::vector<model::specie_name> left_species = {"A", "B", "C"};
    statedef->addCompSpecs("Left", left_species);

    // add reaction1
    const std::vector<model::specie_name> left_r1_lhs = {"A", "B"};
    const std::vector<model::specie_name> left_r1_rhs = {"C"};
    statedef->addCompReac("Left", left_r1_lhs, left_r1_rhs, r1_kcst);

    // add reaction2
    const std::vector<model::specie_name> left_r2_lhs = {"C"};
    const std::vector<model::specie_name> left_r2_rhs = {"A", "B"};
    statedef->addCompReac("Left", left_r2_lhs, left_r2_rhs, r2_kcst);

    statedef->addCompDiff("Left", "A", dcst);
    statedef->addCompDiff("Left", "B", dcst);
    statedef->addCompDiff("Left", "C", dcst);

    // Right compartment

    statedef->addComp("Right");
    const std::vector<model::specie_name> right_species = {"C", "D", "E"};
    statedef->addCompSpecs("Right", right_species);

    // add reaction1
    const std::vector<model::specie_name> right_r1_lhs = {"C", "D"};
    const std::vector<model::specie_name> right_r1_rhs = {"E"};
    statedef->addCompReac("Right", right_r1_lhs, right_r1_rhs, r1_kcst);

    // add reaction2
    const std::vector<model::specie_name> right_r2_lhs = {"E"};
    const std::vector<model::specie_name> right_r2_rhs = {"C", "D"};
    statedef->addCompReac("Right", right_r2_lhs, right_r2_rhs, r2_kcst);
    statedef->addCompDiff("Right", "C", dcst);
    statedef->addCompDiff("Right", "D", dcst);
    statedef->addCompDiff("Right", "E", dcst);

    compartmentCounts.emplace_back("Left", "A", 10);
    compartmentCounts.emplace_back("Left", "B", 10);
    compartmentCounts.emplace_back("Right", "E", 10);
}

SReacUnitTestSimdef::SReacUnitTestSimdef() {
    statedef->addComp("comp_i");
    statedef->addComp("comp_o");
    statedef->addPatch("patch", "comp_i", model::compartment_id("comp_o"));

    statedef->addCompSpecs("comp_i", {"A", "B"});
    statedef->addCompSpecs("comp_o", {"C"});
    statedef->addPatchSpecs("patch", {"D"});

    const std::vector<model::specie_name> R1_lhs = {"A", "B"};
    statedef->addCompReac("comp_i", R1_lhs, {}, 1e8);

    const std::vector<model::specie_name> R2_lhs_p = {"D", "D", "D"};
    const std::vector<model::specie_name> R2_lhs_o = {"C", "C"};
    const std::vector<model::specie_name> R2_rhs_i = {"B", "B", "B", "B"};
    statedef->addSurfReac("patch", {}, R2_lhs_p, R2_lhs_o, R2_rhs_i, {}, {}, 1e21);

    const PetscScalar N0A = 4000.0;
    const PetscScalar N0B = 0.0;
    const PetscScalar N0C = 2000.0;
    const PetscScalar N0D = 6000.0;

    compartmentCounts.emplace_back("comp_i", "A", N0A);
    compartmentCounts.emplace_back("comp_i", "B", N0B);
    compartmentCounts.emplace_back("comp_o", "C", N0C);
    patchCounts.emplace_back("patch", "D", N0D);
}

ValidationSimdef::ValidationSimdef() {
    statedef->addComp("comp1");
    const std::vector<model::specie_name> comp_species = {
        "A_foi",
        "A_for",
        "B_for",
        "A_soA2",
        "C_soA2",
        "A_soAA",
        "B_soAA",
        "C_soAA",
        "A_soAB",
        "B_soAB",
        "C_soAB",
        "A_toA3",
        "C_toA3",
        "A_toA2B",
        "B_toA2B",
        "C_toA2B",
    };
    statedef->addCompSpecs("comp1", comp_species);
    // First order irreversible

    const std::vector<model::specie_name> R1_foi_lhs = {"A_foi"};
    const std::vector<model::specie_name> R1_foi_rhs = {};
    statedef->addCompReac("comp1", R1_foi_lhs, R1_foi_rhs, KCST_foi());

    // First order reversible
    const std::vector<model::specie_name> R1_for_lhs = {"A_for"};
    const std::vector<model::specie_name> R1_for_rhs = {"B_for"};
    statedef->addCompReac("comp1", R1_for_lhs, R1_for_rhs, KCST_f_for());

    const std::vector<model::specie_name> R2_for_lhs = {"B_for"};
    const std::vector<model::specie_name> R2_for_rhs = {"A_for"};
    statedef->addCompReac("comp1", R2_for_lhs, R2_for_rhs, KCST_b_for());

    // Second order irreversible A2
    const std::vector<model::specie_name> R1_soA2_lhs = {"A_soA2", "A_soA2"};
    const std::vector<model::specie_name> R1_soA2_rhs = {"C_soA2"};
    statedef->addCompReac("comp1", R1_soA2_lhs, R1_soA2_rhs, KCST_soA2());

    // Second order irreversible AA
    const std::vector<model::specie_name> R1_soAA_lhs = {"A_soAA", "B_soAA"};
    const std::vector<model::specie_name> R1_soAA_rhs = {"C_soAA"};
    statedef->addCompReac("comp1", R1_soAA_lhs, R1_soAA_rhs, KCST_soAA());

    // Second order irreversible AB
    const std::vector<model::specie_name> R1_soAB_lhs = {"A_soAB", "B_soAB"};
    const std::vector<model::specie_name> R1_soAB_rhs = {"C_soAB"};
    statedef->addCompReac("comp1", R1_soAB_lhs, R1_soAB_rhs, KCST_soAB());

    // Third order irreversible A3
    const std::vector<model::specie_name> R1_toA3_lhs = {"A_toA3", "A_toA3", "A_toA3"};
    const std::vector<model::specie_name> R1_toA3_rhs = {"C_toA3"};
    statedef->addCompReac("comp1", R1_toA3_lhs, R1_toA3_rhs, KCST_toA3());

    // Third order irreversible A2B
    const std::vector<model::specie_name> R1_toA2B_lhs = {"A_toA2B", "A_toA2B", "B_toA2B"};
    const std::vector<model::specie_name> R1_toA2B_rhs = {"C_toA2B"};
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

ValidationSimdef::compartment_counts_t ValidationSimdef::getfoi_CompartementCounts(
    bool zero) const {
    compartment_counts_t counts;

    counts.emplace_back("comp1", "A_foi", zero ? 0 : N_foi());
    return counts;
}

ValidationSimdef::compartment_counts_t ValidationSimdef::getfor_CompartementCounts(
    bool zero) const {
    compartment_counts_t counts;

    counts.emplace_back("comp1", "A_for", zero ? 0 : COUNT_for());
    counts.emplace_back("comp1", "B_for", 0.0);

    return counts;
}

ValidationSimdef::compartment_concs_t ValidationSimdef::getsoA2_CompartementConc(bool zero) const {
    compartment_concs_t counts;

    counts.emplace_back("comp1", "A_soA2", zero ? 0 : CONCA_soA2());

    return counts;
}

ValidationSimdef::compartment_concs_t ValidationSimdef::getsoAA_CompartementConc(bool zero) const {
    compartment_concs_t counts;

    counts.emplace_back("comp1", "A_soAA", zero ? 0 : CONCA_soAA());
    counts.emplace_back("comp1", "B_soAA", zero ? 0 : CONCB_soAA());

    return counts;
}

ValidationSimdef::compartment_concs_t ValidationSimdef::getsoAB_CompartementConc(bool zero) const {
    compartment_concs_t counts;

    counts.emplace_back("comp1", "A_soAB", zero ? 0 : CONCA_soAB());
    counts.emplace_back("comp1", "B_soAB", zero ? 0 : CONCB_soAB());

    return counts;
}

ValidationSimdef::compartment_concs_t ValidationSimdef::gettoA3_CompartementConc(bool zero) const {
    compartment_concs_t counts;

    counts.emplace_back("comp1", "A_toA3", zero ? 0 : CONCA_toA3());

    return counts;
}

ValidationSimdef::compartment_concs_t ValidationSimdef::gettoA2B_CompartementConc(bool zero) const {
    compartment_concs_t counts;

    counts.emplace_back("comp1", "A_toA2B", zero ? 0 : CONCA_toA2B());
    counts.emplace_back("comp1", "B_toA2B", zero ? 0 : CONCB_toA2B());

    return counts;
}


SurfaceReactionsValidationSimdef::SurfaceReactionsValidationSimdef() {
    statedef->addComp("comp1");
    const std::vector<model::specie_name> comp_species = {
        "A_foi",
        "A_for",
        "B_for",
        "A_soA2",
        "C_soA2",
        "A_soAA",
        "B_soAA",
        "C_soAA",
        "A_soAB",
        "B_soAB",
        "C_soAB",
        "A_toA3",
        "C_toA3",
        "A_toA2B",
        "B_toA2B",
        "C_toA2B",
    };
    statedef->addCompSpecs("comp1", comp_species);
    statedef->addPatch("patch1", "comp1");
    const std::vector<model::specie_name> patch_species = {
        "A_so2d",
        "B_so2d",
        "C_so2d",
    };
    statedef->addPatchSpecs("patch1", patch_species);
    // First order irreversible
    const std::vector<model::specie_name> R1_foi_lhs = {"A_foi"};
    const std::vector<model::specie_name> R1_foi_rhs = {};
    statedef->addCompReac("comp1", R1_foi_lhs, R1_foi_rhs, KCST_foi());

    // First order reversible
    const std::vector<model::specie_name> R1_for_lhs = {"A_for"};
    const std::vector<model::specie_name> R1_for_rhs = {"B_for"};
    statedef->addCompReac("comp1", R1_for_lhs, R1_for_rhs, KCST_f_for());

    const std::vector<model::specie_name> R2_for_lhs = {"B_for"};
    const std::vector<model::specie_name> R2_for_rhs = {"A_for"};
    statedef->addCompReac("comp1", R2_for_lhs, R2_for_rhs, KCST_b_for());

    // Second order irreversible A2
    const std::vector<model::specie_name> R1_soA2_lhs = {"A_soA2", "A_soA2"};
    const std::vector<model::specie_name> R1_soA2_rhs = {"C_soA2"};
    statedef->addCompReac("comp1", R1_soA2_lhs, R1_soA2_rhs, KCST_soA2());

    // Second order irreversible AA
    const std::vector<model::specie_name> R1_soAA_lhs = {"A_soAA", "B_soAA"};
    const std::vector<model::specie_name> R1_soAA_rhs = {"C_soAA"};
    statedef->addCompReac("comp1", R1_soAA_lhs, R1_soAA_rhs, KCST_soAA());

    // Second order irreversible AB
    const std::vector<model::specie_name> R1_soAB_lhs = {"A_soAB", "B_soAB"};
    const std::vector<model::specie_name> R1_soAB_rhs = {"C_soAB"};
    statedef->addCompReac("comp1", R1_soAB_lhs, R1_soAB_rhs, KCST_soAB());

    // Third order irreversible A3
    const std::vector<model::specie_name> R1_toA3_lhs = {"A_toA3", "A_toA3", "A_toA3"};
    const std::vector<model::specie_name> R1_toA3_rhs = {"C_toA3"};
    statedef->addCompReac("comp1", R1_toA3_lhs, R1_toA3_rhs, KCST_toA3());

    // Third order irreversible A2B
    const std::vector<model::specie_name> R1_toA2B_lhs = {"A_toA2B", "A_toA2B", "B_toA2B"};
    const std::vector<model::specie_name> R1_toA2B_rhs = {"C_toA2B"};
    statedef->addCompReac("comp1", R1_toA2B_lhs, R1_toA2B_rhs, KCST_toA2B());

    // Second order irreversible 2D
    const std::vector<model::specie_name> R1_so2d_lhs = {"A_so2d", "B_so2d"};
    const std::vector<model::specie_name> R1_so2d_rhs = {"C_so2d"};
    statedef->addSurfReac("patch1", {}, R1_so2d_lhs, {}, {}, R1_so2d_rhs, {}, KCST_so2d());

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
    compartment_counts_t counts;

    counts.emplace_back("comp1", "A_foi", N_foi());

    return counts;
}

SurfaceReactionsValidationSimdef::compartment_counts_t
SurfaceReactionsValidationSimdef::getfor_CompartementCounts() const {
    compartment_counts_t counts;

    counts.emplace_back("comp1", "A_for", COUNT_for());
    counts.emplace_back("comp1", "B_for", 0.0);

    return counts;
}

SurfaceReactionsValidationSimdef::compartment_concs_t
SurfaceReactionsValidationSimdef::getsoA2_CompartementConc() const {
    compartment_concs_t counts;

    counts.emplace_back("comp1", "A_soA2", CONCA_soA2());

    return counts;
}

SurfaceReactionsValidationSimdef::compartment_concs_t
SurfaceReactionsValidationSimdef::getsoAA_CompartementConc() const {
    compartment_concs_t counts;

    counts.emplace_back("comp1", "A_soAA", CONCA_soAA());
    counts.emplace_back("comp1", "B_soAA", CONCB_soAA());

    return counts;
}

SurfaceReactionsValidationSimdef::compartment_concs_t
SurfaceReactionsValidationSimdef::getsoAB_CompartementConc() const {
    compartment_concs_t counts;

    counts.emplace_back("comp1", "A_soAB", CONCA_soAB());
    counts.emplace_back("comp1", "B_soAB", CONCB_soAB());

    return counts;
}

SurfaceReactionsValidationSimdef::compartment_concs_t
SurfaceReactionsValidationSimdef::gettoA3_CompartementConc() const {
    compartment_concs_t counts;

    counts.emplace_back("comp1", "A_toA3", CONCA_toA3());

    return counts;
}

SurfaceReactionsValidationSimdef::compartment_concs_t
SurfaceReactionsValidationSimdef::gettoA2B_CompartementConc() const {
    compartment_concs_t counts;

    counts.emplace_back("comp1", "A_toA2B", CONCA_toA2B());
    counts.emplace_back("comp1", "B_toA2B", CONCB_toA2B());

    return counts;
}

SurfaceReactionsValidationSimdef::patch_counts_t
SurfaceReactionsValidationSimdef::getso2d_PatchCount() const {
    patch_counts_t counts;

    counts.emplace_back("patch1", "A_so2d", COUNTA_so2d());
    counts.emplace_back("patch1", "B_so2d", COUNTB_so2d());

    return counts;
}

CaBurstSimdef::CaBurstSimdef() {
    statedef->addComp("__MESH__");
    const std::vector<model::specie_name> comp_species = {"Ca",
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

    const std::vector<model::specie_name> iCBsf1_f_lhs = {"Ca", "iCBsf"};
    const std::vector<model::specie_name> iCBsf1_f_rhs = {"iCBsCa"};
    statedef->addCompReac("__MESH__", iCBsf1_f_lhs, iCBsf1_f_rhs, iCBsf1_f_kcst());
    statedef->addCompReac("__MESH__", iCBsf1_f_rhs, iCBsf1_f_lhs, iCBsf1_b_kcst());

    const std::vector<model::specie_name> iCBsCa_f_lhs = {"Ca", "iCBsCa"};
    const std::vector<model::specie_name> iCBsCa_f_rhs = {"iCBCaCa"};
    statedef->addCompReac("__MESH__", iCBsCa_f_lhs, iCBsCa_f_rhs, iCBsCa_f_kcst());
    statedef->addCompReac("__MESH__", iCBsCa_f_rhs, iCBsCa_f_lhs, iCBsCa_b_kcst());

    const std::vector<model::specie_name> iCBsf2_f_lhs = {"Ca", "iCBsf"};
    const std::vector<model::specie_name> iCBsf2_f_rhs = {"iCBCaf"};
    statedef->addCompReac("__MESH__", iCBsf2_f_lhs, iCBsf2_f_rhs, iCBsf2_f_kcst());
    statedef->addCompReac("__MESH__", iCBsf2_f_rhs, iCBsf2_f_lhs, iCBsf2_b_kcst());

    const std::vector<model::specie_name> iCBCaf_f_lhs = {"Ca", "iCBCaf"};
    const std::vector<model::specie_name> iCBCaf_f_rhs = {"iCBCaCa"};
    statedef->addCompReac("__MESH__", iCBCaf_f_lhs, iCBCaf_f_rhs, iCBCaf_f_kcst());
    statedef->addCompReac("__MESH__", iCBCaf_f_rhs, iCBCaf_f_lhs, iCBCaf_b_kcst());

    const std::vector<model::specie_name> CBsf1_f_lhs = {"Ca", "CBsf"};
    const std::vector<model::specie_name> CBsf1_f_rhs = {"CBsCa"};
    statedef->addCompReac("__MESH__", CBsf1_f_lhs, CBsf1_f_rhs, CBsf1_f_kcst());
    statedef->addCompReac("__MESH__", CBsf1_f_rhs, CBsf1_f_lhs, CBsf1_b_kcst());

    const std::vector<model::specie_name> CBsCa_f_lhs = {"Ca", "CBsCa"};
    const std::vector<model::specie_name> CBsCa_f_rhs = {"CBCaCa"};
    statedef->addCompReac("__MESH__", CBsCa_f_lhs, CBsCa_f_rhs, CBsCa_f_kcst());
    statedef->addCompReac("__MESH__", CBsCa_f_rhs, CBsCa_f_lhs, CBsCa_b_kcst());

    const std::vector<model::specie_name> CBsf2_f_lhs = {"Ca", "CBsf"};
    const std::vector<model::specie_name> CBsf2_f_rhs = {"CBCaf"};
    statedef->addCompReac("__MESH__", CBsf2_f_lhs, CBsf2_f_rhs, CBsf2_f_kcst());
    statedef->addCompReac("__MESH__", CBsf2_f_rhs, CBsf2_f_lhs, CBsf2_b_kcst());

    const std::vector<model::specie_name> CBCaf_f_lhs = {"Ca", "CBCaf"};
    const std::vector<model::specie_name> CBCaf_f_rhs = {"CBCaCa"};
    statedef->addCompReac("__MESH__", CBCaf_f_lhs, CBCaf_f_rhs, CBCaf_f_kcst());
    statedef->addCompReac("__MESH__", CBCaf_f_rhs, CBCaf_f_lhs, CBCaf_b_kcst());

    const std::vector<model::specie_name> PVca_f_lhs = {"Ca", "PV"};
    const std::vector<model::specie_name> PVca_f_rhs = {"PVCa"};
    statedef->addCompReac("__MESH__", PVca_f_lhs, PVca_f_rhs, PVca_f_kcst());
    statedef->addCompReac("__MESH__", PVca_f_rhs, PVca_f_lhs, PVca_b_kcst());

    const std::vector<model::specie_name> PVmg_f_lhs = {"Mg", "PV"};
    const std::vector<model::specie_name> PVmg_f_rhs = {"PVMg"};
    statedef->addCompReac("__MESH__", PVmg_f_lhs, PVmg_f_rhs, PVmg_f_kcst());
    statedef->addCompReac("__MESH__", PVmg_f_rhs, PVmg_f_lhs, PVmg_b_kcst());

    const std::vector<model::specie_name> CaInflux_lhs = {};
    const std::vector<model::specie_name> CaInflux_rhs = {"Ca"};
    statedef->addCompReac("__MESH__", CaInflux_lhs, CaInflux_rhs, 0.0);

    statedef->addCompDiff("__MESH__", "Ca", DCST());
    statedef->addCompDiff("__MESH__", "CBsf", DCB());
    statedef->addCompDiff("__MESH__", "CBsCa", DCB());
    statedef->addCompDiff("__MESH__", "CBCaf", DCB());
    statedef->addCompDiff("__MESH__", "PV", DPV());
    statedef->addCompDiff("__MESH__", "PVCa", DPV());
    statedef->addCompDiff("__MESH__", "PVMg", DPV());

    statedef->addPatch("__MESH_BOUNDARY__", "__MESH__");
    const std::vector<model::specie_name> patch_species = {
        "Pump",
        "CaPump",
    };
    statedef->addPatchSpecs("__MESH_BOUNDARY__", patch_species);

    const std::vector<model::specie_name> PumpD_f_ilhs = {"Ca"};
    const std::vector<model::specie_name> PumpD_f_slhs = {"Pump"};
    const std::vector<model::specie_name> PumpD_f_srhs = {"CaPump"};
    statedef->addSurfReac(
        "__MESH_BOUNDARY__", PumpD_f_ilhs, PumpD_f_slhs, {}, {}, PumpD_f_srhs, {}, P_f_kcst());

    const std::vector<model::specie_name> PumpD_b_slhs = {"CaPump"};
    const std::vector<model::specie_name> PumpD_b_irhs = {"Ca"};
    const std::vector<model::specie_name> PumpD_b_srhs = {"Pump"};
    statedef->addSurfReac(
        "__MESH_BOUNDARY__", {}, PumpD_b_slhs, {}, PumpD_b_irhs, PumpD_b_srhs, {}, P_b_kcst());

    const std::vector<model::specie_name> PumpD_k_slhs = {"CaPump"};
    const std::vector<model::specie_name> PumpD_k_srhs = {"Pump"};
    statedef->addSurfReac(
        "__MESH_BOUNDARY__", {}, PumpD_k_slhs, {}, {}, PumpD_k_srhs, {}, P_k_kcst());
}

}  // namespace zee
