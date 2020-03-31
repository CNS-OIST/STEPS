#pragma once

#include <string>
#include <vector>

#include "opsplit/statedef.hpp"
#include "opsplit/vocabulary.hpp"

namespace zee {

struct CompartmentCount {
    CompartmentCount(model::compartment_id t_compartment,
                     model::specie_name t_specie,
                     PetscScalar t_num_mols)
        : compartment(std::move(t_compartment))
        , specie(std::move(t_specie))
        , num_mols(t_num_mols) {}
    const model::compartment_id compartment;
    const model::specie_name specie;
    const PetscScalar num_mols;
};

struct CompartmentConc {
    CompartmentConc(model::compartment_id t_compartment,
                    model::specie_name t_specie,
                    PetscScalar t_concentration)
        : compartment(std::move(t_compartment))
        , specie(std::move(t_specie))
        , concentration(t_concentration) {}
    const model::compartment_id compartment;
    const model::specie_name specie;
    const PetscScalar concentration;
};

struct PatchCount {
    PatchCount(model::patch_id t_patch, model::specie_name t_specie, PetscScalar t_num_mols)
        : patch(std::move(t_patch))
        , specie(std::move(t_specie))
        , num_mols(t_num_mols) {}
    const model::patch_id patch;
    const model::specie_name specie;
    const PetscScalar num_mols;
};

class Simdef {
  public:
    using compartment_counts_t = std::vector<CompartmentCount>;
    using compartment_concs_t = std::vector<CompartmentConc>;
    using patch_counts_t = std::vector<PatchCount>;

    inline std::unique_ptr<Statedef>& getStatedef() noexcept {
        return statedef;
    }

    inline const compartment_counts_t& getCompartementCounts() const noexcept {
        return compartmentCounts;
    }

    inline const patch_counts_t& getPatchCounts() const noexcept {
        return patchCounts;
    }


  protected:
    Simdef();
    std::unique_ptr<Statedef> statedef;
    compartment_counts_t compartmentCounts;
    patch_counts_t patchCounts;
};

struct SimpleSimdef: public Simdef {
    SimpleSimdef();
};

struct MultiCompartmentSimdef: public Simdef {
    MultiCompartmentSimdef();
};

struct SReacUnitTestSimdef: public Simdef {
    SReacUnitTestSimdef();
};

struct SurfaceReactionSimdef: public Simdef {
    SurfaceReactionSimdef();
};

struct ValidationSimdef: public Simdef {
    ValidationSimdef();

    // First order irreversible
    static constexpr PetscScalar KCST_foi() {
        return 5.0;
    }  // The reaction constant
    static constexpr PetscScalar N_foi() {
        return 50;
    }  // Can set count or conc
    static constexpr PetscInt NITER_foi() {
        return 100000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_foi() {
        return 2.0 / 100;
    }

    // First order reversible
    static constexpr PetscScalar KCST_f_for() {
        return 10.0;
    }  // The reaction constant
    static constexpr PetscScalar KCST_b_for() {
        return 2.0;
    }
    static constexpr PetscScalar COUNT_for() {
        return 100000;
    }  // Can set count or conc
    static constexpr PetscInt NITER_for() {
        return 10;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_for() {
        return 1.0 / 100;
    }

    // Second order irreversible A2
    static constexpr PetscScalar KCST_soA2() {
        return 10.0e6;
    }  // The reaction constant
    static constexpr PetscScalar CONCA_soA2() {
        return 10.0e-6;
    }
    static constexpr PetscInt NITER_soA2() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_soA2() {
        return 3.0 / 100;
    }

    // Second order irreversible AA
    static constexpr PetscScalar KCST_soAA() {
        return 5.0e6;
    }  // The reaction constant
    static constexpr PetscScalar CONCA_soAA() {
        return 20.0e-6;
    }
    static constexpr PetscScalar CONCB_soAA() {
        return CONCA_soAA();
    }
    static constexpr PetscInt NITER_soAA() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_soAA() {
        return 2.0 / 100;
    }

    // Second order irreversible AB
    static constexpr PetscScalar KCST_soAB() {
        return 5.0e6;
    }  // The reaction constant
    static constexpr PetscScalar CONCA_soAB() {
        return 1.0e-6;
    }
    static constexpr PetscScalar n_soAB() {
        return 2.0;
    }
    static constexpr PetscScalar CONCB_soAB() {
        return CONCA_soAB() / n_soAB();
    }
    static constexpr PetscInt NITER_soAB() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_soAB() {
        return 1.0 / 100;
    }

    // Third order irreversible A3
    static constexpr PetscScalar KCST_toA3() {
        return 1.0e12;
    }  // The reaction constant
    static constexpr PetscScalar CONCA_toA3() {
        return 10.0e-6;
    }
    static constexpr PetscInt NITER_toA3() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_toA3() {
        return 3.0 / 100;
    }

    // Third order irreversible A2B
    static constexpr PetscScalar KCST_toA2B() {
        return 0.1e12;
    }  // The reaction constant
    static constexpr PetscScalar CONCA_toA2B() {
        return 30.0e-6;
    }
    static constexpr PetscScalar CONCB_toA2B() {
        return 20.0e-6;
    }
    static constexpr PetscInt NITER_toA2B() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_toA2B() {
        return 1.0 / 100;
    }

    // Second order irreversible 2D
    static constexpr PetscScalar KCST_so2d() {
        return 10.0e10;
    }  // The reaction constant
    static constexpr PetscInt NITER_so2d() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_so2d() {
        return 2.0 / 100;
    }

    compartment_counts_t getfoi_CompartementCounts(bool zero = false) const;
    compartment_counts_t getfor_CompartementCounts(bool zero = false) const;
    compartment_concs_t getsoA2_CompartementConc(bool zero = false) const;
    compartment_concs_t getsoAA_CompartementConc(bool zero = false) const;
    compartment_concs_t getsoAB_CompartementConc(bool zero = false) const;
    compartment_concs_t gettoA3_CompartementConc(bool zero = false) const;
    compartment_concs_t gettoA2B_CompartementConc(bool zero = false) const;
};


struct SurfaceReactionsValidationSimdef: public Simdef {
    SurfaceReactionsValidationSimdef();

    // First order irreversible
    static constexpr PetscScalar KCST_foi() {
        return 5.0;
    }  // The reaction constant
    static constexpr PetscScalar N_foi() {
        return 50;
    }  // Can set count or conc
    static constexpr PetscInt NITER_foi() {
        return 100000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_foi() {
        return 2.0 / 100;
    }

    // First order reversible
    static constexpr PetscScalar KCST_f_for() {
        return 10.0;
    }  // The reaction constant
    static constexpr PetscScalar KCST_b_for() {
        return 2.0;
    }
    static constexpr PetscScalar COUNT_for() {
        return 100000.0;
    }  // Can set count or conc
    static constexpr PetscInt NITER_for() {
        return 10;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_for() {
        return 1.0 / 100;
    }

    // Second order irreversible A2
    static constexpr PetscScalar KCST_soA2() {
        return 10.0e6;
    }  // The reaction constant
    static constexpr PetscScalar CONCA_soA2() {
        return 10.0e-6;
    }
    static constexpr PetscInt NITER_soA2() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_soA2() {
        return 3.0 / 100;
    }

    // Second order irreversible AA
    static constexpr PetscScalar KCST_soAA() {
        return 5.0e6;
    }  // The reaction constant
    static constexpr PetscScalar CONCA_soAA() {
        return 20.0e-6;
    }
    static constexpr PetscScalar CONCB_soAA() {
        return CONCA_soAA();
    }
    static constexpr PetscInt NITER_soAA() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_soAA() {
        return 2.0 / 100;
    }

    // Second order irreversible AB
    static constexpr PetscScalar KCST_soAB() {
        return 5.0e6;
    }  // The reaction constant
    static constexpr PetscScalar CONCA_soAB() {
        return 1.0e-6;
    }
    static constexpr PetscScalar n_soAB() {
        return 2.0;
    }
    static constexpr PetscScalar CONCB_soAB() {
        return CONCA_soAB() / n_soAB();
    }
    static constexpr PetscInt NITER_soAB() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_soAB() {
        return 1.0 / 100;
    }

    // Third order irreversible A3
    static constexpr PetscScalar KCST_toA3() {
        return 1.0e12;
    }  // The reaction constant
    static constexpr PetscScalar CONCA_toA3() {
        return 10.0e-6;
    }
    static constexpr PetscInt NITER_toA3() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_toA3() {
        return 3.0 / 100;
    }

    // Third order irreversible A2B
    static constexpr PetscScalar KCST_toA2B() {
        return 0.1e12;
    }  // The reaction constant
    static constexpr PetscScalar CONCA_toA2B() {
        return 30.0e-6;
    }
    static constexpr PetscScalar CONCB_toA2B() {
        return 20.0e-6;
    }
    static constexpr PetscInt NITER_toA2B() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_toA2B() {
        return 1.0 / 100;
    }

    // Second order irreversible 2D
    static constexpr PetscScalar KCST_so2d() {
        return 10.0e10;
    }  // The reaction constant
    static constexpr PetscInt NITER_so2d() {
        return 1000;
    }  // The number of iterations
    static constexpr PetscScalar tolerance_so2d() {
        return 2.0 / 100;
    }
    static constexpr PetscScalar COUNTA_so2d() {
        return 100;
    }  // The number of Molecule A
    static constexpr PetscScalar COUNTB_so2d() {
        return 100 / 2;
    }  // The number of Molecule B

    compartment_counts_t getfoi_CompartementCounts() const;
    compartment_counts_t getfor_CompartementCounts() const;
    compartment_concs_t getsoA2_CompartementConc() const;
    compartment_concs_t getsoAA_CompartementConc() const;
    compartment_concs_t getsoAB_CompartementConc() const;
    compartment_concs_t gettoA3_CompartementConc() const;
    compartment_concs_t gettoA2B_CompartementConc() const;
    patch_counts_t getso2d_PatchCount() const;
};

struct CaBurstSimdef: public Simdef {
    CaBurstSimdef();

    static constexpr PetscScalar init_pot() {
        return -60e-3;
    }
    static constexpr PetscScalar TEMPERATURE() {
        return 34.0;
    }
    static constexpr PetscInt Q10() {
        return 3;
    }
    static constexpr PetscScalar FARADAY() {
        return 96485.3365;
    }  // Faraday constant: unit of FARADAY is C/mol
    static constexpr PetscScalar R() {
        return 8.3144621;
    }  // Molar Gas Constant: unit of R is J/mol K
    static constexpr PetscScalar AVOGADRO() {
        return 6.02214129e23;
    }  // Avogadro constant: unit of AVOGADRO is /mol
    static constexpr PetscScalar E_CHARGE() {
        return 1.602176565e-19;
    }  // Elementary charge: unit of E_CHARGE is C
    static PetscScalar Qt() {
        return std::pow(Q10(), ((TEMPERATURE() - 23) / 10));
    }
    static PetscScalar Qt_mslo() {
        return std::pow(Q10(), ((TEMPERATURE() - 25) / 10));
    }
    static constexpr PetscScalar Ra() {
        return 235.7 * 1.0e-2;
    }  // BULK RESISTIVITY
    static constexpr PetscScalar memb_capac() {
        return 1.5e-2;
    }  // MEMBRANE CAPACITANCE

    // CaP channels density & permiability per channel

    static constexpr PetscScalar CaP_P() {
        return 2.5e-20;
    }
    static constexpr PetscScalar CaP_ro() {
        return 3.8e13;
    }

    // CaP channel parameters
    static constexpr PetscScalar vhalfm() {
        return -29.458;
    }
    static constexpr PetscScalar cvm() {
        return 8.429;
    }
    PetscScalar minf_cap(PetscScalar V) {
        PetscScalar vhalfm = -29.458;
        PetscScalar cvm = 8.429;
        PetscScalar vshift = 0.0;
        return (1.0 / (1.0 + std::exp(-(V - vhalfm - vshift) / cvm)));
    }
    PetscScalar tau_cap(PetscScalar V) {
        PetscScalar vshift = 0.0;
        if ((V - vshift) >= -40) {
            return (0.2702 +
                    1.1622 * std::exp(-(V + 26.798 - vshift) * (V + 26.798 - vshift) / 164.19));
        } else {
            return (0.6923 * std::exp((V - vshift) / 1089.372));
        }
    }
    PetscScalar alpha_cap(PetscScalar V) {
        return (minf_cap(V) / tau_cap(V));
    }
    PetscScalar beta_cap(PetscScalar V) {
        return ((1.0 - minf_cap(V)) / tau_cap(V));
    }

    // Intitial conditions
    static constexpr PetscScalar CaP_m0_p() {
        return 0.92402;
    }
    static constexpr PetscScalar CaP_m1_p() {
        return 0.073988;
    }
    static constexpr PetscScalar CaP_m2_p() {
        return 0.0019748;
    }
    static constexpr PetscScalar CaP_m3_p() {
        return 1.7569e-05;
    }

    // CaT channels density & permiability per channel
    static constexpr PetscScalar CaT_P() {
        return 1.65e-20;
    }
    static constexpr PetscScalar CaT_ro() {
        return 3.7576e12;
    }
    PetscScalar minf_cat(PetscScalar V) {
        PetscScalar vhalfm = -52.0;
        PetscScalar cvm = -5.0;
        PetscScalar vshift = 0.0;
        return (1.0 / (1.0 + std::exp((V - vhalfm - vshift) / cvm)));
    }
    PetscScalar taum_cat(PetscScalar V) {
        PetscScalar vshift = 0.0;
        if (V > -90.0) {
            return (1.0 + 1.0 / (std::exp((V + 40.0 - vshift) / 9.0) +
                                 std::exp(-(V + 102.0 - vshift) / 18.0)));
        } else {
            return 1.0;
        }
    }
    PetscScalar hinf_cat(PetscScalar V) {
        PetscScalar vhalfh = -72.0;
        PetscScalar cvh = 7.0;
        PetscScalar vshift = 0.0;
        return (1.0 / (1.0 + std::exp((V - vhalfh - vshift) / cvh)));
    }
    PetscScalar tauh_cat(PetscScalar V) {
        PetscScalar vshift = 0.0;
        return (15.0 + 1.0 / (std::exp((V + 32.0 - vshift) / 7.0)));
    }
    PetscScalar alpham_cat(PetscScalar V) {
        return (minf_cat(V) / taum_cat(V));
    }
    PetscScalar betam_cat(PetscScalar V) {
        return ((1 - minf_cat(V)) / taum_cat(V));
    }
    PetscScalar alphah_cat(PetscScalar V) {
        return (hinf_cat(V) / tauh_cat(V));
    }
    PetscScalar betah_cat(PetscScalar V) {
        return ((1 - hinf_cat(V)) / tauh_cat(V));
    }

    static constexpr PetscScalar CaT_m0h0_p() {
        return 0.58661;
    }
    static constexpr PetscScalar CaT_m1h0_p() {
        return 0.23687;
    }
    static constexpr PetscScalar CaT_m2h0_p() {
        return 0.023912;
    }
    static constexpr PetscScalar CaT_m0h1_p() {
        return 0.10564;
    }
    static constexpr PetscScalar CaT_m1h1_p() {
        return 0.042658;
    }
    static constexpr PetscScalar CaT_m2h1_p() {
        return 0.0043063;
    }

    // BK channels density & conductance per channel
    static constexpr PetscScalar BK_G() {
        return 2.1e-10;
    }
    static constexpr PetscScalar BK_ro() {
        return 2.0238e12;
    }
    static constexpr PetscScalar BK_rev() {
        return -77e-3;
    }
    static constexpr PetscScalar Qo() {
        return 0.73;
    }
    static constexpr PetscScalar Qc() {
        return -0.67;
    }

    static constexpr PetscScalar pf0() {
        return 2.39;
    }
    static constexpr PetscScalar pf1() {
        return 5.4918;
    }
    static constexpr PetscScalar pf2() {
        return 24.6205;
    }
    static constexpr PetscScalar pf3() {
        return 142.4546;
    }
    static constexpr PetscScalar pf4() {
        return 211.0220;
    }

    static constexpr PetscScalar pb0() {
        return 3936.0;
    }
    static constexpr PetscScalar pb1() {
        return 687.3251;
    }
    static constexpr PetscScalar pb2() {
        return 234.5875;
    }
    static constexpr PetscScalar pb3() {
        return 103.2204;
    }
    static constexpr PetscScalar pb4() {
        return 11.6581;
    }

    static constexpr PetscScalar k1() {
        return 1.0e6;
    }
    static constexpr PetscScalar onoffrate() {
        return 1.0e3;
    }
    static constexpr PetscScalar L0() {
        return 1806.0;
    }

    static constexpr PetscScalar Kc() {
        return 8.63e-6;
    }
    static constexpr PetscScalar Ko() {
        return 0.6563e-6;
    }

    PetscScalar c_01() {
        return 4.0 * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar c_12() {
        return 3.0 * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar c_23() {
        return 2.0 * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar c_34() {
        return 1.0 * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar o_01() {
        return 4.0 * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar o_12() {
        return 3.0 * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar o_23() {
        return 2.0 * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar o_34() {
        return 1.0 * k1() * onoffrate() * Qt_mslo();
    }

    PetscScalar c_10() {
        return 1.0 * Kc() * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar c_21() {
        return 2.0 * Kc() * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar c_32() {
        return 3.0 * Kc() * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar c_43() {
        return 4.0 * Kc() * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar o_10() {
        return 1.0 * Ko() * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar o_21() {
        return 2.0 * Ko() * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar o_32() {
        return 3.0 * Ko() * k1() * onoffrate() * Qt_mslo();
    }
    PetscScalar o_43() {
        return 4.0 * Ko() * k1() * onoffrate() * Qt_mslo();
    }

    PetscScalar f_0(PetscScalar mV) {
        return pf0() * Qt_mslo() *
               (std::exp((Qo() * FARADAY() * mV) / (R() * (TEMPERATURE() + 273.15))));
    }
    PetscScalar f_1(PetscScalar mV) {
        return pf1() * Qt_mslo() *
               (std::exp((Qo() * FARADAY() * mV) / (R() * (TEMPERATURE() + 273.15))));
    }
    PetscScalar f_2(PetscScalar mV) {
        return pf2() * Qt_mslo() *
               (std::exp((Qo() * FARADAY() * mV) / (R() * (TEMPERATURE() + 273.15))));
    }
    PetscScalar f_3(PetscScalar mV) {
        return pf3() * Qt_mslo() *
               (std::exp((Qo() * FARADAY() * mV) / (R() * (TEMPERATURE() + 273.15))));
    }
    PetscScalar f_4(PetscScalar mV) {
        return pf4() * Qt_mslo() *
               (std::exp((Qo() * FARADAY() * mV) / (R() * (TEMPERATURE() + 273.15))));
    }

    PetscScalar b_0(PetscScalar mV) {
        return pb0() * Qt_mslo() *
               (std::exp((Qc() * FARADAY() * mV) / (R() * (TEMPERATURE() + 273.15))));
    }
    PetscScalar b_1(PetscScalar mV) {
        return pb1() * Qt_mslo() *
               (std::exp((Qc() * FARADAY() * mV) / (R() * (TEMPERATURE() + 273.15))));
    }
    PetscScalar b_2(PetscScalar mV) {
        return pb2() * Qt_mslo() *
               (std::exp((Qc() * FARADAY() * mV) / (R() * (TEMPERATURE() + 273.15))));
    }
    PetscScalar b_3(PetscScalar mV) {
        return pb3() * Qt_mslo() *
               (std::exp((Qc() * FARADAY() * mV) / (R() * (TEMPERATURE() + 273.15))));
    }
    PetscScalar b_4(PetscScalar mV) {
        return pb4() * Qt_mslo() *
               (std::exp((Qc() * FARADAY() * mV) / (R() * (TEMPERATURE() + 273.15))));
    }

    static constexpr PetscScalar BK_C0_p() {
        return 0.99997;
    }
    static constexpr PetscScalar BK_C1_p() {
        return 4.3619e-07;
    }
    static constexpr PetscScalar BK_C2_p() {
        return 4.1713e-09;
    }
    static constexpr PetscScalar BK_C3_p() {
        return 4.4449e-11;
    }
    static constexpr PetscScalar BK_C4_p() {
        return 6.3132e-14;
    }

    static constexpr PetscScalar BK_O0_p() {
        return 2.5202e-05;
    }
    static constexpr PetscScalar BK_O1_p() {
        return 1.1765e-06;
    }
    static constexpr PetscScalar BK_O2_p() {
        return 6.6148e-08;
    }
    static constexpr PetscScalar BK_O3_p() {
        return 2.4392e-09;
    }
    static constexpr PetscScalar BK_O4_p() {
        return 4.0981e-11;
    }

    static constexpr PetscScalar SK_G() {
        return 1.0e-11;
    }
    static constexpr PetscScalar SK_ro() {
        return 31.0e10;
    }
    static constexpr PetscScalar SK_rev() {
        return -77e-3;
    }

    static constexpr PetscScalar invc1() {
        return 80.0;
    }
    static constexpr PetscScalar invc2() {
        return 80.0;
    }
    static constexpr PetscScalar invc3() {
        return 200.0;
    }

    static constexpr PetscScalar invo1() {
        return 1000.0;
    }
    static constexpr PetscScalar invo2() {
        return 100.0;
    }

    static constexpr PetscScalar diro1() {
        return 160.0;
    }
    static constexpr PetscScalar diro2() {
        return 1200.0;
    }
    static constexpr PetscScalar dirc2() {
        return 200.0e6;
    }
    static constexpr PetscScalar dirc3() {
        return 160.0e6;
    }
    static constexpr PetscScalar dirc4() {
        return 80.0e6;
    }

    PetscScalar invc1_t() {
        return invc1() * Qt();
    }
    PetscScalar invc2_t() {
        return invc2() * Qt();
    }
    PetscScalar invc3_t() {
        return invc3() * Qt();
    }

    PetscScalar invo1_t() {
        return invo1() * Qt();
    }
    PetscScalar invo2_t() {
        return invo2() * Qt();
    }

    PetscScalar diro1_t() {
        return diro1() * Qt();
    }
    PetscScalar diro2_t() {
        return diro2() * Qt();
    }

    PetscScalar dirc2_t() {
        return dirc2() * Qt() / 3.0;
    }
    PetscScalar dirc3_t() {
        return dirc3() * Qt() / 3.0;
    }
    PetscScalar dirc4_t() {
        return dirc4() * Qt() / 3.0;
    }

    static constexpr PetscScalar SK_C1_p() {
        return 0.96256;
    }
    static constexpr PetscScalar SK_C2_p() {
        return 0.036096;
    }
    static constexpr PetscScalar SK_C3_p() {
        return 0.0010829;
    }
    static constexpr PetscScalar SK_C4_p() {
        return 6.4973e-06;
    }

    static constexpr PetscScalar SK_O1_p() {
        return 0.00017326;
    }
    static constexpr PetscScalar SK_O2_p() {
        return 7.7967e-05;
    }

    // leak current channel density & conductance per channel
    static constexpr PetscScalar L_G() {
        return 4.0e-14;
    }
    static constexpr PetscScalar L_ro() {
        return 25.0e10;
    }
    static constexpr PetscScalar L_rev() {
        return -61e-3;
    }

    // Pump parameters
    static constexpr PetscScalar P_f_kcst() {
        return 3e9;
    }
    static constexpr PetscScalar P_b_kcst() {
        return 1.75e4;
    }
    static constexpr PetscScalar P_k_kcst() {
        return 7.255e4;
    }

    // CALCIUM BUFFERING MODEL
    // Ca concentrations
    static constexpr PetscScalar Ca_oconc() {
        return 2e-3;
    }
    static constexpr PetscScalar Ca_iconc() {
        return 45e-9;
    }

    // Mg concentrations
    static constexpr PetscScalar Mg_conc() {
        return 590e-6;
    }

    // Buffer concentrations
    static constexpr PetscScalar iCBsf_conc() {
        return 27.704e-6;
    }
    static constexpr PetscScalar iCBCaf_conc() {
        return 2.6372e-6;
    }
    static constexpr PetscScalar iCBsCa_conc() {
        return 1.5148e-6;
    }
    static constexpr PetscScalar iCBCaCa_conc() {
        return 0.14420e-6;
    }

    static constexpr PetscScalar CBsf_conc() {
        return 110.82e-6;
    }
    static constexpr PetscScalar CBCaf_conc() {
        return 10.549e-6;
    }
    static constexpr PetscScalar CBsCa_conc() {
        return 6.0595e-6;
    }
    static constexpr PetscScalar CBCaCa_conc() {
        return 0.57682e-6;
    }

    static constexpr PetscScalar PV_conc() {
        return 3.2066e-6;
    }
    static constexpr PetscScalar PVCa_conc() {
        return 16.252e-6;
    }
    static constexpr PetscScalar PVMg_conc() {
        return 60.541e-6;
    }

    static constexpr PetscScalar DCST() {
        return 0.223e-9;
    }
    static constexpr PetscScalar DCB() {
        return 0.028e-9;
    }
    static constexpr PetscScalar DPV() {
        return 0.043e-9;
    }

    static constexpr PetscScalar iCBsf1_f_kcst() {
        return 4.35e7;
    }
    static constexpr PetscScalar iCBsf1_b_kcst() {
        return 35.8;
    }

    static constexpr PetscScalar iCBsCa_f_kcst() {
        return 0.55e7;
    }
    static constexpr PetscScalar iCBsCa_b_kcst() {
        return 2.6;
    }

    static constexpr PetscScalar iCBsf2_f_kcst() {
        return 0.55e7;
    }
    static constexpr PetscScalar iCBsf2_b_kcst() {
        return 2.6;
    }

    static constexpr PetscScalar iCBCaf_f_kcst() {
        return 4.35e7;
    }
    static constexpr PetscScalar iCBCaf_b_kcst() {
        return 35.8;
    }

    static constexpr PetscScalar CBsf1_f_kcst() {
        return 4.35e7;
    }
    static constexpr PetscScalar CBsf1_b_kcst() {
        return 35.8;
    }

    static constexpr PetscScalar CBsCa_f_kcst() {
        return 0.55e7;
    }
    static constexpr PetscScalar CBsCa_b_kcst() {
        return 2.6;
    }

    static constexpr PetscScalar CBsf2_f_kcst() {
        return 0.55e7;
    }
    static constexpr PetscScalar CBsf2_b_kcst() {
        return 2.6;
    }

    static constexpr PetscScalar CBCaf_f_kcst() {
        return 4.35e7;
    }
    static constexpr PetscScalar CBCaf_b_kcst() {
        return 35.8;
    }

    static constexpr PetscScalar PVca_f_kcst() {
        return 10.7e7;
    }
    static constexpr PetscScalar PVca_b_kcst() {
        return 0.95;
    }

    static constexpr PetscScalar PVmg_f_kcst() {
        return 0.8e6;
    }
    static constexpr PetscScalar PVmg_b_kcst() {
        return 25.0;
    }
};

}  // namespace zee
