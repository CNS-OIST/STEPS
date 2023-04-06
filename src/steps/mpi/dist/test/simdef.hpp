#pragma once

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include "geom/dist/fwd.hpp"
#include "model/fwd.hpp"
#include "mpi/dist/tetopsplit/definition/fwd.hpp"
#include "util/vocabulary.hpp"

namespace steps {
namespace dist {

struct CompartmentCount {
  CompartmentCount(model::species_name t_species, osh::Real t_num_mols)
      : species(std::move(t_species)), num_mols(t_num_mols) {}
  model::species_name species;
  osh::Real num_mols;
};

struct CompartmentConc {
  CompartmentConc(model::species_name t_species, osh::Real t_concentration)
      : species(std::move(t_species)), concentration(t_concentration) {}
  model::species_name species;
  osh::Real concentration;
};

struct PatchCount {
  PatchCount(model::patch_id t_patch, model::species_name t_species,
             osh::Real t_num_mols)
      : patch(std::move(t_patch)), species(std::move(t_species)),
        num_mols(t_num_mols) {}
  const model::patch_id patch;
  const model::species_name species;
  const osh::Real num_mols;
};

class Simdef {
public:
  Simdef(const steps::model::Model &model, const steps::dist::DistMesh &mesh);
  using compartment_counts_t =
      std::unordered_map<model::compartment_id, std::vector<CompartmentCount>>;
  using compartment_concs_t =
      std::unordered_map<model::compartment_id, std::vector<CompartmentConc>>;
  using patch_counts_t = std::vector<PatchCount>;

  inline std::unique_ptr<Statedef> &getStatedef() noexcept { return statedef; }

  inline compartment_counts_t
  getCompartementCounts(osh::Real factor = 1.0) const noexcept {
    auto comps_count = compartmentCounts;
    for (auto &comp_counts : comps_count) {
      for (auto &count : comp_counts.second) {
        count.num_mols *= factor;
      }
    }
    return comps_count;
  }

  inline const patch_counts_t &getPatchCounts() const noexcept {
    return patchCounts;
  }

protected:
  Simdef();
  std::unique_ptr<Statedef> statedef;
  compartment_counts_t compartmentCounts;
  patch_counts_t patchCounts;
};

struct SimpleSimdef : public Simdef {
  SimpleSimdef(const steps::model::Model &model,
               const steps::dist::DistMesh &mesh);
};

struct SingleCompDistSimdef : public Simdef {
  SingleCompDistSimdef(const steps::model::Model &model,
                       const steps::dist::DistMesh &mesh);
};

struct MultiCompartmentDiffSimdef : public Simdef {
  MultiCompartmentDiffSimdef(const steps::model::Model &model,
                             const steps::dist::DistMesh &mesh);
};

struct MultiCompartmentSimdef : public Simdef {
  MultiCompartmentSimdef(const steps::model::Model &model,
                         const steps::dist::DistMesh &mesh);
};

struct SReacUnitTestSimdef : public Simdef {
  SReacUnitTestSimdef(const steps::model::Model &model,
                      const steps::dist::DistMesh &mesh);
};

struct GHKCurrentUnitTestSimdef : public Simdef {
  GHKCurrentUnitTestSimdef(const steps::model::Model &model,
                           const steps::dist::DistMesh &mesh);
};

struct SurfaceReactionSimdef : public Simdef {
  SurfaceReactionSimdef(const steps::model::Model &model,
                        const steps::dist::DistMesh &mesh);
};

struct ValidationSimdef : public Simdef {
  ValidationSimdef(const steps::model::Model &model,
                   const steps::dist::DistMesh &mesh);

  // First order irreversible
  static constexpr osh::Real KCST_foi() { return 5.0; } // The reaction constant
  static constexpr osh::Real N_foi() { return 50; }     // Can set count or conc
  static constexpr osh::I64 NITER_foi() {
    return 100000;
  } // The number of iterations
  static constexpr osh::Real tolerance_foi() { return 2.0 / 100; }

  // First order reversible
  static constexpr osh::Real KCST_f_for() {
    return 10.0;
  } // The reaction constant
  static constexpr osh::Real KCST_b_for() { return 2.0; }
  static constexpr osh::Real COUNT_for() {
    return 100000;
  } // Can set count or conc
  static constexpr osh::I64 NITER_for() {
    return 10;
  } // The number of iterations
  static constexpr osh::Real tolerance_for() { return 1.0 / 100; }

  // Second order irreversible A2
  static constexpr osh::Real KCST_soA2() {
    return 10.0e6;
  } // The reaction constant
  static constexpr osh::Real CONCA_soA2() { return 10.0e-6; }
  static constexpr osh::I64 NITER_soA2() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_soA2() { return 3.0 / 100; }

  // Second order irreversible AA
  static constexpr osh::Real KCST_soAA() {
    return 5.0e6;
  } // The reaction constant
  static constexpr osh::Real CONCA_soAA() { return 20.0e-6; }
  static constexpr osh::Real CONCB_soAA() { return CONCA_soAA(); }
  static constexpr osh::I64 NITER_soAA() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_soAA() { return 2.0 / 100; }

  // Second order irreversible AB
  static constexpr osh::Real KCST_soAB() {
    return 5.0e6;
  } // The reaction constant
  static constexpr osh::Real CONCA_soAB() { return 1.0e-6; }
  static constexpr osh::Real n_soAB() { return 2.0; }
  static constexpr osh::Real CONCB_soAB() { return CONCA_soAB() / n_soAB(); }
  static constexpr osh::I64 NITER_soAB() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_soAB() { return 1.0 / 100; }

  // Third order irreversible A3
  static constexpr osh::Real KCST_toA3() {
    return 1.0e12;
  } // The reaction constant
  static constexpr osh::Real CONCA_toA3() { return 10.0e-6; }
  static constexpr osh::I64 NITER_toA3() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_toA3() { return 3.0 / 100; }

  // Third order irreversible A2B
  static constexpr osh::Real KCST_toA2B() {
    return 0.1e12;
  } // The reaction constant
  static constexpr osh::Real CONCA_toA2B() { return 30.0e-6; }
  static constexpr osh::Real CONCB_toA2B() { return 20.0e-6; }
  static constexpr osh::I64 NITER_toA2B() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_toA2B() { return 1.0 / 100; }

  // Second order irreversible 2D
  static constexpr osh::Real KCST_so2d() {
    return 10.0e10;
  } // The reaction constant
  static constexpr osh::I64 NITER_so2d() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_so2d() { return 2.0 / 100; }

  compartment_counts_t getfoi_CompartementCounts(bool zero = false) const;
  compartment_counts_t getfor_CompartementCounts(bool zero = false) const;
  compartment_concs_t getsoA2_CompartementConc(bool zero = false) const;
  compartment_concs_t getsoAA_CompartementConc(bool zero = false) const;
  compartment_concs_t getsoAB_CompartementConc(bool zero = false) const;
  compartment_concs_t gettoA3_CompartementConc(bool zero = false) const;
  compartment_concs_t gettoA2B_CompartementConc(bool zero = false) const;
};

struct SurfaceReactionsValidationSimdef : public Simdef {
  SurfaceReactionsValidationSimdef(const steps::model::Model &model,
                                   const steps::dist::DistMesh &mesh);

  // First order irreversible
  static constexpr osh::Real KCST_foi() { return 5.0; } // The reaction constant
  static constexpr osh::Real N_foi() { return 50; }     // Can set count or conc
  static constexpr osh::I64 NITER_foi() {
    return 100000;
  } // The number of iterations
  static constexpr osh::Real tolerance_foi() { return 2.0 / 100; }

  // First order reversible
  static constexpr osh::Real KCST_f_for() {
    return 10.0;
  } // The reaction constant
  static constexpr osh::Real KCST_b_for() { return 2.0; }
  static constexpr osh::Real COUNT_for() {
    return 100000.0;
  } // Can set count or conc
  static constexpr osh::I64 NITER_for() {
    return 10;
  } // The number of iterations
  static constexpr osh::Real tolerance_for() { return 1.0 / 100; }

  // Second order irreversible A2
  static constexpr osh::Real KCST_soA2() {
    return 10.0e6;
  } // The reaction constant
  static constexpr osh::Real CONCA_soA2() { return 10.0e-6; }
  static constexpr osh::I64 NITER_soA2() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_soA2() { return 3.0 / 100; }

  // Second order irreversible AA
  static constexpr osh::Real KCST_soAA() {
    return 5.0e6;
  } // The reaction constant
  static constexpr osh::Real CONCA_soAA() { return 20.0e-6; }
  static constexpr osh::Real CONCB_soAA() { return CONCA_soAA(); }
  static constexpr osh::I64 NITER_soAA() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_soAA() { return 2.0 / 100; }

  // Second order irreversible AB
  static constexpr osh::Real KCST_soAB() {
    return 5.0e6;
  } // The reaction constant
  static constexpr osh::Real CONCA_soAB() { return 1.0e-6; }
  static constexpr osh::Real n_soAB() { return 2.0; }
  static constexpr osh::Real CONCB_soAB() { return CONCA_soAB() / n_soAB(); }
  static constexpr osh::I64 NITER_soAB() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_soAB() { return 1.0 / 100; }

  // Third order irreversible A3
  static constexpr osh::Real KCST_toA3() {
    return 1.0e12;
  } // The reaction constant
  static constexpr osh::Real CONCA_toA3() { return 10.0e-6; }
  static constexpr osh::I64 NITER_toA3() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_toA3() { return 3.0 / 100; }

  // Third order irreversible A2B
  static constexpr osh::Real KCST_toA2B() {
    return 0.1e12;
  } // The reaction constant
  static constexpr osh::Real CONCA_toA2B() { return 30.0e-6; }
  static constexpr osh::Real CONCB_toA2B() { return 20.0e-6; }
  static constexpr osh::I64 NITER_toA2B() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_toA2B() { return 1.0 / 100; }

  // Second order irreversible 2D
  static constexpr osh::Real KCST_so2d() {
    return 10.0e10;
  } // The reaction constant
  static constexpr osh::I64 NITER_so2d() {
    return 1000;
  } // The number of iterations
  static constexpr osh::Real tolerance_so2d() { return 2.0 / 100; }
  static constexpr osh::Real COUNTA_so2d() {
    return 100;
  } // The number of Molecule A
  static constexpr osh::Real COUNTB_so2d() {
    return 100 / 2;
  } // The number of Molecule B

  compartment_counts_t getfoi_CompartementCounts() const;
  compartment_counts_t getfor_CompartementCounts() const;
  compartment_concs_t getsoA2_CompartementConc() const;
  compartment_concs_t getsoAA_CompartementConc() const;
  compartment_concs_t getsoAB_CompartementConc() const;
  compartment_concs_t gettoA3_CompartementConc() const;
  compartment_concs_t gettoA2B_CompartementConc() const;
  patch_counts_t getso2d_PatchCount() const;
};

struct CaBurstSimdef : public Simdef {
  CaBurstSimdef(const steps::model::Model &model,
                const steps::dist::DistMesh &mesh);

  static constexpr osh::Real init_pot() { return -60e-3; }
  static constexpr osh::Real TEMPERATURE() { return 34.0; }
  static constexpr osh::I64 Q10() { return 3; }
  static constexpr osh::Real FARADAY() {
    return 96485.3365;
  } // Faraday constant: unit of FARADAY is C/mol
  static constexpr osh::Real R() {
    return 8.3144621;
  } // Molar Gas Constant: unit of R is J/mol K
  static constexpr osh::Real AVOGADRO() {
    return 6.02214129e23;
  } // Avogadro constant: unit of AVOGADRO is /mol
  static constexpr osh::Real E_CHARGE() {
    return 1.602176565e-19;
  } // Elementary charge: unit of E_CHARGE is C
  static osh::Real Qt() { return std::pow(Q10(), ((TEMPERATURE() - 23) / 10)); }
  static osh::Real Qt_mslo() {
    return std::pow(Q10(), ((TEMPERATURE() - 25) / 10));
  }
  static constexpr osh::Real Ra() { return 235.7 * 1.0e-2; } // BULK RESISTIVITY

  static constexpr osh::Real memb_capac() {
    return 1.5e-2;
  } // MEMBRANE CAPACITANCE

  // MEMBRANE CAPACITANCE for complete model
  static constexpr osh::Real Cm() { return 0.64e-2; }

  static constexpr osh::Real memb_capac_proximal() { return 1.2 * Cm(); }

  static constexpr osh::Real memb_capac_spiny() { return 5.3 * Cm(); }

  // CaP channels density & permiability per channel

  static constexpr osh::Real CaP_P() { return 2.5e-20; }
  static constexpr osh::Real CaP_ro() { return 3.8e13; }

  // CaP channel parameters
  static constexpr osh::Real vhalfm() { return -29.458; }
  static constexpr osh::Real cvm() { return 8.429; }
  osh::Real minf_cap(osh::Real V) {
    osh::Real vhalfm = -29.458;
    osh::Real cvm = 8.429;
    osh::Real vshift = 0.0;
    return (1.0 / (1.0 + std::exp(-(V - vhalfm - vshift) / cvm)));
  }
  osh::Real tau_cap(osh::Real V) {
    osh::Real vshift = 0.0;
    if ((V - vshift) >= -40) {
      return (0.2702 + 1.1622 * std::exp(-(V + 26.798 - vshift) *
                                         (V + 26.798 - vshift) / 164.19));
    } else {
      return (0.6923 * std::exp((V - vshift) / 1089.372));
    }
  }
  osh::Real alpha_cap(osh::Real V) { return (minf_cap(V) / tau_cap(V)); }
  osh::Real beta_cap(osh::Real V) { return ((1.0 - minf_cap(V)) / tau_cap(V)); }

  // Intitial conditions
  static constexpr osh::Real CaP_m0_p() { return 0.92402; }
  static constexpr osh::Real CaP_m1_p() { return 0.073988; }
  static constexpr osh::Real CaP_m2_p() { return 0.0019748; }
  static constexpr osh::Real CaP_m3_p() { return 1.7569e-05; }

  // CaT channels density & permiability per channel
  static constexpr osh::Real CaT_P() { return 1.65e-20; }
  static constexpr osh::Real CaT_ro(bool with_ampa) {
    if (with_ampa)
      return 1.9636e12;
    else
      return 3.7576e12;
  }
  osh::Real minf_cat(osh::Real V) {
    osh::Real vhalfm = -52.0;
    osh::Real cvm = -5.0;
    osh::Real vshift = 0.0;
    return (1.0 / (1.0 + std::exp((V - vhalfm - vshift) / cvm)));
  }
  osh::Real taum_cat(osh::Real V) {
    osh::Real vshift = 0.0;
    if (V > -90.0) {
      return (1.0 + 1.0 / (std::exp((V + 40.0 - vshift) / 9.0) +
                           std::exp(-(V + 102.0 - vshift) / 18.0)));
    } else {
      return 1.0;
    }
  }
  osh::Real hinf_cat(osh::Real V) {
    osh::Real vhalfh = -72.0;
    osh::Real cvh = 7.0;
    osh::Real vshift = 0.0;
    return (1.0 / (1.0 + std::exp((V - vhalfh - vshift) / cvh)));
  }
  osh::Real tauh_cat(osh::Real V) {
    osh::Real vshift = 0.0;
    return (15.0 + 1.0 / (std::exp((V + 32.0 - vshift) / 7.0)));
  }
  osh::Real alpham_cat(osh::Real V) { return (minf_cat(V) / taum_cat(V)); }
  osh::Real betam_cat(osh::Real V) { return ((1 - minf_cat(V)) / taum_cat(V)); }
  osh::Real alphah_cat(osh::Real V) { return (hinf_cat(V) / tauh_cat(V)); }
  osh::Real betah_cat(osh::Real V) { return ((1 - hinf_cat(V)) / tauh_cat(V)); }

  static constexpr osh::Real CaT_m0h0_p() { return 0.58661; }
  static constexpr osh::Real CaT_m1h0_p() { return 0.23687; }
  static constexpr osh::Real CaT_m2h0_p() { return 0.023912; }
  static constexpr osh::Real CaT_m0h1_p() { return 0.10564; }
  static constexpr osh::Real CaT_m1h1_p() { return 0.042658; }
  static constexpr osh::Real CaT_m2h1_p() { return 0.0043063; }

  // BK channels density & conductance per channel
  static constexpr osh::Real BK_G() { return 2.1e-10; }

  // Basically, 2.0238e12 is the experimental/publication value, so it is in the
  // default parameter script, and 1.5 * 2.0238e12 is model adjustment and was
  // added to the simulation script (Weiliang) the 1.5 multiplier is set at the
  // beginning of the python script in steps 3
  static constexpr osh::Real BK_ro() { return 1.5 * 2.0238e12; }
  static constexpr osh::Real BK_rev() { return -77e-3; }
  static constexpr osh::Real Qo() { return 0.73; }
  static constexpr osh::Real Qc() { return -0.67; }

  static constexpr osh::Real pf0() { return 2.39; }
  static constexpr osh::Real pf1() { return 5.4918; }
  static constexpr osh::Real pf2() { return 24.6205; }
  static constexpr osh::Real pf3() { return 142.4546; }
  static constexpr osh::Real pf4() { return 211.0220; }

  static constexpr osh::Real pb0() { return 3936.0; }
  static constexpr osh::Real pb1() { return 687.3251; }
  static constexpr osh::Real pb2() { return 234.5875; }
  static constexpr osh::Real pb3() { return 103.2204; }
  static constexpr osh::Real pb4() { return 11.6581; }

  static constexpr osh::Real k1() { return 1.0e6; }
  static constexpr osh::Real onoffrate() { return 1.0e3; }
  static constexpr osh::Real L0() { return 1806.0; }

  static constexpr osh::Real Kc() { return 8.63e-6; }
  static constexpr osh::Real Ko() { return 0.6563e-6; }

  osh::Real c_01() { return 4.0 * k1() * onoffrate() * Qt_mslo(); }
  osh::Real c_12() { return 3.0 * k1() * onoffrate() * Qt_mslo(); }
  osh::Real c_23() { return 2.0 * k1() * onoffrate() * Qt_mslo(); }
  osh::Real c_34() { return 1.0 * k1() * onoffrate() * Qt_mslo(); }
  osh::Real o_01() { return 4.0 * k1() * onoffrate() * Qt_mslo(); }
  osh::Real o_12() { return 3.0 * k1() * onoffrate() * Qt_mslo(); }
  osh::Real o_23() { return 2.0 * k1() * onoffrate() * Qt_mslo(); }
  osh::Real o_34() { return 1.0 * k1() * onoffrate() * Qt_mslo(); }

  osh::Real c_10() { return 1.0 * Kc() * k1() * onoffrate() * Qt_mslo(); }
  osh::Real c_21() { return 2.0 * Kc() * k1() * onoffrate() * Qt_mslo(); }
  osh::Real c_32() { return 3.0 * Kc() * k1() * onoffrate() * Qt_mslo(); }
  osh::Real c_43() { return 4.0 * Kc() * k1() * onoffrate() * Qt_mslo(); }
  osh::Real o_10() { return 1.0 * Ko() * k1() * onoffrate() * Qt_mslo(); }
  osh::Real o_21() { return 2.0 * Ko() * k1() * onoffrate() * Qt_mslo(); }
  osh::Real o_32() { return 3.0 * Ko() * k1() * onoffrate() * Qt_mslo(); }
  osh::Real o_43() { return 4.0 * Ko() * k1() * onoffrate() * Qt_mslo(); }

  osh::Real f_0(osh::Real mV) {
    return pf0() * Qt_mslo() *
           (std::exp((Qo() * FARADAY() * mV) /
                     (R() * (TEMPERATURE() + 273.15))));
  }
  osh::Real f_1(osh::Real mV) {
    return pf1() * Qt_mslo() *
           (std::exp((Qo() * FARADAY() * mV) /
                     (R() * (TEMPERATURE() + 273.15))));
  }
  osh::Real f_2(osh::Real mV) {
    return pf2() * Qt_mslo() *
           (std::exp((Qo() * FARADAY() * mV) /
                     (R() * (TEMPERATURE() + 273.15))));
  }
  osh::Real f_3(osh::Real mV) {
    return pf3() * Qt_mslo() *
           (std::exp((Qo() * FARADAY() * mV) /
                     (R() * (TEMPERATURE() + 273.15))));
  }
  osh::Real f_4(osh::Real mV) {
    return pf4() * Qt_mslo() *
           (std::exp((Qo() * FARADAY() * mV) /
                     (R() * (TEMPERATURE() + 273.15))));
  }

  osh::Real b_0(osh::Real mV) {
    return pb0() * Qt_mslo() *
           (std::exp((Qc() * FARADAY() * mV) /
                     (R() * (TEMPERATURE() + 273.15))));
  }
  osh::Real b_1(osh::Real mV) {
    return pb1() * Qt_mslo() *
           (std::exp((Qc() * FARADAY() * mV) /
                     (R() * (TEMPERATURE() + 273.15))));
  }
  osh::Real b_2(osh::Real mV) {
    return pb2() * Qt_mslo() *
           (std::exp((Qc() * FARADAY() * mV) /
                     (R() * (TEMPERATURE() + 273.15))));
  }
  osh::Real b_3(osh::Real mV) {
    return pb3() * Qt_mslo() *
           (std::exp((Qc() * FARADAY() * mV) /
                     (R() * (TEMPERATURE() + 273.15))));
  }
  osh::Real b_4(osh::Real mV) {
    return pb4() * Qt_mslo() *
           (std::exp((Qc() * FARADAY() * mV) /
                     (R() * (TEMPERATURE() + 273.15))));
  }

  static constexpr osh::Real BK_C0_p() { return 0.99997; }
  static constexpr osh::Real BK_C1_p() { return 4.3619e-07; }
  static constexpr osh::Real BK_C2_p() { return 4.1713e-09; }
  static constexpr osh::Real BK_C3_p() { return 4.4449e-11; }
  static constexpr osh::Real BK_C4_p() { return 6.3132e-14; }

  static constexpr osh::Real BK_O0_p() { return 2.5202e-05; }
  static constexpr osh::Real BK_O1_p() { return 1.1765e-06; }
  static constexpr osh::Real BK_O2_p() { return 6.6148e-08; }
  static constexpr osh::Real BK_O3_p() { return 2.4392e-09; }
  static constexpr osh::Real BK_O4_p() { return 4.0981e-11; }

  static constexpr osh::Real SK_G() { return 1.0e-11; }
  static constexpr osh::Real SK_ro() { return 31.0e10; }
  static constexpr osh::Real SK_rev() { return -77e-3; }

  static constexpr osh::Real invc1() { return 80.0; }
  static constexpr osh::Real invc2() { return 80.0; }
  static constexpr osh::Real invc3() { return 200.0; }

  static constexpr osh::Real invo1() { return 1000.0; }
  static constexpr osh::Real invo2() { return 100.0; }

  static constexpr osh::Real diro1() { return 160.0; }
  static constexpr osh::Real diro2() { return 1200.0; }
  static constexpr osh::Real dirc2() { return 200.0e6; }
  static constexpr osh::Real dirc3() { return 160.0e6; }
  static constexpr osh::Real dirc4() { return 80.0e6; }

  osh::Real invc1_t() { return invc1() * Qt(); }
  osh::Real invc2_t() { return invc2() * Qt(); }
  osh::Real invc3_t() { return invc3() * Qt(); }

  osh::Real invo1_t() { return invo1() * Qt(); }
  osh::Real invo2_t() { return invo2() * Qt(); }

  osh::Real diro1_t() { return diro1() * Qt(); }
  osh::Real diro2_t() { return diro2() * Qt(); }

  osh::Real dirc2_t() { return dirc2() * Qt() / 3.0; }
  osh::Real dirc3_t() { return dirc3() * Qt() / 3.0; }
  osh::Real dirc4_t() { return dirc4() * Qt() / 3.0; }

  static constexpr osh::Real SK_C1_p() { return 0.96256; }
  static constexpr osh::Real SK_C2_p() { return 0.036096; }
  static constexpr osh::Real SK_C3_p() { return 0.0010829; }
  static constexpr osh::Real SK_C4_p() { return 6.4973e-06; }

  static constexpr osh::Real SK_O1_p() { return 0.00017326; }
  static constexpr osh::Real SK_O2_p() { return 7.7967e-05; }

  // AMPA rate constants
  static constexpr osh::Real AMPA_G() { return 7e-12; }
  static constexpr osh::Real AMPA_TotalG() { return 500e-9; }
  static constexpr osh::Real AMPA_receptors() {
    return AMPA_TotalG() / AMPA_G();
  }
  static constexpr osh::Real AMPA_rev() { return 0.0e3; }
  static constexpr osh::Real rb() { return 13e6; }
  static constexpr osh::Real ru1() { return 0.0059e3; }
  static constexpr osh::Real ru2() { return 86e3; }
  static constexpr osh::Real ro() { return 2.7e3; }
  static constexpr osh::Real rc() { return 0.2e3; }
  static constexpr osh::Real rd() { return 0.9e3; }
  static constexpr osh::Real rr() { return 0.064e3; }

  // leak current channel density & conductance per channel
  static constexpr osh::Real Rm() { return 122e3 / 1.0e4; }
  static constexpr osh::Real g_leak_proximal() { return 1.2 / Rm(); }
  static constexpr osh::Real g_leak_spiny() { return 5.3 / Rm(); }
  static constexpr osh::Real L_ro_proximal() { return 25.0e10; }
  // return different L_G value based on if the model is the full model
  // or the background model
  static constexpr osh::Real L_G(bool is_full_model) {
    if (is_full_model)
      return g_leak_proximal() / L_ro_proximal();
    else
      return 4.0e-14;
  }
  static constexpr osh::Real L_ro_spiny() { return g_leak_spiny() / L_G(true); }
  static constexpr osh::Real L_ro() { return 25.0e10; }
  static constexpr osh::Real L_rev() { return -61e-3; }

  // Pump parameters
  static constexpr osh::Real P_f_kcst() { return 3e9; }
  static constexpr osh::Real P_b_kcst() { return 1.75e4; }
  static constexpr osh::Real P_k_kcst() { return 7.255e4; }

  // CALCIUM BUFFERING MODEL
  // Ca concentrations
  static constexpr osh::Real Ca_oconc() { return 2e-3; }
  static constexpr osh::Real Ca_iconc() { return 45e-9; }

  // Mg concentrations
  static constexpr osh::Real Mg_conc() { return 590e-6; }

  // Buffer concentrations
  static constexpr osh::Real iCBsf_conc() { return 27.704e-6; }
  static constexpr osh::Real iCBCaf_conc() { return 2.6372e-6; }
  static constexpr osh::Real iCBsCa_conc() { return 1.5148e-6; }
  static constexpr osh::Real iCBCaCa_conc() { return 0.14420e-6; }

  static constexpr osh::Real CBsf_conc() { return 110.82e-6; }
  static constexpr osh::Real CBCaf_conc() { return 10.549e-6; }
  static constexpr osh::Real CBsCa_conc() { return 6.0595e-6; }
  static constexpr osh::Real CBCaCa_conc() { return 0.57682e-6; }

  static constexpr osh::Real PV_conc() { return 3.2066e-6; }
  static constexpr osh::Real PVCa_conc() { return 16.252e-6; }
  static constexpr osh::Real PVMg_conc() { return 60.541e-6; }

  static constexpr osh::Real DCST() { return 0.223e-9; }
  static constexpr osh::Real DCB() { return 0.028e-9; }
  static constexpr osh::Real DPV() { return 0.043e-9; }

  static constexpr osh::Real iCBsf1_f_kcst() { return 4.35e7; }
  static constexpr osh::Real iCBsf1_b_kcst() { return 35.8; }

  static constexpr osh::Real iCBsCa_f_kcst() { return 0.55e7; }
  static constexpr osh::Real iCBsCa_b_kcst() { return 2.6; }

  static constexpr osh::Real iCBsf2_f_kcst() { return 0.55e7; }
  static constexpr osh::Real iCBsf2_b_kcst() { return 2.6; }

  static constexpr osh::Real iCBCaf_f_kcst() { return 4.35e7; }
  static constexpr osh::Real iCBCaf_b_kcst() { return 35.8; }

  static constexpr osh::Real CBsf1_f_kcst() { return 4.35e7; }
  static constexpr osh::Real CBsf1_b_kcst() { return 35.8; }

  static constexpr osh::Real CBsCa_f_kcst() { return 0.55e7; }
  static constexpr osh::Real CBsCa_b_kcst() { return 2.6; }

  static constexpr osh::Real CBsf2_f_kcst() { return 0.55e7; }
  static constexpr osh::Real CBsf2_b_kcst() { return 2.6; }

  static constexpr osh::Real CBCaf_f_kcst() { return 4.35e7; }
  static constexpr osh::Real CBCaf_b_kcst() { return 35.8; }

  static constexpr osh::Real PVca_f_kcst() { return 10.7e7; }
  static constexpr osh::Real PVca_b_kcst() { return 0.95; }

  static constexpr osh::Real PVmg_f_kcst() { return 0.8e6; }
  static constexpr osh::Real PVmg_b_kcst() { return 25.0; }
  // #pumpnbs per unit area (im m2) is Total pump times AVOGADRO's NUMBER (1e-11
  // mol/m2 * 6.022e23 /mol )
  static constexpr osh::Real pumpnbs_per_area() { return 6.022141e12; }
};

struct CaBurstFullSimdef : public CaBurstSimdef {
  CaBurstFullSimdef(const steps::model::Model &model,
                    const steps::dist::DistMesh &mesh);
  container::surface_reaction_id ampacc1, ampac1c2;
};

struct CaBurstIntegrationTestSimdef: public Simdef {
    CaBurstIntegrationTestSimdef(const steps::model::Model& model,
                                 const steps::dist::DistMesh& mesh);

    static constexpr double compConc() {
        return 1000 / (1e-18 * 6.022e23);
    }
    static constexpr double potential() {
        return -65.0e-3;
    }
};

struct Rallpack3Simdef : public Simdef {
  Rallpack3Simdef(const steps::model::Model &model,
                  const steps::dist::DistMesh &mesh);

  // Potassium conductance, Siemens/m^2
  static constexpr osh::Real K_G() { return 360.0; }

  // Sodium conductance, Siemens/m^2
  static constexpr osh::Real Na_G() { return 1200.0; }

  // Leak conductance, Siemens/m^2
  static constexpr osh::Real L_G() { return 0.25; }

  // Potassium reversal potential, V
  static constexpr osh::Real K_rev() { return -77e-3; }

  static constexpr osh::Real Na_rev() { return 50e-3; }

  static constexpr osh::Real leak_rev() { return -65.0e-3; }

  static constexpr osh::Real K_ro() { return 18.0e12; }

  static constexpr osh::Real Na_ro() { return 60.0e12; }

  // Total leak conductance for ideal cylinder:
  static constexpr osh::Real L_G_tot() {
    auto surfarea_cyl = 1.0 * 3.141592653589 * 1000 * 1e-12;
    return L_G() * surfarea_cyl;
  }

  // A table of potassium density factors at -65mV, found in getpops. n0, n1,
  // n2, n3, n4
  static std::array<osh::Real, 5> K_FACS() {
    return {0.216750577045, 0.40366011853, 0.281904943772, 0.0874997924409,
            0.0101845682113};
  }

  // A table of sodium density factors. m0h1, m1h1, m2h1, m3h1, m0h0, m1h0,
  // m2h0, m3h0
  static std::array<osh::Real, 8> NA_FACS() {
    return {0.343079175644,    0.0575250437508,  0.00321512825945,
            5.98988373918e-05, 0.506380603793,   0.0849062503811,
            0.00474548939393,  8.84099403236e-05};
  }

  static constexpr osh::Real Ra() { return 1.0; }

  static constexpr osh::Real SIM_END() { return 0.1; }

  static constexpr osh::Real Iinj() { return 0.1e-9; }

  osh::Real _a_m(osh::Real mV) {
    return ((
        ((0.1 * (25 - (mV + 65.)) / (std::exp((25 - (mV + 65.)) / 10.) - 1)))));
  }

  osh::Real _b_m(osh::Real mV) {
    return ((((4. * std::exp(-((mV + 65.) / 18.))))));
  }

  osh::Real _a_h(osh::Real mV) {
    return ((((0.07 * std::exp((-(mV + 65.) / 20.))))));
  }

  osh::Real _b_h(osh::Real mV) {
    return ((((1. / (std::exp((30 - (mV + 65.)) / 10.) + 1)))));
  }

  osh::Real _a_n(osh::Real mV) {
    return (((
        (0.01 * (10 - (mV + 65.)) / (std::exp((10 - (mV + 65.)) / 10.) - 1)))));
  }

  osh::Real _b_n(osh::Real mV) {
    return ((((0.125 * std::exp(-(mV + 65.) / 80.)))));
  }
};

} // namespace dist
} // namespace steps
