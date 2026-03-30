#pragma once

#include <Omega_h_adj.hpp>
#include <type_traits>

#include "geom/dist/distmesh.hpp"
#include "kproc/diffusions.hpp"
#include "math/distributions.hpp"
#include "mol_state.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"
#include "operator/diffusion_operator.hpp"
#include "operator/rleaping_operator.hpp"
#include "operator/rssa_operator.hpp"
#include "operator/ssa_operator.hpp"
#include "rng/rng.hpp"
#include "simulation_data.hpp"
#include "util/mpitools.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

struct MembraneResistivity {
    MembraneResistivity() = default;
    MembraneResistivity(osh::Real t_resistivity, osh::Real t_reversal_potential)
        : resistivity(t_resistivity)
        , reversal_potential(t_reversal_potential) {}
    osh::Real resistivity;
    osh::Real reversal_potential;
};

class Simulation {
  public:
    Simulation(DistMesh& t_mesh, rng::RNG& t_rng);
    virtual ~Simulation() noexcept;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Required methods
    ///////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////
    // General methods //
    /////////////////////

    // Simulation control
    virtual void reset() = 0;
    virtual void run(osh::Real end_time) = 0;

    // Data getting / setting
    virtual osh::Real getTime() const noexcept = 0;
    virtual SSAMethod ssaMethod() const noexcept = 0;
    virtual osh::Real getTemp() const noexcept = 0;
    virtual void setTemp(const osh::Real temp) noexcept = 0;
    virtual osh::Real getDiffusionTolerance() const = 0;
    virtual void setDiffusionTolerance(osh::Real tolerance) = 0;
    virtual osh::Real getDiffusionNormalApproximationThreshold() const = 0;
    virtual void setDiffusionNormalApproximationThreshold(osh::Real threshold) = 0;
    virtual osh::Real getDiffusionCrankNicolsonThreshold() const = 0;
    virtual void setDiffusionCrankNicolsonThreshold(osh::Real threshold) = 0;
    virtual uint getDiffusionLeapThreshold() const = 0;
    virtual void setDiffusionLeapThreshold(uint leap_thresh) = 0;
    virtual uint getDiffusionMaxDtSkips() const = 0;
    virtual void setDiffusionMaxDtSkips(uint max_skips) = 0;
    virtual osh::Real getDiffusionMinDtFactor() const = 0;
    virtual void setDiffusionMinDtFactor(osh::Real factor) = 0;
    virtual uint getReactionSSAThreshold() const = 0;
    virtual void setReactionSSAThreshold(uint thresh) = 0;
    virtual uint getReactionSSASteps() const = 0;
    virtual void setReactionSSASteps(uint steps) = 0;
    virtual uint getReactionLComputePeriod() const = 0;
    virtual void setReactionLComputePeriod(uint period) = 0;
    virtual osh::Real getReactionTolerance() const = 0;
    virtual void setReactionTolerance(osh::Real tolerance) = 0;
    virtual osh::Real getReactionTheta() const = 0;
    virtual void setReactionTheta(osh::Real theta) = 0;
    virtual std::string getSolverName() const = 0;
    virtual void setDiffApplyThreshold(osh::Real threshold) = 0;

#if USE_PETSC
    // E-field specific

    virtual osh::Real getEfieldDT() const = 0;
    virtual void setEfieldDT(const osh::Real dt) const = 0;

    virtual void setPetscOptions(const std::string& s) = 0;

#endif  // USE_PETSC

    // Debugging / monitoring
    virtual void dumpDepGraphToFile(const std::string& path) const = 0;
    virtual std::string createStateReport() const = 0;
    virtual osh::I64 getDiffExtent(bool local = false) const = 0;
    virtual osh::I64 getReacExtent(bool local = false) const = 0;
    virtual double getEFieldTime() const noexcept = 0;
    virtual double getRDTime() const noexcept = 0;
    virtual double getDiffusionTime() const noexcept = 0;
    virtual double getReactionTime() const noexcept = 0;
    virtual std::map<std::string, double> getReactionDebugInfo(bool local = false) const = 0;
    virtual std::map<std::string, double> getDiffusionDebugInfo(bool local = false) const = 0;

    ///////////////////////////
    // Location: Compartment //
    ///////////////////////////

    // Species
    virtual osh::Real getCompSpecCount(const model::compartment_id& compartment,
                                       const model::species_name& species) const = 0;
    virtual void setCompSpecCount(const model::compartment_id& compartment,
                                  const model::species_name& spec,
                                  osh::Real n,
                                  const math::DistributionMethod distribution) = 0;

    virtual osh::Real getCompSpecConc(const model::compartment_id& compartment,
                                      const model::species_name& species) const = 0;
    virtual void setCompSpecConc(const model::compartment_id& compartment,
                                 const model::species_name& spec,
                                 osh::Real conc,
                                 const math::DistributionMethod distribution) = 0;

    virtual bool getCompSpecClamped(const model::compartment_id& compartment,
                                    const model::species_name& spec) const = 0;
    virtual void setCompSpecClamped(const model::compartment_id& compartment,
                                    const model::species_name& spec,
                                    bool clamped) = 0;

    // Reactions
    virtual osh::Real getCompReacK(const model::compartment_id& compartment,
                                   const model::reaction_id& reac) const = 0;
    virtual void setCompReacK(const model::compartment_id& compartment,
                              const model::reaction_id& reac,
                              osh::Real kcst) = 0;

    virtual osh::I64 getCompReacExtent(const model::compartment_id& compartment,
                                       const model::reaction_id& reac) const = 0;

    virtual osh::I64 getCompComplexReacExtent(const model::compartment_id& compartment,
                                              const model::complex_reaction_id& reac) const = 0;

    // Diffusions
    virtual osh::Real getCompDiffD(const model::compartment_id& compartment,
                                   const model::diffusion_id& diff) const = 0;
    virtual void setCompDiffD(const model::compartment_id& compartment,
                              const model::diffusion_id& diff,
                              osh::Real dcst) = 0;

    // Complexes
    virtual osh::Real getCompComplexCount(
        const model::compartment_id& compartment,
        const model::complex_name& complex,
        const std::vector<std::vector<steps::model::SubunitStateFilter>>& f) const = 0;
    virtual void setCompComplexCount(
        const model::compartment_id& compartment,
        const model::complex_name& complex,
        const std::vector<std::vector<steps::model::SubunitStateFilter>>& i,
        osh::Real num_molecules,
        math::DistributionMethod distribution) = 0;

    virtual osh::Real getCompComplexSUSCount(
        const model::compartment_id& compartment,
        const model::complex_name& complex,
        const std::vector<std::vector<steps::model::SubunitStateFilter>>& f,
        model::complex_substate_id m) const = 0;

    /////////////////////
    // Location: Patch //
    /////////////////////

    // Species
    virtual osh::Real getPatchSpecCount(const model::patch_id& compartment,
                                        const model::species_name& species) const = 0;
    virtual void setPatchSpecCount(const model::patch_id& patch_id,
                                   const model::species_name& species,
                                   osh::Real num_molecules,
                                   const math::DistributionMethod) = 0;

    virtual bool getPatchSpecClamped(const model::patch_id& patch,
                                     const model::species_name& spec) const = 0;
    virtual void setPatchSpecClamped(const model::patch_id& patch,
                                     const model::species_name& spec,
                                     bool clamped) = 0;

    // Complexes
    virtual osh::Real getPatchComplexCount(
        const model::patch_id& patch,
        const model::complex_name& complex,
        const std::vector<std::vector<steps::model::SubunitStateFilter>>& f) const = 0;
    virtual void setPatchComplexCount(
        const model::patch_id& patch,
        const model::complex_name& complex,
        const std::vector<std::vector<steps::model::SubunitStateFilter>>& i,
        osh::Real num_molecules,
        math::DistributionMethod distribution) = 0;

    virtual osh::Real getPatchComplexSUSCount(
        const model::patch_id& cidx,
        const model::complex_name& complex,
        const std::vector<std::vector<steps::model::SubunitStateFilter>>& f,
        model::complex_substate_id m) const = 0;

    // Reactions
    virtual osh::Real getPatchSReacK(const model::patch_id& patchId,
                                     const model::surface_reaction_id& reactionId) const = 0;
    virtual void setPatchSReacK(const model::patch_id& patchId,
                                const model::surface_reaction_id& reactionId,
                                osh::Real kCst) = 0;

    virtual osh::I64 getPatchSReacExtent(const model::patch_id& patch,
                                         const model::surface_reaction_id& reac) const = 0;

    virtual osh::I64 getPatchComplexSReacExtent(
        const model::patch_id& patch,
        const model::complex_surface_reaction_id& reac) const = 0;

    virtual osh::I64 getPatchVDepSReacExtent(const model::patch_id& patch,
                                             const model::vdep_surface_reaction_id& reac) const = 0;

    virtual osh::I64 getPatchVDepComplexSReacExtent(
        const model::patch_id& patch,
        const model::vdep_complex_surface_reaction_id& reac) const = 0;

    ////////////////////////
    // Location: Membrane //
    ////////////////////////

#if USE_PETSC
    // E-field values
    virtual MembraneResistivity getMembRes(const model::membrane_id& membrane) const = 0;
    virtual void setMembRes(const model::membrane_id& membrane,
                            osh::Real resistivity,
                            osh::Real reversal_potential) = 0;

    virtual void setMembVolRes(const model::membrane_id& membrane, double ro) = 0;

    virtual void setMembCapac(const model::membrane_id& membrane, double capacitance) = 0;

    virtual void setMembPotential(const model::membrane_id& memb, osh::Real value) = 0;

    virtual void setMembIClamp(const model::membrane_id& membrane, osh::Real current) = 0;

#endif  // USE_PETSC

    //////////////////////////////////
    // Location: Diffusion boundary //
    //////////////////////////////////

    virtual bool getDiffBoundarySpecDiffusionActive(
        const mesh::diffusion_boundary_name& diffusion_boundary_name,
        const model::species_name& spec_id) const = 0;
    virtual void setDiffBoundarySpecDiffusionActive(
        const mesh::diffusion_boundary_name& diffusion_boundary_name,
        const model::species_name& spec_id,
        bool set_active) = 0;

    virtual void setDiffBoundarySpecDcst(const mesh::diffusion_boundary_name& diffb,
                                         const model::species_name& spec,
                                         osh::Real dcst) = 0;

    ///////////////////////////
    // Location: Tetrahedron //
    ///////////////////////////

    // Species
    virtual void getBatchTetSpecCountsNP(const osh::GO* indices,
                                         size_t input_size,
                                         const model::species_name& species,
                                         double* counts,
                                         size_t output_size,
                                         bool local) const = 0;
    virtual void setBatchTetSpecCountsNP(const osh::GO* indices,
                                         size_t input_size,
                                         const model::species_name& s,
                                         double* counts,
                                         size_t output_size,
                                         bool local) = 0;

    virtual void getBatchTetSpecConcsNP(const osh::GO* indices,
                                        size_t input_size,
                                        const model::species_name& species,
                                        double* counts,
                                        size_t output_size,
                                        bool local) const = 0;
    virtual void setBatchTetSpecConcsNP(const osh::GO* indices,
                                        size_t input_size,
                                        const model::species_name& s,
                                        double* concs,
                                        size_t output_size,
                                        bool local) = 0;

    virtual bool getTetSpecClamped(osh::GO tet,
                                   const model::species_name& s,
                                   bool local = false) const = 0;
    virtual void setTetSpecClamped(osh::GO tet,
                                   const model::species_name& s,
                                   bool clamped,
                                   bool local = false) = 0;

    // Reactions
    virtual osh::Real getTetReacK(osh::GO tet,
                                  const model::reaction_id reac,
                                  bool local = false) const = 0;
    virtual void setTetReacK(osh::GO tet,
                             const model::reaction_id reac,
                             osh::Real kcst,
                             bool local = false) = 0;

    virtual osh::Real getTetComplexReacK(osh::GO tet,
                                         const model::complex_reaction_id reac,
                                         bool local = false) const = 0;
    virtual void setTetComplexReacK(osh::GO tet,
                                    const model::complex_reaction_id reac,
                                    osh::Real kcst,
                                    bool local = false) = 0;

    // Diffusions
    virtual osh::Real getTetDiffD(osh::GO tet,
                                  const model::diffusion_id diff,
                                  osh::GO direc_tet,
                                  bool local = false) const = 0;
    virtual void setTetDiffD(osh::GO tet,
                             const model::diffusion_id diff,
                             double dcst,
                             osh::GO direc_tet,
                             bool local = false) = 0;

    // E-field values
#if USE_PETSC
    virtual void getBatchTetVsNP(const osh::GO* indices,
                                 size_t input_size,
                                 osh::Real* voltages,
                                 size_t output_size,
                                 bool local = false) const = 0;
    virtual void setBatchTetVsNP(const osh::GO* indices,
                                 size_t input_size,
                                 osh::Real* voltages,
                                 size_t output_size,
                                 bool local = false) = 0;

    virtual bool getTetVClamped(osh::GO vertex, bool local = false) const = 0;
    virtual void setTetVClamped(osh::GO vertex, bool clamped, bool local = false) = 0;
#endif  // USE_PETSC

    ////////////////////////
    // Location: Triangle //
    ////////////////////////

    // Species
    virtual void getBatchTriSpecCountsNP(const osh::GO* indices,
                                         size_t input_size,
                                         const model::species_name& s,
                                         double* counts,
                                         size_t output_size,
                                         bool local) const = 0;

    virtual void setBatchTriSpecCountsNP(const osh::GO* indices,
                                         size_t input_size,
                                         const model::species_name& s,
                                         double* counts,
                                         size_t output_size,
                                         bool local) = 0;

    virtual bool getTriSpecClamped(osh::GO tri,
                                   const model::species_name& s,
                                   bool local = false) const = 0;
    virtual void setTriSpecClamped(osh::GO tri,
                                   const model::species_name& s,
                                   bool clamped,
                                   bool local = false) = 0;

    // Reactions
    virtual osh::Real getTriSReacK(osh::GO triangle,
                                   const model::surface_reaction_id& reactionId,
                                   bool local) const = 0;
    virtual void setTriSReacK(osh::GO triangle,
                              const model::surface_reaction_id& reactionId,
                              osh::Real kCst,
                              bool local) = 0;

    virtual osh::Real getTriComplexSReacK(osh::GO triangle,
                                          const model::complex_surface_reaction_id& reactionId,
                                          bool local) const = 0;
    virtual void setTriComplexSReacK(osh::GO triangle,
                                     const model::complex_surface_reaction_id& reactionId,
                                     osh::Real kCst,
                                     bool local) = 0;

    // E-field values
#if USE_PETSC
    virtual void getBatchTriVsNP(const osh::GO* indices,
                                 size_t input_size,
                                 osh::Real* voltages,
                                 size_t output_size,
                                 bool local = false) const = 0;
    virtual void setBatchTriVsNP(const osh::GO* indices,
                                 size_t input_size,
                                 osh::Real* voltages,
                                 size_t output_size,
                                 bool local = false) = 0;

    virtual void getBatchTriSReacIsNP(const osh::GO* indices,
                                      size_t input_size,
                                      const model::surface_reaction_id reac,
                                      osh::Real* currents,
                                      size_t output_size,
                                      bool local = false) const = 0;

    virtual void getBatchTriComplexSReacIsNP(const osh::GO* indices,
                                             size_t input_size,
                                             const model::complex_surface_reaction_id reac,
                                             osh::Real* currents,
                                             size_t output_size,
                                             bool local = false) const = 0;

    virtual void getBatchTriVDepSReacIsNP(const osh::GO* indices,
                                          size_t input_size,
                                          const model::vdep_surface_reaction_id reac,
                                          osh::Real* currents,
                                          size_t output_size,
                                          bool local = false) const = 0;

    virtual void getBatchTriVDepComplexSReacIsNP(const osh::GO* indices,
                                                 size_t input_size,
                                                 const model::vdep_complex_surface_reaction_id reac,
                                                 osh::Real* currents,
                                                 size_t output_size,
                                                 bool local = false) const = 0;

    virtual void getBatchTriOhmicIsNP(const osh::GO* indices,
                                      size_t input_size,
                                      const model::ohmic_current_id curr,
                                      osh::Real* currents,
                                      size_t output_size,
                                      bool local = false) const = 0;

    virtual void getBatchTriComplexOhmicIsNP(const osh::GO* indices,
                                             size_t input_size,
                                             const model::complex_ohmic_current_id curr,
                                             osh::Real* currents,
                                             size_t output_size,
                                             bool local = false) const = 0;

    virtual void getBatchTriGHKIsNP(const osh::GO* indices,
                                    size_t input_size,
                                    const model::ghk_current_id curr,
                                    osh::Real* currents,
                                    size_t output_size,
                                    bool local = false) const = 0;

    virtual void getBatchTriComplexGHKIsNP(const osh::GO* indices,
                                           size_t input_size,
                                           const model::complex_ghk_current_id curr,
                                           osh::Real* currents,
                                           size_t output_size,
                                           bool local = false) const = 0;

    virtual void getBatchTriIsNP(const osh::GO* indices,
                                 size_t input_size,
                                 osh::Real* currents,
                                 size_t output_size,
                                 bool local = false) const = 0;

    virtual void getBatchTriOhmicErevsNP(const osh::GO* indices,
                                         size_t input_size,
                                         const model::ohmic_current_id& ohmic_current,
                                         double* rv,
                                         size_t output_size,
                                         bool local) const = 0;

    virtual void getBatchTriComplexOhmicErevsNP(
        const osh::GO* indices,
        size_t input_size,
        const model::complex_ohmic_current_id& ohmic_current,
        double* rv,
        size_t output_size,
        bool local) const = 0;

    virtual void setTriOhmicErev(osh::GO triangle,
                                 const model::ohmic_current_id& ohmic_current,
                                 double reversal_potential,
                                 bool local) = 0;

    virtual void setTriComplexOhmicErev(osh::GO triangle,
                                        const model::complex_ohmic_current_id& ohmic_current,
                                        double reversal_potential,
                                        bool local) = 0;

    virtual bool getTriVClamped(osh::GO vertex, bool local = false) const = 0;
    virtual void setTriVClamped(osh::GO vertex, bool clamped, bool local = false) = 0;

    virtual osh::Real getTriIClamp(osh::GO tri, bool local = false) const = 0;
    virtual void setTriIClamp(osh::GO tri, osh::Real current, bool local = false) = 0;

    virtual MembraneResistivity getTriRes(osh::GO tri, bool local = false) const = 0;
    virtual void setTriRes(const osh::GO tri,
                           osh::Real res,
                           osh::Real erev,
                           bool local = false) = 0;

    virtual osh::Real getTriCapac(osh::GO tri, bool local = false) const = 0;
    virtual void setTriCapac(const osh::GO tri, osh::Real c, bool local = false) = 0;

#endif  // USE_PETSC

    //////////////////////
    // Location: Vertex //
    //////////////////////

#if USE_PETSC
    // E-field value
    virtual void getBatchVertVsNP(const osh::GO* indices,
                                  size_t input_size,
                                  osh::Real* voltages,
                                  size_t output_size,
                                  bool local = false) const = 0;
    virtual void setBatchVertVsNP(const osh::GO* indices,
                                  size_t input_size,
                                  osh::Real* voltages,
                                  size_t output_size,
                                  bool local = false) = 0;

    virtual bool getVertVClamped(osh::GO vertex, bool local = false) const = 0;
    virtual void setVertVClamped(osh::GO vertex, bool clamped, bool local = false) = 0;

    virtual osh::Real getVertIClamp(osh::GO vertex, bool local = false) const = 0;
    virtual void setVertIClamp(osh::GO vertex, osh::Real current, bool local = false) = 0;

#endif  // USE_PETSC

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Convenience methods
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////
    // Location: Tetrahedron //
    ///////////////////////////
    double getTetSpecCount(osh::GO tet, const model::species_name& s, bool local) const;
    void setTetSpecCount(osh::GO tet, const model::species_name& s, double count, bool local);

    double getTetSpecConc(osh::GO tet, const model::species_name& s, bool local) const;
    void setTetSpecConc(osh::GO tet, const model::species_name& s, double conc, bool local);

#if USE_PETSC
    double getTetV(osh::GO tet, bool local) const;
    void setTetV(osh::GO tet, double v, bool local);
#endif  // USE_PETSC

    ////////////////////////
    // Location: Triangle //
    ////////////////////////
    double getTriSpecCount(osh::GO tri, const model::species_name& s, bool local) const;
    void setTriSpecCount(osh::GO tri, const model::species_name& s, double count, bool local);

#if USE_PETSC
    double getTriV(osh::GO tri, bool local) const;
    void setTriV(osh::GO tri, double v, bool local);
    double getTriOhmicErev(osh::GO tri, const model::ohmic_current_id& curr, bool local) const;
    double getTriComplexOhmicErev(osh::GO tri,
                                  const model::complex_ohmic_current_id& curr,
                                  bool local) const;
    double getTriSReacI(osh::GO tri, const model::surface_reaction_id& reac, bool local) const;
    double getTriComplexSReacI(osh::GO tri,
                               const model::complex_surface_reaction_id& reac,
                               bool local) const;
    double getTriVDepSReacI(osh::GO tri,
                            const model::vdep_surface_reaction_id& reac,
                            bool local) const;
    double getTriVDepComplexSReacI(osh::GO tri,
                                   const model::vdep_complex_surface_reaction_id& reac,
                                   bool local) const;
    double getTriOhmicI(osh::GO tri, const model::ohmic_current_id& curr, bool local) const;
    double getTriComplexOhmicI(osh::GO tri,
                               const model::complex_ohmic_current_id& curr,
                               bool local) const;
    double getTriGHKI(osh::GO tri, const model::ghk_current_id& curr, bool local) const;
    double getTriComplexGHKI(osh::GO tri,
                             const model::complex_ghk_current_id& curr,
                             bool local) const;
    double getTriI(osh::GO tri, bool local) const;
#endif  // USE_PETSC

    //////////////////////
    // Location: Vertex //
    //////////////////////
#if USE_PETSC
    double getVertV(osh::GO vert, bool local) const;
    void setVertV(osh::GO vert, double v, bool local);
#endif  // USE_PETSC

  protected:
    inline MPI_Comm comm() const noexcept {
        return mesh.comm_impl();
    }

    const int comm_rank;
    const int comm_size;

    DistMesh& mesh;
    rng::RNG& rng;

    /// total time in seconds spent performing reactions
    double reactions_timer{};
    /// total time in seconds spent performing diffusions
    double diffusions_timer{};
    /// total time in seconds spent performing efield
    double efield_timer{};
};


///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


template <SSAMethod SSA = SSAMethod::SSA,
          NextEventSearchMethod SearchMethod = NextEventSearchMethod::Direct,
          DiffusionMethod DiffMethod = DiffusionMethod::ConstantDiffDt>
class OmegaHSimulation: public Simulation {
  public:
    using super_type = Simulation;
    using mesh_type = DistMesh;

    OmegaHSimulation(steps::model::Model& model,
                     mesh_type& t_mesh,
                     const rng::RNGptr& r,
                     bool t_indepKProcs,
                     bool isEfield);
    virtual ~OmegaHSimulation();

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Required methods
    ///////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////
    // General methods //
    /////////////////////

    // Simulation control
    void reset() override;
    void run(osh::Real end_time) override;

    // Data getting / setting
    inline osh::Real getTime() const noexcept override {
        return state_time;
    }

    SSAMethod ssaMethod() const noexcept override {
        return SSA;
    }

    inline osh::Real getTemp() const noexcept override {
        return statedef->getTemp();
    }
    inline void setTemp(const osh::Real temp) noexcept override {
        statedef->setTemp(temp);
    }

    osh::Real getDiffusionTolerance() const override;
    void setDiffusionTolerance(osh::Real tolerance) override;

    osh::Real getDiffusionNormalApproximationThreshold() const override;
    void setDiffusionNormalApproximationThreshold(osh::Real threshold) override;

    osh::Real getDiffusionCrankNicolsonThreshold() const override;
    void setDiffusionCrankNicolsonThreshold(osh::Real threshold) override;

    uint getDiffusionLeapThreshold() const override;
    void setDiffusionLeapThreshold(uint leap_thresh) override;

    uint getDiffusionMaxDtSkips() const override;
    void setDiffusionMaxDtSkips(uint max_skips) override;

    osh::Real getDiffusionMinDtFactor() const override;
    void setDiffusionMinDtFactor(osh::Real factor) override;

    uint getReactionSSAThreshold() const override;
    void setReactionSSAThreshold(uint thresh) override;

    uint getReactionSSASteps() const override;
    void setReactionSSASteps(uint steps) override;

    uint getReactionLComputePeriod() const override;
    void setReactionLComputePeriod(uint period) override;

    osh::Real getReactionTolerance() const override;
    void setReactionTolerance(osh::Real tolerance) override;

    osh::Real getReactionTheta() const override;
    void setReactionTheta(osh::Real theta) override;

    std::string getSolverName() const override;
    void setDiffApplyThreshold(osh::Real threshold) override;

#if USE_PETSC
    // E-field specific

    osh::Real getEfieldDT() const override;
    void setEfieldDT(const osh::Real dt) const override;

    void setPetscOptions(const std::string& s) override;

#endif  // USE_PETSC

    // Debugging / monitoring

    /// Dump the dependency graph of kproc in a file specified by path
    void dumpDepGraphToFile(const std::string& path) const override;

    std::string createStateReport() const override;

    osh::I64 getDiffExtent(bool local = false) const override;
    osh::I64 getReacExtent(bool local = false) const override;

    double getEFieldTime() const noexcept override {
        return efield_timer;
    }
    double getRDTime() const noexcept override {
        return diffusions_timer + reactions_timer;
    }
    double getDiffusionTime() const noexcept override {
        return diffusions_timer;
    }
    double getReactionTime() const noexcept override {
        return reactions_timer;
    }
    std::map<std::string, double> getReactionDebugInfo(bool local = false) const override;
    std::map<std::string, double> getDiffusionDebugInfo(bool local = false) const override;

    ///////////////////////////
    // Location: Compartment //
    ///////////////////////////

    // Species
    osh::Real getCompSpecCount(const model::compartment_id& compartment,
                               const model::species_name& species) const override;
    void setCompSpecCount(const model::compartment_id& compartment,
                          const model::species_name& spec,
                          osh::Real n,
                          const math::DistributionMethod distribution) override;
    osh::Real getCompSpecConc(const model::compartment_id& compartment,
                              const model::species_name& species) const override;
    void setCompSpecConc(const model::compartment_id& compartment,
                         const model::species_name& spec,
                         osh::Real conc,
                         const math::DistributionMethod distribution) override;

    bool getCompSpecClamped(const model::compartment_id& compartment,
                            const model::species_name& spec) const override;
    void setCompSpecClamped(const model::compartment_id& compartment,
                            const model::species_name& spec,
                            bool clamped) override;

    // Reactions
    osh::Real getCompReacK(const model::compartment_id& compartment,
                           const model::reaction_id& reac) const override;
    void setCompReacK(const model::compartment_id& compartment,
                      const model::reaction_id& reac,
                      osh::Real kcst) override;

    osh::I64 getCompReacExtent(const model::compartment_id& compartment,
                               const model::reaction_id& reac) const override;

    osh::I64 getCompComplexReacExtent(const model::compartment_id& compartment,
                                      const model::complex_reaction_id& reac) const override;

    // Diffusions
    osh::Real getCompDiffD(const model::compartment_id& compartment,
                           const model::diffusion_id& diff) const override;
    void setCompDiffD(const model::compartment_id& compartment,
                      const model::diffusion_id& diff,
                      osh::Real dcst) override;

    // Complexes
    osh::Real getCompComplexCount(
        const model::compartment_id& compartment,
        const model::complex_name& complex,
        const std::vector<std::vector<steps::model::SubunitStateFilter>>& f) const override;
    void setCompComplexCount(const model::compartment_id& compartment,
                             const model::complex_name& complex,
                             const std::vector<std::vector<steps::model::SubunitStateFilter>>& i,
                             osh::Real num_molecules,
                             math::DistributionMethod distribution) override;

    osh::Real getCompComplexSUSCount(
        const model::compartment_id& compartment,
        const model::complex_name& complex,
        const std::vector<std::vector<steps::model::SubunitStateFilter>>& f,
        model::complex_substate_id m) const override;

    /////////////////////
    // Location: Patch //
    /////////////////////

    // Species
    osh::Real getPatchSpecCount(const model::patch_id& patch,
                                const model::species_name& species) const override;
    void setPatchSpecCount(const model::patch_id& patch,
                           const model::species_name& species,
                           osh::Real num_molecules,
                           const math::DistributionMethod distribution) override;

    bool getPatchSpecClamped(const model::patch_id& patch,
                             const model::species_name& spec) const override;
    void setPatchSpecClamped(const model::patch_id& patch,
                             const model::species_name& spec,
                             bool clamped) override;

    // Complexes
    osh::Real getPatchComplexCount(
        const model::patch_id& patch,
        const model::complex_name& complex,
        const std::vector<std::vector<steps::model::SubunitStateFilter>>& f) const override;
    void setPatchComplexCount(const model::patch_id& patch,
                              const model::complex_name& complex,
                              const std::vector<std::vector<steps::model::SubunitStateFilter>>& i,
                              osh::Real num_molecules,
                              math::DistributionMethod distribution) override;

    osh::Real getPatchComplexSUSCount(
        const model::patch_id& cidx,
        const model::complex_name& complex,
        const std::vector<std::vector<steps::model::SubunitStateFilter>>& f,
        model::complex_substate_id m) const override;

    // Reactions
    osh::Real getPatchSReacK(const model::patch_id& patchId,
                             const model::surface_reaction_id& reactionId) const override;
    void setPatchSReacK(const model::patch_id& patchId,
                        const model::surface_reaction_id& reactionId,
                        osh::Real kCst) override;

    osh::I64 getPatchSReacExtent(const model::patch_id& patch,
                                 const model::surface_reaction_id& reac) const override;

    osh::I64 getPatchComplexSReacExtent(
        const model::patch_id& patch,
        const model::complex_surface_reaction_id& reac) const override;

    osh::I64 getPatchVDepSReacExtent(const model::patch_id& patch,
                                     const model::vdep_surface_reaction_id& reac) const override;

    osh::I64 getPatchVDepComplexSReacExtent(
        const model::patch_id& patch,
        const model::vdep_complex_surface_reaction_id& reac) const override;

    ////////////////////////
    // Location: Membrane //
    ////////////////////////

#if USE_PETSC
    // E-field values
    MembraneResistivity getMembRes(const model::membrane_id& membrane) const override;
    void setMembRes(const model::membrane_id& membrane,
                    osh::Real resistivity,
                    osh::Real reversal_potential) override;

    void setMembVolRes(const model::membrane_id& membrane, double ro) override;

    void setMembCapac(const model::membrane_id& membrane, double capacitance) override;

    void setMembPotential(const model::membrane_id& memb, osh::Real value) override;

    void setMembIClamp(const model::membrane_id& membrane, osh::Real current) override;

#endif  // USE_PETSC

    //////////////////////////////////
    // Location: Diffusion boundary //
    //////////////////////////////////

    bool getDiffBoundarySpecDiffusionActive(
        const mesh::diffusion_boundary_name& diffusion_boundary_name,
        const model::species_name& spec_id) const override;
    void setDiffBoundarySpecDiffusionActive(
        const mesh::diffusion_boundary_name& diffusion_boundary_name,
        const model::species_name& spec_id,
        bool set_active) override;

    void setDiffBoundarySpecDcst(const mesh::diffusion_boundary_name& diffb,
                                 const model::species_name& spec,
                                 osh::Real dcst) override;

    ///////////////////////////
    // Location: Tetrahedron //
    ///////////////////////////

    // Species
    void getBatchTetSpecCountsNP(const osh::GO* indices,
                                 size_t input_size,
                                 const model::species_name& species,
                                 double* counts,
                                 size_t output_size,
                                 bool local) const override;
    void setBatchTetSpecCountsNP(const osh::GO* indices,
                                 size_t input_size,
                                 const model::species_name& s,
                                 double* counts,
                                 size_t output_size,
                                 bool local) override;

    void getBatchTetSpecConcsNP(const osh::GO* indices,
                                size_t input_size,
                                const model::species_name& species,
                                double* counts,
                                size_t output_size,
                                bool local) const override;
    void setBatchTetSpecConcsNP(const osh::GO* indices,
                                size_t input_size,
                                const model::species_name& s,
                                double* concs,
                                size_t output_size,
                                bool local) override;

    bool getTetSpecClamped(osh::GO tet,
                           const model::species_name& s,
                           bool local = false) const override;
    void setTetSpecClamped(osh::GO tet,
                           const model::species_name& s,
                           bool clamped,
                           bool local = false) override;

    // Reactions
    osh::Real getTetReacK(osh::GO tet,
                          const model::reaction_id reac,
                          bool local = false) const override;
    void setTetReacK(osh::GO tet,
                     const model::reaction_id reac,
                     osh::Real kcst,
                     bool local = false) override;

    osh::Real getTetComplexReacK(osh::GO tet,
                                 const model::complex_reaction_id reac,
                                 bool local = false) const override;
    void setTetComplexReacK(osh::GO tet,
                            const model::complex_reaction_id reac,
                            osh::Real kcst,
                            bool local = false) override;

    // Diffusions
    osh::Real getTetDiffD(osh::GO tet,
                          const model::diffusion_id diff,
                          osh::GO direc_tet,
                          bool local = false) const override;
    void setTetDiffD(osh::GO tet,
                     const model::diffusion_id diff,
                     double dcst,
                     osh::GO direc_tet,
                     bool local = false) override;

#if USE_PETSC
    // E-field values
    void getBatchTetVsNP(const osh::GO* indices,
                         size_t input_size,
                         osh::Real* voltages,
                         size_t output_size,
                         bool local = false) const override;
    void setBatchTetVsNP(const osh::GO* indices,
                         size_t input_size,
                         osh::Real* voltages,
                         size_t output_size,
                         bool local = false) override;

    bool getTetVClamped(osh::GO vertex, bool local = false) const override;
    void setTetVClamped(osh::GO vertex, bool clamped, bool local = false) override;
#endif  // USE_PETSC

    ////////////////////////
    // Location: Triangle //
    ////////////////////////

    // Species
    void getBatchTriSpecCountsNP(const osh::GO* indices,
                                 size_t input_size,
                                 const model::species_name& s,
                                 double* counts,
                                 size_t output_size,
                                 bool local) const override;

    void setBatchTriSpecCountsNP(const osh::GO* indices,
                                 size_t input_size,
                                 const model::species_name& s,
                                 double* counts,
                                 size_t output_size,
                                 bool local) override;

    bool getTriSpecClamped(osh::GO tri,
                           const model::species_name& s,
                           bool local = false) const override;
    void setTriSpecClamped(osh::GO tri,
                           const model::species_name& s,
                           bool clamped,
                           bool local = false) override;

    // Reactions
    osh::Real getTriSReacK(osh::GO triangle,
                           const model::surface_reaction_id& reactionId,
                           bool local) const override;
    void setTriSReacK(osh::GO triangle,
                      const model::surface_reaction_id& reactionId,
                      osh::Real kCst,
                      bool local) override;

    osh::Real getTriComplexSReacK(osh::GO triangle,
                                  const model::complex_surface_reaction_id& reactionId,
                                  bool local) const override;
    void setTriComplexSReacK(osh::GO triangle,
                             const model::complex_surface_reaction_id& reactionId,
                             osh::Real kCst,
                             bool local) override;

#if USE_PETSC
    // E-field values
    void getBatchTriVsNP(const osh::GO* indices,
                         size_t input_size,
                         osh::Real* voltages,
                         size_t output_size,
                         bool local = false) const override;
    void setBatchTriVsNP(const osh::GO* indices,
                         size_t input_size,
                         osh::Real* voltages,
                         size_t output_size,
                         bool local = false) override;

    void getBatchTriSReacIsNP(const osh::GO* indices,
                              size_t input_size,
                              const model::surface_reaction_id reac,
                              osh::Real* currents,
                              size_t output_size,
                              bool local = false) const override;

    void getBatchTriComplexSReacIsNP(const osh::GO* indices,
                                     size_t input_size,
                                     const model::complex_surface_reaction_id reac,
                                     osh::Real* currents,
                                     size_t output_size,
                                     bool local = false) const override;

    void getBatchTriVDepSReacIsNP(const osh::GO* indices,
                                  size_t input_size,
                                  const model::vdep_surface_reaction_id reac,
                                  osh::Real* currents,
                                  size_t output_size,
                                  bool local = false) const override;

    void getBatchTriVDepComplexSReacIsNP(const osh::GO* indices,
                                         size_t input_size,
                                         const model::vdep_complex_surface_reaction_id reac,
                                         osh::Real* currents,
                                         size_t output_size,
                                         bool local = false) const override;

    void getBatchTriOhmicIsNP(const osh::GO* indices,
                              size_t input_size,
                              const model::ohmic_current_id curr,
                              osh::Real* currents,
                              size_t output_size,
                              bool local = false) const override;

    void getBatchTriComplexOhmicIsNP(const osh::GO* indices,
                                     size_t input_size,
                                     const model::complex_ohmic_current_id curr,
                                     osh::Real* currents,
                                     size_t output_size,
                                     bool local = false) const override;

    void getBatchTriGHKIsNP(const osh::GO* indices,
                            size_t input_size,
                            const model::ghk_current_id curr,
                            osh::Real* currents,
                            size_t output_size,
                            bool local = false) const override;

    void getBatchTriComplexGHKIsNP(const osh::GO* indices,
                                   size_t input_size,
                                   const model::complex_ghk_current_id curr,
                                   osh::Real* currents,
                                   size_t output_size,
                                   bool local = false) const override;

    void getBatchTriIsNP(const osh::GO* indices,
                         size_t input_size,
                         osh::Real* currents,
                         size_t output_size,
                         bool local = false) const override;

    void getBatchTriOhmicErevsNP(const osh::GO* indices,
                                 size_t input_size,
                                 const model::ohmic_current_id& ohmic_current,
                                 double* rv,
                                 size_t output_size,
                                 bool local) const override;

    void getBatchTriComplexOhmicErevsNP(const osh::GO* indices,
                                        size_t input_size,
                                        const model::complex_ohmic_current_id& ohmic_current,
                                        double* rv,
                                        size_t output_size,
                                        bool local) const override;

    void setTriOhmicErev(osh::GO triangle,
                         const model::ohmic_current_id& ohmic_current,
                         double reversal_potential,
                         bool local) override;
    void setTriComplexOhmicErev(osh::GO triangle,
                                const model::complex_ohmic_current_id& ohmic_current,
                                double reversal_potential,
                                bool local) override;

    bool getTriVClamped(osh::GO vertex, bool local = false) const override;
    void setTriVClamped(osh::GO vertex, bool clamped, bool local = false) override;

    osh::Real getTriIClamp(osh::GO tri, bool local = false) const override;
    void setTriIClamp(osh::GO tri, osh::Real current, bool local = false) override;

    MembraneResistivity getTriRes(osh::GO tri, bool local = false) const override;
    void setTriRes(const osh::GO tri, osh::Real res, osh::Real erev, bool local = false) override;

    osh::Real getTriCapac(osh::GO tri, bool local = false) const override;
    void setTriCapac(const osh::GO tri, osh::Real c, bool local = false) override;


#endif  // USE_PETSC

    //////////////////////
    // Location: Vertex //
    //////////////////////

#if USE_PETSC
    // E-field value
    void getBatchVertVsNP(const osh::GO* indices,
                          size_t input_size,
                          osh::Real* voltages,
                          size_t output_size,
                          bool local = false) const override;
    void setBatchVertVsNP(const osh::GO* indices,
                          size_t input_size,
                          osh::Real* voltages,
                          size_t output_size,
                          bool local = false) override;

    bool getVertVClamped(osh::GO vertex, bool local = false) const override;
    void setVertVClamped(osh::GO vertex, bool clamped, bool local = false) override;

    osh::Real getVertIClamp(osh::GO vertex, bool local = false) const override;
    void setVertIClamp(osh::GO vertex, osh::Real current, bool local = false) override;

#endif  // USE_PETSC

  private:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Internal methods
    ///////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////
    // General methods //
    /////////////////////

    /**
     * Fill vectors with the number of species per elements, owned elements, and
     * boundaries
     */
    static void compute_num_species_per_elements(mesh_type& t_mesh,
                                                 const Statedef& statedef,
                                                 osh::LOs& num_species_per_owned_elems,
                                                 osh::LOs& num_species_per_elems,
                                                 std::optional<osh::LOs>& num_species_per_bounds);
    void init(std::unique_ptr<Statedef>&& t_statedef);

    /** Evolve reactions and diffusions by the time step rd_dt
     *
     * Check "simulation.run" for more info on the general framework
     *
     * Note: here we also update state_time
     *
     * @param rd_dt: time step of reactions and diffusions
     **/
    void evolve_rd(const osh::Real rd_dt);

    /**
     * Update reactions and diffusions up to end_time and decide the time-step
     * (rd_dt).
     *
     * Check "simulation.run" for more info on the general framework
     *
     * @param end_time: reaction-diffusion updates up to end_time
     */
    void run_rd(const osh::Real end_time);

    /**
     * Evolve the system by the efield time step ef_dt.
     *
     * Check "simulation.run" for more info on the general framework
     *
     * @param ef_dt: time step of the efield
     */
    void evolve(const osh::Real ef_dt);

    ///////////////////////////
    // Location: Compartment //
    ///////////////////////////

    // Species
    osh::Real getOwnedCompSpecCount(const model::compartment_id& compartment,
                                    const model::species_name& species) const;
    void setOwnedCompSpecCount(const model::compartment_id& compartment,
                               const model::species_name& species,
                               osh::Real num_molecules,
                               const math::DistributionMethod distribution);

    osh::Real getOwnedCompSpecConc(const model::compartment_id& compartment,
                                   const model::species_name& species) const;

    // Complexes
    osh::Real getOwnedCompComplexCount(
        const model::compartment_id& compartment,
        const model::complex_name& complex,
        const std::vector<util::strongid_vector<model::complex_substate_id,
                                                steps::model::SubunitStateFilter>>& f) const;
    void setOwnedCompComplexCount(const model::compartment_id& compartment,
                                  const model::complex_name& complex,
                                  const util::strongid_vector<model::complex_substate_id, uint>& i,
                                  osh::Real num_molecules,
                                  const math::DistributionMethod distribution);

    osh::Real getOwnedCompComplexSUSCount(
        const model::compartment_id& compartment,
        const model::complex_name& complex,
        const std::vector<
            util::strongid_vector<model::complex_substate_id, steps::model::SubunitStateFilter>>& f,
        model::complex_substate_id m) const;


    /////////////////////
    // Location: Patch //
    /////////////////////

    // Species


    // Complexes
    osh::Real getOwnedPatchComplexCount(
        const model::patch_id& patch,
        const model::complex_name& complex,
        const std::vector<util::strongid_vector<model::complex_substate_id,
                                                steps::model::SubunitStateFilter>>& f) const;
    void setOwnedPatchComplexCount(const model::patch_id& patch,
                                   const model::complex_name& complex,
                                   const util::strongid_vector<model::complex_substate_id, uint>& i,
                                   osh::Real num_molecules,
                                   const math::DistributionMethod distribution);

    osh::Real getOwnedPatchComplexSUSCount(
        const model::patch_id& patch,
        const model::complex_name& complex,
        const std::vector<
            util::strongid_vector<model::complex_substate_id, steps::model::SubunitStateFilter>>& f,
        model::complex_substate_id m) const;

    ///////////////////////////
    // Location: Tetrahedron //
    ///////////////////////////

    // Species
    void setOwnedElementSpecCount(const model::compartment_id& compartment,
                                  const mesh::tetrahedron_id_t element,
                                  const model::species_name& species,
                                  osh::Real num_molecules);

    void getBatchElemValsNP(const osh::GO* indices,
                            size_t input_size,
                            const model::species_name& species,
                            osh::Real* counts,
                            bool useConc,
                            bool local) const;
    void setBatchElemValsNP(const osh::GO* indices,
                            size_t input_size,
                            const model::species_name& species,
                            osh::Real* counts,
                            bool useConc,
                            bool local) const;

    // Reactions
    template <typename Reacs, typename MReacID>
    osh::Real getTetReacK(const Reacs& reacs, osh::GO tet, const MReacID reac, bool local) const;
    template <typename Reacs, typename MReacID>
    void setTetReacK(Reacs& reacs, osh::GO tet, const MReacID reac, osh::Real kcst, bool local);

    ////////////////////////
    // Location: Triangle //
    ////////////////////////

    // Species
    void getBatchBoundSpecCountNP(const osh::GO* indices,
                                  size_t input_size,
                                  const model::species_name& species,
                                  osh::Real* counts,
                                  bool local) const;
    void setBatchBoundSpecCountNP(const osh::GO* indices,
                                  size_t input_size,
                                  const model::species_name& species,
                                  osh::Real* counts,
                                  bool local) const;

    // Reactions
    template <typename Reacs, typename MReacID>
    osh::Real getTriSReacK(const Reacs& reacs,
                           osh::GO triangle,
                           const MReacID& reactionId,
                           bool local) const;
    template <typename Reacs, typename MReacID>
    void setTriSReacK(Reacs& reacs,
                      osh::GO triangle,
                      const MReacID& reactionId,
                      osh::Real kCst,
                      bool local);

#if USE_PETSC
    template <typename SReacT, typename SReacID>
    void getBatchTriSReacIsNP(const SReacT& reacs,
                              const osh::GO* indices,
                              size_t input_size,
                              const SReacID reac,
                              osh::Real* currents,
                              bool local) const;

    template <typename CurrdefT, typename MCurrID>
    void getBatchTriOhmicIsNP(const osh::GO* indices,
                              size_t input_size,
                              const MCurrID curr,
                              osh::Real* currents,
                              bool local) const;

    template <typename SReacT, typename MReacID>
    void getBatchTriGHKIsNP(const SReacT& reacs,
                            const osh::GO* indices,
                            size_t input_size,
                            const MReacID curr,
                            osh::Real* currents,
                            bool local) const;

    template <typename CurrT, typename MCurrID>
    void getBatchTriOhmicErevsNP(const gsl::span<const osh::GO>& triangles,
                                 const MCurrID& ohmic_current,
                                 const gsl::span<double>& erev,
                                 bool local) const;

    template <typename CurrdefT, typename MCurrID>
    void setTriOhmicErev(osh::GO triangle,
                         const MCurrID& ohmic_current,
                         double reversal_potential,
                         bool local);
#endif  // USE_PETSC

    /////////////////////
    // Utility methods //
    /////////////////////

    template <typename T>
    T allReduce(const T& v,
                MPI_Op op = MPI_SUM,
                MPI_Datatype mpi_type = util::mpi_get_type<T>()) const {
        T global_v;
        auto err = MPI_Allreduce(&v, &global_v, 1, mpi_type, op, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
        return global_v;
    }

    static util::strongid_vector<model::complex_substate_id, uint> _convertComplexState(
        const std::vector<std::vector<steps::model::SubunitStateFilter>>& f);

    static std::vector<
        util::strongid_vector<model::complex_substate_id, steps::model::SubunitStateFilter>>
    _convertComplexFilters(const std::vector<std::vector<steps::model::SubunitStateFilter>>& f);

    template <typename LT>
    inline LT getLocalInd(osh::GO idx, bool local, bool owned = true) const {
        if (local) {
            if (idx == static_cast<osh::LO>(idx)) {
                return LT(static_cast<osh::LO>(idx));
            } else {
                return {};
            }
        } else {
            using GT = decltype(mesh.getGlobalIndex(LT(static_cast<osh::LO>(idx))));
            return mesh.getLocalIndex(GT(idx), owned);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Attributes
    ///////////////////////////////////////////////////////////////////////////////////////////////

    mesh_type& mesh;
    bool indepKProcs;

    std::unique_ptr<Statedef> statedef;
    std::unique_ptr<SimulationInput> input;
    std::unique_ptr<SimulationData<SSA, SearchMethod, DiffMethod>> data;

    osh::I64 num_iterations{};
    osh::Real state_time{};

    bool outdated_diffusions{false};
};

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<Simulation> GetSimulation(steps::model::Model& model,
                                          DistMesh& t_mesh,
                                          const rng::RNGptr& r,
                                          int ssaMethod,
                                          int searchMethod,
                                          int diffMethod,
                                          bool t_indepKProcs,
                                          bool isEfield);

template <SSAMethod SSA>
std::unique_ptr<Simulation> _getSimulation(steps::model::Model& model,
                                           DistMesh& t_mesh,
                                           const rng::RNGptr& r,
                                           NextEventSearchMethod searchMethod,
                                           DiffusionMethod diffMethod,
                                           bool t_indepKProcs,
                                           bool isEfield);


template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
std::unique_ptr<Simulation> _getSimulation(steps::model::Model& model,
                                           DistMesh& mesh,
                                           const rng::RNGptr& r,
                                           DiffusionMethod diffMethod,
                                           bool indepKProcs,
                                           bool isEfield);

}  // namespace steps::dist
