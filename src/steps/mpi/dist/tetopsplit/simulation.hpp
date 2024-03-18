#pragma once

#include <Omega_h_adj.hpp>

#include "geom/dist/distmesh.hpp"
#include "kproc/diffusions.hpp"
#include "math/distributions.hpp"
#include "mol_state.hpp"
#include "operator/diffusion_operator.hpp"
#include "operator/rssa_operator.hpp"
#include "operator/ssa_operator.hpp"
#include "rng/rng.hpp"
#include "simulation_data.hpp"

namespace steps::dist {

struct CompartmentCount {
    CompartmentCount(model::species_name t_species, osh::Real t_num_mols)
        : species(std::move(t_species))
        , num_mols(t_num_mols) {}
    model::species_name species;
    osh::Real num_mols;
};

struct CompartmentConc {
    CompartmentConc(model::species_name t_species, osh::Real t_concentration)
        : species(std::move(t_species))
        , concentration(t_concentration) {}
    model::species_name species;
    osh::Real concentration;
};

struct PatchCount {
    PatchCount(model::patch_id t_patch, model::species_name t_species, osh::Real t_num_mols)
        : patch(std::move(t_patch))
        , species(std::move(t_species))
        , num_mols(t_num_mols) {}
    const model::patch_id patch;
    const model::species_name species;
    const osh::Real num_mols;
};

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
    using compartment_counts_t =
        std::unordered_map<model::compartment_id, std::vector<CompartmentCount>>;
    using compartment_concs_t =
        std::unordered_map<model::compartment_id, std::vector<CompartmentConc>>;
    using patch_counts_t = std::vector<PatchCount>;

    Simulation(DistMesh& t_mesh, rng::RNG& t_rng, std::ostream& t_outstream);
    virtual ~Simulation() noexcept;

    virtual void setCompSpecCount(
        const model::compartment_id& compartment,
        const model::species_name& species,
        osh::Real num_molecules,
        const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM) = 0;

    void setCompSpecCount(
        const compartment_counts_t& counts,
        const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM);
    virtual void setCompSpecCount(
        const model::compartment_id& compartment,
        const std::vector<CompartmentCount>& counts,
        const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM) = 0;
    void setCompSpecConc(
        const compartment_concs_t& concentrations,
        const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM);
    virtual void setCompSpecConc(
        const model::compartment_id& compartment,
        const std::vector<CompartmentConc>& concs,
        const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM) = 0;

    virtual void setPatchSpecCount(
        const model::patch_id& /*patch_id*/,
        const model::species_name& /*species*/,
        osh::Real /*num_molecules*/,
        const math::DistributionMethod = math::DistributionMethod::DIST_UNIFORM) {
        throw std::logic_error("NOT IMPLEMENTED");
    }
    virtual void setPatchSpecCount(
        const patch_counts_t& /*counts*/,
        const math::DistributionMethod = math::DistributionMethod::DIST_UNIFORM) {
        throw std::logic_error("NOT IMPLEMENTED");
    }
    virtual const std::vector<mesh::triangle_id_t>& getGHKBoundaries() const {
        throw std::logic_error("NOT IMPLEMENTED");
    }
    virtual osh::Reals getGHKCurrents() const {
        throw std::logic_error("NOT IMPLEMENTED");
    }
    virtual osh::Real getTotalGHKCurrent() const {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    virtual void setKCstSReac(const model::patch_id& /*patchId*/,
                              container::surface_reaction_id /*reactionId*/,
                              osh::Real /*kCst*/) {
        throw std::logic_error("NOT IMPLEMENTED");
    }
    virtual void setTriSpecCount(const model::patch_id& /*patch*/,
                                 const model::species_name& /*species*/,
                                 osh::Real /*num_molecules*/) {
        throw std::logic_error("NOT IMPLEMENTED");
    }
    virtual osh::Real getPatchSpecCount(const model::patch_id& /*compartment*/,
                                        const model::species_name& /*species*/) const {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    virtual void setPotential(osh::Real) {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    virtual osh::Reals getPotentialOnVertices(const model::patch_id&) {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    virtual osh::Real getMaxPotentialOnVertices(const model::patch_id&) {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    virtual osh::Real getMinPotentialOnVertices(const model::patch_id&) {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    virtual osh::Reals getPotentialOnTriangles(const model::patch_id&) {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    virtual osh::Real getMaxPotentialOnTriangles(const model::patch_id&) {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    virtual osh::Real getMinPotentialOnTriangles(const model::patch_id&) {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    virtual osh::Real getLocalMaxPotentialOnVertices(const model::patch_id&) {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    virtual osh::Real getLocalMinPotentialOnVertices(const model::patch_id&) {
        throw std::logic_error("NOT_IMPLEMENTED");
    }

#if USE_PETSC
    virtual std::optional<osh::Real> getEfieldDt() const {
        throw std::logic_error("NOT_IMPLEMENTED");
    }

    virtual void setEfieldDt(const osh::Real) const {
        throw std::logic_error("NOT_IMPLEMENTED");
    }
#endif  // USE_PETSC

    virtual osh::Real getCompSpecConc(const model::compartment_id& compartment,
                                      const model::species_name& species) const = 0;
    virtual osh::Real getOwnedCompSpecConc(const model::compartment_id& compartment,
                                           const model::species_name& species) const = 0;
    virtual osh::Real getCompSpecCount(const model::compartment_id& compartment,
                                       const model::species_name& species) const = 0;

#if USE_PETSC
    virtual std::pair<mesh::triangle_ids, osh::Reals> getOhmicCurrents(
        const model::membrane_id& membrane_id,
        const model::channel_id& l_id) const = 0;

    virtual osh::Real getTotalOhmicCurrent(const model::membrane_id& mem_id,
                                           const model::channel_id& chan_id) const = 0;
#endif  // USE_PETSC

    virtual osh::Real getOwnedCompSpecCount(const model::compartment_id& compartment,
                                            const model::species_name& species) const = 0;
    virtual std::pair<std::reference_wrapper<const mesh::tetrahedron_ids>, std::vector<osh::LO>>
    getOwnedElemSpecCount(const model::species_name& species) const = 0;

    virtual std::pair<std::vector<mesh::tetrahedron_global_id_t>, std::vector<osh::LO>>
    getElemSpecCount(const model::species_name& species) const = 0;

    virtual void setOwnedElementSpecCount(const model::compartment_id& compartment,
                                          mesh::tetrahedron_id_t element,
                                          const model::species_name& species,
                                          osh::Real num_mols) = 0;

    virtual void setDiffOpBinomialThreshold(osh::Real threshold) = 0;

    virtual void exportMolStateToVTK(const std::string& filename) = 0;

    virtual osh::I64 getDiffOpExtent(bool local = false) const = 0;
    virtual osh::I64 getSSAOpExtent(bool local = false) const = 0;
    virtual osh::I64 getNIterations() const noexcept = 0;
    virtual osh::Real getIterationTimeStep() const noexcept = 0;

    virtual std::string createStateReport() = 0;

    virtual void init(std::unique_ptr<Statedef>&& statedef) = 0;
    virtual void reset() = 0;
    virtual void run(osh::Real end_time) = 0;
    /// \brief Emit the given message on every rank.
    /// the message is prefixed by "[RANK] "
    void log_all(const std::string& message) const;
    /// Emit the given message on rank 0 only
    void log_once(const std::string& message, bool force_stdout = false) const;
    /**
     * log progress to stdout
     *
     * progress = 100*i/tot %
     *
     * @param i progress index
     * @param tot progress total
     * @param name what is progressing
     */
    void log_progress(const double i, const double tot, const std::string& name = "total") const;

    virtual void setDiffusionBoundaryActive(const mesh::diffusion_boundary_name& boundary_id,
                                            const model::species_name& species,
                                            bool set_active) = 0;
    const std::vector<unsigned int>& get_diffusion_rank_exchanges() const noexcept {
        return diffusion_rank_exchanges;
    }

    virtual SSAMethod ssaMethod() const noexcept = 0;

    void log_diffusion_exchanges() const;

    const DistMesh& getMesh() const noexcept {
        return mesh;
    }

    DistMesh& getMesh() noexcept {
        return mesh;
    }

    double getElapsedSSA() const noexcept {
        return reactions_timer;
    }

    double getElapsedDiff() const noexcept {
        return diffusions_timer;
    }

    double getElapsedEField() const noexcept {
        return efield_timer;
    }

    const int comm_rank;
    const int comm_size;

    inline MPI_Comm comm() const noexcept {
        return getMesh().comm_impl();
    }

  protected:
    /**
     * \name for diffusion debug only
     * \{
     */
    /// matrix N, N with N the number of ranks
    /// M[i, j] provides number of molecules that went from rank i to rank j
    std::vector<unsigned int> diffusion_rank_exchanges;
    /** \} */
    DistMesh& mesh;
    rng::RNG& rng;
    /// total time in seconds spent performing reactions
    double reactions_timer{};
    /// total time in seconds spent performing diffusions
    double diffusions_timer{};
    /// total time in seconds spent performing efield
    double efield_timer{};

    std::ostream& outstream;
};

template <SSAMethod SSA = SSAMethod::SSA,
          NextEventSearchMethod SearchMethod = NextEventSearchMethod::Direct>
class OmegaHSimulation: public Simulation {
  public:
    using super_type = Simulation;
    using mesh_type = DistMesh;

    OmegaHSimulation(mesh_type& t_mesh,
                     rng::RNG& t_rng,
                     std::ostream& t_outstream,
                     bool t_indepKProcs);

    void setCompSpecCount(const model::compartment_id& compartment,
                          const model::species_name& species,
                          osh::Real num_molecules,
                          const math::DistributionMethod distribution =
                              math::DistributionMethod::DIST_UNIFORM) override;
    void setCompSpecCount(const model::compartment_id& compartment,
                          const std::vector<CompartmentCount>& counts,
                          const math::DistributionMethod distribution =
                              math::DistributionMethod::DIST_UNIFORM) override;
    void setCompSpecConc(
        const model::compartment_id& compartment,
        const model::species_name& species,
        osh::Real concentration,
        const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM);
    void setCompSpecConc(const model::compartment_id& compartment,
                         const std::vector<CompartmentConc>& concentrations,
                         const math::DistributionMethod distribution =
                             math::DistributionMethod::DIST_UNIFORM) override;
    void setOwnedElementSpecCount(const model::compartment_id& compartment,
                                  const mesh::tetrahedron_id_t element,
                                  const model::species_name& species,
                                  osh::Real num_molecules) override;

    osh::Real getCompSpecConc(const model::compartment_id& compartment,
                              const model::species_name& species) const override;
    osh::Real getOwnedCompSpecConc(const model::compartment_id& compartment,
                                   const model::species_name& species) const override;
    osh::Real getCompSpecCount(const model::compartment_id& compartment,
                               const model::species_name& species) const override;

#if USE_PETSC
    std::pair<mesh::triangle_ids, osh::Reals> getOhmicCurrents(
        const model::membrane_id& mem_id,
        const model::channel_id& chan_id) const override;
    osh::Real getTotalOhmicCurrent(const model::membrane_id& mem_id,
                                   const model::channel_id& chan_id) const override;
#endif  // USE_PETSC
    /**
     * \brief Distribute a certain number of molecules on a patch.
     *
     * Check setCompCount for more details on distribution methods
     *
     * \param patch patch identifier
     * \param species species name
     * \param num_molecules number of molecules to distribute among the patch
     * \param is_uniform_distribution true: uniform distribution, false:
     * multinomial
     */
    void setPatchSpecCount(const model::patch_id& patch,
                           const model::species_name& species,
                           osh::Real num_molecules,
                           const math::DistributionMethod distribution =
                               math::DistributionMethod::DIST_UNIFORM) override;
    /**
     * \brief Distribute molecules on multiple patches.
     *
     * \param counts list of molecules assignments on patches
     */
    void setPatchSpecCount(const patch_counts_t& counts,
                           const math::DistributionMethod distribution =
                               math::DistributionMethod::DIST_UNIFORM) override;

    /**
     * \brief Get the total number of molecules on a patch.
     *
     * \param patch patch identifier
     * \param species species name
     */
    osh::Real getPatchSpecCount(const model::patch_id& patch,
                                const model::species_name& species) const override;

    osh::Real getOwnedCompSpecCount(const model::compartment_id& compartment,
                                    const model::species_name& species) const override;

    std::pair<std::reference_wrapper<const mesh::tetrahedron_ids>, std::vector<osh::LO>>
    getOwnedElemSpecCount(const model::species_name& species) const override;

    std::pair<std::vector<mesh::tetrahedron_global_id_t>, std::vector<osh::LO>> getElemSpecCount(
        const model::species_name& species) const override;

    void getBatchElemValsNP(const osh::GO* indices,
                            size_t input_size,
                            const model::species_name& species,
                            osh::Real* counts,
                            bool useConc = false,
                            bool local = false) const;

    void setBatchElemValsNP(const osh::GO* indices,
                            size_t input_size,
                            const model::species_name& species,
                            osh::Real* counts,
                            bool useConc = false,
                            bool local = false) const;

    void getBatchBoundSpecCountNP(const osh::GO* indices,
                                  size_t input_size,
                                  const model::species_name& species,
                                  osh::Real* counts,
                                  bool local = false) const;

    void setBatchBoundSpecCountNP(const osh::GO* indices,
                                  size_t input_size,
                                  const model::species_name& species,
                                  osh::Real* counts,
                                  bool local = false) const;

    void getBatchVertVsNP(const osh::GO* indices,
                          size_t input_size,
                          osh::Real* voltages,
                          bool local = false) const;

    void getBatchTriVsNP(const osh::GO* indices,
                         size_t input_size,
                         osh::Real* voltages,
                         bool local = false) const;

    void getBatchTetVsNP(const osh::GO* indices,
                         size_t input_size,
                         osh::Real* voltages,
                         bool local = false) const;

    void getBatchTriOhmicIsNP(const osh::GO* indices,
                              size_t input_size,
                              const model::ohmic_current_id curr,
                              osh::Real* currents,
                              bool local = false) const;

    void getBatchTriGHKIsNP(const osh::GO* indices,
                            size_t input_size,
                            const model::ghk_current_id curr,
                            osh::Real* currents,
                            bool local = false) const;

    void setDiffOpBinomialThreshold(osh::Real threshold) override;
    osh::Real getIterationTimeStep() const noexcept override;

    void exportMolStateToVTK(const std::string& filename) override;

    osh::I64 getDiffOpExtent(bool local = false) const override;
    osh::I64 getSSAOpExtent(bool local = false) const override;
    osh::I64 getNIterations() const noexcept override;

    std::string createStateReport() override;

    void init(std::unique_ptr<Statedef>&& t_statedef) override;
    void reset() override;
    /** Run the simulation up to end_time from state_time.
     *
     * Every time run is called its results are checkpointed at state_time. This
     * function decides the time step for evolve (which moves the simulation one
     * time step further and updates state_time). If end_time-state_time is not
     * a multiple of the time step we reduce the last time step to align it with
     * end_time.
     *
     * end_time-state_time is usually the sampling frequency of  the simulation.
     * In other words, for every recording point at time rt = i*dt we call
     * simulation.run(rt) and we record results. In this sense it is the
     * sampling frequency. Internally the simulation proceeds with its own time
     * steps (efield_dt and reaction-diffusion dt). As explaine,d time steps can
     * be reduced to align the simulation to the recordings.
     *
     * @param end_time : run up to end_time from state_time.
     */
    void run(osh::Real end_time) override;

    inline osh::Real getTime() const noexcept {
        return state_time;
    }

    SSAMethod ssaMethod() const noexcept override {
        return SSA;
    }

    /// Get temperature
    inline osh::Real getTemp() const noexcept {
        return statedef->getTemp();
    }
    /// Set temperature
    inline void setTemp(const osh::Real temp) noexcept {
        statedef->setTemp(temp);
    }

    /**
     * Fill vectors with the number of species per elements, owned elements, and
     * boundaries
     */
    static void compute_num_species_per_elements(mesh_type& t_mesh,
                                                 const Statedef& statedef,
                                                 osh::LOs& num_species_per_owned_elems,
                                                 osh::LOs& num_species_per_elems,
                                                 std::optional<osh::LOs>& num_species_per_bounds);

    virtual void setPotential(osh::Real value) override {
        const auto assign = OMEGA_H_LAMBDA(osh::LO v) {
            input->potential_on_vertices_w[v] = value;
        };
        osh::parallel_for(input->potential_on_vertices_w.size(),
                          assign,
                          "OmegaHSimulation::setPotential");
    }

    /**
     * \brief Set the potential of a membrane and its associated conductor
     * volume
     *
     * \param memb membrane identifier
     * \param value potential in Volts
     */
    void setMembPotential(const model::membrane_id& memb, osh::Real value);

    /**
     * \brief Set the leakage parameters of a membrane:
     * the resistivity and the reversal potential
     *
     * \param membrane membrane identifier
     * \param resistivity resistivity
     * \param reversal_potential reversal potential
     */
    void setMembRes(const model::membrane_id& membrane,
                    osh::Real resistivity,
                    osh::Real reversal_potential);

    /**
     * \brief Get leakage parameters of a membrane:
     * the resistivity and the reversal potential
     *
     * \param membrane membrane identifier
     */
    MembraneResistivity getMembRes(const model::membrane_id& membrane);

    /** Get resistivity and reversal potential of leak on a triangle
     *
     * \param tri triangle identifier
     * \return resistivity and reversal potential
     */
    MembraneResistivity getTriRes(osh::GO tri, bool local = false) const;

    /** Set resistivity and reversal potential of leak on a triangle
     *
     * \param tri triangle identifier
     * \param res resistivity
     * \param erev reversal potential
     */
    void setTriRes(const osh::GO tri, osh::Real res, osh::Real erev, bool local = false);

    /** Get capacitance on a triangle
     *
     * \param tri triangle identifier
     * \return capacitance
     */
    osh::Real getTriCapac(osh::GO tri, bool local = false) const;

    /** Set capacitance on a triangle
     *
     * \param tri triangle identifier
     * \param c capacitance
     */
    void setTriCapac(const osh::GO tri, osh::Real c, bool local = false);

    /** Get current on a vertex
     *
     * \param vertex: global vertex identifier
     * \return current (Amps)
     */
    osh::Real getVertIClamp(osh::GO vertex, bool local = false) const;

    /** Set current on vertex
     *
     * /param vertex: global vertex index
     * /param current: current (Amps)
     */
    void setVertIClamp(osh::GO vertex, osh::Real current, bool local = false);

    double getTriOhmicErev(osh::GO triangle,
                           const model::ohmic_current_id& ohmic_current,
                           bool local) const;
    void getBatchTriOhmicErevsNP(const gsl::span<const osh::GO>& triangles,
                                 const model::ohmic_current_id& ohmic_current,
                                 const gsl::span<double>& erev,
                                 bool local) const;
    void setTriOhmicErev(osh::GO triangle,
                         const model::ohmic_current_id& ohmic_current,
                         double reversal_potential,
                         bool local);

    /** \brief Return GHK boundaries */
    const std::vector<mesh::triangle_id_t>& getGHKBoundaries() const override;
    /** \brief Return GHK current */
    virtual osh::Reals getGHKCurrents() const override;
    /** \brief Return GHK current for all the processes */
    virtual osh::Real getTotalGHKCurrent() const override;
    /** \brief Set a reaction constant in simulation */
    void setKCstSReac(const model::patch_id& patchId,
                      container::surface_reaction_id reactionId,
                      osh::Real kCst) override {
        data->kproc_state.surfaceReactions().updateKCst(patchId, reactionId, kCst, mesh);
    }

    void setPatchSReacK(const model::patch_id& patchId,
                        const model::surface_reaction_id& reactionId,
                        osh::Real kCst);

    osh::Reals getPotentialOnVertices(const model::patch_id& patch) override;

    osh::Real getMaxPotentialOnVertices(const model::patch_id& patch) override {
        return mesh.get_MPI_max(getPotentialOnVertices(patch));
    }

    osh::Real getMinPotentialOnVertices(const model::patch_id& patch) override {
        return mesh.get_MPI_min(getPotentialOnVertices(patch));
    }

    osh::Reals getPotentialOnTriangles(const model::patch_id& patch) override;

    osh::Real getMaxPotentialOnTriangles(const model::patch_id& patch) override {
        return mesh.get_MPI_max(getPotentialOnTriangles(patch));
    }

    osh::Real getMinPotentialOnTriangles(const model::patch_id& patch) override {
        return mesh.get_MPI_min(getPotentialOnTriangles(patch));
    }

    osh::Real getLocalMaxPotentialOnVertices(const model::patch_id& patch) override {
        auto vals = getPotentialOnVertices(patch);
        if (!vals.size()) {
            return std::numeric_limits<osh::Real>::quiet_NaN();
        } else {
            return *std::max_element(vals.begin(), vals.end());
        }
    }

    osh::Real getLocalMinPotentialOnVertices(const model::patch_id& patch) override {
        auto vals = getPotentialOnVertices(patch);
        if (!vals.size()) {
            return std::numeric_limits<osh::Real>::quiet_NaN();
        } else {
            return *std::min_element(vals.begin(), vals.end());
        }
    }

#if USE_PETSC
    std::optional<osh::Real> getEfieldDt() const override {
        if (data->efield)
            return data->efield->getDt();
        return std::nullopt;
    }

    void setEfieldDt(const osh::Real dt) const override {
        if (data->efield)
            data->efield->setDt(dt);
        else
            throw std::logic_error("NO E-FIELD WHERE I CAN ASSIGN DT");
    }

    void setPetscOptions() {
        if (data->efield) {
            data->efield->setPetscOptions();
        } else {
            throw std::logic_error("E-Field is not in use.");
        }
    }
#endif  // USE_PETSC

    void setDiffusionBoundaryActive(const mesh::diffusion_boundary_name& diffusion_boundary_name,
                                    const model::species_name& spec_id,
                                    bool set_active) override;

    bool getDiffusionBoundaryActive(const mesh::diffusion_boundary_name& diffusion_boundary_name,
                                    const model::species_name& spec_id);

    /**
     * Set a current clamp on a membrane
     * \param membrane membrane identifier
     * \param current stimulus to apply (Amps)
     */
    void setMembIClamp(const model::membrane_id& membrane, osh::Real current);

    /// Dump the dependency graph of kproc in a file specified by path
    void dumpDepGraphToFile(const std::string& path) const;

  private:
    /**
     * Update the number of molecules of a certain species in the current rank
     *
     * In case is_uniform is false, the distribution of molecules among the
     * elements is stochastic following the probability: V_elem/V_rank where
     * V_elem is the volume of an element and V_rank is the volume of all the
     * elements owned by the rank
     *
     * In case is_uniform is true molecules are evenly split following the
     * partial volume ratio: V_elem/V_rank with stochastic rounding
     *
     * \param molecules Molecules pool instance
     * \param comp_id compartment identifier
     * \param species species identifier
     * \param num_molecules number of molecules to set
     * \param distribution distribution type. Check distributions.hpp for more
     * information
     */
    void setOwnedCompSpecCount(
        const model::compartment_id& compartment,
        const model::species_name& species,
        osh::Real num_molecules,
        const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM);

    /**
     * Evolve the system by the efield time step ef_dt.
     *
     * Check "simulation.run" for more info on the general framework
     *
     * @param ef_dt: time step of the efield
     */
    void evolve(const osh::Real ef_dt);

    /**
     * Update reactions and diffusions up to end_time and decide the time-step
     * (rd_dt).
     *
     * Check "simulation.run" for more info on the general framework
     *
     * @param end_time: reaction-diffusion updates up to end_time
     */
    void run_rd(const osh::Real end_time);

    /** Evolve reactions and diffusions by the time step rd_dt
     *
     * Check "simulation.run" for more info on the general framework
     *
     * Note: here we also update state_time
     *
     * @param rd_dt: time step of reactions and diffusions
     **/
    void evolve_rd(const osh::Real rd_dt);

    void initialize_discretized_rates();

    mesh_type& mesh;
    const osh::LOs elems2verts;
    const osh::Reals coords;
    osh::I64 num_iterations{};
    osh::Real state_time{};
    std::unique_ptr<Statedef> statedef;
    std::unique_ptr<SimulationInput> input;
    std::unique_ptr<SimulationData<SSA, SearchMethod>> data;
    bool indepKProcs;

    /// provide MPI rank that owned a given element, for diffusion debugging
    /// only
    osh::LOs element_ranks;
};

// explicit template instantiation declarations

extern template class OmegaHSimulation<SSAMethod::SSA, NextEventSearchMethod::Direct>;
extern template class OmegaHSimulation<SSAMethod::RSSA, NextEventSearchMethod::Direct>;
extern template class OmegaHSimulation<SSAMethod::SSA, NextEventSearchMethod::GibsonBruck>;

}  // namespace steps::dist
