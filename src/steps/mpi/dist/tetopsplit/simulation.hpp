#pragma once

#include <Omega_h_adj.hpp>

#include "geom/dist/distmesh.hpp"
#include "kproc/diffusions.hpp"
#include "math/distributions.hpp"
#include "mol_state.hpp"
#include "mpi/dist/test/simulation.hpp"
#include "operator/diffusion_operator.hpp"
#include "operator/rssa_operator.hpp"
#include "operator/ssa_operator.hpp"
#include "rng/rng.hpp"
#include "simulation_data.hpp"

namespace steps {
namespace dist {

template <SSAMethod SSA = SSAMethod::SSA, typename RNG = std::mt19937,
          typename NumMolecules = osh::LO,
          NextEventSearchMethod SearchMethod = NextEventSearchMethod::Direct>
class OmegaHSimulation : public Simulation<RNG> {
public:
  using super_type = Simulation<RNG>;
  using mesh_type = DistMesh;
  using molecules_type = NumMolecules;

  OmegaHSimulation(const ScenarioInput &t_scenario, mesh_type &t_mesh,
                   RNG &t_rng, std::ostream &t_outstream);

  void setCompCount(const model::compartment_id& compartment,
                    const model::species_name& species,
                    osh::Real num_molecules,
                    const math::DistributionMethod distribution =
                        math::DistributionMethod::DIST_UNIFORM) override;
  void setCompCount(const model::compartment_id& compartment,
                    const std::vector<CompartmentCount>& counts,
                    const math::DistributionMethod distribution =
                        math::DistributionMethod::DIST_UNIFORM) override;
  void setCompConc(
      const model::compartment_id& compartment,
      const model::species_name& species,
      osh::Real concentration,
      const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM);
  void setCompConc(const model::compartment_id& compartment,
                   const std::vector<CompartmentConc>& concentrations,
                   const math::DistributionMethod distribution =
                       math::DistributionMethod::DIST_UNIFORM) override;
  void setOwnedElementCount(const model::compartment_id &compartment,
                            const mesh::tetrahedron_id_t element,
                            const model::species_name &species,
                            osh::Real num_molecules) override;

  osh::Real getCompConc(const model::compartment_id &compartment,
                        const model::species_name &species) const override;
  osh::Real getOwnedCompConc(const model::compartment_id &compartment,
                             const model::species_name &species) const override;
  osh::Real getCompCount(const model::compartment_id &compartment,
                         const model::species_name &species) const override;

#if USE_PETSC
  std::pair<mesh::triangle_ids, osh::Reals>
  getOhmicCurrents(const model::membrane_id &mem_id,
                   const model::channel_id &chan_id) const override;
  osh::Real
  getTotalOhmicCurrent(const model::membrane_id &mem_id,
                       const model::channel_id &chan_id) const override;
#endif // USE_PETSC
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
  void setPatchCount(const model::patch_id& patch,
                     const model::species_name& species,
                     osh::Real num_molecules,
                     const math::DistributionMethod distribution =
                         math::DistributionMethod::DIST_UNIFORM) override;
  /**
   * \brief Distribute molecules on multiple patches.
   *
   * \param counts list of molecules assignments on patches
   */
  void setPatchCount(const Simdef::patch_counts_t& counts,
                     const math::DistributionMethod distribution =
                         math::DistributionMethod::DIST_UNIFORM) override;

  /**
   * \brief Get the total number of molecules on a patch.
   *
   * \param patch patch identifier
   * \param species species name
   */
  osh::Real getPatchCount(const model::patch_id &patch,
                          const model::species_name &species) const override;

  osh::Real
  getOwnedCompCount(const model::compartment_id &compartment,
                    const model::species_name &species) const override;

  std::pair<std::reference_wrapper<const mesh::tetrahedron_ids>,
            std::vector<osh::LO>>
  getOwnedElemCount(const model::species_name &species) const override;

  std::pair<std::vector<mesh::tetrahedron_global_id_t>, std::vector<osh::LO>>
  getElemCount(const model::species_name &species) const override;

  void getBatchElemValsNP(const osh::GO* indices,
                          size_t input_size,
                          const model::species_name& species,
                          osh::Real* counts,
                          bool useConc = false) const;

  void setBatchElemValsNP(const osh::GO* indices,
                          size_t input_size,
                          const model::species_name& species,
                          osh::Real* counts,
                          bool useConc = false) const;

  void getBatchBoundCountNP(const osh::GO* indices,
                            size_t input_size,
                            const model::species_name& species,
                            osh::Real* counts) const;

  void setBatchBoundCountNP(const osh::GO* indices,
                            size_t input_size,
                            const model::species_name& species,
                            osh::Real* counts) const;

  void getBatchVertVsNP(const osh::GO* indices, size_t input_size, osh::Real* voltages) const;

  void getBatchTriVsNP(const osh::GO* indices, size_t input_size, osh::Real* voltages) const;

  void getBatchTetVsNP(const osh::GO* indices, size_t input_size, osh::Real* voltages) const;

  void getBatchTriOhmicIsNP(const osh::GO* indices,
                            size_t input_size,
                            const model::ohmic_current_id curr,
                            osh::Real* currents) const;

  void getBatchTriGHKIsNP(const osh::GO* indices,
                          size_t input_size,
                          const model::ghk_current_id curr,
                          osh::Real* currents) const;

  void setDiffOpBinomialThreshold(osh::Real threshold) override;
  osh::Real getIterationTimeStep() const noexcept override;

  void exportMolStateToVTK(const std::string &filename) override;

  osh::I64 getDiffOpExtent(bool local = false) const override;
  osh::I64 getSSAOpExtent(bool local = false) const override;
  osh::I64 getNIterations() const noexcept override;

  std::string createStateReport() override;

  void init(std::unique_ptr<Statedef> &&t_statedef) override;
  void reset() override;
  /** Run the simulation up to end_time from state_time.
   *
   * Every time run is called its results are checkpointed at state_time. This
   * function decides the time step for evolve (which moves the simulation one
   * time step further and updates state_time). If end_time-state_time is not a
   * multiple of the time step we reduce the last time step to align it with
   * end_time.
   *
   * end_time-state_time is usually the sampling frequency of  the simulation.
   * In other words, for every recording point at time rt = i*dt we call
   * simulation.run(rt) and we record results. In this sense it is the sampling
   * frequency. Internally the simulation proceeds with its own time steps
   * (efield_dt and reaction-diffusion dt). As explaine,d time steps can be
   * reduced to align the simulation to the recordings.
   *
   * @param end_time : run up to end_time from state_time.
   */
  void run(osh::Real end_time) override;

  inline osh::Real getTime() const noexcept { return state_time; }

  SSAMethod ssaMethod() const noexcept override { return SSA; }

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
  static void compute_num_species_per_elements(
      mesh_type &t_mesh, const Statedef &statedef,
      osh::LOs &num_species_per_owned_elems, osh::LOs &num_species_per_elems,
      boost::optional<osh::LOs> &num_species_per_bounds);
  mesh_type &mesh_impl() noexcept { return mesh; }

  virtual void setPotential(osh::Real value) override {
    const auto assign = OMEGA_H_LAMBDA(osh::LO v) {
      input->potential_on_vertices_w[v] = value;
    };
    osh::parallel_for(input->potential_on_vertices_w.size(), assign,
                      "OmegaHSimulation::setPotential");
  }

  /**
   * \brief Set the potential of a membrane and its associated conductor volume
   *
   * \param memb membrane identifier
   * \param value potential in Volts
   */
  void setMembPotential(const model::membrane_id &memb, osh::Real value);

  /** Get current on a vertex
   *
   * \param vertex: global vertex identifier
   * \return current (Amps)
   */
  osh::Real getVertIClamp(mesh::vertex_global_id_t vertex) const;

  /** Set current on vertex
   *
   * /param vertex: global vertex index
   * /param current: current (Amps)
   */
  void setVertIClamp(mesh::vertex_global_id_t vertex, osh::Real current);

  /** \brief Return GHK boundaries */
  const std::vector<mesh::triangle_id_t> &getGHKBoundaries() const override;
  /** \brief Return GHK current */
  virtual osh::Reals getGHKCurrents() const override;
  /** \brief Return GHK current for all the processes */
  virtual osh::Real getTotalGHKCurrent() const override;
  /** \brief Set a reaction constant in simulation */
  void setKCstSReac(const model::patch_id &patchId,
                    container::surface_reaction_id reactionId,
                    osh::Real kCst) override {
    data->kproc_state.surfaceReactions().updateKCst(patchId, reactionId, kCst,
                                                    mesh);
  }

  void setPatchSReacK(const model::patch_id &patchId,
                      const model::surface_reaction_id &reactionId,
                      osh::Real kCst);

  osh::Reals getPotentialOnVertices(const model::patch_id &patch) override;

  osh::Real getMaxPotentialOnVertices(const model::patch_id &patch) override {
    return mesh.get_MPI_max(getPotentialOnVertices(patch));
  }

  osh::Real getMinPotentialOnVertices(const model::patch_id &patch) override {
    return mesh.get_MPI_min(getPotentialOnVertices(patch));
  }

  osh::Reals getPotentialOnTriangles(const model::patch_id &patch) override;

  osh::Real getMaxPotentialOnTriangles(const model::patch_id &patch) override {
    return mesh.get_MPI_max(getPotentialOnTriangles(patch));
  }

  osh::Real getMinPotentialOnTriangles(const model::patch_id &patch) override {
    return mesh.get_MPI_min(getPotentialOnTriangles(patch));
  }

  osh::Real
  getLocalMaxPotentialOnVertices(const model::patch_id &patch) override {
    auto vals = getPotentialOnVertices(patch);
    if (!vals.size()) {
      return std::numeric_limits<osh::Real>::quiet_NaN();
    }
    else {
      return *std::max_element(vals.begin(), vals.end());
    }
  }

  osh::Real
  getLocalMinPotentialOnVertices(const model::patch_id &patch) override {
    auto vals = getPotentialOnVertices(patch);
    if (!vals.size()) {
      return std::numeric_limits<osh::Real>::quiet_NaN();
    }
    else {
      return *std::min_element(vals.begin(), vals.end());
    }
  }

#if USE_PETSC
  boost::optional<osh::Real> getEfieldDt() const override {
    if (data->efield)
      return data->efield->getDt();
    return boost::none;
  }

  void setEfieldDt(const osh::Real dt) const override {
    if (data->efield)
      data->efield->setDt(dt);
    else
      throw std::logic_error("NO E-FIELD WHERE I CAN ASSIGN DT");
  }

  void setEfieldTolerances(double atol, double rtol, KSPNormType norm_type) {
    if (data->efield) {
      data->efield->setTolerances(atol, rtol, norm_type);
    } else {
      throw std::logic_error("E-Field is not in use.");
    }
  }
#endif // USE_PETSC

  void setDiffusionBoundaryActive(
      const mesh::diffusion_boundary_name &diffusion_boundary_name,
      const model::species_name &spec_id, bool set_active) override;

  bool getDiffusionBoundaryActive(
      const mesh::diffusion_boundary_name &diffusion_boundary_name,
      const model::species_name &spec_id);

  /**
   * Set a current clamp on a membrane
   * \param membrane membrane identifier
   * \param current stimulus to apply (Amps)
   */
  void setMembIClamp(const model::membrane_id& membrane, osh::Real current);


private:
  /**
   * Update the number of molecules of a certain species in the current rank
   *
   * In case is_uniform is false, the distribution of molecules among the
   * elements is stochastic following the probability: V_elem/V_rank where
   * V_elem is the volume of an element and V_rank is the volume of all the
   * elements owned by the rank
   *
   * In case is_uniform is true molecules are evenly split following the partial
   * volume ratio: V_elem/V_rank with stochastic rounding
   *
   * \param molecules Molecules pool instance
   * \param comp_id compartment identifier
   * \param species species identifier
   * \param num_molecules number of molecules to set
   * \param distribution distribution type. Check distributions.hpp for more information
   */
  void setOwnedCompCount(
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

  mesh_type &mesh;
  const osh::LOs elems2verts;
  const osh::Reals coords;
  osh::I64 num_iterations{};
  osh::Real state_time{};
  std::unique_ptr<Statedef> statedef;
  std::unique_ptr<SimulationInput<RNG, NumMolecules>> input;
  std::unique_ptr<SimulationData<SSA, RNG, NumMolecules, SearchMethod>> data;

  /// provide MPI rank that owned a given element, for diffusion debugging only
  osh::LOs element_ranks;
};

// explicit template instantiation declarations
extern template class OmegaHSimulation<SSAMethod::SSA, std::mt19937, osh::I32,
                                       NextEventSearchMethod::Direct>;
extern template class OmegaHSimulation<SSAMethod::RSSA, std::mt19937, osh::I32,
                                       NextEventSearchMethod::Direct>;
extern template class OmegaHSimulation<SSAMethod::SSA, std::mt19937, osh::I32,
                                       NextEventSearchMethod::GibsonBruck>;

extern template class OmegaHSimulation<SSAMethod::SSA, std::mt19937, osh::I64,
                                       NextEventSearchMethod::Direct>;
extern template class OmegaHSimulation<SSAMethod::RSSA, std::mt19937, osh::I64,
                                       NextEventSearchMethod::Direct>;
extern template class OmegaHSimulation<SSAMethod::SSA, std::mt19937, osh::I64,
                                       NextEventSearchMethod::GibsonBruck>;

extern template class OmegaHSimulation<SSAMethod::SSA, steps::rng::RNG, osh::I32,
                                       NextEventSearchMethod::Direct>;
extern template class OmegaHSimulation<SSAMethod::RSSA, steps::rng::RNG, osh::I32,
                                       NextEventSearchMethod::Direct>;
extern template class OmegaHSimulation<SSAMethod::SSA, steps::rng::RNG, osh::I32,
                                       NextEventSearchMethod::GibsonBruck>;

extern template class OmegaHSimulation<SSAMethod::SSA, steps::rng::RNG, osh::I64,
                                       NextEventSearchMethod::Direct>;
extern template class OmegaHSimulation<SSAMethod::RSSA, steps::rng::RNG, osh::I64,
                                       NextEventSearchMethod::Direct>;
extern template class OmegaHSimulation<SSAMethod::SSA, steps::rng::RNG, osh::I64,
                                       NextEventSearchMethod::GibsonBruck>;

extern template class Simulation<steps::rng::RNG>;

} // namespace dist
} // namespace steps
