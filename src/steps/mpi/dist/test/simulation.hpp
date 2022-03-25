#pragma once

#include <fstream>
#include <memory>
#include <random>

#include <mpi.h>

#include "math/distributions.hpp"
#include "mpi/dist/test/scenario.hpp"
#include "mpi/dist/test/simdef.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"
#include "rng/rng.hpp"


namespace steps {
namespace dist {

template <typename RNG> class Simulation {
public:
  Simulation(const ScenarioInput &t_scenario, DistMesh &t_mesh, RNG &t_rng,
             std::ostream &t_outstream);
  virtual ~Simulation() noexcept;

  virtual void setCompCount(
      const model::compartment_id& compartment,
      const model::species_name& species,
      osh::Real num_molecules,
      const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM) = 0;

  void setCompCount(
      const Simdef::compartment_counts_t& counts,
      const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM);
  virtual void setCompCount(
      const model::compartment_id& compartment,
      const std::vector<CompartmentCount>& counts,
      const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM) = 0;
  void setCompConc(
      const Simdef::compartment_concs_t& concentrations,
      const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM);
  virtual void setCompConc(
      const model::compartment_id& compartment,
      const std::vector<CompartmentConc>& concs,
      const math::DistributionMethod distribution = math::DistributionMethod::DIST_UNIFORM) = 0;

  virtual void setPatchCount(
      const model::patch_id& /*patch_id*/,
      const model::species_name& /*species*/,
      osh::Real /*num_molecules*/,
      const math::DistributionMethod = math::DistributionMethod::DIST_UNIFORM) {
      throw std::logic_error("NOT IMPLEMENTED");
  }
  virtual void setPatchCount(
      const Simdef::patch_counts_t& /*counts*/,
      const math::DistributionMethod = math::DistributionMethod::DIST_UNIFORM) {
      throw std::logic_error("NOT IMPLEMENTED");
  }
  virtual const std::vector<mesh::triangle_id_t> &getGHKBoundaries() const {
    throw std::logic_error("NOT IMPLEMENTED");
  }
  virtual osh::Reals getGHKCurrents() const {
    throw std::logic_error("NOT IMPLEMENTED");
  }
  virtual osh::Real getTotalGHKCurrent() const {
    throw std::logic_error("NOT IMPLEMENTED");
  }

  virtual void setKCstSReac(const model::patch_id & /*patchId*/,
                            container::surface_reaction_id /*reactionId*/,
                            osh::Real /*kCst*/) {
    throw std::logic_error("NOT IMPLEMENTED");
  }
  virtual void setTriCount(const model::patch_id & /*patch*/,
                           const model::species_name & /*species*/,
                           osh::Real /*num_molecules*/) {
    throw std::logic_error("NOT IMPLEMENTED");
  }
  virtual osh::Real
  getPatchCount(const model::patch_id & /*compartment*/,
                const model::species_name & /*species*/) const {
    throw std::logic_error("NOT IMPLEMENTED");
  }

  virtual void setPotential(osh::Real) {
    throw std::logic_error("NOT IMPLEMENTED");
  }

  virtual osh::Reals getPotentialOnVertices(const model::patch_id &) {
    throw std::logic_error("NOT IMPLEMENTED");
  }

  virtual osh::Real getMaxPotentialOnVertices(const model::patch_id &) {
    throw std::logic_error("NOT IMPLEMENTED");
  }

  virtual osh::Real getMinPotentialOnVertices(const model::patch_id &) {
    throw std::logic_error("NOT IMPLEMENTED");
  }

  virtual osh::Reals getPotentialOnTriangles(const model::patch_id &) {
    throw std::logic_error("NOT IMPLEMENTED");
  }

  virtual osh::Real getMaxPotentialOnTriangles(const model::patch_id &) {
    throw std::logic_error("NOT IMPLEMENTED");
  }

  virtual osh::Real getMinPotentialOnTriangles(const model::patch_id &) {
    throw std::logic_error("NOT IMPLEMENTED");
  }

  virtual osh::Real getLocalMaxPotentialOnVertices(const model::patch_id &) {
    throw std::logic_error("NOT IMPLEMENTED");
  }

  virtual osh::Real getLocalMinPotentialOnVertices(const model::patch_id &) {
    throw std::logic_error("NOT_IMPLEMENTED");
  }

#if USE_PETSC
  virtual boost::optional<osh::Real> getEfieldDt() const {
    throw std::logic_error("NOT_IMPLEMENTED");
  }

  virtual void setEfieldDt(const osh::Real) const {
    throw std::logic_error("NOT_IMPLEMENTED");
  }
#endif // USE_PETSC

  virtual osh::Real getCompConc(const model::compartment_id &compartment,
                                const model::species_name &species) const = 0;
  virtual osh::Real
  getOwnedCompConc(const model::compartment_id &compartment,
                   const model::species_name &species) const = 0;
  virtual osh::Real getCompCount(const model::compartment_id &compartment,
                                 const model::species_name &species) const = 0;

#if USE_PETSC
  virtual std::pair<mesh::triangle_ids, osh::Reals>
  getOhmicCurrents(const model::membrane_id &membrane_id,
                   const model::channel_id &l_id) const = 0;

  virtual osh::Real
  getTotalOhmicCurrent(const model::membrane_id &mem_id,
                       const model::channel_id &chan_id) const = 0;
#endif // USE_PETSC

  virtual osh::Real
  getOwnedCompCount(const model::compartment_id &compartment,
                    const model::species_name &species) const = 0;
  virtual std::pair<std::reference_wrapper<const mesh::tetrahedron_ids>,
                    std::vector<osh::LO>>
  getOwnedElemCount(const model::species_name &species) const = 0;

  virtual std::pair<std::vector<mesh::tetrahedron_global_id_t>,
                    std::vector<osh::LO>>
  getElemCount(const model::species_name &species) const = 0;

  virtual void setOwnedElementCount(const model::compartment_id &compartment,
                                    mesh::tetrahedron_id_t element,
                                    const model::species_name &species,
                                    osh::Real num_mols) = 0;

  virtual void setDiffOpBinomialThreshold(osh::Real threshold) = 0;

  virtual void exportMolStateToVTK(const std::string &filename) = 0;

  virtual osh::I64 getDiffOpExtent(bool local = false) const = 0;
  virtual osh::I64 getSSAOpExtent(bool local = false) const = 0;
  virtual osh::I64 getNIterations() const noexcept = 0;
  virtual osh::Real getIterationTimeStep() const noexcept = 0;

  virtual std::string createStateReport() = 0;

  virtual void init(std::unique_ptr<Statedef> &&statedef) = 0;
  virtual void reset() = 0;
  virtual void run(osh::Real end_time) = 0;
  /// \brief Emit the given message on every rank.
  /// the message is prefixed by "[RANK] "
  void log_all(const std::string &message) const;
  /// Emit the given message on rank 0 only
  void log_once(const std::string &message, bool force_stdout = false) const;
  /**
   * log progress to stdout
   *
   * progress = 100*i/tot %
   *
   * @param i progress index
   * @param tot progress total
   * @param name what is progressing
   */
  void log_progress(const double i, const double tot,
                    const std::string &name = "total") const;

  virtual void
  setDiffusionBoundaryActive(const mesh::diffusion_boundary_name &boundary_id,
                             const model::species_name &species,
                             bool set_active) = 0;
  const std::vector<unsigned int> &get_diffusion_rank_exchanges() const
      noexcept {
    return diffusion_rank_exchanges;
  }

  virtual SSAMethod ssaMethod() const noexcept = 0;

  void log_diffusion_exchanges() const;

  const DistMesh &getMesh() const noexcept { return mesh; }

  const ScenarioInput &getScenario() const noexcept { return scenario; }

  DistMesh &getMesh() noexcept { return mesh; }

  double getElapsedSSA() const noexcept { return reactions_timer; }

  double getElapsedDiff() const noexcept {
    return diffusions_timer;
  }

  double getElapsedEField() const noexcept { return efield_timer; }

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
  /// tell whether diffusion debugging is enabled or not
  bool debug_diffusion_{false};
  /// matrix N, N with N the number of ranks
  /// M[i, j] provides number of molecules that went from rank i to rank j
  std::vector<unsigned int> diffusion_rank_exchanges;
  /** \} */
  const ScenarioInput &scenario;
  DistMesh &mesh;
  RNG &rng;
  /// total time in seconds spent performing reactions
  double reactions_timer{};
  /// total time in seconds spent performing diffusions
  double diffusions_timer{};
  /// total time in seconds spent performing efield
  double efield_timer{};

  std::ostream &outstream;
};

// explicit template instantiation declarations
extern template class Simulation<std::mt19937>;

extern template class Simulation<steps::rng::RNG>;

} // namespace dist
} // namespace steps
