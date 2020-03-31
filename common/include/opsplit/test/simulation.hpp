#pragma once

#include <memory>
#include <random>

#include <boost/timer/timer.hpp>

#include "opsplit/fwd.hpp"
#include "opsplit/test/scenario.hpp"
#include "opsplit/test/simdef.hpp"

namespace zee {

using cpu_times = boost::timer::cpu_times;

template <typename RNG>
class Simulation {
  public:
    Simulation(const ScenarioInput& t_scenario, DistMesh& t_mesh, RNG& t_rng);
    virtual ~Simulation() noexcept;

    virtual void setCompCount(const model::compartment_id& compartment,
                              const model::specie_name& specie,
                              PetscScalar num_molecules) = 0;

    virtual void setCompCount(const Simdef::compartment_counts_t& counts) = 0;
    virtual void setCompConc(const Simdef::compartment_concs_t& concs) = 0;

    virtual void setPatchCount(const model::patch_id& /*patch_id*/,
                               const model::specie_name& /*specie*/,
                               PetscScalar /*num_molecules*/) {
        throw std::logic_error("NOT IMPLEMENTED");
    }
    virtual void setPatchCount(const Simdef::patch_counts_t& /*counts*/) {
        throw std::logic_error("NOT IMPLEMENTED");
    }
    virtual PetscScalar getPatchCount(const model::patch_id& /*compartment*/,
                                      const model::specie_name& /*specie*/) const {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    virtual PetscScalar getCompConc(const model::compartment_id& compartment,
                                    const model::specie_name& specie) const = 0;
    virtual PetscScalar getOwnedCompConc(const model::compartment_id& compartment,
                                         const model::specie_name& specie) const = 0;
    virtual PetscScalar getCompCount(const model::compartment_id& compartment,
                                     const model::specie_name& specie) const = 0;
    virtual PetscScalar getOwnedCompCount(const model::compartment_id& compartment,
                                          const model::specie_name& specie) const = 0;
    virtual void setOwnedElementCount(const model::compartment_id& compartment,
                                      PetscInt element,
                                      const model::specie_name& specie,
                                      PetscScalar num_mols) = 0;

    virtual void setDiffOpBinomialThreshold(PetscScalar threshold) = 0;

    virtual void exportMolStateToVTK(const std::string& filename) = 0;

    virtual PetscInt64 getDiffOpExtent() const = 0;
    virtual PetscInt64 getSSAOpExtent() const = 0;
    virtual PetscInt64 getNIterations() const noexcept = 0;
    virtual PetscScalar getIterationTimeStep() const noexcept = 0;

    virtual std::string createStateReport() = 0;

    virtual void init(std::unique_ptr<Statedef>&& statedef) = 0;
    virtual void reset() = 0;
    virtual void run(PetscScalar end_time) = 0;
    /// \brief Emit the given message on every rank.
    /// the message is prefixed by "[RANK] "
    virtual void log_all(const std::string& message) const = 0;
    /// Emit the given message on rank 0 only
    virtual void log_once(const std::string& message) const = 0;

    const std::vector<unsigned int>& get_diffusion_rank_exchanges() const noexcept {
        return diffusion_rank_exchanges;
    }

    void log_diffusion_exchanges() const;

    const DistMesh& getMesh() const noexcept {
        return mesh;
    }

    DistMesh& getMesh() noexcept {
        return mesh;
    }

    cpu_times getElapsedSSA() const noexcept {
        return reactions_timer.elapsed();
    }

    cpu_times getElapsedDiff() const noexcept {
        return diffusion_timer.elapsed();
    }

    const int comm_rank;
    const int comm_size;

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
    const ScenarioInput& scenario;
    DistMesh& mesh;
    RNG rng;
    boost::timer::cpu_timer reactions_timer;
    boost::timer::cpu_timer diffusion_timer;
};

// explicit template instantiation declarations
extern template class Simulation<std::mt19937>;

}  // namespace zee
