#pragma once

#include "distributions.hpp"
#include "kproc/diffusions.hpp"
#include "mesh.hpp"
#include "mol_state.hpp"
#include "operator/rssa_operator.hpp"
#include "operator/ssa_operator.hpp"
#include "opsplit/test/simulation.hpp"
#include "simulation_data.hpp"

namespace zee {

template <osh::Int Dim, SSAMethod SSA, typename RNG>
class OmegaHSimulation: public Simulation<RNG> {
  public:
    using super_type = Simulation<RNG>;
    OmegaHSimulation(const ScenarioInput& t_scenario, OmegaHMesh<Dim>& t_mesh, RNG& t_rng);

    void setCompCount(const model::compartment_id& compartment,
                      const model::specie_name& specie,
                      PetscScalar num_molecules) override;
    void setOwnedCompCount(const model::compartment_id& compartment,
                           const model::specie_name& specie,
                           PetscScalar num_molecules);
    void setCompCount(const Simdef::compartment_counts_t& counts) override;
    void setCompConc(const model::compartment_id& compartment,
                     const model::specie_name& specie,
                     PetscScalar concentration);
    void setCompConc(const Simdef::compartment_concs_t& concs) override;
    void setOwnedElementCount(const model::compartment_id& compartment,
                              PetscInt element,
                              const model::specie_name& specie,
                              PetscScalar num_molecules) override;

    PetscScalar getCompConc(const model::compartment_id& compartment,
                            const model::specie_name& specie) const override;
    PetscScalar getOwnedCompConc(const model::compartment_id& compartment,
                                 const model::specie_name& specie) const override;
    PetscScalar getCompCount(const model::compartment_id& compartment,
                             const model::specie_name& specie) const override;

    void setPatchCount(const model::patch_id& patch_id,
                       const model::specie_name& specie,
                       PetscScalar num_molecules) override;

    void setPatchCount(const Simdef::patch_counts_t& counts) override;

    PetscScalar getPatchCount(const model::patch_id& compartment,
                              const model::specie_name& specie) const override;

    PetscScalar getOwnedCompCount(const model::compartment_id& compartment,
                                  const model::specie_name& specie) const override;

    void setDiffOpBinomialThreshold(PetscScalar threshold) override;
    PetscScalar getIterationTimeStep() const noexcept override;

    void exportMolStateToVTK(const std::string& filename) override;

    PetscInt64 getDiffOpExtent() const override;
    PetscInt64 getSSAOpExtent() const override;
    PetscInt64 getNIterations() const noexcept override;

    std::string createStateReport() override;

    void init(std::unique_ptr<Statedef>&& t_statedef) override;
    void reset() override;
    void run(PetscScalar end_time) override;

    void log_all(const std::string& message) const override;
    void log_once(const std::string& message) const override;

  private:
    /**
     * Given a number of molecules to distribute over the entire mesh,
     * distribute them among the available ranks in regards of their measure.
     * \param comp_id compartment identifier
     * \param num_molecules Total number of molecules to spread on the mesh
     * \param result output container where number of molecules per rank is stored
     */
    void distribute_num_molecules_on_ranks(const model::compartment_id& comp_id,
                                           PetscScalar num_molecules,
                                           std::vector<PetscScalar>& result);

    /**
     * Update the number of molecules of a certain specie in the current rank
     * among the elements in regarding of their measure.
     * \param molecules Molecules pool instance
     * \param comp_id compartement identifier
     * \param specie specie identifier
     * \param num_molecules number of molecules to set
     */
    void setOwnedCompCount(MolState& molecules,
                           const model::compartment_id& comp_id,
                           const model::specie_name& specie,
                           PetscScalar num_molecules);

    void initialize_discretized_rates();

    OmegaHMesh<Dim>& mesh;
    const osh::LOs elems2verts;
    const osh::Reals coords;
    PetscInt64 num_iterations{};
    PetscScalar state_time{};
    std::unique_ptr<Statedef> statedef;
    std::unique_ptr<SimulationInput<RNG>> input;
    std::unique_ptr<SimulationData<Dim, SSA, RNG>> data;

    /// provide MPI rank that owned a given element, for diffusion debugging only
    osh::Read<osh::LO> element_ranks;
};

// explicit template instantiation declarations
extern template class OmegaHSimulation<2, SSAMethod::SSA, std::mt19937>;
extern template class OmegaHSimulation<2, SSAMethod::RSSA, std::mt19937>;
extern template class OmegaHSimulation<3, SSAMethod::SSA, std::mt19937>;
extern template class OmegaHSimulation<3, SSAMethod::RSSA, std::mt19937>;

}  // namespace zee
