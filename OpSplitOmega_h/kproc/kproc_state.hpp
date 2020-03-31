#pragma once

#include <ostream>
#include <set>
#include <vector>

#include "../common.hpp"
#include "../mesh.hpp"
#include "reactions.hpp"
#include "surface_reactions.hpp"

namespace zee {

class Statedef;

/**
 * \brief Class encapsulates the Kinetic process state of a STEPS simulation.
 *
 * This class encapsulates all Kinetic Processes (KProc) such as Reac
 * in the simulation.
 *
 * This class has no corresponding class in STEPS as KProc state is stored in
 * TetOpSplit/Tetexact solver class without encapsulation.
 */

using KProcID = std::pair<KProcType, size_t>;

/**
 * Hold all dependencies of a given kinetic process
 */
struct KProcDeps {
    using Propensities = std::vector<KProcID>;
    using Occupancy = std::vector<MolState::ElementID>;

    /**
     * \brief the kinetic processes that depend on a change of propensity
     */
    Propensities propensities;

    /**
     * \brief the pairs (mesh element, specie) that depend on a change of occupancy
     */
    Occupancy occupancy;
};


class KProcState {
  public:
    /** KProcState constructor
     * \param statedef model definition
     * \param mesh distributed mesh
     * \param discovery if true, dry run KProcState to record dependencies.
     */
    template <osh::Int Dim>
    KProcState(const Statedef& statedef, OmegaHMesh<Dim>& mesh, bool discovery = false);

    void updateAllPropensities(const MolState& mol_state);
    void updatePropensities(const MolState& mol_state, const KProcID& event);
    void updateOccupancy(DiffusionVariables& diffusionVariables,
                         MolState& mol_state,
                         PetscScalar simulation_time,
                         const KProcID& event) const;
    void updateMolState(MolState& mol_state, const KProcID& event) const;

    inline PetscScalar getSSA_A0() const noexcept {
        return ssa_a0;
    }

    inline const kproc::Reactions& reactions() const noexcept;

    inline const kproc::SurfaceReactions& surfaceReactions() const noexcept;

    template <class Generator>
    KProcID getNextEvent(Generator& rng) const;

    std::pair<osh::LOs, boost::optional<osh::LOs>> getNumberOfSpeciesPerOwnedElement() const;

    const osh::LOs& getNumberOfSpeciesPerElement() const noexcept {
        return reactions().getNumberOfSpeciesPerElement();
    }

    const KProcDeps& extractDependenciesFromEvent(const KProcID& event) const;


    void report(std::ostream& report_stream) const;

  private:
    using DependenciesMap = std::map<MolState::ElementID, std::vector<KProcID>>;

    template <typename KineticProcesses>
    void cacheDependencies(const KineticProcesses& processes,
                           const DependenciesMap& dependency_map,
                           std::vector<KProcDeps>& dependencies);

    template <typename KineticProcesses>
    void collateDependencies(const KineticProcesses& processes,
                             DependenciesMap& dependency_map) const;


    void setupDependencies();

    PetscScalar updateSSA_A0();

    [[noreturn]] void updateOccupancy(mesh::boundary_id element,
                                      container::specie_id specie,
                                      DiffusionVariables& diffusionVariables,
                                      MolState& mol_state,
                                      PetscScalar simulation_time) const;
    void updateOccupancy(mesh::element_id element,
                         container::specie_id specie,
                         DiffusionVariables& diffusionVariables,
                         MolState& mol_state,
                         PetscScalar simulation_time) const;


    std::vector<KProcDeps> reactions_dependencies_;
    std::vector<KProcDeps> surface_reactions_dependencies_;

    const Statedef& pStatedef;

    PetscScalar ssa_a0{};

    kproc::Reactions reactions_;
    kproc::SurfaceReactions surface_reactions_;
};

inline const kproc::Reactions& KProcState::reactions() const noexcept {
    return reactions_;
}

inline const kproc::SurfaceReactions& KProcState::surfaceReactions() const noexcept {
    return surface_reactions_;
}

extern template KProcState::KProcState(const Statedef& statedef,
                                       OmegaHMesh<2>& mesh,
                                       bool discovery);
extern template KProcState::KProcState(const Statedef& statedef,
                                       OmegaHMesh<3>& mesh,
                                       bool discovery);

extern template KProcID KProcState::getNextEvent(std::mt19937& rng) const;

}  // namespace zee
