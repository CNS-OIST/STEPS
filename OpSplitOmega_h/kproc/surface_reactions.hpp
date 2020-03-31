#pragma once
/**
 * \file surface_reactions.hpp
 * Provide the \a SurfaceReactions class
 */

#include <cmath>
#include <petscsys.h>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../common.hpp"
#include "../mesh.hpp"
#include "../mol_state.hpp"
#include "reactions.hpp"


#include <opsplit/compdef.hpp>
#include <opsplit/diffdef.hpp>
#include <opsplit/patchdef.hpp>
#include <opsplit/reacdef.hpp>
#include <opsplit/sreacdef.hpp>
#include <opsplit/statedef.hpp>
#include <opsplit/vocabulary.hpp>

namespace zee {

namespace kproc {

class SurfaceReactions {
  public:
    /**
     * Surface reactions constructor
     * \param statedef model definition
     * \param mesh distributed mesh
     * \param discovery if true, dry run KProcState to record dependencies.
     */
    template <osh::Int Dim>
    SurfaceReactions(const Statedef& statedef, OmegaHMesh<Dim>& mesh, bool discovery = false);

    inline size_t size() const noexcept;

    inline const SReacdef& getReacDef(size_t index) const noexcept {
        return reacdefs_[index];
    }

    /**
     * @brief Returns a list of molecular state elements that effect the propensity of the reaction
     * 'index'.
     */
    const std::vector<MolState::ElementID>& getPropensityDependency(size_t index) const {
        return reaction_lhs_[index];
    }

    /**
     * @brief Returns a list of molecular state elements updated in the
     * event of the reaction identified by the index occuring.
     */
    std::pair<std::reference_wrapper<const std::vector<model::region_id>>,
              std::reference_wrapper<const std::vector<MolState::ElementID>>>
    getMolStateElementsUpdates(size_t index) const {
        return std::make_pair(std::cref(region_ids_[index]), std::cref(reaction_upd_[index]));
    }

    inline mesh::element_id getInnerCompartmentElementId(size_t index) const noexcept {
        return inner_compartment_element_id_[index];
    }

    inline const boost::optional<mesh::element_id>& getOuterCompartmentElementId(size_t index) const
        noexcept {
        return outer_compartment_element_id_[index];
    }

    inline PetscScalar getRate(size_t index) const noexcept;

    PetscScalar updateRate(const MolState& mol_state, size_t index);

    PetscScalar computeRate(const MolState& mol_state, size_t index) const;

    void apply(MolState& mol_state, size_t index) const;

    inline const osh::LOs& getNumberOfSpeciesPerBoundaryOwnedElement() const noexcept {
        return num_species_per_boundary_;
    }

    inline KProcType getKProcType() const {
        return KProcType::SReac;
    }

  private:
    using Stoichiometry = std::vector<PetscInt>;

    template <SReacdef::PoolChangeType PoolChange>
    std::tuple<std::vector<MolState::ElementID>,
               SurfaceReactions::Stoichiometry,
               std::vector<model::region_id>>
    reactionMolStateDependencyAndStoichiometry(
        const SReacdef& reacdef,
        mesh::boundary_id patch_element_id,
        mesh::element_id inner_compartment_element_id,
        const boost::optional<mesh::element_id>& outer_compartment_element) const;

    /**
     *  \brief Reaction rate constant in the particular mesh boundary.
     * Multiplier of the propensity of the reaction in the mesh boundary.
     *
     * \tparam Dim mesh dimension
     * \param mesh distributed mesh
     * \param reacdef definition of the reaction
     * \param boundary index of a mesh boundary
     * \param inner_compartment_element mesh element in the inner compartment
     * \param outer_compartment_element optional mesh element in the outer compartment
     * \return the reaction rate constant.
     */
    template <osh::Int Dim>
    PetscScalar compute_ccst(OmegaHMesh<Dim>& mesh,
                             const SReacdef& reacdef,
                             mesh::boundary_id boundary,
                             mesh::element_id inner_compartment_element,
                             const boost::optional<mesh::element_id>& outer_compartment_element);

    /// local_index of element in inner compartment for ith surface reaction
    std::vector<mesh::element_id> inner_compartment_element_id_;
    /// local_index of element in outer compartment for ith surface reaction
    std::vector<boost::optional<mesh::element_id>> outer_compartment_element_id_;

    /**
     * \name
     * size of dim 1: number of surface reactions
     * size of dim 2: number of reactants involved
     * \{
     */
    /// specie-element id of each reactant in the ith surface reaction
    std::vector<std::vector<MolState::ElementID>> reaction_lhs_;
    /// stoichiometry coefficient of each reactant in the ith surface reaction
    std::vector<Stoichiometry> stoichiometry_lhs_;
    /** \} */

    /**
     * \name
     * size of dim 1: number of surface reactions
     * size of dim 2: number of species involved
     * \{
     */
    /// specie-element id of each specie in the ith surface reaction
    std::vector<std::vector<MolState::ElementID>> reaction_upd_;
    /// stoichiometry difference of each specie in the ith surface reaction
    /// \TODO TCL: might be able to reduce memory footprint by accessing directly in reacdefs
    std::vector<Stoichiometry> stoichiometry_upd_;
    /// patch/compartment for each specie involved in the ith surface reaction
    std::vector<std::vector<model::region_id>> region_ids_;
    /** \} */

    /// number of species in the ith boundary
    osh::LOs num_species_per_boundary_;
    /// reaction definition for the ith surface reaction
    std::vector<std::reference_wrapper<SReacdef>> reacdefs_;
    /// propensity rate constant for ith surface reaction
    std::vector<PetscScalar> ccsts_;
    /// rate for ith surface reaction
    std::vector<PetscScalar> rates_;
};

inline size_t SurfaceReactions::size() const noexcept {
    return reacdefs_.size();
}

inline PetscScalar SurfaceReactions::getRate(size_t index) const noexcept {
    return rates_[index];
}

// explicit template instantiation declarations
extern template SurfaceReactions::SurfaceReactions(const zee::Statedef& statedef,
                                                   zee::OmegaHMesh<2>& mesh,
                                                   bool discovery);
extern template SurfaceReactions::SurfaceReactions(const zee::Statedef& statedef,
                                                   zee::OmegaHMesh<3>& mesh,
                                                   bool discovery);

}  // namespace kproc

}  // namespace zee
