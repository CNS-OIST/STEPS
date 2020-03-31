#pragma once
/**
 * \file reactions.hpp
 * Provide the \a Reactions class
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

#include "opsplit/compdef.hpp"
#include "opsplit/diffdef.hpp"
#include "opsplit/reacdef.hpp"
#include "opsplit/statedef.hpp"
#include "opsplit/vocabulary.hpp"

namespace zee {


namespace kproc {

class Reactions {
  public:
    using index_type = size_t;
    /**
     *
     * \param statedef model definition
     * \param mesh distributed mesh object
     * \param discovery if true, dry run reactions construction to record mol state dependencies
     */
    template <osh::Int Dim>
    Reactions(const Statedef& statedef, OmegaHMesh<Dim>& mesh, bool discovery = false);

    inline size_t size() const noexcept;
    inline mesh::element_id getOwnerPoint(size_t index) const noexcept;
    inline PetscScalar getRate(size_t index) const noexcept;
    inline const Reacdef& getReacDef(size_t index) const noexcept;

    void report(std::ostream& report_stream, size_t index) const;

    /**
     * \brief Update the rate of the KProc and/or occupancy.
     *
     * Update the rate of the KProc, if period is not 0.0,
     * and the KProc is diff/sdiff, also update the occupancy.
     *
     * Note that this design is different from STEPS 3.x,
     * where the occupancy data is stored in Tet/Tri.
     */
    PetscScalar updateRate(const MolState& mol_state, size_t index);

    /**
     * \brief Compute the rate of the KProc and/or occupancy.
     */
    PetscScalar computeRate(const MolState& mol_state, size_t index) const;

    void apply(MolState& mol_state, size_t index) const;

    /**
     * @brief Returns a list of molecular state elements that effect the propensity of the reaction
     * index.
     */
    std::vector<MolState::ElementID> getPropensityDependency(size_t index) const;

    /**
     * @brief Returns a list of molecular state elements updated in the
     * event of the reaction identified by the index occuring.
     */
    std::pair<std::vector<model::region_id>, std::vector<MolState::ElementID>>
    getMolStateElementsUpdates(size_t index) const;

    const osh::LOs& getNumberOfSpeciesPerOwnedElement() const noexcept {
        return num_species_per_owned_element_;
    }

    const osh::LOs& getNumberOfSpeciesPerElement() const noexcept {
        return num_species_per_element_;
    }

    inline KProcType getKProcType() const noexcept {
        return KProcType::Reac;
    }

  private:
    /**
     * \brief Reaction rate constant in the particular mesh element.
     * Multiplier of the propensity of the reaction in the mesh element.
     *
     * \param reacdef definition of the reaction
     * \param element index of a mesh element
     * \return the reaction rate constant.
     */
    PetscScalar compute_ccst(const Reacdef& reacdef, mesh::element_id element) const;

    std::vector<Reacdef> reacdefs_;
    std::vector<mesh::element_id> ownerPoints_;
    std::vector<PetscScalar> ccsts_;
    std::vector<PetscScalar> rates_;
    osh::LOs num_species_per_element_;
    osh::LOs num_species_per_owned_element_;
    const MeasureInfo& measureInfo;
};

inline size_t Reactions::size() const noexcept {
    return ownerPoints_.size();
}

inline mesh::element_id Reactions::getOwnerPoint(size_t index) const noexcept {
    return ownerPoints_[index];
}

inline PetscScalar Reactions::getRate(size_t index) const noexcept {
    return rates_[index];
}

inline const Reacdef& Reactions::getReacDef(size_t index) const noexcept {
    return reacdefs_[index];
}

// explicit template instantiation declarations
extern template Reactions::Reactions(const zee::Statedef& statedef,
                                     zee::OmegaHMesh<2>& mesh,
                                     bool discovery);
extern template Reactions::Reactions(const zee::Statedef& statedef,
                                     zee::OmegaHMesh<3>& mesh,
                                     bool discovery);

}  // namespace kproc
}  // namespace zee
