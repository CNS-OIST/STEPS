#pragma once
/**
 * \file reactions.hpp
 * Provide the \a Reactions class
 */

#include <cmath>
#include <functional>
#include <random>
#include <vector>

#include "../mol_state.hpp"

#include "fwd.hpp"
#include "geom/dist/fwd.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/fwd.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "reactions_iterator.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist::kproc {

class Reactions {
  public:
    using iterator_type = reactions_iterator<Reactions>;
    using const_iterator_type = reactions_iterator<const Reactions>;
    /**
     *
     * \param statedef model definition
     * \param mesh distributed mesh object
     */
    Reactions(const Statedef& statedef, DistMesh& mesh, MolState& mol_state);

    inline size_t size() const noexcept;
    inline mesh::tetrahedron_id_t getOwnerPoint(size_t index) const noexcept;
    inline const Reacdef& getReacDef(size_t index) const noexcept;

    void report(std::ostream& report_stream, size_t index) const;

    /**
     * \brief Compute the rate of the KProc.
     */
    osh::Real computeRate(const MolState& mol_state, size_t index) const;

    const std::vector<MolStateElementID>& updateMolStateAndOccupancy(MolState& mol_state,
                                                                     size_t index,
                                                                     osh::Real event_time) const;

    /**
     * \brief Returns a list of molecular state elements that effect the
     * propensity of the reaction index.
     */
    inline const std::vector<MolStateElementID>& getPropensityDependency(
        size_t index) const noexcept {
        return reactions_lhs_[index];
    }

    /**
     * \brief Returns a list of molecular state elements updated in the
     * event of the reaction identified by the index occuring.
     */
    inline const std::vector<MolStateElementID>& getMolStateElementsUpdates(
        size_t index) const noexcept {
        return reactions_upd_[index];
    }

    inline KProcType getKProcType() const noexcept {
        return KProcType::Reac;
    }

    /// \return an iterator to the beginning
    inline iterator_type begin() noexcept;

    /// \return an iterator to the end
    inline iterator_type end() noexcept;

    /// \return an iterator to the beginning
    inline const_iterator_type begin() const noexcept;

    /// \return an iterator to the end
    inline const_iterator_type end() const noexcept;

  private:
    /**
     * \brief Reaction rate constant in the particular mesh element.
     * Multiplier of the propensity of the reaction in the mesh element.
     *
     * \param reacdef definition of the reaction
     * \param element index of a mesh element
     * \return the reaction rate constant.
     */
    osh::Real compute_ccst(const Reacdef& reacdef, mesh::tetrahedron_id_t element) const;

    std::vector<std::reference_wrapper<Reacdef>> reacdefs_;
    std::vector<mesh::tetrahedron_id_t> ownerPoints_;
    std::vector<osh::Real> ccsts_;
    std::vector<std::vector<MolStateElementID>> reactions_upd_;
    /// species-element id of each reactant in the ith surface reaction
    std::vector<std::vector<MolStateElementID>> reactions_lhs_;
    std::vector<std::vector<osh::I64>> stoichiometry_change_;

    const Measure& measureInfo;
};

inline size_t Reactions::size() const noexcept {
    return ownerPoints_.size();
}

inline mesh::tetrahedron_id_t Reactions::getOwnerPoint(size_t index) const noexcept {
    return ownerPoints_[index];
}

inline const Reacdef& Reactions::getReacDef(size_t index) const noexcept {
    return reacdefs_[index];
}

inline Reactions::iterator_type Reactions::begin() noexcept {
    return {*this, size_t{0}};
}

inline Reactions::iterator_type Reactions::end() noexcept {
    return {*this, this->size()};
}

inline Reactions::const_iterator_type Reactions::begin() const noexcept {
    return {*this};
}

inline Reactions::const_iterator_type Reactions::end() const noexcept {
    return {*this, this->size()};
}

}  // namespace steps::dist::kproc
