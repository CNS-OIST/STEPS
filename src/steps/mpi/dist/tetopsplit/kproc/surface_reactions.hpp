#pragma once
/**
 * \file surface_reactions.hpp
 * Provide the \a SurfaceReactions class
 */

#include <cmath>
#include <random>
#include <vector>

#include <Omega_h_adj.hpp>

#include "../mol_state.hpp"
#include "geom/dist/distmesh.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "mpi/dist/tetopsplit/kproc/surface_reactions.hpp"
#include "reactions.hpp"
#include "reactions_iterator.hpp"
#include "util/vocabulary.hpp"

namespace steps {
namespace dist {
namespace kproc {

/**
 * \brief SurfaceReactionBase is a container for regular surface reactions and
 * voltage dependent surface reactions. It factors common functionalities
 * between the two types of surface reactions.
 *
 * \tparam PropensityType either a scalar or a function that takes as input a
 * potential
 */
template <typename PropensityType> class SurfaceReactionsBase {
public:
  using iterator_type =
      reactions_iterator<SurfaceReactionsBase<PropensityType>>;
  using const_iterator_type =
      reactions_iterator<const SurfaceReactionsBase<PropensityType>>;
  /**
   * Surface reactions constructor
   * \param statedef model definition
   * \param mesh distributed mesh
   */
  template <typename NumMolecules>
  SurfaceReactionsBase(const Statedef& statedef, DistMesh& mesh, MolState<NumMolecules>& mol_state);

  inline size_t size() const noexcept;

  inline const SReacdefBase<PropensityType> &getReacDef(size_t index) const
      noexcept;

  /**
   * \brief Returns a list of molecular state elements that effect the
   * propensity of the reaction 'index'.
   */
  inline const std::vector<MolStateElementID> &
  getPropensityDependency(size_t index) const noexcept;

  /**
   * \return kind of surface reaction
   */
  static std::string name();

  /**
   * \brief Returns a list of molecular state elements updated in the
   * event of the reaction identified by the index occurring.
   */
  inline const std::vector<MolStateElementID> &
  getMolStateElementsUpdates(size_t index) const noexcept;

  inline mesh::tetrahedron_id_t getInnerCompartmentElementId(size_t index) const
      noexcept;

  inline const boost::optional<mesh::tetrahedron_id_t> &
  getOuterCompartmentElementId(size_t index) const noexcept;


  /** Computes the exchange rate r in moles/s
   *
   * The link with the current I is: r = I/(z*F)
   * where:
   *
   * - z is the valence
   * - F is the Faraday constant = electron_charge * N_avogadro
   *
   * @tparam NumMolecules
   */
  template <typename NumMolecules>
  osh::Real computeRate(const MolState<NumMolecules> &mol_state,
                        size_t index) const;

  template <typename MoleculeType>
  void apply(MolState<MoleculeType> &mol_state, size_t index) const;

  template <typename NumMolecules>
  const std::vector<MolStateElementID>& updateMolStateAndOccupancy(
      MolState<NumMolecules>& mol_state,
      size_t index,
      const osh::Real event_time) const;

  /**
   * \name Iterators
   * \{
   */

  /// \return an iterator to the beginning
  inline iterator_type begin() noexcept;

  /// \return an iterator to the end
  inline iterator_type end() noexcept;

  /// \return an iterator to the beginning
  inline const_iterator_type begin() const noexcept;

  /// \return an iterator to the end
  inline const_iterator_type end() const noexcept;

  /** \} */

  /**
   * \brief A report for surf. reac. index
   *
   * \param report_stream stream
   * \param index index
   */
  void report(std::ostream &report_stream, size_t index) const;

  const std::vector<mesh::triangle_id_t> &boundaries() const noexcept {
    return boundary_id_;
  }

protected:
  /**
   * \brief Factor of the kinetic constant.
   * \param index kproc index
   * \return the geometric part of the kin constant.
   */
  osh::Real kinConstantGeomFactor(const DistMesh &mesh, size_t index) const;

  using Stoichiometry = std::vector<osh::I64>;

  template <typename SReacdefBase<PropensityType>::PoolChangeType PoolChange>
  std::tuple<std::vector<MolStateElementID>, Stoichiometry,
             std::vector<model::region_id>>
  reactionMolStateDependencyAndStoichiometry(
      const SReacdefBase<PropensityType> &reacdef,
      mesh::triangle_id_t patch_element_id,
      mesh::tetrahedron_id_t inner_compartment_element_id,
      const boost::optional<mesh::tetrahedron_id_t> &outer_compartment_element)
      const;

  /// local_index of element on the patch for ith surface reaction
  std::vector<mesh::triangle_id_t> boundary_id_;
  /// local_index of element in inner compartment for ith surface reaction
  std::vector<mesh::tetrahedron_id_t> inner_compartment_element_id_;
  /// local_index of element in outer compartment for ith surface reaction
  std::vector<boost::optional<mesh::tetrahedron_id_t>>
      outer_compartment_element_id_;

  /**
   * \name
   * size of dim 1: number of surface reactions
   * size of dim 2: number of reactants involved
   * \{
   */
  /// species-element id of each reactant in the ith surface reaction
  std::vector<std::vector<MolStateElementID>> reaction_lhs_;
  /// stoichiometry coefficient of each reactant in the ith surface reaction
  std::vector<Stoichiometry> stoichiometry_lhs_;
  /** \} */

  /**
   * \name
   * size of dim 1: number of surface reactions
   * size of dim 2: number of species involved
   * \{
   */
  /// species-element id of each species in the ith surface reaction
  std::vector<std::vector<MolStateElementID>> reaction_upd_;
  /// stoichiometry difference of each species in the ith surface reaction
  /// \TODO TCL: might be able to reduce memory footprint by accessing directly
  /// in reacdefs
  std::vector<Stoichiometry> stoichiometry_upd_;
  /** \} */

  /// reaction definition for the ith surface reaction
  std::vector<std::reference_wrapper<SReacdefBase<PropensityType>>> reacdefs_;
  /// propensity rate constant for ith surface reaction
  std::vector<osh::Real> ccsts_;

  const Statedef &state_def_;
};

template <typename PropensityType>
inline size_t SurfaceReactionsBase<PropensityType>::size() const noexcept {
  return reacdefs_.size();
}

template <typename PropensityType>
inline const SReacdefBase<PropensityType> &
SurfaceReactionsBase<PropensityType>::getReacDef(size_t index) const noexcept {
  return reacdefs_[index];
}

template <typename PropensityType>
const std::vector<MolStateElementID> &
SurfaceReactionsBase<PropensityType>::getPropensityDependency(
    size_t index) const noexcept {
  return reaction_lhs_[index];
}

template <typename PropensityType>
inline const std::vector<MolStateElementID> &
SurfaceReactionsBase<PropensityType>::getMolStateElementsUpdates(
    size_t index) const noexcept {
  return reaction_upd_[index];
}

template <typename PropensityType>
inline mesh::tetrahedron_id_t
SurfaceReactionsBase<PropensityType>::getInnerCompartmentElementId(
    size_t index) const noexcept {
  return inner_compartment_element_id_[index];
}

template <typename PropensityType>
inline const boost::optional<mesh::tetrahedron_id_t> &
SurfaceReactionsBase<PropensityType>::getOuterCompartmentElementId(
    size_t index) const noexcept {
  return outer_compartment_element_id_[index];
}

template <typename PropensityType>
inline typename SurfaceReactionsBase<PropensityType>::iterator_type
SurfaceReactionsBase<PropensityType>::begin() noexcept {
  return {*this};
}

template <typename PropensityType>
inline typename SurfaceReactionsBase<PropensityType>::iterator_type
SurfaceReactionsBase<PropensityType>::end() noexcept {
  return {*this, this->size()};
}

template <typename PropensityType>
inline typename SurfaceReactionsBase<PropensityType>::const_iterator_type
SurfaceReactionsBase<PropensityType>::begin() const noexcept {
  return {*this};
}

template <typename PropensityType>
inline typename SurfaceReactionsBase<PropensityType>::const_iterator_type
SurfaceReactionsBase<PropensityType>::end() const noexcept {
  return {*this, this->size()};
}

//-------------------------------------------------------

// explicit template instantiation declarations


extern template osh::Real
SurfaceReactionsBase<SReacInfo>::computeRate(const MolState<osh::LO> &mol_state,
                                             size_t index) const;
extern template osh::Real
SurfaceReactionsBase<VDepInfo>::computeRate(const MolState<osh::LO> &mol_state,
                                            size_t index) const;
extern template osh::Real
SurfaceReactionsBase<SReacInfo>::computeRate(const MolState<osh::GO> &mol_state,
                                             size_t index) const;
extern template osh::Real
SurfaceReactionsBase<VDepInfo>::computeRate(const MolState<osh::GO> &mol_state,
                                            size_t index) const;
extern template osh::Real
SurfaceReactionsBase<GHKInfo>::computeRate(const MolState<osh::LO> &mol_state,
                                           size_t index) const;
extern template osh::Real
SurfaceReactionsBase<GHKInfo>::computeRate(const MolState<osh::GO> &mol_state,
                                           size_t index) const;

extern template osh::Real
SurfaceReactionsBase<SReacInfo>::kinConstantGeomFactor(const DistMesh &mesh,
                                                       size_t index) const;

extern template const std::vector<MolStateElementID>&
SurfaceReactionsBase<SReacInfo>::updateMolStateAndOccupancy(MolState<osh::LO>& mol_state,
                                                            size_t index,
                                                            const osh::Real event_time) const;
extern template const std::vector<MolStateElementID>&
SurfaceReactionsBase<VDepInfo>::updateMolStateAndOccupancy(MolState<osh::LO>& mol_state,
                                                           size_t index,
                                                           const osh::Real event_time) const;
extern template const std::vector<MolStateElementID>&
SurfaceReactionsBase<SReacInfo>::updateMolStateAndOccupancy(MolState<osh::GO>& mol_state,
                                                            size_t index,
                                                            const osh::Real event_time) const;
extern template const std::vector<MolStateElementID>&
SurfaceReactionsBase<VDepInfo>::updateMolStateAndOccupancy(MolState<osh::GO>& mol_state,
                                                           size_t index,
                                                           const osh::Real event_time) const;
extern template const std::vector<MolStateElementID>&
SurfaceReactionsBase<GHKInfo>::updateMolStateAndOccupancy(MolState<osh::LO>& mol_state,
                                                          size_t index,
                                                          const osh::Real event_time) const;
extern template const std::vector<MolStateElementID>&
SurfaceReactionsBase<GHKInfo>::updateMolStateAndOccupancy(MolState<osh::GO>& mol_state,
                                                          size_t index,
                                                          const osh::Real event_time) const;
extern template void
SurfaceReactionsBase<VDepInfo>::report(std::ostream &report_stream,
                                       size_t index) const;
extern template void
SurfaceReactionsBase<SReacInfo>::report(std::ostream &report_stream,
                                        size_t index) const;
extern template void
SurfaceReactionsBase<GHKInfo>::report(std::ostream &report_stream,
                                      size_t index) const;

extern template SurfaceReactionsBase<SReacInfo>::SurfaceReactionsBase(const Statedef& statedef,
                                                                      DistMesh& mesh,
                                                                      MolState<osh::LO>& mol_state);
extern template SurfaceReactionsBase<SReacInfo>::SurfaceReactionsBase(const Statedef& statedef,
                                                                      DistMesh& mesh,
                                                                      MolState<osh::GO>& mol_state);
//-------------------------------------------------------

class SurfaceReactions : public SurfaceReactionsBase<SReacInfo> {
public:
  template <typename NumMolecules>
  SurfaceReactions(const Statedef& statedef, DistMesh& mesh, MolState<NumMolecules>& mol_state);

  inline constexpr KProcType getKProcType() const { return KProcType::SReac; }

  void updateKCst(const model::patch_id &patch_id,
                  container::surface_reaction_id id, osh::Real kCst,
                  DistMesh &mesh) {
    for (size_t k = 0; k < size(); k++) {
      if (patch_id == reacdefs_[k].get().patchdef().getID() &&
          reacdefs_[k].get().surfaceReactionID() == id) {
        ccsts_[k] = kinConstantGeomFactor(mesh, k) * kCst;
      }
    }
  }

private:
  void updateCcst(DistMesh &mesh);
};

//-------------------------------------------------------

extern template SurfaceReactions::SurfaceReactions(const Statedef& statedef,
                                                   DistMesh& mesh,
                                                   MolState<osh::LO>& mol_state);
extern template SurfaceReactions::SurfaceReactions(const Statedef& statedef,
                                                   DistMesh& mesh,
                                                   MolState<osh::GO>& mol_state);

class VDepSurfaceReactions : public SurfaceReactionsBase<VDepInfo> {
public:
  template <typename NumMolecules>
  VDepSurfaceReactions(const Statedef& statedef, DistMesh& mesh, MolState<NumMolecules>& mol_state);

  inline constexpr KProcType getKProcType() const {
    return KProcType::VDepSReac;
  }

  void kCstUpdate(osh::Reals &potential_on_verts);

private:
  std::vector<osh::Real> kinConstantGeomFactor_;
  const osh::LOs tri2verts_;
};

extern template VDepSurfaceReactions::VDepSurfaceReactions(const Statedef& statedef,
                                                           DistMesh& mesh,
                                                           MolState<osh::LO>& mol_state);
extern template VDepSurfaceReactions::VDepSurfaceReactions(const Statedef& statedef,
                                                           DistMesh& mesh,
                                                           MolState<osh::GO>& mol_state);

//-------------------------------------------------------

class GHKSurfaceReactions : public SurfaceReactionsBase<GHKInfo> {
public:
  template <typename NumMolecules>
  GHKSurfaceReactions(const Statedef& statedef, DistMesh& mesh, MolState<NumMolecules>& mol_state);

  inline constexpr KProcType getKProcType() const {
    return KProcType::GHKSReac;
  }
  /** Get currents
   *
   * /Note: This must be called after finalizeCurrents. Otherwise we get the charge flow not the
   * currents
   */
  inline osh::Reals currents() const noexcept {
      return (currents_);
  }

    const osh::Write<osh::GO> &getTri2Curr(const model::ghk_current_id &curr_id) const;

    void updatePotential(osh::Reals &potential_on_verts);

    /** Compute the GHK-voltage dependent reaction rate
     *
     * the concentrations are in n_molecules/m^3
     *
     * The GHK current is split in 2 reactions (as it is the case in general for
     * the rest of the code, not steps 3). Thus, the return value is a pair with
     * sign that depends on which direction the GHK current is flowing. Baudouin
     * verified the signs very carefully (since he solved a bug on these
     * currents where the sign was swapped).
     *
     * the rate is: permeability [m/s] * nuFoRT [1] * (conc_i - conc_o *
     * eNuFoRT) [n_molecules/m^3] * n_channels [1] = J [A/m^2]/(valence [1] *
     * Q_charge [C])
     *
     * This comes from a comparison with the GHK current formula on wikipedia.
     * In fact, confronting the 2 equations we have: rate = J [A/m^2] *
     * N_avogadro [mol] /(valence [1] * FARADAY [A/mol]).
     *
     * Since Q_charge [C] * N_avogadro [mol] = FARADAY we have the
     * aforementioned dimensional formula for the rate: J [A/m^2]/(valence [1] *
     * Q_charge [C])
     *
     * @tparam NumMolecules: Number of molecules. Strong_id
     * @param mol_state: number of molecules on the membrane (and neighbouring
     * compartments)
     * @param index: reaction index (to identify the reaction)
     * @return rate: J [A/m^2]/(valence [1] * Q_charge [C])
     */
    template <typename NumMolecules>
    osh::Real computeRate(const MolState<NumMolecules> &mol_state,
                          size_t index) const;

    void resetCurrents();

    /** Update currents_ with the net charge flow for the particular reaction
     *
     * /Note: when this function is called currents_ is actually the count of charge exchanged
     * during the time step. Only after finalizeCurrents we get the real currents
     *
     * /param reaction_index
     */
    void updateChargeFlow(size_t reaction_index);

    /**  Rescale currents by the period
     *
     * Before this function is called the currents are actually measured in
     * Coulomb and represent the total charge flow for a particular reaction. In formulae:
     *
     * i (before rescaling) = n_reactions_in_ef_dt * valence * Q_charge.
     *
     * where the valence is the net exchange of charges per reaction. It is usually the same as the
     * ion valence except if multiple ions are involved.
     *
     * @param period it is the efield_dt
     */
    void finalizeCurrents(double period) {
        const auto rescaleByThePeriod = OMEGA_H_LAMBDA(osh::LO reac_idx) {
            currents_[reac_idx] /= period;
        };
        osh::parallel_for(currents_.size(), rescaleByThePeriod);
    }

    // Number of specific ghk reactions per triangle (inward and outward)
    constexpr uint rpt() const { return 2; }

private:
    /// a mapping from a triangle id to its vertices ids
    const osh::LOs tri2verts_;
    /// a mapping from a reaction index to a ghk current
    osh::Write<osh::Real> currents_;
    /// a mapping from a reaction index to the inner element volume
    osh::Reals inner_element_vol_;
    /// a mapping from a reaction index to the outer element volume
    osh::Reals outer_element_vol_;
    /// a mapping from a reaction index to the potential on the triangle 
    osh::Write<osh::Real> potential_on_boundary_;

    /// a mapping from a ghk current identifier to a tri to reaction idx mapping
    /// the second mapping associates a local tri id to 2 reaction indices
    /// (corresponding to the inward and outward GHK reactions)
    // TODO use strong id for the value
    std::map<model::ghk_current_id, osh::Write<osh::GO>> curr2tri2reac_;
};

extern template GHKSurfaceReactions::GHKSurfaceReactions(const Statedef& statedef,
                                                         DistMesh& mesh,
                                                         MolState<osh::LO>& mol_state);
extern template GHKSurfaceReactions::GHKSurfaceReactions(const Statedef& statedef,
                                                         DistMesh& mesh,
                                                         MolState<osh::GO>& mol_state);


} // namespace kproc
} // namespace dist
} // namespace steps
