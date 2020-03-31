#pragma once

#include <random>

#include <boost/optional/optional.hpp>

#include "../common.hpp"
#include "../kproc/kproc_state.hpp"
#include "../kproc/reactions.hpp"
#include "../mesh.hpp"
#include "collection.hpp"

namespace zee {

/**
 * \brief RSSA operator.
 *
 * This class implements the RSSA operator in the OpSplit solution.
 *
 * The core function of this class is run(), which partially corresponds to the SSA component
 * the simulation core loop tetexact::Tetexact::run() / tetopsplit::TetOpSplit::run() in STEPS.
 */

template <typename RNG>
class RSSAOperator {
  public:
    /**
     * Constructor
     * \param mol_state molecular state
     * \param dv molecules transferred to neighbours
     * \param kproc_state kinetic processes state
     * \param mesh distributed mesh
     */
    template <osh::Int Dim>
    RSSAOperator(MolState& mol_state,
                 DiffusionVariables& dv,
                 KProcState& kproc_state,
                 const OmegaHMesh<Dim>& mesh,
                 RNG& t_rng);


    PetscInt64 getExtent() const;


    void run(PetscScalar period);

  private:
    /**
     * \brief Check all reaction rates bounds and update those as necessary.
     *
     * \param mol_state molecular state
     */
    void checkAndUpdateReactionRatesBounds(const MolState& mol_state);
    /**
     * \brief Check and update the reaction rates as necessary after a reaction that happened in
     * elementId.
     *
     * \param mol_state molecular state
     * \param elementId element identifier
     * \param reaction_indices Reaction stoichiometric coefficients.
     */
    void checkAndUpdateReactionRatesBounds(
        const MolState& mol_state,
        const mesh::element_id elementId,
        const boost::optional<std::vector<PetscInt>>& reaction_indices = {});
    /**
     * \brief Bounds generator for molecules.
     *
     * Bounds are defined based on the static member delta_rel_ as (nc * (1 - delta_rel_), nc * (1 +
     * delta_rel_)) where nc is the number of molecules. A special treatment is applied for a small
     * number of molecules.
     * \param nc number of molecules
     * \param pPoolLB Lower bound on the number of molecules (updated by the procedure).
     * \param pPoolUB Upper bound on the number of molecules (updated by the procedure).
     */
    static inline void applyBounds(osh::LO nc, osh::LO& pPoolLB, osh::LO& pPoolUB);
    /**
     * \brief Simulate a candidate reaction
     */
    kproc::Reactions::index_type singleOutReaction();

    MolState& pMolState;
    DiffusionVariables& pDiffusionVariables;
    KProcState& pKProcState;

    const mesh::element_ids& owned_elements_;

    // bounds on molecule populations
    MolState mol_state_lower_bound_;
    MolState mol_state_upper_bound_;

    // dependent reactions
    std::map<std::pair<mesh::element_id, container::specie_id>,
             std::set<kproc::Reactions::index_type>>
        dependent_reactions_;

    // reaction propensity rates
    std::vector<PetscScalar> a_lower_bound_;
    VectorSumUpdater a_upper_bound_;

    // uniform distribution
    std::uniform_real_distribution<double> uniform_;

    RNG& rng_;

    // counter
    PetscInt64 extent{};

    // approximate relative distance of the molecules lower/upper bounds from the actual molecules
    // count.
    static constexpr PetscScalar delta_rel_ = 0.05;
};


// explicit instantiation declarations
extern template class RSSAOperator<std::mt19937>;
extern template RSSAOperator<std::mt19937>::RSSAOperator(MolState& mol_state,
                                                         DiffusionVariables& dv,
                                                         KProcState& kproc_state,
                                                         const OmegaHMesh<2>& mesh,
                                                         std::mt19937& t_rng);
extern template RSSAOperator<std::mt19937>::RSSAOperator(MolState& mol_state,
                                                         DiffusionVariables& dv,
                                                         KProcState& kproc_state,
                                                         const OmegaHMesh<3>& mesh,
                                                         std::mt19937& t_rng);


}  // namespace zee
