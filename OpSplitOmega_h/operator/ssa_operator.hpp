#pragma once

#include <random>

#include "../common.hpp"

#include "../kproc/kproc_state.hpp"
#include "../mesh.hpp"

namespace zee {

/**
 * \brief SSA operator.
 *
 * This class implements the SSA operator in the OpSplit solution.
 *
 * The core function of this class is run(), which partially corresponds to the SSA component
 * the simulation core loop tetexact::Tetexact::run() / tetopsplit::TetOpSplit::run() in STEPS.
 */
template <typename RNG>
class SSAOperator {
  public:
    template <osh::Int Dim>
    SSAOperator(MolState& mol_state,
                DiffusionVariables& dv,
                KProcState& kproc_state,
                const OmegaHMesh<Dim>& mesh,
                RNG& t_rng);

    PetscInt64 getExtent() const;

    void run(PetscScalar period);

  private:
    MolState& pMolState;
    DiffusionVariables& pDiffusionVariables;
    KProcState& pKProcState;
    RNG& rng_;

    const mesh::element_ids& owned_elements_;

    PetscInt64 extent{};
};

// explicit instantiation declarations
extern template class SSAOperator<std::mt19937>;
extern template SSAOperator<std::mt19937>::SSAOperator(MolState& mol_state,
                                                       DiffusionVariables& dv,
                                                       KProcState& kproc_state,
                                                       const OmegaHMesh<2>& mesh,
                                                       std::mt19937& t_rng);
extern template SSAOperator<std::mt19937>::SSAOperator(MolState& mol_state,
                                                       DiffusionVariables& dv,
                                                       KProcState& kproc_state,
                                                       const OmegaHMesh<3>& mesh,
                                                       std::mt19937& t_rng);

}  // namespace zee
