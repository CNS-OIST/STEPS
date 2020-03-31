#include "ssa_operator.hpp"

#include "../kproc/diffusions.hpp"
#include "mpitools.hpp"

#undef MPI_Allreduce

namespace zee {

template <typename RNG>
template <osh::Int Dim>
SSAOperator<RNG>::SSAOperator(MolState& mol_state,
                              DiffusionVariables& dv,
                              KProcState& kproc_state,
                              const OmegaHMesh<Dim>& mesh,
                              RNG& t_rng)
    : pMolState(mol_state)
    , pDiffusionVariables(dv)
    , pKProcState(kproc_state)
    , rng_(t_rng)
    , owned_elements_(mesh.getOwnedElems()) {}

template <typename RNG>
PetscInt64 SSAOperator<RNG>::getExtent() const {
    PetscInt64 global_ext;
    auto ierr = MPI_Allreduce(&extent, &global_ext, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    if (ierr != 0) {
        MPI_Abort(MPI_COMM_WORLD, ierr);
    }
    return global_ext;
}

template <typename RNG>
void SSAOperator<RNG>::run(PetscScalar period) {
    PetscScalar cumulative_dt = 0.0;
    pKProcState.updateAllPropensities(pMolState);
    while (true) {
        PetscScalar a0 = pKProcState.getSSA_A0();
        if (a0 <= 0.0) {
            break;
        }
        std::exponential_distribution<double> exponential(a0);
        cumulative_dt += exponential(rng_);
        if (cumulative_dt > period) {
            break;
        }
        auto kinetic_process_id = pKProcState.getNextEvent<RNG>(rng_);
        // TODO(BDM): improve occupancy computation.
        pKProcState.updateOccupancy(pDiffusionVariables,
                                    pMolState,
                                    cumulative_dt,
                                    kinetic_process_id);
        pKProcState.updateMolState(pMolState, kinetic_process_id);
        pKProcState.updatePropensities(pMolState, kinetic_process_id);
        extent += 1;
    }

    osh::for_each(owned_elements_.data().begin(), owned_elements_.data().end(), [&](osh::LO e) {
        mesh::element_id elem(e);
        for (container::specie_id s(0); s < pMolState.numSpecies(elem); ++s) {
            const auto update_time_interval = period -
                                              pDiffusionVariables.last_update_time(elem, s);
            pDiffusionVariables.occupancy(elem, s) +=
                update_time_interval * static_cast<Omega_h::Real>(pMolState(elem, s));
        }
    });
}

// explicit instantiation definitions
template class SSAOperator<std::mt19937>;
template SSAOperator<std::mt19937>::SSAOperator(MolState& mol_state,
                                                DiffusionVariables& dv,
                                                KProcState& kproc_state,
                                                const OmegaHMesh<2>& mesh,
                                                std::mt19937& t_rng);
template SSAOperator<std::mt19937>::SSAOperator(MolState& mol_state,
                                                DiffusionVariables& dv,
                                                KProcState& kproc_state,
                                                const OmegaHMesh<3>& mesh,
                                                std::mt19937& t_rng);

}  // namespace zee
