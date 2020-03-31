#pragma once

#include "common.hpp"

#include "kproc/diffusions.hpp"
#include "kproc/kproc_state.hpp"
#include "operator/diffusion_operator.hpp"

namespace zee {

/**
 * Input data required by the diffusion simulation
 */
template <typename RNG>
class SimulationInput {
  public:
    SimulationInput(const osh::LOs& t_species_per_elements, osh::LO t_num_iterations, RNG& t_rng)
        : pools(t_species_per_elements)
        , num_iterations(t_num_iterations)
        , molecules_leaving(t_rng)
        , species_per_element(t_species_per_elements) {}

    SimulationInput(const osh::LOs& t_species_per_owned_elements,
                    const boost::optional<osh::LOs>& t_species_per_owned_element_boundaries,
                    const osh::LOs& t_species_per_element,
                    const osh::LO t_num_iterations,
                    RNG& t_rng)
        : pools(t_species_per_owned_elements, t_species_per_owned_element_boundaries)
        , num_iterations(t_num_iterations)
        , molecules_leaving(t_rng)
        , species_per_element(t_species_per_element) {}

    /// number of molecules per specie per triangle/tetrahedron
    MolState pools;
    /// number of diffusion iterations
    const osh::LO num_iterations;
    /// functor to determine number of molecules leaving an element (triangle/tetrahedron)
    const LeavingMolecules<RNG> molecules_leaving;
    /// Species per all elements (owned or not)
    const osh::LOs species_per_element;
};


namespace {

template <SSAMethod>
struct ssa_traits {};

template <>
struct ssa_traits<SSAMethod::SSA> {
    template <typename RNG>
    using ssa_operator_type = SSAOperator<RNG>;
};

template <>
struct ssa_traits<SSAMethod::RSSA> {
    template <typename RNG>
    using ssa_operator_type = RSSAOperator<RNG>;
};

}  // namespace

/**
 * Internal data used by simulation
 */
template <osh::Int Dim, SSAMethod SSA, typename RNG>
class SimulationData {
  public:
    using ssa_operator_type = typename ssa_traits<SSA>::template ssa_operator_type<RNG>;

    SimulationData(OmegaHMesh<Dim>& mesh,
                   const Statedef& statedef,
                   SimulationInput<RNG>& input,
                   bool molecules_pools_force_dist_for_variable_sized,
                   RNG& t_rng)
        : pools(input.pools)
        , diffusions(mesh, input, molecules_pools_force_dist_for_variable_sized)
        , kproc_state(statedef, mesh)
        , ssaOp(pools, diffusions.dv(), kproc_state, mesh, t_rng)
        , diffOp(mesh, t_rng, pools, diffusions, time_delta) {}

    SimulationData(const SimulationData&) = delete;

    MolState& pools;
    osh::Real time_delta{};
    Diffusions<Dim, RNG> diffusions;
    KProcState kproc_state;
    ssa_operator_type ssaOp;
    DiffusionOperator<Dim, RNG> diffOp;

    void reset() {
        diffusions.reset();
    }
    osh::Real updateIterationTimeStep() {
        time_delta = 1. / diffusions.rates().rates_max_sum();
        assert(time_delta > 0.0);
        return time_delta;
    }
};

}  // namespace zee
