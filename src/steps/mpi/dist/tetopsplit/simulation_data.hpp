#pragma once

#include <fstream>

#include "kproc/diffusions.hpp"
#include "kproc/kproc_state.hpp"
#include "operator/diffusion_operator.hpp"
#include "operator/fwd.hpp"
#if USE_PETSC
#include "operator/efield_operator.hpp"
#endif  // USE_PETSC

namespace steps::dist {

/**
 * Input data required by the diffusion simulation
 */
class SimulationInput {
  public:
    SimulationInput(const osh::LOs& t_species_per_elements,
                    osh::LO t_num_iterations,
                    rng::RNG& t_rng)
        : pools(t_species_per_elements)
        , num_iterations(t_num_iterations)
        , molecules_leaving(t_rng)
        , species_per_element(t_species_per_elements) {}

    SimulationInput(const osh::LOs& t_species_per_owned_elements,
                    const std::optional<osh::LOs>& t_species_per_owned_element_boundaries,
                    const osh::LOs& t_species_per_element,
                    const osh::LO t_num_iterations,
                    rng::RNG& t_rng,
                    osh::LO num_vertices)
        : pools(t_species_per_owned_elements, true, t_species_per_owned_element_boundaries)
        , num_iterations(t_num_iterations)
        , molecules_leaving(t_rng)
        , species_per_element(t_species_per_element)
        , potential_on_vertices_w(num_vertices, DEFAULT_MEMB_POT)
        , current_on_vertices_w(num_vertices, 0) {}

    /// number of molecules per species per triangle/tetrahedron
    MolState pools;
    /// number of diffusion iterations
    const osh::LO num_iterations;
    /// functor to determine number of molecules leaving an element
    /// (triangle/tetrahedron)
    const kproc::LeavingMolecules molecules_leaving;
    /// Species per all elements (owned or not)
    const osh::LOs species_per_element;
    /// potential on vertices for E-Field, r-w
    osh::Write<osh::Real> potential_on_vertices_w;
    /// current on vertices for E-Field, r-w
    osh::Write<osh::Real> current_on_vertices_w;
};

namespace {

template <SSAMethod>
struct ssa_traits {};

template <>
struct ssa_traits<SSAMethod::SSA> {
    template <NextEventSearchMethod SearchMethod>
    using ssa_operator_type = SSAOperator<SearchMethod>;
};

template <>
struct ssa_traits<SSAMethod::RSSA> {
    template <NextEventSearchMethod /* SearchMethod */>
    using ssa_operator_type = RSSAOperator;
};

}  // namespace

/**
 * Internal data used by simulation
 */
template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
class SimulationData {
  public:
    using ssa_operator_type = typename ssa_traits<SSA>::template ssa_operator_type<SearchMethod>;

    SimulationData(DistMesh& mesh,
                   const Statedef& statedef,
                   SimulationInput& input,
                   rng::RNG& t_rng,
                   bool indepKProcs)
        : pools(input.pools)
        , diffusions(mesh, input)
        , kproc_state(statedef, mesh, pools, indepKProcs)
        , ssaOp(pools, kproc_state, t_rng, osh::Reals(input.potential_on_vertices_w))
        , diffOp(mesh, t_rng, pools, diffusions, kproc_state) {
#if USE_PETSC
        if (statedef.is_efield_enabled()) {
            efield.emplace(mesh, statedef, kproc_state.ghkCurrentsBoundaries(), pools);
        }
#endif  // USE_PETSC
    }

    SimulationData(const SimulationData&) = delete;

    MolState& pools;
    osh::Real time_delta{};
    kproc::Diffusions diffusions;
    kproc::KProcState kproc_state;
    ssa_operator_type ssaOp;
    DiffusionOperator diffOp;
#if USE_PETSC
    std::optional<EFieldOperator> efield;
#endif  // USE_PETSC
    bool active_diffusions{};

    void reset(const osh::Real state_time) {
        diffusions.reset();
        pools.reset(state_time);
        ssaOp.reset();
        kproc_state.resetCurrents();
    }

    osh::Real updateIterationTimeStep() {
        auto global_max_sums = diffusions.global_rates_max_sum();
        active_diffusions = (global_max_sums > std::numeric_limits<osh::Real>::epsilon());

        // if no diffusions are present. time delta can be set to infinity.
        time_delta = active_diffusions ? 1.0 / global_max_sums
                                       : std::numeric_limits<osh::Real>::infinity();

        return time_delta;
    }
};

}  // namespace steps::dist
