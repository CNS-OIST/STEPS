#pragma once

#include <fstream>

#include "geom/dist/distmesh.hpp"
#include "kproc/diffusions.hpp"
#include "kproc/kproc_state.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"
#include "operator/diffusion_operator.hpp"
#include "operator/fwd.hpp"
#include "operator/rleaping_operator.hpp"
#include "util/common.hpp"
#if USE_PETSC
#include "operator/efield_operator.hpp"
#endif  // USE_PETSC

namespace steps::dist {

/**
 * Input data required by the diffusion simulation
 */
class SimulationInput {
  public:
    SimulationInput(const osh::LOs& t_species_per_owned_elements,
                    const std::optional<osh::LOs>& t_species_per_owned_element_boundaries,
                    const osh::LOs& t_species_per_element,
                    const osh::LOs& substates_per_complexes,
                    const osh::LO t_num_iterations,
                    rng::RNG& t_rng,
                    osh::LO num_vertices,
                    const DistMesh& mesh,
                    const Statedef& statedef_)
        : pools(mesh,
                statedef_,
                t_species_per_owned_elements,
                substates_per_complexes,
                {mesh.comm_rank(), mesh.comm_size()},
                true,
                t_species_per_owned_element_boundaries)
        , num_iterations(t_num_iterations)
        , molecules_leaving(t_rng)
        , species_per_element(t_species_per_element)
        , potential_on_vertices_w(num_vertices, DEFAULT_MEMB_POT)
        , current_on_vertices_w(num_vertices, 0)
        , current_on_triangles_w(mesh.owned_bounds_mask().size(), 0)
        , capacitance_on_triangles_w(t_species_per_owned_element_boundaries.has_value()
                                         ? t_species_per_owned_element_boundaries->size()
                                         : 0,
                                     0)
        , conductivity_on_triangles_w(t_species_per_owned_element_boundaries.has_value()
                                          ? t_species_per_owned_element_boundaries->size()
                                          : 0,
                                      0)
        , reversal_potential_on_triangles_w(t_species_per_owned_element_boundaries.has_value()
                                                ? t_species_per_owned_element_boundaries->size()
                                                : 0,
                                            0) {}

    void reset(const Statedef& statedef, DistMesh& mesh) const {
        std::fill(potential_on_vertices_w.begin(), potential_on_vertices_w.end(), DEFAULT_MEMB_POT);
        std::fill(current_on_vertices_w.begin(), current_on_vertices_w.end(), 0);

        for (auto& memb: statedef.membranes()) {
            auto capac = memb->capacitance();
            for (const Patchdef& patchdef: memb->getPatchdefs()) {
                for (const auto tri: patchdef.patch().getTris(true)) {
                    capacitance_on_triangles_w[tri.get()] = capac;
                }
            }
        }
    }

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
    /// current on triangles for E-Field, r-w
    osh::Write<osh::Real> current_on_triangles_w;
    /// capacitance on triangles for E-Field, r-w
    osh::Write<osh::Real> capacitance_on_triangles_w;
    /// conductance on triangles for E-Field, r-w
    osh::Write<osh::Real> conductivity_on_triangles_w;
    /// reversal potential on triangles for E-Field, r-w
    osh::Write<osh::Real> reversal_potential_on_triangles_w;
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
struct ssa_traits<SSAMethod::RLeaping> {
    template <NextEventSearchMethod /* SearchMethod */>
    using ssa_operator_type = RLeapingOperator;
};

template <>
struct ssa_traits<SSAMethod::RSSA> {
    template <NextEventSearchMethod /* SearchMethod */>
    using ssa_operator_type = RSSAOperator;
};

template <DiffusionMethod>
struct diff_traits {};

template <>
struct diff_traits<DiffusionMethod::ConstantDiffDt> {
    using diff_operator_type = DiffusionOperator;
};

template <>
struct diff_traits<DiffusionMethod::TauLeapingDiffDt> {
    using diff_operator_type = TauLeapingDiffusionOperator;
};

}  // namespace

/**
 * Internal data used by simulation
 */
template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
class SimulationData {
  public:
    using ssa_operator_type = typename ssa_traits<SSA>::template ssa_operator_type<SearchMethod>;
    using diff_operator_type = typename diff_traits<DiffMethod>::diff_operator_type;

    SimulationData(DistMesh& mesh,
                   const Statedef& statedef,
                   SimulationInput& input,
                   rng::RNG& t_rng,
                   bool indepKProcs)
        : pools(input.pools)
        , diffusions(mesh, statedef, input)
        , kproc_state(statedef, mesh, pools, indepKProcs)
        , ssaOp(pools, kproc_state, t_rng)
        , diffOp(mesh, t_rng, pools, diffusions, kproc_state) {
#if USE_PETSC
        if (statedef.is_efield_enabled()) {
            efield.emplace(mesh, statedef, pools);
        }
#endif  // USE_PETSC
        pools.finalize_complex_occupancy();
        initialize_diffusions();
    }

    SimulationData(const SimulationData&) = delete;

    MolState& pools;
    kproc::Diffusions diffusions;
    kproc::KProcState kproc_state;
    ssa_operator_type ssaOp;
    diff_operator_type diffOp;
#if USE_PETSC
    std::optional<EFieldOperator> efield;
#endif  // USE_PETSC

    void reset(const osh::Real state_time) {
        diffusions.reset();
        pools.reset(state_time);
        ssaOp.reset();
        diffOp.reset();
        kproc_state.reset();
        initialize_diffusions();
#if USE_PETSC
        // efield->resetStiffnessMatrix();
#endif  // USE_PETSC
    }

    void initialize_diffusions() {
        diffusions.initialize_discretized_rates();
        diffOp.initialize();
    }
};

}  // namespace steps::dist
