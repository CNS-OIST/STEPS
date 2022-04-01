#include "ca_burst_integration_test.hpp"

#include <iostream>

#include "geom/dist/distmesh.hpp"
#include "math/constants.hpp"
#include "model/model.hpp"
#include "mpi/dist/test/simulation.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "mpi/dist/tetopsplit/kproc/surface_reactions.hpp"

namespace steps {
namespace dist {

CaBurstIntegrationTest::CaBurstIntegrationTest(const ScenarioInput& t_input)
    : Scenario("CaBurstIntegrationTest", "Unit test for GHK current and E-field", t_input) {}

std::unique_ptr<Statedef> CaBurstIntegrationTest::createStatedef(
    const simulation_t& simulation) const {
    steps::model::Model model;
    CaBurstIntegrationTestSimdef simdef(model, simulation.getMesh());
    return std::move(simdef.getStatedef());
}

void CaBurstIntegrationTest::register_compartments(DistMesh& mesh) const {
    mesh.addComp("comp1", model::compartment_label(1));
    mesh.addComp("comp2", model::compartment_label(2));
}

void CaBurstIntegrationTest::fill_compartments(simulation_t& simulation) const {
    steps::model::Model model;
    CaBurstIntegrationTestSimdef simdef(model, simulation.getMesh());
    simulation.setCompConc("comp1", {{"C", CaBurstIntegrationTestSimdef::compConc()}});
    simulation.setCompConc("comp2", {{"C", CaBurstIntegrationTestSimdef::compConc() / 2}});
}

void CaBurstIntegrationTest::fill_patches(simulation_t& simulation) const {
    steps::model::Model model;
    CaBurstIntegrationTestSimdef simdef(model, simulation.getMesh());
    for (const auto& initializer: simdef.getPatchCounts()) {
        simulation.setPatchCount(initializer.patch,
                                 initializer.species,
                                 initializer.num_mols * input.num_mols_factor);
    }
}

void CaBurstIntegrationTest::run_simulation_impl(simulation_t& simulation) {
    const osh::Real time_step = 1e-4;

    const auto n_time_steps = static_cast<size_t>(input.end_time / time_step);

    C_count.reserve(n_time_steps);
    A_count.reserve(n_time_steps);
    B_count.reserve(n_time_steps);
    C_surf_count.reserve(n_time_steps);
    V_min.reserve(n_time_steps);
    V_max.reserve(n_time_steps);
    i_ghk.reserve(n_time_steps);
    t.reserve(n_time_steps);

    simulation.reset();
    simulation.setPotential(CaBurstIntegrationTestSimdef::potential());
    fill_compartments(simulation);
    fill_patches(simulation);

    for (size_t k = 0; k <= n_time_steps; ++k) {
        simulation.run(k * time_step);
        C_count.emplace_back(simulation.getCompCount("comp1", "C"));
        A_count.emplace_back(simulation.getCompCount("comp1", "A"));
        B_count.emplace_back(simulation.getCompCount("comp1", "B"));
        C_surf_count.emplace_back(simulation.getPatchCount("patch1", "C"));
        V_min.emplace_back(simulation.getMinPotentialOnVertices("patch1"));
        V_max.emplace_back(simulation.getMaxPotentialOnVertices("patch1"));
        i_ghk.emplace_back(simulation.getTotalGHKCurrent());

        t.emplace_back(k * time_step);
        // uncomment this to get progress
        // simulation.log_progress(k, n_time_steps);
    }
}

int CaBurstIntegrationTest::check_and_log_results_impl(simulation_t& simulation) const {
    // uncomment this to get the traces
    simulation.log_once(
        "step, time_step, C_count, A_count, B_count, Csurf_count, V_min, V_max, i_ghk_1\n");
    for (size_t i = 0; i < t.size(); ++i) {
        std::stringstream ss;
        ss.precision(13);
        ss << i << ' ' << t[i] << ' ' << C_count[i] << ' ' << A_count[i] << ' ' << B_count[i] << ' '
           << C_surf_count[i] << ' ' << V_min[i] << ' ' << V_max[i] << ' ' << i_ghk[i] << std::endl;
        simulation.log_once(ss.str());
    }

    // This test actually does not test anything in particular. It verifies the interactions of the
    // various mechanisms that appear in ca_burst. These values are something that I have seen it
    // converges too
    const double expected_v_max = -0.0639663;
    if (std::abs(V_max.back() - expected_v_max) > 1e-4) {
        std::stringstream ss;
        ss << "The interactions of the mechanisms have changed. Also the order in which the operations are performed "
              "can change this test. This test is mostly to see if it runs up until this point. Further "
              "investigation is suggested. New final V_max: "
           << V_max.back()
           << " "
              "different from the expected value: "
           << expected_v_max;
        simulation.log_once(ss.str());
        return 1;
    }

    return 0;
}

}  // namespace dist
}  // namespace steps
