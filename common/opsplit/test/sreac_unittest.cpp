#include <hadoken/format/format.hpp>

#include "opsplit/compdef.hpp"
#include "opsplit/diffdef.hpp"
#include "opsplit/mesh.hpp"
#include "opsplit/patchdef.hpp"
#include "opsplit/reacdef.hpp"
#include "opsplit/sreacdef.hpp"
#include "opsplit/test/simdef.hpp"
#include "opsplit/test/simulation.hpp"
#include "sreac_unittest.hpp"

namespace zee {

using hadoken::scat;

SReacUnitTest::SReacUnitTest(const ScenarioInput& t_input)
    : Scenario("SReacUnitTest", "Unit test for surface reactions", t_input) {}

std::unique_ptr<Statedef> SReacUnitTest::createStatedef() const {
    SReacUnitTestSimdef simdef;
    return std::move(simdef.getStatedef());
}

void SReacUnitTest::register_compartments(DistMesh& mesh) const {
    mesh.addComp("comp_i", model::compartment_label(1));
    mesh.addComp("comp_o", model::compartment_label(2));
}

void SReacUnitTest::fill_compartments(simulation_t& simulation) const {
    SReacUnitTestSimdef simdef;
    for (const auto& compartment: simdef.getCompartementCounts()) {
        simulation.setCompCount(compartment.compartment,
                                compartment.specie,
                                compartment.num_mols * input.num_mols_factor);
    }
    for (const auto& initializer: simdef.getPatchCounts()) {
        simulation.setPatchCount(initializer.patch,
                                 initializer.specie,
                                 initializer.num_mols * input.num_mols_factor);
    }
}

void SReacUnitTest::run_simulation_impl(simulation_t& simulation) {
    simulation.run(input.end_time);
}

int SReacUnitTest::check_and_log_results_impl(const simulation_t& simulation) const {
    PetscScalar A = simulation.getCompCount("comp_i", "A");
    PetscScalar B = simulation.getCompCount("comp_i", "B");
    PetscScalar C = simulation.getCompCount("comp_o", "C");
    PetscScalar D = simulation.getPatchCount("patch", "D");
    simulation.log_once(scat("A: ", A, ", B: ", B, ", C: ", C, ", D: ", D));
    bool status = EXIT_SUCCESS;
    if (static_cast<PetscInt>(A) != 0) {
        simulation.log_once(scat("Error: expected concentration of A to be 0 but got ", A));
        status = EXIT_FAILURE;
    }
    if (static_cast<PetscInt>(B) != 0) {
        simulation.log_once(scat("Error: expected concentration of B to be 0 but got ", B));
        status = EXIT_FAILURE;
    }
    if (static_cast<PetscInt>(C) != 0) {
        simulation.log_once(scat("Error: expected concentration of C to be 0 but got ", C));
        status = EXIT_FAILURE;
    }
    if (static_cast<PetscInt>(D) != 3000) {
        simulation.log_once(scat("Error: expected 3000 mols D but got ", D));
        status = EXIT_FAILURE;
    }
    return status;
}

}  // namespace zee
