#include <map>

#include <hadoken/format/format.hpp>

#include "opsplit/compdef.hpp"
#include "opsplit/diffdef.hpp"
#include "opsplit/mesh.hpp"
#include "opsplit/patchdef.hpp"
#include "opsplit/reacdef.hpp"
#include "opsplit/sreacdef.hpp"
#include "opsplit/test/simdef.hpp"
#include "opsplit/test/simulation.hpp"
#include "simple_model.hpp"

namespace zee {

using hadoken::scat;

SimpleModel2::SimpleModel2(const ScenarioInput& t_input)
    : Scenario("SimpleModel", "Simple container in Chen_FNeuroinf_2017", t_input) {}

std::unique_ptr<Statedef> SimpleModel2::createStatedef() const {
    SimpleSimdef simdef;
    return std::move(simdef.getStatedef());
}

void SimpleModel2::register_compartments(DistMesh& mesh) const {
    mesh.addComp("comp", input.comp_label);
}

void SimpleModel2::fill_compartments(simulation_t& simulation) const {
    SimpleSimdef simdef;
    for (const auto& compartment: simdef.getCompartementCounts()) {
        simulation.setCompCount(compartment.compartment,
                                compartment.specie,
                                compartment.num_mols * input.num_mols_factor);
    }
}

void SimpleModel2::run_simulation_impl(simulation_t& simulation) {
    simulation.run(input.end_time);
}

int SimpleModel2::check_and_log_results_impl(const simulation_t& simulation) const {
    PetscScalar A = simulation.getCompCount("comp", "A");
    PetscScalar B = simulation.getCompCount("comp", "B");
    PetscScalar C = simulation.getCompCount("comp", "C");
    PetscScalar D = simulation.getCompCount("comp", "D");
    PetscScalar E = simulation.getCompCount("comp", "E");
    PetscScalar F = simulation.getCompCount("comp", "F");
    PetscScalar G = simulation.getCompCount("comp", "G");
    PetscScalar H = simulation.getCompCount("comp", "H");
    PetscScalar I = simulation.getCompCount("comp", "I");
    PetscScalar J = simulation.getCompCount("comp", "J");
    simulation.log_once(scat("A: ",
                             A,
                             ", B: ",
                             B,
                             ", C: ",
                             C,
                             ", D: ",
                             D,
                             ", E: ",
                             E,
                             ", F: ",
                             F,
                             ", G: ",
                             G,
                             ", H: ",
                             H,
                             ", I: ",
                             I,
                             ", J: ",
                             J));
    return EXIT_SUCCESS;
}

}  // namespace zee
