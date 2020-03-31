#include <map>
#include <vector>

#include <hadoken/format/format.hpp>

#include "diff_sim.hpp"
#include "opsplit/compdef.hpp"
#include "opsplit/diffdef.hpp"
#include "opsplit/mesh.hpp"
#include "opsplit/patchdef.hpp"
#include "opsplit/reacdef.hpp"
#include "opsplit/sreacdef.hpp"
#include "opsplit/statedef.hpp"

#include "opsplit/test/simulation.hpp"

namespace zee {

using hadoken::scat;

DiffusionOnly::DiffusionOnly(const ScenarioInput& t_input)
    : Scenario("diffusionOnly", "Only diffuse one specie", t_input) {}

// FIXME(TCL)
/*
void DiffusionOnly::parse_options() {
    PetscOptionsGetInt(nullptr, nullptr, "-comp_tag", &comp_tag, nullptr);
    PetscOptionsGetScalar(nullptr, nullptr, "-t", &end_time, nullptr);
    PetscOptionsGetScalar(nullptr, nullptr, "-i", &interval, nullptr);
    PetscOptionsGetScalar(nullptr, nullptr, "-d", &dcst, nullptr);
    PetscOptionsGetScalar(nullptr, nullptr, "-n", &NMols, nullptr);
}
*/

std::unique_ptr<Statedef> DiffusionOnly::createStatedef() const {
    auto statedef = std::make_unique<Statedef>();
    statedef->addComp("comp");
    std::vector<model::specie_name> species = {"A"};
    statedef->addCompSpecs("comp", species);

    statedef->addCompDiff("comp", "A", dcst);

    return statedef;
}

void DiffusionOnly::register_compartments(DistMesh& mesh) const {
    mesh.addComp("comp", comp_tag);
}

void DiffusionOnly::fill_compartments(simulation_t& simulation) const {
    if (input.do_inject_in_elt == -1) {
        simulation.setCompCount("comp", "A", NMols);
    } else {
        simulation.setOwnedElementCount("comp", input.do_inject_in_elt, "A", NMols);
    }
}

void DiffusionOnly::run_simulation_impl(simulation_t& simulation) {
    int iter = 0;

    for (PetscScalar t = 0.0; t < input.end_time; t += input.do_interval) {  // NOLINT
        simulation.log_once(scat("time_step=", t));
        simulation.exportMolStateToVTK(scat("count_", iter, ".vtk"));
        simulation.run(t);
        iter++;
    }
}

}  // namespace zee
