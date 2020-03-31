#include <map>

#include <hadoken/format/format.hpp>

#include "mpitools.hpp"
#include "multi_comp.hpp"
#include "opsplit/compdef.hpp"
#include "opsplit/diffdef.hpp"
#include "opsplit/mesh.hpp"
#include "opsplit/patchdef.hpp"
#include "opsplit/reacdef.hpp"
#include "opsplit/sreacdef.hpp"
#include "opsplit/test/simdef.hpp"
#include "opsplit/test/simulation.hpp"

namespace zee {

using hadoken::scat;

MultipleCompartment::MultipleCompartment(const ScenarioInput& t_input)
    : Scenario("MultipleCompartment", "Simulation occurs in 2 compartments", t_input) {}

std::unique_ptr<Statedef> MultipleCompartment::createStatedef() const {
    MultiCompartmentSimdef simdef;
    return std::move(simdef.getStatedef());
}

void MultipleCompartment::register_compartments(DistMesh& mesh) const {
    mesh.addComp("Left", model::compartment_label(1));
    mesh.addComp("Right", model::compartment_label(2));
}

void MultipleCompartment::fill_compartments(simulation_t& simulation) const {
    MultiCompartmentSimdef simdef;
    for (const auto& compartment_count: simdef.getCompartementCounts()) {
        simulation.setCompCount(compartment_count.compartment,
                                compartment_count.specie,
                                compartment_count.num_mols * input.num_mols_factor);
    }
}

void MultipleCompartment::log_num_molecules(const simulation_t& simulation) const {
    PetscScalar A_Left = simulation.getCompCount("Left", "A");
    PetscScalar B_Left = simulation.getCompCount("Left", "B");
    PetscScalar C_Left = simulation.getCompCount("Left", "C");
    PetscScalar C_Right = simulation.getCompCount("Right", "C");
    PetscScalar D_Right = simulation.getCompCount("Right", "D");
    PetscScalar E_Right = simulation.getCompCount("Right", "E");

    simulation.log_once(scat("Comp Left has ", PetscRealPart(A_Left), " of spec A"));
    simulation.log_once(scat("Comp Left has ", PetscRealPart(B_Left), "of spec B"));
    simulation.log_once(scat("Comp Left has ", PetscRealPart(C_Left), "of spec C"));
    simulation.log_once(scat("Comp Right has ", PetscRealPart(C_Right), "of spec C"));
    simulation.log_once(scat("Comp Right has ", PetscRealPart(D_Right), "of spec D"));
    simulation.log_once(scat("Comp Right has ", PetscRealPart(E_Right), "of spec E"));
}

void MultipleCompartment::run_simulation_impl(simulation_t& simulation) {
    PetscScalar interval = 0.00001;
    simulation.log_once(scat("Binomial Threshold: ", input.threshold));
    log_num_molecules(simulation);
    for (auto i = 1; i < 100; i++) {
        const std::string output_name = "conc_" + std::to_string(i) + ".vtk";
        simulation.exportMolStateToVTK(output_name);
        simulation.run(interval * i);
        log_num_molecules(simulation);
        simulation.log_once(scat("Number of iterations: ", simulation.getNIterations()));
    }
}

}  // namespace zee
