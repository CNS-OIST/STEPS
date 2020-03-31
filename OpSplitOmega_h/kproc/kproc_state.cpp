#include "kproc_state.hpp"

#include <stdexcept>

#include <Omega_h_for.hpp>
#include <hadoken/format/format.hpp>

#include "diffusions.hpp"
#include "opsplit/compdef.hpp"
#include "opsplit/statedef.hpp"
#include "opsplit/vocabulary.hpp"
#include "reactions.hpp"


namespace zee {

//------------------------------------------------------------------

template <osh::Int Dim>
KProcState::KProcState(const Statedef& statedef, OmegaHMesh<Dim>& mesh, bool discovery)
    : pStatedef(statedef)
    , reactions_(statedef, mesh, discovery)
    , surface_reactions_(statedef, mesh, discovery) {
    if (!discovery) {
        setupDependencies();
    }
}

//------------------------------------------------------------------

template <typename KineticProcesses>
void KProcState::collateDependencies(const KineticProcesses& processes,
                                     DependenciesMap& dependency_map) const {
    for (size_t k = 0; k < processes.size(); k++) {
        auto dependencies = processes.getPropensityDependency(k);
        for (const auto& d: dependencies) {
            dependency_map[d].push_back({processes.getKProcType(), k});
        }
    }
}

//------------------------------------------------------------------

template <typename KineticProcesses>
void KProcState::cacheDependencies(const KineticProcesses& processes,
                                   const DependenciesMap& dependency_map,
                                   std::vector<KProcDeps>& dependencies) {
    std::vector<std::pair<std::set<KProcID>, std::set<MolState::ElementID>>> unique_deps(
        processes.size());
    // compute dependencies in a temporary datastructure that avoid taking the same dependency
    // more than once
    for (size_t k = 0; k < processes.size(); ++k) {
        std::vector<model::region_id> regions;
        std::vector<MolState::ElementID> mol_state_elements_to_update;
        std::tie(regions, mol_state_elements_to_update) = processes.getMolStateElementsUpdates(k);
        for (size_t l = 0; l < mol_state_elements_to_update.size(); l++) {
            const MolState::ElementID mol_state_element = mol_state_elements_to_update[l];
            container::specie_id spec = std::get<1>(mol_state_element);
            auto it = dependency_map.find(mol_state_element);
            if (it != dependency_map.end()) {
                const std::vector<KProcID>& deps = it->second;
                unique_deps[k].first.insert(deps.begin(), deps.end());
            }
            auto is_specie_diffused = boost::apply_visitor(
                [&spec, this](auto v) { return pStatedef.getDefinition(v).isDiffused(spec); },
                regions[l]);
            if (is_specie_diffused) {
                unique_deps[k].second.insert(mol_state_element);
            }
        }
    }

    // update the datastructure given in parameter accordingly to the previous computation
    dependencies.resize(processes.size());
    for (auto pid = 0u; pid < processes.size(); ++pid) {
        {
            auto& deps = dependencies[pid].propensities;
            const auto& udeps = unique_deps[pid].first;
            std::copy(udeps.begin(), udeps.end(), std::back_inserter(deps));
        }
        {
            auto& deps = dependencies[pid].occupancy;
            const auto& udeps = unique_deps[pid].second;
            std::copy(udeps.begin(), udeps.end(), std::back_inserter(deps));
        }
    }
}

//------------------------------------------------------------------

void KProcState::setupDependencies() {
    DependenciesMap dependency_map;

    // collate kinetic processes dependencies
    collateDependencies(reactions(), dependency_map);
    collateDependencies(surfaceReactions(), dependency_map);
    // caching the kinetic processes dependencies
    reactions_dependencies_.resize(reactions().size());
    surface_reactions_dependencies_.resize(surfaceReactions().size());
    cacheDependencies(reactions(), dependency_map, reactions_dependencies_);
    cacheDependencies(surfaceReactions(), dependency_map, surface_reactions_dependencies_);
}

//------------------------------------------------------------------

void KProcState::updateMolState(MolState& mol_state, const KProcID& event) const {
    switch (event.first) {
    case KProcType ::Reac:
        reactions().apply(mol_state, event.second);
        break;
    case KProcType::SReac:
        surfaceReactions().apply(mol_state, event.second);
        break;
    case KProcType::SDiff:
    case KProcType::VDepSReac:
    case KProcType::Diff:
        throw std::logic_error(
            hadoken::scat("Unhandled kinetic process ", static_cast<int>(event.first)));
    }
}

//------------------------------------------------------------------

const KProcDeps& KProcState::extractDependenciesFromEvent(const KProcID& event) const {
    const KProcDeps* dependencies{};
    switch (event.first) {
    case KProcType::Reac:
        dependencies = &reactions_dependencies_[event.second];
        break;
    case KProcType::SReac:
        dependencies = &surface_reactions_dependencies_[event.second];
        break;
    case KProcType::SDiff:
    case KProcType::VDepSReac:
    case KProcType::Diff:
        throw std::invalid_argument(
            hadoken::scat("Unhandled kinetic process ", static_cast<int>(event.first)));
    }
    return *dependencies;
}

//------------------------------------------------------------------

void KProcState::updateOccupancy(DiffusionVariables& diffusionVariables,
                                 MolState& mol_state,
                                 const PetscScalar simulation_time,
                                 const KProcID& event) const {
    const auto& dependencies = extractDependenciesFromEvent(event);
    const KProcDeps::Occupancy& occupancy_updates = dependencies.occupancy;
    for (const auto& d: occupancy_updates) {
        const auto& element_id = std::get<0>(d);
        const auto& specie_id = std::get<1>(d);

        // Careful here: SDiff are not yet handled
        boost::apply_visitor(
            [&](auto entity) {
                this->updateOccupancy(
                    entity, specie_id, diffusionVariables, mol_state, simulation_time);
            },
            element_id);
    }
}

void KProcState::updateOccupancy(mesh::boundary_id,
                                 container::specie_id,
                                 DiffusionVariables&,
                                 MolState&,
                                 PetscScalar) const {
    throw std::logic_error("Occupancy updates are not yet handled for surfaces.");
}

void KProcState::updateOccupancy(mesh::element_id element,
                                 container::specie_id specie,
                                 DiffusionVariables& diffusionVariables,
                                 MolState& mol_state,
                                 const PetscScalar simulation_time) const {
    const auto update_time_interval = simulation_time -
                                      diffusionVariables.last_update_time(element, specie);
    diffusionVariables.occupancy(element, specie) +=
        update_time_interval * static_cast<Omega_h::Real>(mol_state(element, specie));
    diffusionVariables.last_update_time(element, specie) = simulation_time;
}

//------------------------------------------------------------------

void KProcState::updateAllPropensities(const MolState& mol_state) {
    ssa_a0 = 0.0;
    for (size_t index = 0; index < reactions_.size(); ++index) {
        ssa_a0 += reactions_.updateRate(mol_state, index);
    }
    for (size_t index = 0; index < surface_reactions_.size(); ++index) {
        ssa_a0 += surface_reactions_.updateRate(mol_state, index);
    }
}

//------------------------------------------------------------------

void KProcState::updatePropensities(const MolState& mol_state, const KProcID& event) {
    const auto& dependencies = extractDependenciesFromEvent(event);
    for (const auto& k: dependencies.propensities) {
        switch (k.first) {
        case KProcType::Reac:
            reactions_.updateRate(mol_state, k.second);
            break;
        case KProcType::SReac:
            surface_reactions_.updateRate(mol_state, k.second);
            break;
        case KProcType::SDiff:
        case KProcType::VDepSReac:
        case KProcType::Diff:
            throw std::invalid_argument(
                hadoken::scat("Unhandled kinetic process ", static_cast<int>(k.first)));
        }
    }
    updateSSA_A0();
}

//------------------------------------------------------------------

PetscScalar KProcState::updateSSA_A0() {
    ssa_a0 = 0.0;
    for (size_t index = 0; index < reactions().size(); ++index) {
        ssa_a0 += reactions().getRate(index);
    }
    for (size_t index = 0; index < surfaceReactions().size(); ++index) {
        ssa_a0 += surfaceReactions().getRate(index);
    }
    return ssa_a0;
}

//------------------------------------------------------------------

void KProcState::report(std::ostream& report_stream) const {
    report_stream << "KProcState Report" << '\n';
}

//------------------------------------------------------------------

std::pair<osh::LOs, boost::optional<osh::LOs>> KProcState::getNumberOfSpeciesPerOwnedElement()
    const {
    const auto& per_boundary_owned = surfaceReactions().getNumberOfSpeciesPerBoundaryOwnedElement();
    const boost::optional<osh::LOs> boundary_elements =
        per_boundary_owned.size() > 0 ? boost::optional<osh::LOs>(per_boundary_owned) : boost::none;
    return {reactions().getNumberOfSpeciesPerOwnedElement(), boundary_elements};
}

//------------------------------------------------------------------

template <class Generator>
KProcID KProcState::getNextEvent(Generator& rng) const {
    PetscScalar a0(getSSA_A0());
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const PetscScalar selector = a0 * uniform(rng);
    PetscScalar partial_sums = 0.0;
    for (size_t index = 0; index < reactions_.size(); ++index) {
        const PetscScalar rate = reactions_.getRate(index);
        partial_sums += rate;
        if (selector < partial_sums) {
            return {KProcType::Reac, index};
        }
    }
    for (size_t index = 0; index < surface_reactions_.size(); ++index) {
        const PetscScalar rate = surface_reactions_.getRate(index);
        partial_sums += rate;
        if (selector < partial_sums) {
            return {KProcType::SReac, index};
        }
    }
    throw std::invalid_argument("SSA engine is unable to select a KProc.");
}

//------------------------------------------------------------------

// explicit template instantiation definitions
template KProcState::KProcState(const Statedef& statedef, OmegaHMesh<2>& mesh, bool discovery);
template KProcState::KProcState(const Statedef& statedef, OmegaHMesh<3>& mesh, bool discovery);
template KProcID KProcState::getNextEvent(std::mt19937& rng) const;

}  // namespace zee
