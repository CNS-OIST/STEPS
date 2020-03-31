#include <iostream>

#include <Omega_h_map.hpp>
#include <Omega_h_mark.hpp>
#include <hadoken/format/format.hpp>

#include "kproc_state.hpp"
#include "opsplit/vocabulary.hpp"
#include "reactions.hpp"


namespace zee {
namespace kproc {

template <osh::Int Dim>
Reactions::Reactions(const zee::Statedef& statedef, zee::OmegaHMesh<Dim>& mesh, bool discovery)
    : measureInfo(mesh.getMeasureInfo()) {
    const auto& owned_elems_mask = mesh.getOwnedElemsMask();

    osh::Write<osh::LO> num_species_per_element(owned_elems_mask.size(), 0);
    osh::Write<osh::LO> num_species_per_owned_element(owned_elems_mask.size(), 0);

    for (const auto& compartment: statedef.compdefs()) {
        const auto& elements = mesh.getEntities(compartment->getID());
        for (auto k: elements) {
            if (owned_elems_mask[k.get()]) {
                num_species_per_owned_element[k.get()] = compartment->getNSpecs();
            }
            num_species_per_element[k.get()] = compartment->getNSpecs();
            if (discovery) {
                continue;
            }
            if (owned_elems_mask[k.get()]) {
                for (const auto& reacdef: compartment->reacdefs()) {
                    reacdefs_.push_back(*reacdef);
                    ownerPoints_.push_back(k);
                    ccsts_.push_back(compute_ccst(*reacdef, k));
                    rates_.push_back(0.);
                }
            }
        }
    }
    num_species_per_element_ = num_species_per_element;
    num_species_per_owned_element_ = num_species_per_owned_element;
}

//------------------------------------------------------------------

void Reactions::report(std::ostream& report_stream, size_t index) const {
    report_stream << "Type: "
                  << static_cast<int>(KProcType::Reac)  // << " Container Idx: " << getKProcLIdx()
                  << " Rate: " << rates_[index];
    report_stream << "{";
    report_stream << "}" << '\n';
    report_stream << '\n';
}

//------------------------------------------------------------------

PetscScalar Reactions::computeRate(const MolState& mol_state, size_t index) const {
    const auto& lhs = reacdefs_[index].getPoolChangeArray(Reacdef::PoolChangeArrayType::LHS);
    PetscScalar h_mu = 1.0;
    const container::specie_id num_species(
        static_cast<container::specie_id::value_type>(lhs.size()));
    for (container::specie_id specie(0); specie < num_species; specie++) {
        PetscInt lhs_s = -lhs[static_cast<size_t>(specie.get())];
        if (lhs_s == 0) {
            continue;
        }
        auto pool_s = mol_state(this->getOwnerPoint(index), specie);
        if (lhs_s > pool_s) {
            h_mu = 0.0;
            break;
        }
        switch (lhs_s) {
        case 4: {
            h_mu *= (pool_s - 3);
            OMEGA_H_FALLTHROUGH;
        }
        case 3: {
            h_mu *= (pool_s - 2);
            OMEGA_H_FALLTHROUGH;
        }
        case 2: {
            h_mu *= (pool_s - 1);
            OMEGA_H_FALLTHROUGH;
        }
        case 1: {
            h_mu *= pool_s;
            break;
        }
        default: {
            throw std::runtime_error("Reaction rate computation error");
        }
        }
    }
    return h_mu * ccsts_[index];
}

//------------------------------------------------------------------

PetscScalar Reactions::updateRate(const MolState& mol_state, size_t index) {
    rates_[index] = computeRate(mol_state, index);
    return rates_[index];
}

//------------------------------------------------------------------

PetscScalar Reactions::compute_ccst(const Reacdef& reacdef, mesh::element_id element) const {
    const auto measure = measureInfo.element_measure(element);
    PetscScalar scale = 1.0e3 * measure * AVOGADRO;
    PetscInt o1 = reacdef.getOrder() - 1;
    PetscScalar ccst = reacdef.getKcst() * std::pow(scale, static_cast<PetscScalar>(-o1));
    return ccst;
}

//------------------------------------------------------------------

void Reactions::apply(MolState& mol_state, size_t index) const {
    const std::vector<PetscInt>& upd_array = reacdefs_[index].getPoolChangeArray(
        Reacdef::PoolChangeArrayType::UPD);
    const container::specie_id num_species(
        static_cast<container::specie_id::value_type>(upd_array.size()));
    const auto element = this->getOwnerPoint(index);
    for (container::specie_id specie(0); specie < num_species; specie++) {
        assert(mol_state(element, specie) >= -upd_array[static_cast<size_t>(specie.get())] &&
               mol_state(element, specie) <=
                   std::numeric_limits<MolState::elem_type>::max() -
                       std::max(upd_array[static_cast<size_t>(specie.get())], {}));

        mol_state(element,
                  specie) += static_cast<Omega_h::LO>(upd_array[static_cast<size_t>(specie.get())]);
    }
}

//------------------------------------------------------------------

std::vector<MolState::ElementID> Reactions::getPropensityDependency(size_t index) const {
    const std::vector<PetscInt>& lhs = getReacDef(index).getPoolChangeArray(
        Reacdef::PoolChangeArrayType::LHS);
    std::vector<MolState::ElementID> propensity_dependency;
    for (size_t l = 0; l < lhs.size(); ++l) {
        if (lhs[l] != 0) {
            container::specie_id spec_id(static_cast<int>(l));
            propensity_dependency.push_back(
                MolState::mkVolumeElement(getOwnerPoint(index), spec_id));
        }
    }
    return propensity_dependency;
}

//------------------------------------------------------------------

std::pair<std::vector<model::region_id>, std::vector<MolState::ElementID>>
Reactions::getMolStateElementsUpdates(size_t index) const {
    std::vector<MolState::ElementID> mol_state_elements_updates;
    std::vector<model::region_id> region_ids;
    const auto& reac = getReacDef(index);
    const std::vector<PetscInt>& change = reac.getPoolChangeArray(
        Reacdef::PoolChangeArrayType::UPD);
    const model::compartment_id comp_id = reac.compdef().getID();
    for (size_t l = 0; l < change.size(); ++l) {
        if (change[l] != 0) {
            container::specie_id spec(static_cast<int>(l));
            auto element_id = getOwnerPoint(index);
            region_ids.push_back(comp_id);
            mol_state_elements_updates.push_back(MolState::mkVolumeElement(element_id, spec));
        }
    }
    return std::make_pair(region_ids, mol_state_elements_updates);
}

//------------------------------------------------------------------

// explicit template instantiation declarations
template Reactions::Reactions(const zee::Statedef& statedef,
                              zee::OmegaHMesh<2>& mesh,
                              bool discovery);
template Reactions::Reactions(const zee::Statedef& statedef,
                              zee::OmegaHMesh<3>& mesh,
                              bool discovery);
}  // namespace kproc
}  // namespace zee
