#include "reactions.hpp"

#include <Omega_h_map.hpp>
#include <Omega_h_mark.hpp>

#include "kproc_state.hpp"
#include "geom/dist/distmesh.hpp"
#include "math/constants.hpp"
#include "util/vocabulary.hpp"

namespace steps {
namespace dist {
namespace kproc {

template <typename NumMolecules>
Reactions::Reactions(const Statedef& statedef, DistMesh& mesh, MolState<NumMolecules>& mol_state)
    : measureInfo(mesh.getMeasure()) {
    const auto& owned_elems_mask = mesh.owned_elems_mask();

    for (const auto& compartment: statedef.compdefs()) {
        const auto& elements = mesh.getEntities(compartment->getID());
        for (auto k: elements) {
            if (owned_elems_mask[k.get()]) {
                for (const auto& reacdef: compartment->reacdefs()) {
                    reacdefs_.push_back(*reacdef);
                    ownerPoints_.push_back(k);
                    ccsts_.push_back(compute_ccst(*reacdef, k));
                    std::vector<osh::I64> stoichiometry_change;
                    std::vector<MolStateElementID> reaction_upd;
                    std::vector<MolStateElementID> reaction_lhs;
                    const auto& upd_array = reacdef->getPoolChangeUPD();
                    for (size_t spec = 0; spec < upd_array.size(); spec++) {
                        if (upd_array[spec] != 0) {
                            container::species_id speciesId(static_cast<int>(spec));
                            reaction_upd.push_back(mkVolumeElement(k, speciesId));
                            stoichiometry_change.push_back(upd_array[spec]);

                            // track occupancy if the molecule can diffuse (here we have only
                            // molecules, not channel states)
                            if (compartment->isDiffused(speciesId)) {
                                mol_state.track_occupancy_rd(k, speciesId);
                            }
                        }
                    }
                    const auto& lhs_array = reacdef->getPoolChangeLHS();
                    for (size_t spec = 0; spec < lhs_array.size(); ++spec) {
                        if (lhs_array[spec] != 0) {
                            const container::species_id spec_id(static_cast<int>(spec));
                            reaction_lhs.emplace_back(k, spec_id);
                        }
                    }
                    reactions_upd_.push_back(reaction_upd);
                    reactions_lhs_.push_back(reaction_lhs);
                    stoichiometry_change_.push_back(stoichiometry_change);
                }
            }
        }
  }
}

//------------------------------------------------------------------

void Reactions::report(std::ostream &report_stream, size_t index) const {
    getReacDef(index).report(report_stream, ownerPoints_[index]);
}

//------------------------------------------------------------------

template <typename NumMolecules>
osh::Real Reactions::computeRate(const MolState<NumMolecules> &mol_state,
                                 size_t index) const {
  const auto &lhs = reacdefs_[index].get().getPoolChangeLHS();
  osh::Real h_mu = 1.0;
  const container::species_id num_species(
      static_cast<container::species_id::value_type>(lhs.size()));
  for (container::species_id species(0); species < num_species; species++) {
    osh::I64 lhs_s = -lhs[static_cast<size_t>(species.get())];
    if (lhs_s == 0) {
      continue;
    }
    auto pool_s = mol_state(this->getOwnerPoint(index), species);
    if (lhs_s > pool_s) {
      h_mu = 0.0;
      break;
    }
    switch (lhs_s) {
    case 4: {
      h_mu *= static_cast<osh::Real>(pool_s - 3);
      OMEGA_H_FALLTHROUGH;
    }
    case 3: {
      h_mu *= static_cast<osh::Real>(pool_s - 2);
      OMEGA_H_FALLTHROUGH;
    }
    case 2: {
      h_mu *= static_cast<osh::Real>(pool_s - 1);
      OMEGA_H_FALLTHROUGH;
    }
    case 1: {
      h_mu *= static_cast<osh::Real>(pool_s);
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

osh::Real Reactions::compute_ccst(const Reacdef &reacdef,
                                  mesh::tetrahedron_id_t element) const {
  const auto measure = measureInfo.element_measure(element);
  osh::Real scale = 1.0e3 * measure * math::AVOGADRO;
  osh::I64 o1 = reacdef.getOrder() - 1;
  osh::Real ccst =
      reacdef.getKcst() * std::pow(scale, static_cast<osh::Real>(-o1));
  return ccst;
}

//------------------------------------------------------------------

template <typename NumMolecules>
const std::vector<MolStateElementID>& Reactions::updateMolStateAndOccupancy(
    MolState<NumMolecules>& mol_state,
    size_t index,
    const osh::Real event_time) const {
    const auto& upd = reactions_upd_[index];
    const auto& stoichoimetry = stoichiometry_change_[index];
    for (size_t k = 0; k < upd.size(); k++) {
        const auto& elmt = upd[k];
        const auto& s = stoichoimetry[k];

        assert(mol_state(elmt) >= -s);
        assert(mol_state(elmt) <= std::numeric_limits<NumMolecules>::max() - std::max(s, {}));
        mol_state.add_and_update_occupancy(elmt, static_cast<osh::LO>(s), event_time);
    }
    return upd;
}

//------------------------------------------------------------------

// explicit template instantiation definitions
template osh::Real Reactions::computeRate(const MolState<osh::I32> &mol_state,
                                          size_t index) const;
template osh::Real Reactions::computeRate(const MolState<osh::I64> &mol_state,
                                          size_t index) const;
template const std::vector<MolStateElementID>& Reactions::updateMolStateAndOccupancy(
    MolState<osh::I32>& mol_state,
    size_t index,
    const osh::Real event_time) const;
template const std::vector<MolStateElementID>& Reactions::updateMolStateAndOccupancy(
    MolState<osh::I64>& mol_state,
    size_t index,
    const osh::Real event_time) const;

template Reactions::Reactions(const Statedef& statedef,
                              DistMesh& mesh,
                              MolState<osh::LO>& mol_state);
template Reactions::Reactions(const Statedef& statedef,
                              DistMesh& mesh,
                              MolState<osh::GO>& mol_state);

} // namespace kproc
} // namespace dist
} // namespace steps
