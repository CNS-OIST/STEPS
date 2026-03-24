#include "efield.hpp"

#include "../mol_state.hpp"
#include "geom/dist/distmemb.hpp"
#include "geom/dist/distmesh.hpp"
#include "geom/patch.hpp"
#include "math/constants.hpp"
#include "model/chan.hpp"
#include "model/chanstate.hpp"
#include "model/ghkcurr.hpp"
#include "model/ohmiccurr.hpp"
#include "model/surfsys.hpp"
#include "mpi/dist/tetopsplit/definition/complexeventsdef.hpp"
#include "patchdef.hpp"
#include "statedef.hpp"
#include "util/vocabulary.hpp"
#include <stdexcept>
#include <string>

namespace steps::dist {

osh::Real OhmicCurrdefBase::getReversalPotential(mesh::triangle_id_t triangle) const {
    auto it = reversal_potentials.find(triangle);
    if (it != reversal_potentials.end()) {
        return it->second;
    }
    return reversal_potential;
}

void OhmicCurrdefBase::setReversalPotential(mesh::triangle_id_t triangle, osh::Real value) {
    reversal_potentials.emplace(triangle, value);
}

void OhmicCurrdefBase::reset() {
    reversal_potentials.clear();
}

OhmicCurrdef::OhmicCurrdef(const Patchdef& patchdef, const steps::model::OhmicCurr& curr)
    : OhmicCurrdefBase(curr)
    , channel_state(patchdef.getSpecContainerIdx(curr.getChanState())) {}

ComplexOhmicCurrdef::ComplexOhmicCurrdef(const Patchdef& patchdef,
                                         const steps::model::ComplexOhmicCurr& curr)
    : OhmicCurrdefBase(curr)
    , channel_state(curr.getChanState(), patchdef.statedef()) {}

GHKCurrdefBase::GHKCurrdefBase(const Patchdef& patchdef, const steps::model::GHKcurrBase& curr)
    : ion_id(curr.getIon().getID())
    , valence(curr._valence()) {
    if (not curr._infosupplied()) {
        std::ostringstream msg;
        msg << "GHK current " << curr.getID() << ": Undefined permeability.";
        throw std::invalid_argument(msg.str());
    }
    if (not curr._realflux()) {
        throw std::invalid_argument(
            "GHK currents in distributed STEPS do not support "
            "the computeflux=False argument.");
    }
    if (curr._vshift() != 0.0) {
        throw std::invalid_argument(
            "GHK currents in distributed STEPS do not support "
            "the vshift argument.");
    }
}

GHKCurrdef::GHKCurrdef(const Patchdef& patchdef, const steps::model::GHKcurr& curr)
    : GHKCurrdefBase(patchdef, curr)
    , channel_state(patchdef.getSpecContainerIdx(curr.getChanState())) {}

ComplexGHKCurrdef::ComplexGHKCurrdef(const Patchdef& patchdef,
                                     const steps::model::ComplexGHKcurr& curr)
    : GHKCurrdefBase(patchdef, curr)
    , channel_state(curr.getChanState(), patchdef.statedef()) {}

#ifdef USE_PETSC

PetscReal OhmicCurrdefBase::getTriCurrentOnVertex(const osh::Real potential_on_vertex,
                                                  const mesh::triangle_id_t& b_id,
                                                  const MolState& mol_state,
                                                  const DistMesh& mesh,
                                                  const osh::Real sim_time) const {
    // A tri split among the vertexes
    const double Avert = mesh.getTri(b_id).area / 3.0;
    const PetscReal tri_oc_bc = getTriBConVertex(b_id, mol_state, Avert, sim_time);

    return tri_oc_bc * (potential_on_vertex - getReversalPotential(b_id));
}

PetscReal OhmicCurrdef::getTriBConVertex(const mesh::triangle_id_t& b_id,
                                         const MolState& mol_state,
                                         const double Avert,
                                         const osh::Real sim_time) const {
    const auto avg_open_channels = mol_state.get_occupancy_ef(b_id, channel_state, sim_time) / 3.0;

    return avg_open_channels * conductance;
}

PetscReal ComplexOhmicCurrdef::getTriBConVertex(const mesh::triangle_id_t& b_id,
                                                const MolState& mol_state,
                                                const double /*Avert*/,
                                                const osh::Real sim_time) const {
    const auto avg_open_channels = mol_state.get_occupancy_ef(b_id, occupancy_id, sim_time) / 3.0;

    return avg_open_channels * conductance;
}

#endif  // USE_PETSC

Membrane::Membrane(Statedef& statedef_, const DistMemb& membrane)
    : memb_(membrane)
    , statedef(statedef_)
    , patches_ids(membrane.patches().begin(), membrane.patches().end())
    , capacitance_(membrane.getCapacitance())
    , current_([](auto) { return 0.0; }) {
    for (const auto& patchId: patches_ids) {
        const auto& patch = statedef.mesh().getPatch(patchId);
        auto& patchdef = statedef.getPatchdef(patchId);
        patchdefs.emplace_back(patchdef);
        for (auto& ssysName: patch.getSurfsys()) {
            auto& ssys = statedef.model().getSurfsys(ssysName);
            // Ohmic currents
            for (auto& [currId, curr]: ssys._getAllOhmicCurrs()) {
                patchdef.addCurrent(*curr);
            }
            for (auto& [currId, curr]: ssys._getAllComplexOhmicCurrs()) {
                patchdef.addCurrent(*curr);
            }
            // GHK currents
            for (auto& [currId, curr]: ssys._getAllGHKcurrs()) {
                patchdef.addCurrent(*curr);
                patchdef.addGHKReacs(*curr);
            }
            for (auto& [currId, curr]: ssys._getAllComplexGHKcurrs()) {
                patchdef.addCurrent(*curr);
                patchdef.addGHKReacs(*curr);
            }
        }
    }
}

void Membrane::reset() {
    capacitance_ = memb_.getCapacitance();
    // Currents are owned by statedef so they are reset there
    conductivity_ = 0;
    reversal_potential_ = 0;
    current_ = [](auto) { return 0.0; };
}

std::ostream& operator<<(std::ostream& os, OhmicCurrdef const& m) {
    os << "OhmicCurrent.conductance: " << m.conductance
       << "\nOhmicCurrent.reversal_potential: " << m.reversal_potential
       << "\nOhmicCurrent.channel_state: " << std::to_string(m.channel_state);
    return os << "\n";
}

#ifdef USE_PETSC

std::ostream& operator<<(std::ostream& os, const TriMatAndVecs& obj) {
    os << "vert_idxs:\n";
    for (const auto i: obj.face_bf2vertsPETSc) {
        os << i << ' ';
    }
    os << '\n';
    os << "triStiffnessMat:\n";
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            os << obj.triStiffnessPETSc[3 * j + i] << ' ';
        }
        os << '\n';
    }
    os << "triBC:\n";
    for (const auto i: obj.triBC) {
        os << i << ' ';
    }
    os << '\n';
    os << "triI:\n";
    for (const auto i: obj.triI) {
        os << i << ' ';
    }
    os << '\n';

    return os;
}

#endif  // USE_PETSC

}  // namespace steps::dist
