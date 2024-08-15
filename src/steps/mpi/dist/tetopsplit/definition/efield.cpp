#include "efield.hpp"

#include "../mol_state.hpp"
#include "geom/dist/distmesh.hpp"

namespace steps::dist {

osh::Real OhmicCurrent::getReversalPotential(mesh::triangle_id_t triangle) const {
    auto it = reversal_potentials.find(triangle);
    if (it != reversal_potentials.end()) {
        return it->second;
    }
    return reversal_potential;
}

void OhmicCurrent::setReversalPotential(mesh::triangle_id_t triangle, osh::Real value) {
    reversal_potentials.emplace(triangle, value);
}

void OhmicCurrent::reset() {
    reversal_potentials.clear();
}

#ifdef USE_PETSC

PetscReal OhmicCurrent::getTriCurrentOnVertex(const osh::Real potential_on_vertex,
                                              const mesh::triangle_id_t& b_id,
                                              const MolState& mol_state,
                                              const DistMesh& mesh,
                                              const osh::Real sim_time) const {
    // A tri split among the vertexes
    const double Avert = mesh.getTri(b_id).area / 3.0;
    const PetscReal tri_oc_bc = getTriBConVertex(b_id, mol_state, Avert, sim_time);

    return tri_oc_bc * (potential_on_vertex - getReversalPotential(b_id));
}

PetscReal OhmicCurrent::getTriBConVertex(const mesh::triangle_id_t& b_id,
                                         const MolState& mol_state,
                                         const double Avert,
                                         const osh::Real sim_time) const {
    const auto avg_open_channels =
        channel_state ? mol_state.get_occupancy_ef(b_id, *channel_state, sim_time) / 3.0 : Avert;

    return avg_open_channels * conductance;
}

#endif  // USE_PETSC

std::ostream& operator<<(std::ostream& os, OhmicCurrent const& m) {
    return os << "OhmicCurrent.conductance: " << m.conductance
              << "\nOhmicCurrent.reversal_potential: " << m.reversal_potential
              << "\nOhmicCurrent.channel_state: "
              << (m.channel_state ? std::to_string(*m.channel_state) : "not assigned") << '\n';
}

std::ostream& operator<<(std::ostream& os, GHKCurrent const& m) {
    return os << "GHKCurrent.ion_channel_state: " << m.ion_channel_state
              << "\nGHKCurrent.ion_id: " << m.ion_id << "\nGHKCurrent.valence: " << m.valence
              << '\n';
}

std::ostream& operator<<(std::ostream& os, Channel const& m) {
    os << "Channel.channel_states: ";
    for (const auto& i: m.channel_states) {
        os << i << ' ';
    }
    os << "\nChannel.ohmic_currents.size(): " << m.ohmic_currents.size() << '\n';
    for (const auto& i: m.ohmic_currents) {
        os << i;
    }
    os << "Channel.ghk_currents.size(): " << m.ghk_currents.size() << '\n';
    for (const auto& i: m.ghk_currents) {
        os << i;
    }
    return os;
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
