#include "opsplit/mesh.hpp"

#include <hadoken/format/format.hpp>

namespace zee {

using hadoken::scat;

DistMesh::DistMesh(std::string t_filename, PetscScalar t_scale, bool t_label_elems)
    : filename(std::move(t_filename))
    , scale(t_scale)
    , label_elems(t_label_elems) {}

DistMesh::~DistMesh() = default;

PetscInt DistMesh::getTotalNumElements() const {
    auto num_elements = getOwnedNumElements();
    PetscInt total_num_elements{};
    auto err =
        MPI_Allreduce(&num_elements, &total_num_elements, 1, MPIU_INT, MPI_SUM, MPI_COMM_WORLD);
    if (err != 0) {
        MPI_Abort(MPI_COMM_WORLD, err);
    }
    return total_num_elements;
}

const mesh::element_ids& DistMesh::getCompElements(model::compartment_label comp_label) const {
    auto result = complabel2elems.find(comp_label);
    if (result != complabel2elems.end()) {
        return result->second;
    }
    throw std::invalid_argument("Undefined compartment label.");
}

const mesh::element_ids& DistMesh::getCompElements(const model::compartment_id& comp_id) const {
    const auto label = compIdtoLabel.find(comp_id);
    if (label != compIdtoLabel.end()) {
        return getCompElements(label->second);
    }
    throw std::invalid_argument("Undefined compartment id.");
}

void DistMesh::addComp(const model::compartment_id& compartment,
                       model::compartment_label cell_set_label) {
    compIdtoLabel[compartment] = cell_set_label;
    compLabelToId[cell_set_label] = compartment;
    addCompImpl(compartment, cell_set_label);
}

model::compartment_label DistMesh::getCompLabel(const model::compartment_id& compartment) const {
    auto result = compIdtoLabel.find(compartment);
    if (result != compIdtoLabel.end()) {
        return result->second;
    }
    throw std::invalid_argument(scat("Undefined compartment id: ", compartment));
}

const model::compartment_id& DistMesh::getCompID(model::compartment_label comp_label) const {
    auto result = compLabelToId.find(comp_label);
    if (result != compLabelToId.end()) {
        return result->second;
    }
    throw std::invalid_argument("Undefined compartment label.");
}

PetscScalar DistMesh::getTotalCompVol(model::compartment_label comp_label) const {
    PetscScalar own_vol = getOwnedCompVol(comp_label);
    PetscScalar total_vol = 0.0;
    MPI_Allreduce(&own_vol, &total_vol, 1, MPIU_SCALAR, MPI_SUM, MPI_COMM_WORLD);
    return total_vol;
}

PetscScalar DistMesh::getTotalPatchArea(model::patch_id patchId) {
    PetscScalar own_area = getOwnedPatchArea(patchId);
    PetscScalar total_area = 0.0;
    MPI_Allreduce(&own_area, &total_area, 1, MPIU_SCALAR, MPI_SUM, MPI_COMM_WORLD);
    return total_area;
}

}  // namespace zee
