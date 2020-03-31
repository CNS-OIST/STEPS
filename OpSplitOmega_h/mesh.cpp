#include <algorithm>
#include <cassert>
#include <iostream>

#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mark.hpp>
#include <Omega_h_shape.hpp>
#include <hadoken/format/format.hpp>
#include <petsctime.h>

#include "mesh.hpp"
#include "mesh_utils.hpp"
#include "mpitools.hpp"
#include "opsplit/type_traits.hpp"
#include "opsplit/vocabulary.hpp"

#include <iterator>

namespace zee {

using hadoken::scat;

#define INITIAL_COMPARTMENT_ID -1

/// \brief provide information about the neighbors of a tetrahedron
/// Intermediate structure used to fill the flat_multimap
struct TetNeighborStruct {
    /// indices of the tetrahedron elements
    std::array<osh::LO, 4> indices{{-1, -1, -1, -1}};
    /// face identifiers for the given tetrahedron
    std::array<osh::LO, 4> faces{};
    /// distances from the tetrahedron barycenters
    std::array<PetscScalar, 4> distances{};
    /// surface of the faces
    std::array<PetscScalar, 4> areas{};
};

MeasureInfo::MeasureInfo(osh::Int num_compartments,
                         const element_measure_func& t_element_measure_func)
    : rank_(mpi_comm_rank())
    , num_measures_per_rank_(1 + num_compartments)
    , element_measure_func_(t_element_measure_func) {}

void MeasureInfo::init(const mesh::element_ids& t_owned_elements, const osh::LOs& t_elem2compid) {
    const auto comm_size = mpi_comm_size();
    osh::Write<osh::Real> rank2measures(num_measures_per_rank_ * comm_size, 0);
    auto rankmeasures = rank2measures.data() + rank_ * num_measures_per_rank_;
    for (auto element: t_owned_elements) {
        const auto measure = element_measure_func_(element);
        const auto compid = t_elem2compid[element.get()];
        rankmeasures[0] += measure;
        if (compid != INITIAL_COMPARTMENT_ID) {
            rankmeasures[1 + compid] += measure;
        }
    }

    {
        int err = MPI_Allgather(MPI_IN_PLACE,
                                num_measures_per_rank_,
                                MPIU_REAL,
                                rank2measures.data(),
                                num_measures_per_rank_,
                                MPIU_REAL,
                                MPI_COMM_WORLD);
        if (err != 0) {
            MPI_Abort(MPI_COMM_WORLD, err);
        }
        rank2measures_ = rank2measures;
    }

    {
        osh::Write<osh::Real> mesh_measures(num_measures_per_rank_, 0);
        for (auto rank = 0; rank < comm_size; ++rank) {
            for (auto i = 0; i < num_measures_per_rank_; ++i) {
                mesh_measures[i] += rank2measures[rank * num_measures_per_rank_ + i];
            }
        }
        mesh_measures_ = mesh_measures;
    }
}

template <osh::Int Dim>
OmegaHMesh<Dim>::OmegaHMesh(osh::Mesh& t_mesh,
                            const std::string& t_filename,
                            PetscScalar t_scale,
                            bool t_label_elems)
    : DistMesh(t_filename, t_scale, t_label_elems)
    , mesh(t_mesh) {
    importFromFile();
}

template <osh::Int Dim>
void OmegaHMesh<Dim>::init() {
    {
        // ensure all elements have a dedicated compartment
        std::vector<osh::LO> bad_elements;
        for (osh::LO element{}; element < elem2compid.size(); ++element) {
            if (elem2compid[element] == INITIAL_COMPARTMENT_ID) {
                bad_elements.push_back(element);
            }
        }
        if (!bad_elements.empty()) {
            std::cerr << "Error: the following mesh elements miss an associated compartment: ";
            std::copy(bad_elements.begin(),
                      bad_elements.end(),
                      std::ostream_iterator<osh::LO>(std::cerr, " "));
            std::cerr << "\nAbort\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    this->measureInfo->init(this->owned_elems_, this->elem2compid);
}

template <osh::Int Dim>
std::string OmegaHMesh<Dim>::getBackend() const {
    return "Omega_h";
}

template <osh::Int Dim>
void OmegaHMesh<Dim>::importFromFile() {
    // create a distributed mesh object, calls mesh.balance() and then returns it
    {
        auto coords = osh::deep_copy(mesh.coords());
        osh::parallel_for(coords.size(), [=](osh::LO index) { coords[index] *= scale; });
        mesh.set_coords(coords);
    }
    this->mesh.set_parting(OMEGA_H_GHOSTED);
    ownedElemsMask = this->mesh.owned(this->mesh.dim());

    {
        std::vector<osh::LO> owned_elems;
        owned_elems.reserve(static_cast<size_t>(mesh.nelems()));
        for (osh::LO elem = 0; elem < mesh.nelems(); ++elem) {
            if (ownedElemsMask[elem] != 0) {
                owned_elems.push_back(elem);
            }
        }
        osh::Write<osh::LO> owned_elems_array(static_cast<osh::LO>(owned_elems.size()));
        std::copy(owned_elems.begin(), owned_elems.end(), owned_elems_array.data());
        owned_elems_ = osh::Read<osh::LO>(owned_elems_array);
    }

    const auto& coords = mesh.coords();
    const auto& areas = osh::measure_ents_real(&mesh,
                                               Omega_h::FACE,
                                               Omega_h::LOs(mesh.nents(Omega_h::FACE), 0, 1),
                                               coords);
    this->fill_triInfo(coords, areas);

    if (mesh.dim() == 3) {
        this->fill_tetInfo(coords, areas);
        measureFunc = [&](mesh::element_id element) {
            return tetInfo[static_cast<size_t>(element.get())].vol;
        };
    } else {
        osh::Write<osh::LO> neighbors(mesh.nelems());
        osh::parallel_for(
            neighbors.size(), OMEGA_H_LAMBDA(osh::LO elem_id) {
                neighbors[elem_id] = triInfo[static_cast<size_t>(elem_id)].num_neighbors;
            });
        neighbors_per_element_ = neighbors;
        measureFunc = [&](mesh::element_id element) {
            return triInfo[static_cast<size_t>(element.get())].area;
        };
    }
    elem2compid = osh::Write<osh::LO>(mesh.nelems(), INITIAL_COMPARTMENT_ID);
    measureInfo = std::make_unique<MeasureInfo>(num_compartments(), measureFunc);
}

template <osh::Int Dim>
void OmegaHMesh<Dim>::fill_triInfo(const Omega_h::Reals& coords, const Omega_h::Reals& areas) {
    const auto& verts2tris = mesh.ask_up(1, 2);
    const auto& tris2verts = mesh.ask_verts_of(Omega_h::FACE);

    auto fill_triInfo = OMEGA_H_LAMBDA(osh::LO triangle_id) {
        auto& info = triInfo[static_cast<size_t>(triangle_id)];
        const auto tri2verts = osh::gather_verts<3>(tris2verts, triangle_id);
        const auto tri2x = osh::gather_vectors<3, 2>(coords, tri2verts);
        info.centroid = (tri2x[0] + tri2x[1] + tri2x[2]) / 3.;
        info.area = areas[triangle_id];
        for (auto v = 0; v < 3; ++v) {
            auto vid = tris2verts[v];
            info.num_neighbors += static_cast<osh::LO>(
                verts2tris.a2ab[vid + 1] - verts2tris.a2ab[vid] > 1);
        }
    };
    triInfo.resize(static_cast<size_t>(areas.size()));
    osh::parallel_for(mesh.nents(Omega_h::FACE), fill_triInfo);
}

template <osh::Int Dim>
void OmegaHMesh<Dim>::fill_tetInfo(const Omega_h::Reals& coords, const Omega_h::Reals& areas) {
    // return the topology from each tetrahedron to its 4 faces as Omega_h::LOs
    const auto graph = mesh.ask_down(mesh.dim(), mesh.dim() - 1);
    const auto& tets2tris = graph.ab2b;
    // return for each triangle a data structure to get its neighboring tetrahedrons
    const auto tris2tets = mesh.ask_up(mesh.dim() - 1, mesh.dim());
    // for each a = triangle, a2ab gives the offset in ab2b
    const auto& f2fc = tris2tets.a2ab;
    // for each ab = tetrahedron neighboring a triangle, gives its local index
    const auto& fc2c = tris2tets.ab2b;

    const auto volumes = measure_elements_real(&mesh);
    const auto& tets2verts = mesh.ask_elem_verts();  // ask_verts_of(Omega_h::REGION);
    tetInfo.resize(static_cast<size_t>(volumes.size()));
    std::vector<TetNeighborStruct> neighborsInfo(static_cast<size_t>(volumes.size()));
    osh::Write<osh::LO> neighbors_per_element(volumes.size());
    const auto fill_tetInfo = [&](osh::LO osh_tet_id) {
        const auto tetrahedron_id = static_cast<size_t>(osh_tet_id);
        const auto tet2verts = osh::gather_verts<4>(tets2verts, osh_tet_id);
        const auto tet2x = osh::gather_vectors<4, 3>(coords, tet2verts);
        auto& tetinfo = tetInfo[tetrahedron_id];
        auto& neighborinfo = neighborsInfo[tetrahedron_id];
        tetinfo.centroid = barycenter(tet2x);
        tetinfo.vol = volumes[osh_tet_id];
        // returns the local indexes of the four faces of tet
        const auto tet2tris = Omega_h::gather_down<4>(tets2tris, osh_tet_id);
        auto neighbor_index = 0;
        for (auto f = 0; f < 4; ++f) {
            // local ordering of third and fourth neighbour is switched to compare with
            // STEPS
            auto triangle_id = tet2tris[f];
            // the local index of any tetrahedron neighboring tri is stored in fc2c[offset +
            // {0 or 1}]
            auto offset = f2fc[triangle_id];
            if (f2fc[triangle_id + 1] - offset > 1) {
                const auto neighbor_id = fc2c[offset + (fc2c[offset] == osh_tet_id ? 1 : 0)];
                const auto neighbor2verts = osh::gather_verts<4>(tets2verts, neighbor_id);
                const auto neighbor2x = osh::gather_vectors<4, 3>(coords, neighbor2verts);
                const auto neighbor_centroid = barycenter(neighbor2x);
                neighborinfo.indices[static_cast<size_t>(neighbor_index)] = neighbor_id;
                neighborinfo.areas[static_cast<size_t>(neighbor_index)] = areas[triangle_id];
                neighborinfo.distances[static_cast<size_t>(neighbor_index)] = norm(
                    tetinfo.centroid - neighbor_centroid);
                {
                    // find face of `neighbord_id` shared with `osh_tet_id`
                    const auto neighbour2neighbours =
                        gather_neighbours<Dim>(tris2tets, tets2tris, neighbor_id);
                    int n2 = 0;
                    for (const auto n2_id: neighbour2neighbours) {
                        if (n2_id == -1) {
                            continue;
                        }
                        if (osh_tet_id == n2_id) {
                            neighborinfo.faces[static_cast<size_t>(neighbor_index)] = n2;
                            break;
                        }
                        ++n2;
                    }
                }
                neighbor_index++;
            }
        }
        neighbors_per_element[osh_tet_id] = neighbor_index;
    };
    neighbors_per_element_ = neighbors_per_element;
    osh::parallel_for(mesh.nents(Omega_h::REGION) /*mesh.nelems()*/, fill_tetInfo);
    tet_neighbors_real_data_.reshape(neighbors_per_element_);
    tet_neighbors_int_data_.reshape(neighbors_per_element_);

    osh::parallel_for(mesh.nelems(), [&](osh::LO elem) {
        const auto& neighbors_info = neighborsInfo[static_cast<size_t>(elem)];
        for (auto i = 0; i < tet_neighbors_real_data_.size(elem); ++i) {
            tet_neighbors_int_data_(elem, i)[0] = neighbors_info.indices[static_cast<size_t>(i)];
            tet_neighbors_int_data_(elem, i)[1] = neighbors_info.faces[static_cast<size_t>(i)];
            tet_neighbors_real_data_(elem, i)[0] = neighbors_info.distances[static_cast<size_t>(i)];
            tet_neighbors_real_data_(elem, i)[1] = neighbors_info.areas[static_cast<size_t>(i)];
        }
    });
    osh::Write<osh::LO> neighbors_per_owned_element_idx(owned_elems_.size());
    osh::parallel_for(neighbors_per_owned_element_idx.size(), [&](osh::LO ownedElemIdx) {
        neighbors_per_owned_element_idx[ownedElemIdx] =
            neighbors_per_element_[owned_elems_[ownedElemIdx].get()];
    });
    neighbors_per_owned_element_idx_ = neighbors_per_owned_element_idx;
}

template <osh::Int Dim>
const TetStruct& OmegaHMesh<Dim>::getTet(PetscInt index) const {
    return tetInfo[static_cast<size_t>(index)];
}

template <osh::Int Dim>
const TriStruct& OmegaHMesh<Dim>::getTri(PetscInt index) const {
    return triInfo[static_cast<size_t>(index)];
}

template <osh::Int Dim>
void OmegaHMesh<Dim>::addCompImpl(const model::compartment_id& compartment,
                                  model::compartment_label cell_set_label) {
    const mesh::compartment_id comp_id = getCompID(compartment);
    bool register_all_elems{false};
    auto class_set = mesh.class_sets.find(compartment);
    if (class_set == mesh.class_sets.end()) {
        class_set = mesh.class_sets.find(std::to_string(cell_set_label));
        if (class_set == mesh.class_sets.end()) {
            register_all_elems = true;
        }
    }

    mesh::element_ids marked;
    PetscScalar p_ownedCompVol{};

    std::vector<osh::LO> v_owned_elems;
    if (register_all_elems) {
        osh::Write<osh::LO> w_marked(mesh.nelems());
        osh::parallel_for(
            mesh.nelems(), OMEGA_H_LAMBDA(auto elem) { w_marked[elem] = elem; });
        marked = w_marked;
    } else {
        const auto& ret = osh::mark_class_closures(&mesh, Dim, class_set->second);
        marked = osh::collect_marked(ret);
    }

    {
        std::vector<osh::LO> elems;
        for (const auto element: marked) {
            elems.push_back(element.get());
            if (isOwned(element)) {
                p_ownedCompVol += measureFunc(element);
                v_owned_elems.push_back(element.get());
            }
            elem2compid[element.get()] = comp_id.get();
        }
        {
            osh::Write<osh::LO> w_elems(static_cast<osh::LO>(elems.size()));
            // parallel std::copy(elems.begin(), elems.end(), w_elems.begin())
            osh::parallel_for(
                static_cast<osh::LO>(elems.size()), OMEGA_H_LAMBDA(auto index) {
                    const auto element_index = elems[static_cast<size_t>(index)];
                    w_elems[index] = element_index;
                });
            complabel2elems.emplace(cell_set_label, w_elems);
        }
    }

    compid2ownedvol.resize(compid2ownedvol.size() + 1, p_ownedCompVol);

    osh::Write<osh::LO> owned_elems(static_cast<osh::Int>(v_owned_elems.size()));
    std::copy(v_owned_elems.begin(), v_owned_elems.end(), owned_elems.begin());
    comp2owned_elems_.resize(comp2owned_elems_.size() + 1, mesh::element_ids(owned_elems));
    measureInfo = std::make_unique<MeasureInfo>(num_compartments(), measureFunc);
}

template <osh::Int Dim>
PetscScalar OmegaHMesh<Dim>::getOwnedCompVol(model::compartment_label comp_label) const {
    auto apiidIt = compLabelToId.find(comp_label);
    if (apiidIt == compLabelToId.end()) {
        throw std::invalid_argument(scat("Unknown compartment label: ", comp_label));
    }
    auto meshId = apicompid2meshcompid.find(apiidIt->second);
    assert(meshId != apicompid2meshcompid.end());
    return compid2ownedvol[static_cast<size_t>(meshId->second.get())];
}

template <osh::Int Dim>
template <class Tag>
std::vector<osh::ClassPair> OmegaHMesh<Dim>::getClassPairs(const strong_string<Tag>& label) const {
    std::vector<osh::ClassPair> class_pairs;
    auto it = mesh.class_sets.find(label);
    if (it != mesh.class_sets.end()) {
        return it->second;
    } else {
        throw std::logic_error("No compartment/patch named " + label);
    }
}

template <osh::Int Dim>
mesh::element_ids OmegaHMesh<Dim>::getEntities(const model::compartment_id& compartmentId) {
    return getEntitiesImpl(compartmentId, false);
}

template <osh::Int Dim>
mesh::boundary_ids OmegaHMesh<Dim>::getEntities(const model::patch_id& patchId) {
    return getEntitiesImpl(patchId, false);
}

template <osh::Int Dim>
mesh::boundary_ids OmegaHMesh<Dim>::getOwnedEntities(const model::patch_id& patch) {
    return getEntitiesImpl(patch, true);
}

template <osh::Int Dim>
mesh::element_ids OmegaHMesh<Dim>::getOwnedEntities(const model::compartment_id& compartment) {
    return getEntitiesImpl(compartment, true);
}

template <osh::Int Dim>
template <typename Tag>
osh::LOs OmegaHMesh<Dim>::getEntitiesImpl(const strong_string<Tag>& region, bool owned) {
    if (region == "__MESH_BOUNDARY__") {
        if (Dim - 1 != entity_dimension<Dim, strong_string<Tag>>::value) {
            throw std::logic_error("Wrong entity dimension for __MESH_BOUNDARY__.");
        }
        auto tri_owned = getMesh().owned(Dim - 1);
        auto tri2tets = getMesh().ask_up(Dim - 1, Dim);
        std::vector<osh::LO> boundary_elems;
        for (osh::LO boundary = 0; boundary < tri_owned.size(); boundary++) {
            if (owned && !tri_owned[boundary]) {
                continue;
            }
            osh::LO noffsets = tri2tets.a2ab[boundary + 1] - tri2tets.a2ab[boundary];
            if (noffsets == 1) {
                boundary_elems.push_back(boundary);
            }
        }
        osh::Write<osh::LO> tri_indices(static_cast<osh::LO>(boundary_elems.size()));
        std::copy(boundary_elems.begin(), boundary_elems.end(), tri_indices.begin());
        return osh::LOs(tri_indices);
    } else if (region == "__MESH__") {
        if (Dim != entity_dimension<Dim, strong_string<Tag>>::value) {
            throw std::logic_error("Wrong entity dimension for __MESH__.");
        }
        if (owned) {
            return getOwnedElems().data();
        } else {
            osh::Write<osh::LO> elements(getMesh().nents(Dim));
            std::iota(elements.begin(), elements.end(), 0);
            return osh::LOs(elements);
        }
    } else {
        auto class_pairs = getClassPairs(region);
        std::vector<osh::Int> dims(class_pairs.size());
        std::transform(class_pairs.begin(), class_pairs.end(), dims.begin(), [](const auto& v) {
            return v.dim;
        });
        std::set<osh::Int> dim_set(dims.begin(), dims.end());
        if (dim_set.empty()) {
            throw std::logic_error("Empty class sets.");
        }
        if (dim_set.size() != 1) {
            throw std::logic_error("Class sets of different dimensions can't be measured");
        }
        auto dim = *dim_set.begin();
        if (dim != entity_dimension<Dim, strong_string<Tag>>::value) {
            throw std::logic_error(
                hadoken::scat("Expecting the dimension of the region to be ", dim));
        }
        auto marks = osh::mark_class_closures(&getMesh(), dim, class_pairs);
        auto owned_elems_mask = getMesh().owned(dim);
        osh::Write<osh::I8> marked_and_owned(marks.size());
        if (owned) {
            for (osh::LO k = 0; k < marks.size(); k++) {
                marked_and_owned[k] = marks[k] & owned_elems_mask[k];
            }
            return osh::collect_marked(marked_and_owned);
        } else {
            return osh::collect_marked(marks);
        }
    }
}

template <osh::Int Dim>
std::tuple<osh::LOs, osh::Read<osh::Real>, osh::Real> OmegaHMesh<Dim>::measure(
    const model::region_id& region) {
    const auto& elems = boost::apply_visitor(
        [this](auto label) -> osh::LOs { return this->getOwnedEntities(label).data(); }, region);
    const auto dim = boost::apply_visitor(
        [](auto label) -> osh::Int { return entity_dimension<Dim, decltype(label)>::value; },
        region);
    osh::Read<osh::Real> ents_measures =
        osh::measure_ents_real(&getMesh(), dim, elems, getMesh().coords());
    return {elems, ents_measures, get_sum(ents_measures)};
}

template <osh::Int Dim>
PetscInt OmegaHMesh<Dim>::getOwnedNumElements() const {
    return static_cast<PetscInt>(owned_elems_.size());
}

namespace {
template <osh::Int Dim>
struct DimToElementType {};

template <>
struct DimToElementType<2> {
    static const DistMesh::ElementType value = DistMesh::ElementType::Tri;
};

template <>
struct DimToElementType<3> {
    static const DistMesh::ElementType value = DistMesh::ElementType::Tet;
};
}  // namespace

template <osh::Int Dim>
DistMesh::ElementType OmegaHMesh<Dim>::getElementType(PetscInt /*point_idx*/) const {
    return DimToElementType<Dim>::value;
}

template <osh::Int Dim>
std::string OmegaHMesh<Dim>::createReport() const {
    return {"Distribute Mesh Partition Report\n"};
}

// explicit instantiation definitions
template class OmegaHMesh<2>;
template std::vector<osh::ClassPair> OmegaHMesh<2>::getClassPairs(
    const model::patch_id& label) const;
template std::vector<osh::ClassPair> OmegaHMesh<2>::getClassPairs(
    const model::compartment_id& label) const;

template class OmegaHMesh<3>;
template std::vector<osh::ClassPair> OmegaHMesh<3>::getClassPairs(
    const model::patch_id& label) const;
template std::vector<osh::ClassPair> OmegaHMesh<3>::getClassPairs(
    const model::compartment_id& label) const;

}  // namespace zee
