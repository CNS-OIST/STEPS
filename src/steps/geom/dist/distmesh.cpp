#include "distmesh.hpp"

#include <limits>

#include <Omega_h_array_ops.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_mark.hpp>
#include <Omega_h_shape.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <easylogging++.h>
#if USE_PETSC
#include <petscsys.h>
#endif // USE_PETSC

#include "distcomp.hpp"
#include "distpatch.hpp"
#include "math/tetrahedron.hpp"
#include "util/mesh.hpp"
#include "util/mpitools.hpp"

namespace steps {
namespace dist {

DistMesh::DistMesh(osh::Mesh mesh, const std::string& path, osh::Real scale)
    : mesh_(mesh)
    , path_(path)
    , scale_(scale)
    , total_num_elems_(mesh.nglobal_ents(dim()))
    , total_num_bounds_(mesh.nglobal_ents(dim() - 1))
    , total_num_verts_(mesh.nglobal_ents(osh::VERT)) {
    if (mesh_.dim() != dim()) {
        throw std::domain_error("Unsupported mesh dimension : " + std::to_string(mesh_.dim()));
    }
#if USE_PETSC
  const auto petsc_max_int = std::numeric_limits<PetscInt>::max();
  if (total_num_elems() > petsc_max_int) {
      std::ostringstream oss;
      oss << "Number of tetrahedrons (" << total_num_elems()
          << ") is greater than what PETSc can support (" << petsc_max_int
          << ").\nYou may recompile PETSc with option --with-64-bit-indices=1";
      throw std::overflow_error(oss.str());
  }
#endif // USE_PETSC
  if (scale != 0) {
    auto coords = osh::deep_copy(mesh_.coords());
    osh::parallel_for(
        coords.size(),
        OMEGA_H_LAMBDA(osh::LO index) { coords[index] *= scale; });
    mesh_.set_coords(coords);
  }
  this->mesh_.set_parting(OMEGA_H_GHOSTED);
  owned_elems_mask_ = this->mesh_.owned(dim());
  owned_bounds_mask_ = this->mesh_.owned(dim() - 1);
  owned_verts_mask_ = this->mesh_.owned(osh::VERT);

  {
    std::vector<osh::LO> owned_elems;
    owned_elems.reserve(static_cast<size_t>(mesh_.nelems()));
    this->elemLocal2Global = mesh_.globals(dim());
    for (osh::LO tet = 0; tet < mesh_.nelems(); ++tet) {
        mesh::tetrahedron_global_id_t global_index{elemLocal2Global[tet]};
        elemGlobal2Local[global_index] = mesh::tetrahedron_local_id_t(tet);
        if (owned_elems_mask_[tet] != 0) {
            owned_elems.push_back(tet);
        }
    }
    osh::Write<osh::LO> owned_elems_array(
        static_cast<osh::LO>(owned_elems.size()));
    std::copy(owned_elems.begin(), owned_elems.end(), owned_elems_array.data());
    owned_elems_ = osh::Read<osh::LO>(owned_elems_array);
  }

  {
    std::vector<osh::LO> owned_bounds;
    owned_bounds.reserve(static_cast<size_t>(mesh_.nfaces()));
    this->boundLocal2Global = mesh_.globals(dim() - 1);
    for (osh::LO tri = 0; tri < mesh_.nfaces(); ++tri) {
        mesh::triangle_global_id_t global_index{boundLocal2Global[tri]};
        boundGlobal2Local[global_index] = mesh::triangle_local_id_t(tri);
        if (owned_bounds_mask_[tri] != 0) {
            owned_bounds.push_back(tri);
        }
    }
    osh::Write<osh::LO> owned_bounds_array(
        static_cast<osh::LO>(owned_bounds.size()));
    std::copy(owned_bounds.begin(), owned_bounds.end(),
              owned_bounds_array.begin());
    owned_bounds_ = owned_bounds_array;
  }

  {
      std::vector<osh::LO> owned_verts;
      owned_verts.reserve(static_cast<size_t>(mesh_.nverts()));
      this->vertLocal2Global = mesh_.globals(osh::VERT);
      for (osh::LO ver = 0; ver < mesh_.nverts(); ++ver) {
          mesh::vertex_global_id_t global_index{vertLocal2Global[ver]};
          vertGlobal2Local[global_index] = mesh::vertex_local_id_t(ver);
          if (owned_verts_mask_[ver] != 0) {
              owned_verts.push_back(ver);
          }
      }
      osh::Write<osh::LO> owned_verts_array(static_cast<osh::LO>(owned_verts.size()));
      std::copy(owned_verts.begin(), owned_verts.end(), owned_verts_array.begin());
      owned_verts_ = owned_verts_array;
  }

  {
      // Compute BBox
      ownedBBoxMin.fill(std::numeric_limits<osh::Real>::max());
      ownedBBoxMax.fill(std::numeric_limits<osh::Real>::min());
      const auto &tets2verts = mesh_.ask_elem_verts();
      for (const auto elem : owned_elems()) {
          const auto tet2verts = osh::gather_verts<4>(tets2verts, elem.get());
          const auto tet2x = osh::gather_vectors<4, dim()>(coords(), tet2verts);
          for (const auto &p : tet2x) {
              for (int i = 0; i < dim() ; ++i) {
                  ownedBBoxMin[i] = std::min(ownedBBoxMin[i], p[i]);
                  ownedBBoxMax[i] = std::max(ownedBBoxMax[i], p[i]);
              }
          }
      }
  }

  const auto &coords = mesh_.coords();
  const auto &areas = osh::measure_ents_real(
      &mesh_, Omega_h::FACE, Omega_h::LOs(mesh_.nents(Omega_h::FACE), 0, 1),
      coords);
  this->fill_triInfo(coords, areas);

  if (mesh_.dim() == 3) {
    this->fill_tetInfo(coords, areas);
    measureFunc_ = [&](mesh::tetrahedron_local_id_t tet) {
      return tetInfo_[static_cast<size_t>(tet.get())].vol;
    };
  } else {
    osh::Write<osh::LO> neighbors(mesh_.nelems());
    osh::parallel_for(
        neighbors.size(), OMEGA_H_LAMBDA(osh::LO elem_id) {
          neighbors[elem_id] =
              triInfo_[static_cast<size_t>(elem_id)].num_neighbors;
        });
    neighbors_per_element_ = neighbors;
    measureFunc_ = [&](mesh::tetrahedron_local_id_t element) {
      return triInfo_[static_cast<size_t>(element.get())].area;
    };
  }
  elem2compid_ = osh::Write<osh::LO>(mesh_.nelems(), INITIAL_COMPARTMENT_ID);
  measure_ = std::make_unique<Measure>(comm_impl(), num_compartments(), measureFunc_);
  diffusion_boundary_ids_.resize(mesh_.nents(dim() - 1), util::nothing);
}

DistMesh::DistMesh(osh::Library &library, const std::string &path,
                   osh::Real scale)
    : DistMesh(DistMesh::load_mesh(library, path), path, scale) {}

void DistMesh::init() {
  {
    // ensure all elements have a dedicated compartment
    std::vector<osh::LO> bad_elements;
    for (osh::LO element{}; element < elem2compid_.size(); ++element) {
      if (elem2compid_[element] == INITIAL_COMPARTMENT_ID) {
        bad_elements.push_back(element);
      }
    }
    if (!bad_elements.empty()) {
      std::cerr << "Error: the following mesh elements miss an associated "
                   "compartment: ";
      std::copy(bad_elements.begin(), bad_elements.end(),
                std::ostream_iterator<osh::LO>(std::cerr, " "));
      std::cerr << "\nAbort\n";
      MPI_Abort(comm_impl(), 1);
    }
  }
  {
    // fill the \a tet_neighbors_in_comp_index_ member variable

    // compute the number of neighbors that belong to the same compartment for
    // every element
    osh::Write<osh::LO> num_neighbors_per_elem_in_comp(mesh_.nelems());
    osh::parallel_for(
        mesh_.nelems(), OMEGA_H_LAMBDA(osh::LO elem) {
          int num_neighbors_in_comp = 0;
          for (const auto &neighbor : tet_neighbors_int_data_[elem]) {
            if (elem2compid_[elem] == elem2compid_[neighbor[0]]) {
              ++num_neighbors_in_comp;
            }
          }
          num_neighbors_per_elem_in_comp[elem] = num_neighbors_in_comp;
        });
    // now resize the container accordingly and fill it
    tet_neighbors_in_comp_index_.reshape(num_neighbors_per_elem_in_comp);
    osh::parallel_for(
        mesh_.nelems(), OMEGA_H_LAMBDA(osh::LO elem) {
          int idx = 0;
          for (int i = 0; i < tet_neighbors_int_data_.size(elem); ++i) {
            const auto neighbor_elem = tet_neighbors_int_data_(elem, i)[0];
            if (elem2compid_[elem] == elem2compid_[neighbor_elem]) {
              tet_neighbors_in_comp_index_(elem, idx++) += i;
            }
          }
        });
  }
  this->measure_->init(this->owned_elems_, this->elem2compid_);
}

int DistMesh::comm_rank() const noexcept {
    return util::mpi_comm_rank(comm_impl());
}

int DistMesh::comm_size() const noexcept {
    return util::mpi_comm_size(comm_impl());
}

osh::Mesh DistMesh::load_mesh(osh::Library &library, const std::string &path) {
  /**
   *  Create Omega_h mesh object from a GMSH mesh
   *  \param filename use parallel import if the prefix (without leading '_')
   *  of a multi-part mesh is given, let Omega_h do the partitioning otherwise.
   *  \return an Omega_h mesh object
   */
  const auto rank0 = library.world()->rank() == 0;
#ifdef OMEGA_H_USE_GMSH
  {
      using namespace Omega_h::filesystem;
      std::ostringstream part_name;
      part_name << path << '_' << library.world()->rank() + 1 << ".msh";
      if (exists(part_name.str())) {
          if (rank0) {
              std::clog << "Creating Omega_h mesh from GMSH mesh partition "
                        << part_name.str() << '\n';
          }
          return Omega_h::gmsh::read_parallel(path, library.world());
      }
  }
#endif
  if (rank0) {
    CLOG(INFO, "general_log")
        << "Creating Omega_h mesh from GMSH mesh " << path << '\n';
  }
  return Omega_h::gmsh::read(path, library.world());
}

void DistMesh::fill_triInfo(const Omega_h::Reals &coords,
                            const Omega_h::Reals &areas) {
  const auto &verts2tris = mesh_.ask_up(1, 2);
  const auto &tris2verts = mesh_.ask_verts_of(Omega_h::FACE);

  auto fill_triInfo = OMEGA_H_LAMBDA(mesh::triangle_local_id_t tri) {
    auto &info = triInfo_[static_cast<size_t>(tri.get())];
    const auto tri2verts = osh::gather_verts<3>(tris2verts, tri.get());
    const auto tri2x = osh::gather_vectors<3, 3>(coords, tri2verts);
    info.centroid = (tri2x[0] + tri2x[1] + tri2x[2]) / 3.;
    info.area = areas[tri.get()];
    for (auto v = 0; v < 3; ++v) {
      auto vid = tris2verts[v];
      info.num_neighbors += static_cast<osh::LO>(
          verts2tris.a2ab[vid + 1] - verts2tris.a2ab[vid] > 1);
    }
  };
  triInfo_.resize(static_cast<size_t>(areas.size()));
  osh::parallel_for(mesh_.nents(Omega_h::FACE), fill_triInfo);

}

/// \brief provide information about the neighbors of a tetrahedron
/// Intermediate structure used to fill the flat_multimap
struct TetNeighborStruct {
    /// indices of the tetrahedron elements
    std::array<osh::LO, 4> indices{{-1, -1, -1, -1}};
    /// face identifiers for the given tetrahedron
    std::array<osh::LO, 4> faces{};
    /// triangles ids
    std::array<osh::LO, 4> triangle_ids{};
    /// distances from the tetrahedron barycenters
    std::array<osh::Real, 4> distances{};
    /// surface of the faces
    std::array<osh::Real, 4> areas{};
};

void DistMesh::fill_tetInfo(const Omega_h::Reals &coords,
                            const Omega_h::Reals &areas) {
  // return the topology from each tetrahedron to its 4 faces as Omega_h::LOs
  const auto graph = mesh_.ask_down(dim(), dim() - 1);
  const auto &tets2tris = graph.ab2b;
  // return for each triangle a data structure to get its neighboring
  // tetrahedrons
  const auto tris2tets = mesh_.ask_up(dim() - 1, dim());
  // for each a = triangle, a2ab gives the offset in ab2b
  const auto &f2fc = tris2tets.a2ab;
  // for each ab = tetrahedron neighboring a triangle, gives its local index
  const auto &fc2c = tris2tets.ab2b;

  const auto volumes = measure_elements_real(&mesh_);
  const auto &tets2verts =
      mesh_.ask_elem_verts(); // ask_verts_of(Omega_h::REGION);
  tetInfo_.resize(static_cast<size_t>(volumes.size()));
  std::vector<TetNeighborStruct> neighborsInfo(
      static_cast<size_t>(volumes.size()));
  osh::Write<osh::LO> neighbors_per_element(volumes.size());
  const auto fill_tetInfo = [&](mesh::tetrahedron_local_id_t tet) {
    const auto tetrahedron_id = static_cast<size_t>(tet.get());
    const auto tet2verts = osh::gather_verts<4>(tets2verts, tet.get());
    const auto tet2x = osh::gather_vectors<4, 3>(coords, tet2verts);
    auto &tetinfo = tetInfo_[tetrahedron_id];
    auto &neighborinfo = neighborsInfo[tetrahedron_id];
    tetinfo.centroid = barycenter(tet2x);
    tetinfo.vol = volumes[tet.get()];
    // returns the local indexes of the four faces of tet
    const auto tet2tris = Omega_h::gather_down<4>(tets2tris, tet.get());
    auto neighbor_index = 0;
    for (auto f = 0; f < 4; ++f) {
      // local ordering of third and fourth neighbour is switched to compare
      // with STEPS
      auto triangle_id = tet2tris[f];
      // the local index of any tetrahedron neighboring tri is stored in
      // fc2c[offset + {0 or 1}]
      auto offset = f2fc[triangle_id];
      if (f2fc[triangle_id + 1] - offset > 1) {
        const auto neighbor_id =
            fc2c[offset + (fc2c[offset] == tet.get() ? 1 : 0)];
        const auto neighbor2verts =
            osh::gather_verts<4>(tets2verts, neighbor_id);
        const auto neighbor2x =
            osh::gather_vectors<4, 3>(coords, neighbor2verts);
        const auto neighbor_centroid = barycenter(neighbor2x);
        neighborinfo.indices[static_cast<size_t>(neighbor_index)] = neighbor_id;
        neighborinfo.triangle_ids[static_cast<size_t>(neighbor_index)] =
            triangle_id;
        neighborinfo.areas[static_cast<size_t>(neighbor_index)] =
            areas[triangle_id];
        neighborinfo.distances[static_cast<size_t>(neighbor_index)] =
            norm(tetinfo.centroid - neighbor_centroid);
        {
          // find face of `neighbord_id` shared with `tet`
          const auto neighbour2neighbours =
              gather_neighbours<dim()>(tris2tets, tets2tris, neighbor_id);
          int n2 = 0;
          for (const auto n2_id : neighbour2neighbours) {
            if (n2_id == -1) {
              continue;
            }
            if (tet.get() == n2_id) {
              neighborinfo.faces[static_cast<size_t>(neighbor_index)] = n2;
              break;
            }
            ++n2;
          }
        }
        neighbor_index++;
      }
    }
    neighbors_per_element[tet.get()] = neighbor_index;
  };
  neighbors_per_element_ = neighbors_per_element;
  osh::parallel_for(mesh_.nents(Omega_h::REGION) /*mesh.nelems()*/,
                    fill_tetInfo);
  tet_neighbors_real_data_.reshape(neighbors_per_element_);
  tet_neighbors_int_data_.reshape(neighbors_per_element_);

  osh::parallel_for(mesh_.nelems(), [&](osh::LO tet) {
    const auto &neighbors_info = neighborsInfo[static_cast<size_t>(tet)];
    for (auto i = 0; i < tet_neighbors_real_data_.size(tet); ++i) {
      tet_neighbors_int_data_(tet, i)[0] =
          neighbors_info.indices[static_cast<size_t>(i)];
      tet_neighbors_int_data_(tet, i)[1] =
          neighbors_info.faces[static_cast<size_t>(i)];
      tet_neighbors_int_data_(tet, i)[2] =
          neighbors_info.triangle_ids[static_cast<size_t>(i)];
      tet_neighbors_real_data_(tet, i)[0] =
          neighbors_info.distances[static_cast<size_t>(i)];
      tet_neighbors_real_data_(tet, i)[1] =
          neighbors_info.areas[static_cast<size_t>(i)];
    }
  });
  osh::Write<osh::LO> neighbors_per_owned_element_idx(owned_elems_.size());
  osh::parallel_for(
      neighbors_per_owned_element_idx.size(), [&](osh::LO ownedElemIdx) {
        neighbors_per_owned_element_idx[ownedElemIdx] =
            neighbors_per_element_[owned_elems_[ownedElemIdx].get()];
      });
  neighbors_per_owned_element_idx_ = neighbors_per_owned_element_idx;
}

void DistMesh::setTetComp(mesh::tetrahedron_local_id_t tet_index,
                          DistComp *compartment) {
  tetInfo_[tet_index.get()].compPtr = compartment;
}

std::vector<DistComp*> DistMesh::getAllComps() const {
    return distcomps;
}

std::vector<DistPatch*> DistMesh::getAllPatches() const {
  return distpatches;
}

DistComp *DistMesh::getTetComp(mesh::tetrahedron_local_id_t tet_index) const {
    assert(tet_index.valid());
    return tetInfo_[tet_index.get()].compPtr;
}

DistComp *DistMesh::getTetComp(mesh::tetrahedron_global_id_t tet_index) const {
    mesh::compartment_id meshCompId(boost::none);
    auto localInd = getLocalIndex(tet_index);
    if (localInd.valid()) {
        auto *comp = getTetComp(localInd);
        if (comp != nullptr) {
            auto it =
                apicompid2meshcompid.find(model::compartment_id(comp->getID()));
            if (it != apicompid2meshcompid.end()) {
                meshCompId = it->second;
            }
        }
    }
    syncData(&meshCompId, 1, MPI_INT32_T, localInd.valid());
    if (meshCompId.valid()) {
        return distcomps[meshCompId.get()];
    } else {
        return nullptr;
    }
}

double DistMesh::getTetVol(mesh::tetrahedron_global_id_t tet_index) const {
    double vol;
    auto localInd = getLocalIndex(tet_index);
    if (localInd.valid()) {
        vol = getTetVol(localInd);
    }
    syncData(&vol, 1, MPI_DOUBLE, localInd.valid());
    return vol;
}

double DistMesh::getTetVol(mesh::tetrahedron_local_id_t tet_index) const {
    assert(tet_index.valid());
    return tetInfo_[tet_index.get()].vol;
}

std::vector<mesh::tetrahedron_global_id_t>
DistMesh::getTetTetNeighb(mesh::tetrahedron_global_id_t tet_index) const {
    std::vector<osh::GO> neighb_inds;
    int nbNeighbs{};
    auto localInd = getLocalIndex(tet_index);
    if (localInd.valid()) {
        for (const auto& tet : getTetTetNeighb(localInd)) {
            neighb_inds.emplace_back(getGlobalIndex(tet));
        }
        nbNeighbs = neighb_inds.size();
    }
    syncData(&nbNeighbs, 1, MPI_INT, localInd.valid());
    neighb_inds.resize(nbNeighbs);
    syncData(neighb_inds.data(), neighb_inds.size(), MPI_INT64_T, localInd.valid());
    return {neighb_inds.begin(), neighb_inds.end()};
}

std::vector<mesh::tetrahedron_local_id_t>
DistMesh::getTetTetNeighb(mesh::tetrahedron_local_id_t tet_index) const {
    assert(tet_index.valid());
    std::vector<mesh::tetrahedron_local_id_t> neighb_inds;
    for (const auto &neighbor : tet_neighbors_int_data_[tet_index.get()]) {
        neighb_inds.emplace_back(neighbor[0]);
    }
    return neighb_inds;
}

std::vector<mesh::triangle_global_id_t>
DistMesh::getTetTriNeighb(mesh::tetrahedron_global_id_t tet_index) {
    auto localInd = getLocalIndex(tet_index);
    std::vector<mesh::triangle_global_id_t> faces;
    if (localInd.valid()) {
        for (auto &tri : getTetTriNeighb(localInd)) {
            faces.emplace_back(getGlobalIndex(tri));
        }
    } else {
        faces.resize(4);
    }
    syncData(faces.data(), 4, MPI_INT64_T, localInd.valid());
    return faces;
}

std::vector<mesh::triangle_local_id_t>
DistMesh::getTetTriNeighb(mesh::tetrahedron_local_id_t tet_index) {
    assert(tet_index.valid());
    const auto graph = mesh_.ask_down(dim(), dim() - 1);
    const auto tet2tris = Omega_h::gather_down<4>(graph.ab2b, tet_index.get());
    return {tet2tris.begin(), tet2tris.end()};
}

std::vector<double>
DistMesh::getTetBarycenter(mesh::tetrahedron_global_id_t tet_index) const {
    auto localInd = getLocalIndex(tet_index);
    std::vector<double> center;
    center.reserve(dim());
    if (localInd.valid()) {
        center = getTetBarycenter(localInd);
    } else {
        center.resize(dim());
    }
    syncData(center.data(), dim(), MPI_DOUBLE, localInd.valid());
    return center;
}

std::vector<double>
DistMesh::getTetBarycenter(mesh::tetrahedron_local_id_t tet_index) const {
    assert(tet_index.valid());
    auto &center = getTet(tet_index).centroid;
    return {center.begin(), center.end()};
}

mesh::tetrahedron_global_id_t
DistMesh::findTetByPoint(const std::vector<double> &position, bool local) {
    steps::math::point3d pos{position[0], position[1], position[2]};
    mesh::tetrahedron_global_id_t tet(boost::none);

    bool inside = true;
    for (int i = 0; inside and i < dim(); ++i) {
        inside = ownedBBoxMin[i] <= pos[i] and pos[i] <= ownedBBoxMax[i];
    }
    if (inside) {
        const auto &tets2verts = mesh_.ask_elem_verts();
        for (const auto elem : owned_elems()) {
            const auto tet2verts = osh::gather_verts<4>(tets2verts, elem.get());
            const auto tet2x =
                osh::gather_vectors<4, dim()>(coords(), tet2verts);
            std::vector<steps::math::point3d> verts;
            for (const auto &p : tet2x) {
                verts.emplace_back(p[0], p[1], p[2]);
            }
            if (steps::math::tet_inside(verts[0], verts[1], verts[2], verts[3],
                                        pos)) {
                tet = mesh::tetrahedron_global_id_t(getGlobalIndex(elem));
                break;
            }
        }
    }

    if (not local) {
        syncData(&tet, 1, MPI_INT64_T, tet.valid());
    }
    return tet;
}

std::vector<mesh::tetrahedron_global_id_t> DistMesh::getAllTetIndices() {
    return getAllEntities<model::compartment_id, mesh::tetrahedron_global_id_t,
                          mesh::tetrahedron_local_id_t>(
        model::compartment_id("__MESH__"));
}

std::vector<mesh::tetrahedron_local_id_t>
DistMesh::getLocalTetIndices(bool owned) {
    auto local_entities =
        getEntitiesImpl(model::compartment_id("__MESH__"), owned);
    return {local_entities.begin(), local_entities.end()};
}

std::vector<mesh::vertex_global_id_t>
DistMesh::getTet_(mesh::tetrahedron_global_id_t tet_index) {
    std::vector<mesh::vertex_global_id_t> verts;
    verts.reserve(4);
    auto localInd = getLocalIndex(tet_index);
    if (localInd.valid()) {
        for (auto vert : getTet_(localInd)) {
            verts.emplace_back(getGlobalIndex(vert));
        }
    } else {
        verts.resize(4, boost::none);
    }
    syncData(verts.data(), verts.size(), MPI_INT64_T, localInd.valid());
    return verts;
}

std::vector<mesh::vertex_local_id_t>
DistMesh::getTet_(mesh::tetrahedron_local_id_t tet_index) {
    const auto tet2verts = osh::gather_verts<4>(mesh_.ask_elem_verts(), tet_index.get());
    return {tet2verts.begin(), tet2verts.end()};
}

double DistMesh::getTriArea(mesh::triangle_global_id_t tri_index) const {
    double area;
    auto localInd = getLocalIndex(tri_index);
    if (localInd.valid()) {
        area = getTriArea(localInd);
    }
    syncData(&area, 1, MPI_DOUBLE, localInd.valid());
    return area;
}

double DistMesh::getTriArea(mesh::triangle_local_id_t tri_index) const {
    assert(tri_index.valid());
    return triInfo_[tri_index.get()].area;
}

std::vector<mesh::vertex_global_id_t>
DistMesh::getTri_(mesh::triangle_global_id_t tri_index) {
    std::vector<mesh::vertex_global_id_t> vertices(3);
    auto localInd = getLocalIndex(tri_index);
    if (localInd.valid()) {
        auto verts = getTri_(localInd);
        for (uint i = 0 ; i < verts.size() ; ++i) {
            vertices[i] = getGlobalIndex(verts[i]);
        }
    }
    syncData(vertices.data(), 3, MPI_INT64_T, localInd.valid());
    return vertices;
}

std::vector<mesh::vertex_local_id_t>
DistMesh::getTri_(mesh::triangle_local_id_t tri_index) {
    std::vector<mesh::vertex_local_id_t> vertices(3);
    const auto& tris2verts = mesh_.ask_verts_of(Omega_h::FACE);
    const auto verts = osh::gather_verts<3>(tris2verts, tri_index.get());
    for (auto i = 0; i < verts.size(); ++i) {
        vertices[i] = verts[i];
    }
    return vertices;
}

std::vector<mesh::triangle_global_id_t> DistMesh::getSurfTris() {
    return getAllEntities<model::patch_id, mesh::triangle_global_id_t,
                          mesh::triangle_local_id_t>(
        model::patch_id("__MESH_BOUNDARY__"));
}

std::vector<mesh::triangle_local_id_t> DistMesh::getSurfLocalTris() {
    auto local_entities =
        getEntitiesImpl(model::patch_id("__MESH_BOUNDARY__"), true);
    return {local_entities.begin(), local_entities.end()};
}

template <typename Tag, typename Global, typename Local>
std::vector<Global> DistMesh::getAllEntities(const Tag &tag) {
    auto local_entities = getEntitiesImpl(tag, true);

    std::vector<Global> global_entities;
    global_entities.reserve(local_entities.size());
    for (auto &ind : local_entities) {
        global_entities.emplace_back(getGlobalIndex(Local(ind)));
    }

    return allGatherEntities(global_entities, MPI_INT64_T);
}

std::vector<mesh::tetrahedron_global_id_t>
DistMesh::getTriTetNeighb(mesh::triangle_global_id_t tri_index) {
    auto localInd = getLocalIndex(tri_index);
    int nbNeighbs{0};
    std::vector<mesh::tetrahedron_global_id_t> neighbs;
    if (localInd.valid()) {
        for (auto &tri: getTriTetNeighb(localInd)) {
            neighbs.emplace_back(getGlobalIndex(tri));
        }
        nbNeighbs = neighbs.size();
    }
    syncData(&nbNeighbs, 1, MPI_INT, localInd.valid());
    neighbs.resize(nbNeighbs);
    syncData(neighbs.data(), neighbs.size(), MPI_INT64_T, localInd.valid());
    return neighbs;
}

std::vector<mesh::tetrahedron_local_id_t>
DistMesh::getTriTetNeighb(mesh::triangle_local_id_t tri_index) {
    assert(tri_index.valid());
    std::vector<mesh::tetrahedron_local_id_t> neighbs;
    const auto tris2tets = mesh_.ask_up(dim() - 1, dim());
    auto offset = tris2tets.a2ab[tri_index.get()];
    auto nind1 = tris2tets.ab2b[offset];
    auto nind2 = tris2tets.ab2b[offset + 1];
    neighbs.emplace_back(nind1);
    if (nind2 != nind1) {
        neighbs.emplace_back(nind2);
    }
    return neighbs;
}

std::vector<double>
DistMesh::getTriBarycenter(mesh::triangle_global_id_t tri_index) const {
    auto localInd = getLocalIndex(tri_index);
    std::vector<double> center;
    if (localInd.valid()) {
        center = getTriBarycenter(localInd);
    } else {
        center.resize(dim());
    }
    syncData(center.data(), dim(), MPI_DOUBLE, localInd.valid());
    return center;
}

std::vector<double>
DistMesh::getTriBarycenter(mesh::triangle_local_id_t tri_index) const {
    auto &centroid = getTri(tri_index).centroid;
    return {centroid.begin(), centroid.end()};
}

std::vector<mesh::triangle_global_id_t> DistMesh::getAllTriIndices() {
    std::vector<mesh::triangle_global_id_t> global_triangles;
    auto triangles = getLocalTriIndices(true);
    global_triangles.reserve(triangles.size());
    for (auto tri: triangles) {
        global_triangles.emplace_back(getGlobalIndex(tri));
    }
    return allGatherEntities(global_triangles, MPI_INT64_T);
}

std::vector<mesh::triangle_local_id_t>
DistMesh::getLocalTriIndices(bool owned) {
    auto bound_owned = owned_bounds_mask_;
    std::vector<mesh::triangle_local_id_t> boundaries_v;
    for (osh::LO boundary = 0; boundary < bound_owned.size(); boundary++) {
        if (bound_owned[boundary] or not owned) {
            boundaries_v.push_back(boundary);
        }
    }
    return boundaries_v;
}

std::vector<double>
DistMesh::getVertex(mesh::vertex_global_id_t vert_index) const {
    std::vector<double> vert(dim());
    auto localInd = getLocalIndex(vert_index);
    if (localInd.valid()) {
        vert = getVertex(localInd);
    }
    syncData(vert.data(), dim(), MPI_DOUBLE, localInd.valid());

    return vert;
}

std::vector<double>
DistMesh::getVertex(mesh::vertex_local_id_t vert_index) const {
    std::vector<double> vert(dim());
    const auto v = osh::gather_vectors<1, dim()>(
        coords(), osh::Few<osh::LO, 1>{vert_index.get()});
    for (auto i = 0; i < v[0].size(); ++i) {
        vert[i] = v[0][i];
    }
    return vert;
}

std::vector<mesh::vertex_global_id_t> DistMesh::getAllVertIndices() {
    std::vector<mesh::vertex_global_id_t> global_vertices;
    auto vertices = getLocalVertIndices(true);
    global_vertices.reserve(vertices.size());
    for (auto vert: vertices) {
        global_vertices.emplace_back(getGlobalIndex(vert));
    }
    return allGatherEntities(global_vertices, MPI_INT64_T);
}

std::vector<mesh::vertex_local_id_t>
DistMesh::getLocalVertIndices(bool owned) {
    auto vert_owned = owned_verts_mask_;
    std::vector<mesh::vertex_local_id_t> vertices_v;
    for (osh::LO vert = 0; vert < vert_owned.size(); vert++) {
        if (vert_owned[vert] or not owned) {
            vertices_v.push_back(vert);
        }
    }
    return vertices_v;
}

std::vector<double> DistMesh::getBoundMin(bool local) const {
    if (local) {
        return {ownedBBoxMin.begin(), ownedBBoxMin.end()};
    } else {
        std::vector<double> minBound(dim());
        MPI_Allreduce(ownedBBoxMin.data(), minBound.data(), dim(), MPI_DOUBLE,
                      MPI_MIN, comm_impl());
        return minBound;
    }
}

std::vector<double> DistMesh::getBoundMax(bool local) const {
    if (local) {
        return {ownedBBoxMax.begin(), ownedBBoxMax.end()};
    } else {
        std::vector<double> maxBound(dim());
        MPI_Allreduce(ownedBBoxMax.data(), maxBound.data(), dim(), MPI_DOUBLE,
                      MPI_MAX, comm_impl());
        return maxBound;
    }
}

void DistMesh::setTriPatch(mesh::triangle_local_id_t tri_index,
                           DistPatch *patch) {
  triInfo_[tri_index.get()].patchPtr = patch;
}

DistPatch *DistMesh::getTriPatch(mesh::triangle_local_id_t tri_index) const {
    return triInfo_[tri_index.get()].patchPtr;
}

DistPatch *DistMesh::getTriPatch(mesh::triangle_global_id_t tri_index) const {
    mesh::patch_id meshPatchId(boost::none);
    auto localInd = getLocalIndex(tri_index);
    if (localInd.valid()) {
        auto *patch = getTriPatch(localInd);
        if (patch != nullptr) {
            auto it =
                apipatchid2meshpatchid.find(model::patch_id(patch->getID()));
            if (it != apipatchid2meshpatchid.end()) {
                meshPatchId = it->second;
            }
        }
    }
    syncData(&meshPatchId, 1, MPI_INT32_T, localInd.valid());
    if (meshPatchId.valid()) {
        return distpatches[meshPatchId.get()];
    } else {
        return nullptr;
    }
}

mesh::tetrahedron_local_id_t
DistMesh::getLocalIndex(mesh::tetrahedron_global_id_t tet, bool owned) const {
  auto result = elemGlobal2Local.find(tet);
  if (result != elemGlobal2Local.end() and
      (not owned or owned_elems_mask_[result->second.get()] != 0)) {
    return result->second;
  }
  return boost::none;
}

osh::GO DistMesh::getGlobalIndex(mesh::tetrahedron_local_id_t tet) const {
  return elemLocal2Global[tet.get()];
}

mesh::triangle_local_id_t
DistMesh::getLocalIndex(mesh::triangle_global_id_t tri, bool owned) const {
  auto result = boundGlobal2Local.find(tri);
  if (result != boundGlobal2Local.end() and
      (not owned or owned_bounds_mask_[result->second.get()] != 0)) {
    return result->second;
  }
  return boost::none;
}

osh::GO DistMesh::getGlobalIndex(mesh::triangle_local_id_t tri) const {
  return boundLocal2Global[tri.get()];
}

mesh::vertex_local_id_t DistMesh::getLocalIndex(mesh::vertex_global_id_t vert,
                                                bool owned) const {
  auto result = vertGlobal2Local.find(vert);
  if (result != vertGlobal2Local.end() and
      (not owned or owned_verts_mask_[result->second.get()] != 0)) {
    return result->second;
  }
  return boost::none;
}

osh::GO DistMesh::getGlobalIndex(mesh::vertex_local_id_t vert) const {
  return vertLocal2Global[vert.get()];
}

mesh::tetrahedron_ids
DistMesh::getOwnedEntities(const model::compartment_id &compartment) {
  const auto mesh_comp_it = apicompid2meshcompid.find(compartment);
  if (mesh_comp_it != apicompid2meshcompid.end()) {
    return comp2owned_elems_[mesh_comp_it->second.get()];
  }
  return getEntitiesImpl(compartment, true);
}

mesh::triangle_ids DistMesh::getOwnedEntities(const model::patch_id &patch) {
    const auto mesh_patch_it = apipatchid2meshpatchid.find(patch);
    if (mesh_patch_it != apipatchid2meshpatchid.end() and
        static_cast<size_t>(mesh_patch_it->second.get()) < patch2owned_bounds_.size()) {
        return patch2owned_bounds_[mesh_patch_it->second.get()];
    }
    return getEntitiesImpl(patch, true);
}

mesh::tetrahedron_ids
DistMesh::getEntities(const model::compartment_id &compartment) {
  const auto mesh_comp_it = apicompid2meshcompid.find(compartment);
  if (mesh_comp_it != apicompid2meshcompid.end()) {
    return compid2elems_[mesh_comp_it->second];
  }
  return getEntitiesImpl(compartment, false);
}

mesh::triangle_ids DistMesh::getEntities(const model::patch_id &patch) {
  const auto mesh_patch_it = apipatchid2meshpatchid.find(patch);
  if (mesh_patch_it != apipatchid2meshpatchid.end()) {
    return patchid2bounds_[mesh_patch_it->second];
  }
  return getEntitiesImpl(patch, false);
}

std::vector<mesh::tetrahedron_global_id_t>
DistMesh::getTaggedTetrahedrons(const model::compartment_id &comp) {
    return
        getAllEntities<model::compartment_id, mesh::tetrahedron_global_id_t,
                       mesh::tetrahedron_local_id_t>(comp);
}

std::vector<mesh::tetrahedron_local_id_t>
DistMesh::getTaggedLocalTetrahedrons(const model::compartment_id &comp, bool owned) {
    auto elems = getEntitiesImpl(comp, owned);
    return {elems.begin(), elems.end()};
}

std::vector<mesh::triangle_global_id_t>
DistMesh::getTaggedTriangles(const model::patch_id &patch) {
    return getAllEntities<model::patch_id, mesh::triangle_global_id_t,
                          mesh::triangle_local_id_t>(patch);
}

std::vector<mesh::triangle_local_id_t>
DistMesh::getTaggedLocalTriangles(const model::patch_id &patch, bool owned) {
    auto elems = getEntitiesImpl(patch, owned);
    return {elems.begin(), elems.end()};
}

std::vector<mesh::vertex_global_id_t>
DistMesh::getTaggedVertices(const model::vertgroup_id &verts) {
    return getAllEntities<model::vertgroup_id, mesh::vertex_global_id_t,
                          mesh::vertex_local_id_t>(verts);
}

std::vector<mesh::vertex_local_id_t>
DistMesh::getTaggedLocalVertices(const model::vertgroup_id &verts, bool owned) {
    auto elems = getEntitiesImpl(verts, owned);
    return {elems.begin(), elems.end()};
}

std::vector<std::string> DistMesh::getTags(osh::Int dim) const {
    std::vector<std::string> tags;
    for (auto cs : mesh_.class_sets) {
        auto &class_pairs = cs.second;
        std::vector<osh::Int> dims(class_pairs.size());
        std::transform(class_pairs.begin(), class_pairs.end(), dims.begin(),
                       [](const auto &v) { return v.dim; });
        std::set<osh::Int> dim_set(dims.begin(), dims.end());
        if (dim_set.size() == 1 and *dim_set.begin() == dim) {
            tags.emplace_back(cs.first);
        }
    }
    return tags;
}

template <typename Tag>
osh::LOs DistMesh::getEntitiesImpl(const util::strong_string<Tag> &region,
                                   bool owned) {
  // for string delimiter
  auto split = [](const auto &str, const std::string &delimiter) {
    std::vector<std::string> res;
    boost::split(res, str, boost::is_any_of(delimiter));
    return res;
  };
  // WARNING!!! Here we strip the .__BOUNDARY__ part (i.e. "smooth" or "spiny")
  // and look for that tag in the mesh
  auto reg_split = split(region, ".");
  if (reg_split.size() == 2 &&  reg_split[1] == "__BOUNDARY__") {
    if (dim() - 1 != entity_dimension<dim(), util::strong_string<Tag>>::value) {
      throw std::logic_error("Wrong entity dimension for " +  region);
    }
    // WARNING!!! notice that we look for the tetrahedrons marked as "smooth" or
    // "spiny", not the surfaces! Such a hack! this "auto ents =
    // getEntities(model::patch_id(reg_split[0]));" would have asked for the
    // surfaces with this tag
    auto ents = getEntities(model::compartment_id(reg_split[0]));
    std::unordered_set<osh::LO> entss;
    for (auto e : ents) {
     entss.insert(e.get());
    }
    auto bound_owned = owned_bounds_mask_;
    auto bound2elems = mesh_.ask_up(dim() - 1, dim());
    std::vector<osh::LO> boundaries_v;
    for (osh::LO boundary = 0; boundary < bound_owned.size(); boundary++) {
      if (owned && !bound_owned[boundary]) {
        continue;
      }
      const auto num_elems =
          bound2elems.a2ab[boundary + 1] - bound2elems.a2ab[boundary];
      if (num_elems == 1) { 
        if (entss.find(bound2elems.ab2b[bound2elems.a2ab[boundary]]) != entss.end()) {
          boundaries_v.push_back(boundary);
        }
      }
    }
    // convert std::vector to Omega_h array
    osh::Write<osh::LO> boundaries(static_cast<osh::LO>(boundaries_v.size()));
    std::copy(boundaries_v.begin(), boundaries_v.end(), boundaries.begin());
    return boundaries;
  } else if (region == "__MESH_BOUNDARY__") {
    if (dim() - 1 != entity_dimension<dim(), util::strong_string<Tag>>::value) {
      throw std::logic_error("Wrong entity dimension for __MESH_BOUNDARY__.");
    }
    auto bound_owned = owned_bounds_mask_;
    auto bound2elems = mesh_.ask_up(dim() - 1, dim());
    std::vector<osh::LO> boundaries_v;
    for (osh::LO boundary = 0; boundary < bound_owned.size(); boundary++) {
      if (owned && !bound_owned[boundary]) {
        continue;
      }
      const auto num_elems =
          bound2elems.a2ab[boundary + 1] - bound2elems.a2ab[boundary];
      if (num_elems == 1) {
        boundaries_v.push_back(boundary);
      }
    }
    // convert std::vector to Omega_h array
    osh::Write<osh::LO> boundaries(static_cast<osh::LO>(boundaries_v.size()));
    std::copy(boundaries_v.begin(), boundaries_v.end(), boundaries.begin());
    return {boundaries};
  } else if (region == "__MESH__") {
    if (dim() != entity_dimension<dim(), util::strong_string<Tag>>::value) {
      throw std::logic_error("Wrong entity dimension for __MESH__.");
    }
    if (owned) {
        return owned_elems().data();
    } else {
        osh::Write<osh::LO> elements(mesh_.nents(dim()));
        std::iota(elements.begin(), elements.end(), 0);
        return osh::LOs(elements);
    }
  } else {
    auto class_pairs = getClassPairs(region);
    std::vector<osh::Int> dims(class_pairs.size());
    std::transform(class_pairs.begin(), class_pairs.end(), dims.begin(),
                   [](const auto &v) { return v.dim; });
    std::set<osh::Int> dim_set(dims.begin(), dims.end());
    if (dim_set.empty()) {
      throw std::logic_error("Empty class sets.");
    }
    if (dim_set.size() != 1) {
      throw std::logic_error(
          "Class sets of different dimensions can't be measured");
    }
    auto dim = *dim_set.begin();
    if (dim !=
        entity_dimension<DistMesh::dim(), util::strong_string<Tag>>::value) {
      std::ostringstream oss;
      oss << "Expecting the dimension of the region to be " << dim;
      throw std::logic_error(oss.str());
    }
    auto marks = osh::mark_class_closures(&mesh_, dim, class_pairs);
    osh::Write<osh::I8> marked_and_owned(marks.size());
    if (owned) {
        const auto& owned_elems_mask = mesh_.owned(dim);
        for (osh::LO k = 0; k < marks.size(); k++) {
            marked_and_owned[k] = marks[k] & owned_elems_mask[k];
      }
      return osh::collect_marked(marked_and_owned);
    } else {
      return osh::collect_marked(marks);
    }
  }
}

std::vector<osh::ClassPair>
DistMesh::getClassPairs(const model::vertgroup_id &verts) const {
  const auto it = mesh_.class_sets.find(verts);
  if (it != mesh_.class_sets.end()) {
    return it->second;
  } else {
    throw std::logic_error("No such vertex group: " + verts);
  }
}

std::vector<osh::ClassPair>
DistMesh::getClassPairs(const model::patch_id &patch) const {
  const auto it = mesh_.class_sets.find(patch);
  if (it != mesh_.class_sets.end()) {
    return it->second;
  } else {
    throw std::logic_error("No such patch: " + patch);
  }
}

std::vector<osh::ClassPair>
DistMesh::getClassPairs(const model::compartment_id &compartment) const {
  auto it = mesh_.class_sets.find(compartment);
  if (it != mesh_.class_sets.end()) {
    return it->second;
  } else {
    const auto label_it = compIdtoLabel.find(compartment);
    if (label_it != compIdtoLabel.end()) {
      it = mesh_.class_sets.find(std::to_string(label_it->second));
      if (it != mesh_.class_sets.end()) {
        return it->second;
      }
    }
    throw std::logic_error("No such compartment: " + compartment);
  }
}

model::compartment_id
DistMesh::getCompartment(mesh::tetrahedron_id_t element) const noexcept {
  const auto mesh_comp_id(elem2compid_[element.get()]);
  return meshcompid2apicompid[static_cast<size_t>(mesh_comp_id)];
}

std::tuple<osh::LOs, osh::Reals, osh::Real>
DistMesh::measure(const model::region_id &region) {
    const auto& entities =
        std::visit([this](auto label) -> osh::LOs { return this->getOwnedEntities(label).data(); },
                   region);
    const auto dim = std::visit(
        [](auto label) -> osh::Int {
            return entity_dimension<DistMesh::dim(), decltype(label)>::value;
        },
        region);
    auto ents_measures = osh::measure_ents_real(&mesh_, dim, entities, mesh_.coords());
    return {entities, ents_measures, get_sum(ents_measures)};
}

mesh::compartment_id
DistMesh::getCompID(const model::compartment_id &compartment) noexcept {
  const auto nextid =
      mesh::compartment_id(static_cast<mesh::compartment_id::value_type>(
          apicompid2meshcompid.size()));
  const auto &status = apicompid2meshcompid.insert({compartment, nextid});
  if (status.second) {
    meshcompid2apicompid.push_back(compartment);
  }
  return status.first->second;
}

mesh::patch_id
DistMesh::getPatchID(const model::patch_id &patch) noexcept {
  const auto nextid =
      mesh::patch_id(static_cast<mesh::patch_id::value_type>(
          apipatchid2meshpatchid.size()));
  const auto &status = apipatchid2meshpatchid.insert({patch, nextid});
  if (status.second) {
    meshpatchid2apipatchid.push_back(patch);
  }
  return status.first->second;
}

osh::Real DistMesh::total_measure(const model::region_id &region) {
  const auto owned_measure = std::get<2>(measure(region));
  osh::Real total_measure{};
  auto err = MPI_Allreduce(&owned_measure, &total_measure, 1, MPI_DOUBLE, MPI_SUM, comm_impl());
  if (err != MPI_SUCCESS) {
      MPI_Abort(comm_impl(), err);
  }
  return total_measure;
}

osh::Real DistMesh::local_measure(const model::region_id &region) {
    return std::get<2>(measure(region));
}

void DistMesh::addComp(const model::compartment_id &compartment,
                       model::compartment_label cell_set_label,
                       DistComp* comp) {
  compIdtoLabel.emplace(compartment, cell_set_label);
  compLabelToId.emplace(cell_set_label, compartment);
  const mesh::compartment_id comp_id = getCompID(compartment);
  // TODO We test comp here because C++ tests do not create DistComps
  if (comp != nullptr) {
      distcomps.push_back(comp);
  }

  bool register_all_elems{false};
  auto class_set = mesh_.class_sets.find(compartment);
  if (class_set == mesh_.class_sets.end()) {
    class_set = mesh_.class_sets.find(std::to_string(cell_set_label));
    if (class_set == mesh_.class_sets.end()) {
      register_all_elems = true;
    }
  }

  mesh::tetrahedron_ids marked;
  osh::Real p_ownedCompVol{};

  std::vector<osh::LO> v_owned_elems;
  if (register_all_elems) {
    osh::Write<osh::LO> w_marked(mesh_.nelems());
    osh::parallel_for(
        mesh_.nelems(), OMEGA_H_LAMBDA(auto elem) { w_marked[elem] = elem; });
    marked = w_marked;
  } else {
    const auto &ret =
        osh::mark_class_closures(&mesh_, dim(), class_set->second);
    // this sync is needed because classSets is a list of the LOCAL entities that need to be
    // registered
    const auto &synced_ret = mesh_.sync_array(dim(), ret, 1);
    marked = osh::collect_marked(synced_ret);
  }

  {
    std::vector<osh::LO> elems;
    for (const auto element : marked) {
      elems.push_back(element.get());
      if (isOwned(element)) {
        p_ownedCompVol += measureFunc_(element);
        v_owned_elems.push_back(element.get());
      }
      elem2compid_[element.get()] = comp_id.get();
    }
    {
      osh::Write<osh::LO> w_elems(static_cast<osh::LO>(elems.size()));
      // parallel std::copy(elems.begin(), elems.end(), w_elems.begin())
      osh::parallel_for(
          static_cast<osh::LO>(elems.size()), OMEGA_H_LAMBDA(auto index) {
            const auto element_index = elems[static_cast<size_t>(index)];
            w_elems[index] = element_index;
          });
      compid2elems_.emplace(comp_id, w_elems);
    }
  }

  compid2ownedvol.resize(compid2ownedvol.size() + 1, p_ownedCompVol);

  osh::Write<osh::LO> owned_elems(static_cast<osh::Int>(v_owned_elems.size()));
  std::copy(v_owned_elems.begin(), v_owned_elems.end(), owned_elems.begin());
  comp2owned_elems_.resize(comp2owned_elems_.size() + 1,
                           mesh::tetrahedron_ids(owned_elems));
  measure_ = std::make_unique<Measure>(comm_impl(), num_compartments(), measureFunc_);
}

void DistMesh::addComp(const model::compartment_id &compartment,
                       const std::vector<mesh::tetrahedron_global_id_t> &tets,
                       DistComp *comp) {
    const mesh::compartment_id comp_id = getCompID(compartment);
    if (comp != nullptr) {
        distcomps.push_back(comp);
    }
    osh::Real p_ownedCompVol{};
    std::vector<osh::LO> v_owned_elems;
    std::vector<mesh::tetrahedron_local_id_t> elems;
    for (const auto &tet : tets) {
        auto localInd = getLocalIndex(tet, false);
        if (localInd.valid()) {
            elems.push_back(localInd);
            mesh::tetrahedron_local_id_t element{localInd};
            if (isOwned(element)) {
                p_ownedCompVol += measureFunc_(element);
                v_owned_elems.push_back(localInd.get());
            }
            elem2compid_[element.get()] = comp_id.get();
        }
    }
    {
      osh::Write<osh::LO> w_elems(static_cast<osh::LO>(elems.size()));
      // parallel std::copy(elems.begin(), elems.end(), w_elems.begin())
      osh::parallel_for(
          static_cast<osh::LO>(elems.size()), OMEGA_H_LAMBDA(auto index) {
            const auto element_index = elems[static_cast<size_t>(index)].get();
            w_elems[index] = element_index;
          });
      compid2elems_.emplace(comp_id, w_elems);
    }

    compid2ownedvol.resize(compid2ownedvol.size() + 1, p_ownedCompVol);

    osh::Write<osh::LO> owned_elems(
        static_cast<osh::Int>(v_owned_elems.size()));
    std::copy(v_owned_elems.begin(), v_owned_elems.end(), owned_elems.begin());
    comp2owned_elems_.resize(comp2owned_elems_.size() + 1,
                             mesh::tetrahedron_ids(owned_elems));
    measure_ = std::make_unique<Measure>(comm_impl(), num_compartments(),
                                         measureFunc_);
}

void DistMesh::addPatch(const model::patch_id &name, DistPatch *patch) {
    const mesh::patch_id patch_id = getPatchID(name);
    apipatchid2meshpatchid[name] = patch_id;
    distpatches.push_back(patch);
}

void DistMesh::addPatch(const model::patch_id &name,
                        const std::vector<mesh::triangle_global_id_t> &tris,
                        DistPatch *patch) {
    const mesh::patch_id patch_id = getPatchID(name);
    distpatches.push_back(patch);
    std::vector<osh::LO> v_owned_bounds;
    std::vector<mesh::triangle_local_id_t> bounds;
    for (const auto &tri : tris) {
        auto localInd = getLocalIndex(tri, false);
        if (localInd.valid()) {
            bounds.push_back(localInd);
            if (isOwned(localInd)) {
                v_owned_bounds.push_back(localInd.get());
            }
        }
    }
    {
      osh::Write<osh::LO> w_bounds(static_cast<osh::LO>(bounds.size()));
      // parallel std::copy(bounds.begin(), bounds.end(), w_bounds.begin())
      osh::parallel_for(
          static_cast<osh::LO>(bounds.size()), OMEGA_H_LAMBDA(auto index) {
            const auto bound_index = bounds[static_cast<size_t>(index)].get();
            w_bounds[index] = bound_index;
          });
      patchid2bounds_.emplace(patch_id, w_bounds);
    }
    osh::Write<osh::LO> owned_bounds(
        static_cast<osh::Int>(v_owned_bounds.size()));
    std::copy(v_owned_bounds.begin(), v_owned_bounds.end(),
              owned_bounds.begin());
    patch2owned_bounds_.resize(patch2owned_bounds_.size() + 1,
                               mesh::triangle_ids(owned_bounds));
}

void DistMesh::addDiffusionBoundary(
    const mesh::diffusion_boundary_name &name,
    const model::compartment_id &comp1, const model::compartment_id &comp2,
    boost::optional<std::set<mesh::triangle_global_id_t>> triangles) {
  auto it = diff_bound_name_2_index_.find(name);
  if (it != diff_bound_name_2_index_.end()) {
    throw std::invalid_argument("A diffusion boundary named " + name +
                                std::string(" already exists"));
  }
  diff_bound_name_2_index_[name] = diffusion_boundaries_.size();
  diffusion_boundaries_.resize(diffusion_boundaries_.size() + 1);
  auto &db = diffusion_boundaries_.back();
  db.mdl_comp1 = comp1;
  db.mdl_comp2 = comp2;
  db.msh_comp1 = getCompID(comp1);
  db.msh_comp2 = getCompID(comp2);
  auto elems1 = compid2elems_[db.msh_comp1];
  auto elems2 = compid2elems_[db.msh_comp2];
  auto g = mesh_.ask_down(dim(), dim() - 1);
  auto bounding_triangles = [&g](const mesh::tetrahedron_ids &elems,
                                 std::set<mesh::triangle_id_t> &s) {
    for (const auto elem_idx : elems) {
      auto tris = osh::gather_down<dim() + 1>(g.ab2b, elem_idx.get());
      s.insert(tris.begin(), tris.end());
    }
  };
  std::set<mesh::triangle_id_t> tris_comp_1, tris_comp_2;
  bounding_triangles(elems1, tris_comp_1);
  bounding_triangles(elems2, tris_comp_2);
  std::vector<mesh::triangle_id_t> intersect_comp_12;
  std::set_intersection(tris_comp_1.begin(), tris_comp_1.end(),
                        tris_comp_2.begin(), tris_comp_2.end(),
                        std::back_inserter(intersect_comp_12));
  if (triangles) {
    db.triangles.reserve(triangles->size());
    for (auto t : *triangles) {
      auto li = getLocalIndex(t, false);
      if (li.valid()) {
          db.triangles.push_back(li);
      }
    }
    std::sort(db.triangles.begin(), db.triangles.end());
    std::vector<mesh::triangle_id_t> intersect_triangles_comp_12;
    intersect_triangles_comp_12.reserve(db.triangles.size());
    std::set_intersection(intersect_comp_12.begin(), intersect_comp_12.end(),
                          db.triangles.begin(), db.triangles.end(),
                          std::back_inserter(intersect_triangles_comp_12));
    if (db.triangles.size() != intersect_triangles_comp_12.size()) {
      throw std::logic_error("Diffusion boundary: some triangles are not part "
                             "of the boundary between the two compartments");
    }
  } else {
    std::swap(db.triangles, intersect_comp_12);
  }
  size_t diffusion_boundary_idx = diffusion_boundaries_.size() - 1;
  for (const auto &t : db.triangles) {
    diffusion_boundary_ids_[static_cast<size_t>(t.get())] =
        diffusion_boundary_idx;
  }
}

void DistMesh::addMembrane(
    const model::membrane_id name,
    DistMemb *memb) {
  if (membranes_.count(name) > 0) {
    throw std::invalid_argument("A membrane named " + name +
                                std::string(" already exists"));
  }
  membranes_[name] = memb;
}

void DistMesh::syncData(void *buff, int count,
                        MPI_Datatype datatype, bool isRoot, bool unknownCount) const {
    int root;
    int local_root = isRoot ? util::mpi_comm_rank(comm_impl()) : 0;
    auto err = MPI_Allreduce(&local_root, &root, 1, MPI_INT, MPI_MAX, comm_impl());
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm_impl(), err);
    }

    if (unknownCount) {
        err = MPI_Bcast(&count, 1, MPI_INT, root, comm_impl());
        if (err != MPI_SUCCESS) {
            MPI_Abort(comm_impl(), err);
        }
    }

    err = MPI_Bcast(buff, count, datatype, root, comm_impl());
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm_impl(), err);
    }
}

template <typename Entity>
std::vector<Entity>
DistMesh::allGatherEntities(const std::vector<Entity> &entities,
                            MPI_Datatype datatype) {
    int local_size = entities.size();
    std::vector<int> sizes(steps::util::mpi_comm_size(comm_impl()));

    auto err = MPI_Allgather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT,
                             comm_impl());
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm_impl(), err);
    }

    std::vector<int> offsets(sizes.size() + 1);
    std::partial_sum(sizes.begin(), sizes.end(), offsets.begin() + 1);

    std::vector<Entity> all_entities(offsets.back());

    err = MPI_Allgatherv(entities.data(), entities.size(), datatype,
                         all_entities.data(), sizes.data(), offsets.data(),
                         datatype, comm_impl());
    if (err != MPI_SUCCESS) {
        MPI_Abort(comm_impl(), err);
    }

    return all_entities;
}

}  // namespace dist
}  // namespace steps
