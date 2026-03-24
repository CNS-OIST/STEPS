#pragma once

#include <string>
#include <unordered_set>

#include "geom/dist/fwd.hpp"
#include "geom/patch.hpp"
#include "util/error.hpp"
#include "util/mesh.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

class DistPatch: public wm::Patch {
  public:
    /**
     * \brief Create a patch to a distributed mesh using the physical tag defined in the mesh file.
     *
     * Add a patch which contains all triangles with the same physical tag in the Gmsh file.
     * If multiple physical tags are provided for the tetrahedron, the value of the first physical
     * tag will be used. When adding a patch, the implementation automatically assigns an internal
     * id for the patch, which starts from 0 and increases by 1 for every added patch.
     *
     * \attention Parallelism: Collective
     *
     * \param patch Name id of the patch.
     * \param mesh The distributed mesh instance
     * \param physical_tag Physical tag of the patch triangles.
     * \param icomp Reference to the inner compartment.
     * \param ocomp Pointer to the outer compartment. Default: nullptr
     */
    DistPatch(const mesh::patch_name& patch,
              DistMesh& mesh,
              DistComp& icomp,
              DistComp* ocomp = nullptr);

    /**
     * \brief Create a patch to a distributed mesh using the physical tag defined in the mesh file.
     *
     * Add a patch which contains all triangles with the same physical tag in the Gmsh file.
     * If multiple physical tags are provided for the tetrahedron, the value of the first physical
     * tag will be used. When adding a patch, the implementation automatically assigns an internal
     * id for the patch, which starts from 0 and increases by 1 for every added patch.
     *
     * \attention Parallelism: Collective
     *
     * \param patch Name id of the patch.
     * \param mesh The distributed mesh instance
     * \param physical_tag Physical tag of the patch triangles.
     * \param icomp Pointer to the inner compartment.
     * \param ocomp Pointer to the outer compartment. Default: nullptr
     */
    DistPatch(const mesh::patch_name& patch,
              DistMesh& mesh,
              mesh::patch_physical_tag physical_tag,
              DistComp& icomp,
              DistComp* ocomp = nullptr);

    DistPatch(const mesh::patch_name& patch,
              DistMesh& mesh,
              std::string tag,
              DistComp& icomp,
              DistComp* ocomp = nullptr);

    /**
     * \brief Create a patch to a distributed mesh using global indices of the triangles.
     *
     * Add a patch to a distributed mesh using the global indices of the triangles.
     * When adding a patch, the implementation automatically assigns an internal id for the
     * patch, which starts from 0 and increases by 1 for every added patch.
     *
     * Note:
     * 1. Each rank has individual global index list containing both triangles owned
     * by the rank, as well as ghost triangles.
     * 2. The global index of a triangle may be different from its index in the mesh file.
     * It may also be different when running the simulation with different number of cores.
     * But the implementation should guarantee its consistancy within the same simluation.
     *
     * \attention Parallelism: Collective
     *
     * \param patch Name id of the patch.
     * \param mesh the distributed mesh instance
     * \param tri_gidxs Vector of global indices of the patch triangles.
     * \param icomp Pointer to the inner compartment.
     * \param ocomp Pointer to the outer compartment. Default: nullptr
     */
    DistPatch(const mesh::patch_name& patch,
              DistMesh& mesh,
              const std::vector<mesh::triangle_global_id_t>& tri_gidxs,
              DistComp& icomp,
              DistComp* ocomp = nullptr);

    DistPatch(const mesh::patch_name& patch,
              DistMesh& mesh,
              const std::vector<mesh::triangle_local_id_t>& tri_lidxs,
              DistComp& icomp,
              DistComp* ocomp = nullptr);

    void init(const mesh::triangle_local_ids& localInds);

    /**
     * \brief Synchronize patch triangles across ranks
     *
     * In case the user created the patch with triangle lists that only include owned
     * triangles, this method returns a list of local triangle indices that are part of
     * the patch, including ghost triangles.
     *
     * \attention Parallelism: Collective
     *
     * \return Vector of local triangle indices
     */
    mesh::triangle_local_ids syncLocalInds(const mesh::triangle_local_ids& localInds) const;

    /**
     * \brief Get the list of all triangle indices of the patch.
     *
     * \attention Parallelism: Collective
     *
     * \return Vector of global triangle indices.
     */
    std::vector<mesh::triangle_global_id_t> getAllTriIndices() const;

    /**
     * \brief Get the list of local triangle indices of the patch.
     *
     * \attention Parallelism: Local
     *
     * \param owned Whether the tetrahedron are owned by the process.
     * \return Vector of global tetrahedron indices.
     */
    std::vector<mesh::triangle_local_id_t> getLocalTriIndices(bool owned = true) const;

    const mesh::triangle_local_ids& getTris(bool owned = true) const noexcept {
        if (owned) {
            return ownedTriLocalIndices;
        } else {
            return triLocalIndices;
        }
    }

    /**
     * \brief Get the area of patch segment owned by the process.
     *
     * Return the sum area of triangles that belong to the patch
     * and owned by the process.
     *
     * \attention Parallelism: Local
     *
     * \return Area of the patch segment owned by the process.
     */
    inline osh::Real getOwnedArea() const noexcept {
        return ownedArea;
    }

    /**
     * \brief Get the area of patch segment across all processes.
     *
     * Return the sum area of triangles that belong to the patch
     * across all processes. This function returns the prestored value,
     * so it can be called locally.
     *
     * \attention Parallelism: Local
     *
     * \return Total area of the patch segment.
     */
    inline osh::Real getTotalArea() const noexcept {
        return getArea();
    }

    /**
     * \brief Get the lower point of the compartment bounding box.
     *
     * If local is true, returns the lower point of the owned bounding box.
     *
     * \attention Parallelism: Collective / Local
     *
     * \return 3D vector of the coordinates of the lower point of the bounding
     * box.
     */
    std::vector<double> getBoundMin(bool local = false) const;

    /**
     * \brief Get the upper point of the compartment bounding box.
     *
     * If local is true, returns the upper point of the owned bounding box.
     *
     * \attention Parallelism: Collective / Local
     *
     * \return 3D vector of the coordinates of the upper point of the bounding
     * box.
     */
    std::vector<double> getBoundMax(bool local = false) const;

    void setMeshID(mesh::patch_id id) {
        assert(not id_.valid());
        id_ = id;
    }

    mesh::patch_id getMeshID() const noexcept {
        assert(id_.valid());
        return id_;
    }

  private:
    mesh::patch_id id_;
    DistMesh& meshRef;
    osh::Real ownedArea{};
    mesh::triangle_local_ids triLocalIndices;
    mesh::triangle_local_ids ownedTriLocalIndices;

    // Bounding box of the owned elements
    std::array<osh::Real, mesh_dimensions()> ownedBBoxMin{};
    std::array<osh::Real, mesh_dimensions()> ownedBBoxMax{};
};

}  // namespace steps::dist
