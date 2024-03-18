#pragma once

#include <string>
#include <unordered_set>

#include "geom/comp.hpp"
#include "geom/dist/distmesh.hpp"
#include "util/mesh.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

class DistComp: public wm::Comp {
  public:
    /**
     * \brief Add a compartment to a distributed mesh using its name as the physical tag defined
     * in the mesh file.
     *
     * Add a compartment which contains all tetrahedrons with the same physical tag in the Gmsh
     * file. If multiple physical tags are provided for the tetrahedron, the value of the first
     * physical tag will be used. When adding a compartment, the implementation automatically
     * assigns an internal id for the compartment, which starts from 0 and increases by 1 for every
     * added compartment.
     *
     * \attention Parallelism: Collective
     *
     * \param compartment Name id of the compartment.
     * \param mesh the distributed mesh instance
     */
    DistComp(const mesh::compartment_name& compartment, DistMesh& mesh, double cond = 0.0);

    /**
     * \brief Add a compartment to a distributed mesh using the physical tag defined in the mesh
     * file.
     *
     * Add a compartment which contains all tetrahedrons with the same physical tag in the Gmsh
     * file. If multiple physical tags are provided for the tetrahedron, the value of the first
     * physical tag will be used. When adding a compartment, the implementation automatically
     * assigns an internal id for the compartment, which starts from 0 and increases by 1 for every
     * added compartment.
     *
     * \attention Parallelism: Collective
     *
     * \param compartment Name id of the compartment.
     * \param mesh the distributed mesh instance
     * \param physical_tag Physical tag of the compartment tetrahedrons.
     */
    DistComp(const mesh::compartment_name& compartment,
             DistMesh& mesh,
             mesh::compartment_physical_tag physical_tag,
             double cond = 0.0);

    DistComp(const mesh::compartment_name& compartment,
             DistMesh& mesh,
             std::string tag,
             double cond);

    /**
     * \brief Add a compartment to a distributed mesh using global indices of the tetrahedrons.
     *
     * Add a compartment to a distributed mesh using the global indices of the tetrahedrons.
     * When adding a compartment, the implementation automatically assigns an internal id for the
     * compartment, which starts from 0 and increases by 1 for every added compartment.
     * Note:
     * 1. Each rank has individual global index list containing both tetrahedrons owned
     * by the rank, as well as ghost tetrahedrons.
     * 2. The global index of a tetrahedron may be different from its index in the mesh file.
     * It may also be different when running the simulation with different number of cores.
     * But the implementation should guarantee its consistency within the same simulation.
     *
     * \attention Parallelism: Collective
     *
     * \param compartment Name id of the compartment.
     * \param mesh the distributed mesh instance
     * \param global_indices Vector of global indices of the compartment tetrahedrons.
     * \param cond Volume conductance of the tetrahedrons in the compartment
     */
    DistComp(const mesh::compartment_name& compartment,
             DistMesh& mesh,
             const std::vector<mesh::tetrahedron_global_id_t>& global_indices,
             double cond = 0.0);

    DistComp(const mesh::compartment_name& compartment,
             DistMesh& mesh,
             const std::vector<mesh::tetrahedron_local_id_t>& local_indices,
             double cond = 0.0);

    /**
     * \brief Get the list of all tetrahedron indices of the compartment.
     *
     * \attention Parallelism: Collective
     *
     * \return Vector of global tetrahedron indices.
     */
    std::vector<mesh::tetrahedron_global_id_t> getAllTetIndices() const;

    /**
     * \brief Get the list of local tetrahedron indices of the compartment.
     *
     * \attention Parallelism: Local
     *
     * \param owned Whether the tetrahedron are owned by the process.
     * \return Vector of global tetrahedron indices.
     */
    const std::vector<mesh::tetrahedron_local_id_t>& getLocalTetIndices(bool owned = true) const;

    /**
     * \brief Get all the triangles on the surface of the compartment.
     *
     * \attention Parallelism: Collective
     *
     * \return Vector of triangle global indices.
     */
    std::vector<mesh::triangle_global_id_t> getSurfTris();

    /**
     * \brief Get the triangles on the surface of the compartment.
     *
     * \attention Parallelism: Local
     *
     * \return Vector of triangle local indices.
     */
    std::vector<mesh::triangle_local_id_t> getSurfLocalTris();

    /**
     * \brief Get the volume of compartment segment owned by the process.
     *
     * Return the sum of volumes of tetrahedrons that belong to the compartment
     * and owned by the process.
     *
     * \attention Parallelism: Local
     *
     * \return Volume of the compartment segment owned by the process.
     */
    inline osh::Real getOwnedVol() const {
        return ownedVol;
    }

    /**
     * \brief Get the volume of compartment segment across all processes.
     *
     * Return the sum of volumes of tetrahedrons that belong to the compartment
     * across all processes. This function only returns the prestored value,
     * so it can be called locally.
     *
     * \attention Parallelism: Local
     *
     * \return Total volume of the compartment segment.
     */
    inline osh::Real getTotalVol() const {
        return pVol;
    }

    /**
     * \brief Return the compartment conductivity.
     *
     * \attention Parallelism: Collective
     */
    inline double getConductivity() const noexcept {
        return pConductivity;
    }

    /**
     * \brief Set the compartment conductivity.
     *
     * \attention Parallelism: Collective
     */
    inline void setConductivity(double cond) noexcept {
        pConductivity = cond;
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

  private:
    /**
     * \brief Add a tetrahedron to the compartment.
     *
     * Add a tetrahedron to the compartment using its local index.
     * This is an internal method.
     *
     * \attention Parallelism: Local
     *
     * \param local_index Local index of the tetrahedron.
     */
    void _addTet(mesh::tetrahedron_local_id_t local_index);

    /**
     * \brief Compute and update the total volume of the compartment.
     *
     * Compute and update the total volume of the compartment across
     * the whole compmunicator.
     *
     * \attention Parallelism: Collective
     *
     */
    void _computeTotalVol();

    /**
     * \brief Compute the owned bounding box of the compartment.
     *
     * \attention Parallelism: Local
     *
     */
    void _computeBBox();

    DistMesh& meshRef;
    osh::Real ownedVol;
    std::vector<mesh::tetrahedron_local_id_t> tetLocalIndices;
    std::vector<mesh::tetrahedron_local_id_t> ownedTetLocalIndices;

    // Bounding box of the owned elements
    std::array<osh::Real, mesh_dimensions()> ownedBBoxMin{};
    std::array<osh::Real, mesh_dimensions()> ownedBBoxMax{};

    double pConductivity;
};

}  // namespace steps::dist
