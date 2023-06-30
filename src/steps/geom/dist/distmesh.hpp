#pragma once

#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <Omega_h_array_ops.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_owners.hpp>

#include "measure.hpp"

#include "geom/dist/fwd.hpp"
#include "geom/geom.hpp"
#include "math/point.hpp"
#include "util/flat_multimap.hpp"
#include "util/mesh.hpp"
#include "util/optional_num.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

using optional_id_t = util::OptionalNum<size_t>;

class DistMesh: public wm::Geom {
  public:
    using point3d = math::point3d;

    DistMesh(osh::Mesh mesh, const std::string& path, osh::Real scale = 0);
    DistMesh(osh::Library& library, const std::string& path, osh::Real scale = 0);

    /// to be called once all compartments, patches, ... have been declared
    void init();

    static constexpr int dim() noexcept {
        return mesh_dimensions();
    }

    static osh::Mesh load_mesh(osh::Library& library, const std::string& path);

    inline std::string getBackend() const {
        return "Omega_h";
    }

    static constexpr bool use_gmsh() noexcept {
#ifdef OMEGA_H_USE_GMSH
        return true;
#endif
        return false;
    }

    /******************************* Compartment *********************************/

    /**
     * \brief Get the volume of compartment segment owned by the process.
     *
     * Return the sum of volumes of tetrahedrons that belong to the compartment
     * and owned by the process.
     *
     * \attention Parallelism: Local
     *
     * \param compartment Name id of the compartment.
     * \return Volume of the compartment segment owned by the process.
     */
    osh::Real getCompOwnedVol(const mesh::compartment_name& compartment) const;

    /**
     * \brief Get the volume of compartment segment across all processes.
     *
     * Return the sum of volumes of tetrahedrons that belong to the compartment
     * across all processes.
     *
     * \attention Parallelism: Collective
     *
     * \param compartment Name id of the compartment.
     * \return Total volume of the compartment segment.
     */
    osh::Real getCompTotalVol(const mesh::compartment_name& compartment) const;

    /**
     * \brief Get a vector of pointers to all compartments.
     *
     * Get a vector of pointers to all compartments.
     *
     * \attention Parallelism: Local
     *
     * \return A vector of pointers to the compartments.
     */
    std::vector<DistComp*> getAllComps() const;

    /********************************** Tet ************************************/

    /**
     * \brief Set the compartment of a tetrahedron.
     *
     * Set the compartment of a tetrahedron with its local index.
     * This function should only be called by the DistComp class.
     *
     * \attention Parallelism: Local
     *
     * \param tet_index Local index of the tetrahedron.
     * \param compartment Pointer to the DistComp object.
     */
    void setTetComp(mesh::tetrahedron_local_id_t tet_index, DistComp* compartment);

    /**
     * \brief Get the compartment of a tetrahedron.
     *
     * Get the compartment of a tetrahedron with its local index.
     * This function returns nullptr if the tetrahedron has not been assigned to
     * any compartment.
     * This function should only be called by the DistComp or the DistPatch class.
     *
     * \attention Parallelism: Local
     *
     * \param tet_index Local index of the tetrahedron.
     * \return Pointer to the DistComp object in which the tetrahedron is.
     */
    DistComp* getTetComp(mesh::tetrahedron_local_id_t tet_index) const;

    /**
     * \brief Get the compartment of a tetrahedron.
     *
     * Get the compartment of a tetrahedron with its global index.
     * This function returns nullptr if the tetrahedron has not been assigned to
     * any compartment.
     *
     * \attention Parallelism: Collective
     *
     * \param tet_index Global index of the tetrahedron.
     * \return Pointer to the DistComp object in which the tetrahedron is.
     */
    DistComp* getTetComp(mesh::tetrahedron_global_id_t tet_index) const;

    /**
     * \brief Get the volume of a tetrahedron.
     *
     * \attention Parallelism: Collective
     *
     * \param tet_index Global index of the tetrahedron.
     * \return Volume of the tetrahedron.
     */
    double getTetVol(mesh::tetrahedron_global_id_t tet_index) const;

    /**
     * \brief Get the volume of a tetrahedron.
     *
     * \attention Parallelism: Local
     *
     * \param tet_index Local index of the tetrahedron.
     * \return Volume of the tetrahedron.
     */
    double getTetVol(mesh::tetrahedron_local_id_t tet_index) const;

    /**
     * \brief Get the tetrahedron neighbors of a tetrahedron.
     *
     * \attention Parallelism: Collective
     *
     * \param tet_index Global index of the tetrahedron.
     * \return Vector of tetrahedron global indexes of the neighbors.
     */
    std::vector<mesh::tetrahedron_global_id_t> getTetTetNeighb(
        mesh::tetrahedron_global_id_t tet_index) const;

    /**
     * \brief Get the tetrahedron neighbors of a tetrahedron.
     *
     * \attention Parallelism: Local
     *
     * \param tet_index Local index of the tetrahedron.
     * \param owned If true, only owned tetrahedrons will be returned
     * \return Vector of tetrahedron local indexes of the neighbors.
     */
    std::vector<mesh::tetrahedron_local_id_t> getTetTetNeighb(
        mesh::tetrahedron_local_id_t tet_index,
        bool owned) const;

    /**
     * \brief Get the triangle neighbors (faces) of a tetrahedron.
     *
     * \attention Parallelism: Collective
     *
     * \param tet_index Global index of the tetrahedron.
     * \return Vector of triangle global indexes.
     */
    std::vector<mesh::triangle_global_id_t> getTetTriNeighb(
        mesh::tetrahedron_global_id_t tet_index);

    /**
     * \brief Get the triangle neighbors (faces) of a tetrahedron.
     *
     * \attention Parallelism: Local
     *
     * \param tet_index Local index of the tetrahedron.
     * \return Vector of triangle local indexes, potentially non-owned.
     */
    std::vector<mesh::triangle_local_id_t> getTetTriNeighb(mesh::tetrahedron_local_id_t tet_index);


    /**
     * \brief Get the barycenter of a tetrahedron.
     *
     * \attention Parallelism: Collective
     *
     * \param tet_index Global index of the tetrahedron.
     * \return 3D vector of barycenter coordinates.
     */
    std::vector<double> getTetBarycenter(mesh::tetrahedron_global_id_t tet_index) const;

    /**
     * \brief Get the barycenter of a tetrahedron.
     *
     * \attention Parallelism: Local
     *
     * \param tet_index Local index of the tetrahedron.
     * \return 3D vector of barycenter coordinates.
     */
    std::vector<double> getTetBarycenter(mesh::tetrahedron_local_id_t tet_index) const;

    /**
     * \brief Get the tetrahedron that contains a given point.
     *
     * \attention Parallelism: Collective
     *
     * \param position 3D position of the point.
     * \return Global tetrahedron id, invalid if the tetrahedron does not exist.
     */
    mesh::tetrahedron_global_id_t findTetByPoint(const std::vector<double>& position);

    /**
     * \brief Get the tetrahedron that contains a given point.
     *
     * \attention Parallelism: Local
     *
     * \param position 3D position of the point.
     * \return Local tetrahedron id, invalid if no local tetrahedron contains the point
     */
    mesh::tetrahedron_local_id_t findLocalTetByPoint(const point3d& position);
    mesh::tetrahedron_local_id_t findLocalTetByPoint(const std::vector<double>& position);

    /**
     * \brief Check if point belongs to a tetrahedron or not
     *
     * \attention Parallelism: Collective
     *
     * \param pos 3D position of the point.
     * \param tet_id Global tetrahedron id.
     * \return True if point in tet, otherwise False.
     */
    bool isPointInTet(const std::vector<double>& pos, mesh::tetrahedron_global_id_t tet_id);

    /**
     * \brief Check if point belongs to a local tetrahedron or not
     *
     * \attention Parallelism: Local
     *
     * \param pos 3D position of the point.
     * \param tet_id Local tetrahedron id.
     * \return True if point in local tet, otherwise False.
     */
    bool isPointInTet(const std::vector<double>& pos, mesh::tetrahedron_local_id_t tet_id);

    /**
     * \brief Get the list of all tetrahedron indices.
     *
     * \attention Parallelism: Collective
     *
     * \return Vector of global tetrahedron indices.
     */
    std::vector<mesh::tetrahedron_global_id_t> getAllTetIndices();

    /**
     * \brief Get the list of local tetrahedron indices.
     *
     * \attention Parallelism: Collective
     *
     * \param owned Whether the tetrahedron are owned by the process.
     * \return Vector of global tetrahedron indices.
     */
    std::vector<mesh::tetrahedron_local_id_t> getLocalTetIndices(bool owned);

    /**
     * \brief Get the list of vertex indices for the given tetrahedron.
     *
     * \attention Parallelism: Collective
     *
     * \param tet_index Global tetrahedron index.
     * \return Vector of global vertex indices.
     */
    std::vector<mesh::vertex_global_id_t> getTet_(mesh::tetrahedron_global_id_t tet_index);

    /**
     * \brief Get the list of vertex indices for the given tetrahedron.
     *
     * \attention Parallelism: Local
     *
     * \param tet_index Local tetrahedron index.
     * \return Vector of local vertex indices.
     */
    std::vector<mesh::vertex_local_id_t> getTet_(mesh::tetrahedron_local_id_t tet_index);

    void addComp(const model::compartment_id& comp_id,
                 model::compartment_label cell_set_label,
                 DistComp* comp = nullptr);

    void addComp(const model::compartment_id& comp_id,
                 const std::vector<mesh::tetrahedron_global_id_t>& tets,
                 DistComp* comp = nullptr);

    void addComp(const model::compartment_id& comp_id,
                 const std::vector<mesh::tetrahedron_local_id_t>& tets,
                 DistComp* comp = nullptr);

    /********************************** Patch ************************************/

    /**
     * \brief Set the patch of a triangle.
     *
     * Set the patch of a triangle with its local index.
     * This function should only be called by the DistPatch class.
     *
     * \attention Parallelism: Local
     *
     * \param tri_index Local index of the triangle.
     * \param patch Pointer to the DistPatch object.
     */
    void setTriPatch(mesh::triangle_local_id_t tri_index, DistPatch* patch);

    /**
     * \brief Get the patch of a triangle.
     *
     * Get the patch of a triangle with its local index.
     * This function returns nullptr if the triangle has not been assigned to
     * any patch.
     * This function should only be called by the DistComp or the DistPatch class.
     *
     * \attention Parallelism: Local
     *
     * \param tri_index Local index of the triangle.
     * \return Pointer to the DistPatch object in which the triangle is.
     */
    DistPatch* getTriPatch(mesh::triangle_local_id_t tri_index) const;

    /**
     * \brief Get the patch of a triangle.
     *
     * Get the patch of a triangle with its global index.
     * This function returns nullptr if the triangle has not been assigned to
     * any patch.
     *
     * \attention Parallelism: Collective
     *
     * \param tri_index Global index of the triangle.
     * \return Pointer to the DistPatch object in which the triangle is.
     */
    DistPatch* getTriPatch(mesh::triangle_global_id_t tri_index) const;

    /**
     * \brief Get the area of patch segment owned by the process.
     *
     * Return the sum area of triangles that belong to the patch
     * and owned by the process.
     *
     * \attention Parallelism: Local
     *
     * \param patch Name id of the patch.
     * \return Area of the patch segment owned by the process.
     */
    osh::Real getPatchOwnedArea(const mesh::patch_name& patch);

    /**
     * \brief Get the area of patch segment across all processes.
     *
     * Return the sum area of triangles that belong to the patch
     * across all processes.
     *
     * \attention Parallelism: Collective
     *
     * \param patch Name id of the patch.
     * \return Total area of the patch segment.
     */
    osh::Real getPatchTotalArea(const mesh::patch_name& patch) const;

    /**
     * \brief Get a vector of pointers to all patches.
     *
     * Get a vector of pointers to all patches.
     *
     * \attention Parallelism: Local
     *
     * \return A vector of points to the patches.
     */
    std::vector<DistPatch*> getAllPatches() const;

    /********************************** Tri ************************************/

    /**
     * \brief Get the area of a triangle.
     *
     * \attention Parallelism: Collective
     *
     * \param tri_index Global index of the triangle.
     * \return Area of the triangle.
     */
    double getTriArea(mesh::triangle_global_id_t tri_index) const;

    /**
     * \brief Get the area of a triangle.
     *
     * \attention Parallelism: Local
     *
     * \param tri_index Local index of the triangle.
     * \return Area of the triangle.
     */
    double getTriArea(mesh::triangle_local_id_t tri_index) const;

    // Need to temporarily rename to getTri_ because there is another getTri method
    // that does not return a list of vertices.
    // TODO Rename the other one to something else
    std::vector<mesh::vertex_global_id_t> getTri_(mesh::triangle_global_id_t tri_index);
    std::vector<mesh::vertex_local_id_t> getTri_(mesh::triangle_local_id_t tri_index);

    /**
     * \brief Get all the triangles on the surface of the mesh.
     *
     * \attention Parallelism: Collective
     *
     * \return Vector of triangle global indices.
     */
    std::vector<mesh::triangle_global_id_t> getSurfTris();

    /**
     * \brief Get the triangles on the surface of the mesh.
     *
     * \attention Parallelism: Local
     *
     * \param owned Whether the triangles should be owned.
     * \return Vector of triangle local indices.
     */
    std::vector<mesh::triangle_local_id_t> getSurfLocalTris(bool owned);

    /**
     * \brief Get the tetrahedron neighbors of a triangle.
     *
     * \attention Parallelism: Collective
     *
     * \param tet_index Global index of the triangle.
     * \return Vector of tetraheron global indexes.
     */
    std::vector<mesh::tetrahedron_global_id_t> getTriTetNeighb(
        mesh::triangle_global_id_t tri_index);

    /**
     * \brief Get the tetrahedron neighbors of a triangle.
     *
     * \attention Parallelism: Local
     *
     * \param tet_index Local index of the triangle.
     * \param owned If true, only owned tetrahedrons will be returned
     * \return Vector of tetraheron local indexes.
     */
    std::vector<mesh::tetrahedron_local_id_t> getTriTetNeighb(mesh::triangle_local_id_t tri_index,
                                                              bool owned);

    /**
     * \brief Get the barycenter of a triangle.
     *
     * \attention Parallelism: Collective
     *
     * \param tet_index Global index of the triangle.
     * \return 3D vector of barycenter coordinates.
     */
    std::vector<double> getTriBarycenter(mesh::triangle_global_id_t tri_index) const;

    /**
     * \brief Get the barycenter of a triangle.
     *
     * \attention Parallelism: Local
     *
     * \param tet_index Local index of the triangle.
     * \return 3D vector of barycenter coordinates.
     */
    std::vector<double> getTriBarycenter(mesh::triangle_local_id_t tri_index) const;

    /**
     * \brief Get the list of all mesh triangle indices.
     *
     * \attention Parallelism: Collective
     *
     * \return Vector of global triangle indices.
     */
    std::vector<mesh::triangle_global_id_t> getAllTriIndices();

    /**
     * \brief Get the list of local triangle indices.
     *
     * \attention Parallelism: Collective
     *
     * \param owned Whether the triangles should be owned by the process.
     * \return Vector of local triangle indices.
     */
    std::vector<mesh::triangle_local_id_t> getLocalTriIndices(bool owned);

    void addPatch(const model::patch_id& name, DistPatch* patch);

    void addPatch(const model::patch_id& name,
                  const std::vector<mesh::triangle_global_id_t>& tris,
                  DistPatch* patch);

    void addPatch(const model::patch_id& name,
                  const std::vector<mesh::triangle_local_id_t>& tris,
                  DistPatch* patch);

    /********************************** Vert **********************************/

    std::vector<double> getVertex(mesh::vertex_global_id_t vert_index) const;
    std::vector<double> getVertex(mesh::vertex_local_id_t vert_index) const;

    /**
     * \brief Get the list of all mesh vertex indices.
     *
     * \attention Parallelism: Collective
     *
     * \return Vector of global vertex indices.
     */
    std::vector<mesh::vertex_global_id_t> getAllVertIndices();

    /**
     * \brief Get the list of local vertex indices.
     *
     * \attention Parallelism: Collective
     *
     * \param owned Whether the vertices should be owned by the process.
     * \return Vector of local vertex indices.
     */
    std::vector<mesh::vertex_local_id_t> getLocalVertIndices(bool owned);

    /********************************** Shape *********************************/

    /**
     * \brief Get the lower point of the mesh bounding box.
     *
     * If local is true, returns the lower point of the owned bounding box.
     *
     * \attention Parallelism: Collective / Local
     *
     * \return 3D vector of the coordinates of the lower point of the bounding box.
     */
    std::vector<double> getBoundMin(bool local) const;

    /**
     * \brief Get the upper point of the mesh bounding box.
     *
     * If local is true, returns the upper point of the owned bounding box.
     *
     * \attention Parallelism: Collective / Local
     *
     * \return 3D vector of the coordinates of the upper point of the bounding box.
     */
    std::vector<double> getBoundMax(bool local) const;

    /********************************** Universal *****************************/
    /**
     * \brief Get the global index of a mesh element from its local index.
     *
     * Get the global index of a mesh element from its local index.
     * A mesh element can be a vertex/bar/tri/tet.
     *
     * \attention Parallelism: Local
     *
     * \param index Local index of the mesh element.
     * \return Global index of the mesh element.
     */
    mesh::tetrahedron_global_id_t getGlobalIndex(mesh::tetrahedron_local_id_t index) const;

    /**
     * \brief Get the local index of a mesh element from its global index.
     *
     * Get the local index of a tetrahedron from its global index.
     * A mesh element can be a vertex/bar/tri/tet.
     *
     * \attention Parallelism: Local
     *
     * \param index Global index of the mesh element.
     * \return Local index of the mesh element.
     */
    mesh::tetrahedron_local_id_t getLocalIndex(mesh::tetrahedron_global_id_t index,
                                               bool owned) const;

    /**
     * \brief Get the global index of a mesh element from its local index.
     *
     * Get the global index of a mesh element from its local index.
     * A mesh element can be a vertex/bar/tri/tet.
     *
     * \attention Parallelism: Local
     *
     * \param index Local index of the mesh element.
     * \return Global index of the mesh element.
     */
    mesh::triangle_global_id_t getGlobalIndex(mesh::triangle_local_id_t index) const;

    /**
     * \brief Get the local index of a mesh element from its global index.
     *
     * Get the local index of a tetrahedron from its global index.
     * A mesh element can be a vertex/bar/tri/tet.
     *
     * \attention Parallelism: Local
     *
     * \param index Global index of the mesh element.
     * \return Local index of the mesh element.
     */
    mesh::triangle_local_id_t getLocalIndex(mesh::triangle_global_id_t index, bool owned) const;

    mesh::vertex_global_id_t getGlobalIndex(mesh::vertex_local_id_t index) const;
    mesh::vertex_local_id_t getLocalIndex(mesh::vertex_global_id_t index, bool owned) const;

    /// Identifiers of the elements owned by this process
    inline const mesh::tetrahedron_ids& owned_elems() const noexcept {
        return owned_elems_;
    }

    /// Owned and ghost elements managed by this process
    inline const osh::Bytes& owned_elems_mask() const noexcept {
        return owned_elems_mask_;
    }

    /// Number of elements in the entire mesh
    osh::GO total_num_elems() const {
        return total_num_elems_;
    };

    /// Number of elements owned by this process
    inline osh::LO num_elems() const noexcept {
        return owned_elems_.size();
    }

    /// Identifiers of the boundaries owned by this process
    inline const osh::LOs& owned_bounds() const noexcept {
        return owned_bounds_;
    }

    /// Owned and ghost boundaries managed by this process
    inline const osh::Bytes& owned_bounds_mask() const noexcept {
        return owned_bounds_mask_;
    }

    /// Number of boundaries in the entire mesh
    inline osh::GO total_num_bounds() const noexcept {
        return total_num_bounds_;
    };

    /// Number of boundaries owned by this process
    inline osh::LO num_bounds() const noexcept {
        return owned_bounds_.size();
    }

    /// Identifiers of the vertices owned by this process
    inline const osh::LOs& owned_verts() const noexcept {
        return owned_verts_;
    }

    /// Owned and ghost vertices managed by this process
    inline const osh::Bytes& owned_verts_mask() const noexcept {
        return owned_verts_mask_;
    }

    /// Number of vertices in the entire mesh
    osh::GO total_num_verts() const noexcept {
        return total_num_verts_;
    };

    /// Number of vertices owned by this process
    inline osh::LO num_verts() const noexcept {
        return owned_verts_.size();
    }

    /// Create new dist for variable sized
    inline osh::Dist create_dist_for_variable_sized(osh::Int dims, osh::LOs copies2data) {
        return osh::create_dist_for_variable_sized(mesh_.ask_dist(dims), copies2data);
    }

    /// Syn the array
    template <typename T>
    osh::Read<T> sync_array(osh::Int ent_dim, osh::Read<T> a, osh::Int width) {
        return mesh_.sync_array(ent_dim, a, width);
    }

    /// Mapping a2ab from objects of dimension \p from to objects \p to (\f$ from < to \f$)
    inline osh::LOs bounds2elems_a2ab(osh::Int from, osh::Int to) {
        return mesh_.ask_up(from, to).a2ab;
    }

    /// Mapping ab2b from objects of dimension \p from to objects \p to (\f$ from < to \f$)
    inline osh::LOs bounds2elems_ab2b(osh::Int from, osh::Int to) {
        return mesh_.ask_up(from, to).ab2b;
    }

    /// ids of vertices of object of dimension \p dim
    inline osh::LOs ask_verts_of(osh::Int dim) {
        return mesh_.ask_verts_of(dim);
    }

    /// ask elem -> vert mapping
    inline osh::LOs ask_elem_verts() {
        return mesh_.ask_elem_verts();
    }

    /// Global indices
    inline const osh::GOs global_indices(osh::Int dim) {
        return osh::globals_from_owners(&mesh_, dim);
    }

    /// Vertices coordinates
    inline osh::Reals coords() const {
        return mesh_.coords();
    }

    /// Communicator
    inline MPI_Comm comm_impl() const noexcept {
        return mesh_.comm()->get_impl();
    }

    /// Communicator rank
    int comm_rank() const noexcept;

    /// Communicator size
    int comm_size() const noexcept;

    /// Print communicator stats for debugging
    std::string inline comm_stats() const {
        std::stringstream ss;
        ss << "rank/size: " << comm_rank() << '/' << comm_size();
        return ss.str();
    }

    /// Max over the various MPI processes
    template <typename T>
    inline T get_MPI_max(const osh::Read<T>& a) const {
        return osh::get_max(mesh_.comm(), a);
    }

    /// Min over the various MPI processes
    template <typename T>
    inline T get_MPI_min(const osh::Read<T>& a) const {
        return osh::get_min(mesh_.comm(), a);
    }

    /// Sum over the various MPI processes
    template <typename T>
    inline osh::promoted_t<T> get_MPI_sum(const osh::Read<T>& a) const {
        return osh::get_sum(mesh_.comm(), a);
    }

    /**
     * \brief Expose floating point data associated to tetrahedrons neighbors.
     * \FIXME(TCL) this should not be exposed directly
     */
    inline const util::flat_multimap<osh::Real, 2>& tet_neighbors_real_data() const noexcept {
        return tet_neighbors_real_data_;
    }

    /**
     * \brief Expose integral data associated to tetrahedrons neighbors.
     * \FIXME(TCL) this should not be exposed directly
     */
    inline const util::flat_multimap<osh::LO, 3>& tet_neighbors_int_data() const noexcept {
        return tet_neighbors_int_data_;
    }

    /**
     * Index of the tetrahedrons that belong to the same compartment of a given
     * element in \a tet_neighbors_real_data and \a tet_neighbors_int_data
     */
    inline const util::flat_multimap<osh::LO, 1>& tet_neighbors_in_comp_index() const noexcept {
        return tet_neighbors_in_comp_index_;
    }

    /**
     * \brief Get the sum of the measure (volume for compartments and area
     * for patches) across all processes.
     * \param region patch or compartment
     */
    osh::Real total_measure(const model::region_id& region);

    /**
     * \brief Get the the measure (volume for compartments and area
     * for patches) of a region.
     * \param region patch or compartment
     */
    osh::Real local_measure(const model::region_id& region);

    /**
     * \brief Check if a point is owned by this process.
     *
     * Giving the point index of a mesh element, this function returns
     * if the element is owned by the process, i.e. not ghost point.
     * \param point_idx Point index of the element.
     * \return True or 't' if the element is owned by the process. False or 'f' if it is not.
     */
    inline bool isOwned(mesh::tetrahedron_local_id_t point_idx) const noexcept {
        assert(point_idx.get() < owned_elems_mask_.size());
        return owned_elems_mask_[point_idx.get()] != 0;
    }
    inline bool isOwned(mesh::triangle_local_id_t point_idx) const noexcept {
        assert(point_idx.get() < owned_bounds_mask_.size());
        return owned_bounds_mask_[point_idx.get()] != 0;
    }

    inline bool isOwned(mesh::vertex_local_id_t point_idx) const noexcept {
        assert(point_idx.get() < owned_verts_mask_.size());
        return owned_verts_mask_[point_idx.get()] != 0;
    }

    inline osh::Int num_compartments() const noexcept {
        return compid2elems_.size();
    }

    inline const Measure& getMeasure() const noexcept {
        return *measure_;
    }

    /**
     *
     * @param region patch or compartment identifier
     * @return tuple containing all elements the measure of the region
     * \TODO THIS METHOD SHOULD BE CONST
     */
    std::tuple<osh::LOs, osh::Reals, osh::Real> measure(const model::region_id& region);

    struct TetStruct {
        DistComp* compPtr{nullptr};
        osh::Real vol;
        osh::Vector<3> centroid;
    };

    struct TriStruct {
        DistPatch* patchPtr{nullptr};
        osh::LO num_neighbors;
        osh::Real area;
        osh::Vector<3> centroid;
    };

    inline const std::vector<TetStruct>& getTetInfo() const noexcept {
        return tetInfo_;
    }

    inline const std::vector<TriStruct>& getTriInfo() const noexcept {
        return triInfo_;
    }

    inline const TriStruct& getTri(mesh::triangle_id_t id) const noexcept {
        return triInfo_[static_cast<size_t>(id.get())];
    }

    inline const TetStruct& getTet(mesh::tetrahedron_id_t id) const noexcept {
        return tetInfo_[static_cast<size_t>(id.get())];
    }

    const auto& getClassSets() const noexcept {
        return mesh_.class_sets;
    }

    mesh::compartment_id getCompID(const model::compartment_id& compartment) noexcept;

    mesh::patch_id getPatchID(const model::patch_id& patch) noexcept;

    model::compartment_id getCompartment(mesh::tetrahedron_id_t element) const noexcept;

    mesh::compartment_id getCompartmentMeshID(mesh::tetrahedron_id_t element) const noexcept {
        return mesh::compartment_id(elem2compid_[element.get()]);
    }
    /**
     * \param patch patch label
     * \return boundaries owned by this process that belong to the given patch
     */
    mesh::triangle_ids getOwnedEntities(const model::patch_id& patch);

    //    /**
    //     * \param patch patch label
    //     * \return boundaries owned by this process that belong to the given
    //     patch
    //     */
    //    mesh::triangle_ids getOwnedEntities(const mesh::patch_id &patch);
    //
    /**
     * \param compartment compartment label
     * \return elements owned by this process that belong to the given
     * compartment
     */
    mesh::tetrahedron_ids getOwnedEntities(const model::compartment_id& compartment);

    /**
     * \brief Get all elements in this partition that belong to a compartment.
     *
     * Return a vector in which all tetrahedral point indices in the mesh
     * partition is stored, if the tetrahedrons belong to the compartment. Note
     * that this includes ghost points.
     */
    mesh::tetrahedron_ids getEntities(const model::compartment_id& compartment);

    /**
     * \param compartmentId compartment label
     * \return elements in this process (owned or not) that belong to the given
     * compartment
     */
    mesh::triangle_ids getEntities(const model::patch_id& compartmentId);

    /**
     * \brief Get tetrahedrons tagged with the given tag.
     *
     * \attention Parallelism: Collective
     *
     * \return vector of the global tetrahedron indices.
     */
    std::vector<mesh::tetrahedron_global_id_t> getTaggedTetrahedrons(
        const model::compartment_id& comp);

    /**
     * \brief Get tetrahedrons tagged with the given tag in this partition.
     *
     * \attention Parallelism: Local
     *
     * \return vector of the local tetrahedron indices.
     */
    std::vector<mesh::tetrahedron_local_id_t> getTaggedLocalTetrahedrons(
        const model::compartment_id& comp,
        bool owned);

    /**
     * \brief Get triangles tagged with the given tag.
     *
     * \attention Parallelism: Collective
     *
     * \return vector of the global triangle indices.
     */
    std::vector<mesh::triangle_global_id_t> getTaggedTriangles(const model::patch_id& patch);

    /**
     * \brief Get triangles tagged with the given tag in this partition.
     *
     * \attention Parallelism: Local
     *
     * \return vector of the local triangle indices.
     */
    std::vector<mesh::triangle_local_id_t> getTaggedLocalTriangles(const model::patch_id& patch,
                                                                   bool owned);

    /**
     * \brief Get vertices tagged with the given tag.
     *
     * \attention Parallelism: Collective
     *
     * \return vector of the global vertex indices.
     */
    std::vector<mesh::vertex_global_id_t> getTaggedVertices(const model::vertgroup_id& verts);

    /**
     * \brief Get vertices tagged with the given tag in this partition.
     *
     * \attention Parallelism: Local
     *
     * \return vector of the local vertex indices.
     */
    std::vector<mesh::vertex_local_id_t> getTaggedLocalVertices(const model::vertgroup_id& verts,
                                                                bool owned);

    /**
     * \brief Get the names of tagged regions of dimension dim.
     *
     * \attention Parallelism: Local
     *
     * \return vector of names of tagged regions
     */
    std::vector<std::string> getTags(osh::Int dim) const;

    /// provide the number of neighbors of a given element
    inline const osh::Read<osh::LO>& neighbors_per_element() const noexcept {
        return neighbors_per_element_;
    }

    const std::unordered_map<model::compartment_id, model::compartment_label>& compartment_labels()
        const noexcept {
        return compIdtoLabel;
    }
    /**
     * \brief Add to the mesh geometrical information about a diffusion
     * boundary. This is done by providing ids of neighbouring compartment and
     * potentially ids of a custom subset of triangles on the boundary of such
     * compartments.
     */
    void addDiffusionBoundary(
        const mesh::diffusion_boundary_name& name,
        const model::compartment_id& comp1,
        const model::compartment_id& comp2,
        std::optional<std::set<mesh::triangle_global_id_t>> triangles = std::nullopt);

    struct DiffusionBoundary {
        /// Test whether a species with container 1 index is diffusing
        std::vector<bool> comp1_diffusing_species;
        /// Test whether a species with container 2 index is diffusing
        std::vector<bool> comp2_diffusing_species;
        /// Convert a container 1 index into a container 2 index
        std::vector<container::species_id> conv_12;
        /// Convert a container 2 index into a container 1 index
        std::vector<container::species_id> conv_21;
        /// triangles of a diffusing boundary
        std::vector<mesh::triangle_id_t> triangles;
        /// Adjacent compartments id
        mesh::compartment_id msh_comp1, msh_comp2;
        model::compartment_id mdl_comp1, mdl_comp2;
    };

    std::vector<DiffusionBoundary>& diffusionBoundaries() noexcept {
        return diffusion_boundaries_;
    }

    void addMembrane(const model::membrane_id name, DistMemb* memb);

    const std::map<model::membrane_id, DistMemb*>& membranes() const noexcept {
        return membranes_;
    }

    /**
     * \brief Convert a species id in the from_comp_id compartment into a
     * species id in the other compartment
     */
    inline container::species_id convertSpeciesID(mesh::triangle_id_t triangle,
                                                  mesh::compartment_id from_comp_id,
                                                  container::species_id spec_id) const {
        const auto& db = diffusion_boundaries_[*diffusion_boundary_ids_[triangle.get()]];
        return (db.msh_comp1 == from_comp_id)
                   ? diffusion_boundaries_[*diffusion_boundary_ids_[triangle.get()]]
                         .conv_12[spec_id.get()]
                   : diffusion_boundaries_[*diffusion_boundary_ids_[triangle.get()]]
                         .conv_21[spec_id.get()];
    }

    /**
     * \brief  Determine whether a triangle allows diffusion of a species id
     * spec_id defined in the compartment comp_id
     */
    inline bool isActiveDiffusionBoundary(mesh::triangle_id_t triangle,
                                          mesh::compartment_id comp_id,
                                          container::species_id spec_id) const {
        if (diffusion_boundary_ids_[triangle.get()]) {
            return (diffusion_boundaries_[*diffusion_boundary_ids_[triangle.get()]].msh_comp1 ==
                    comp_id)
                       ? diffusion_boundaries_[*diffusion_boundary_ids_[triangle.get()]]
                             .comp1_diffusing_species[spec_id.get()]
                       : diffusion_boundaries_[*diffusion_boundary_ids_[triangle.get()]]
                             .comp2_diffusing_species[spec_id.get()];
        } else {
            return false;
        }
    }

    /**
     * \brief Extract the boundary index from the nickname
     */
    size_t getDiffusionBoundaryIndex(const mesh::diffusion_boundary_name& name) const {
        auto it = diff_bound_name_2_index_.find(name);
        if (it != diff_bound_name_2_index_.end()) {
            return it->second;
        } else {
            throw std::invalid_argument(std::string("Unknown diffusion boundary ") + name);
        }
    }

    // public alias type for segment intersections
    using intersection_list_t = std::vector<std::pair<mesh::tetrahedron_local_id_t, double>>;

    /**
     * \brief Computes the percentage of intersection of a segment with the mesh tets
     *
     * \param init_seg_length Length of segment as called by intersect() to avoid endless recursion
     * \return A vector where each position contains pairs <tet, intersection ratio>
     */
    intersection_list_t intersectDeterministic(const point3d& p_start,
                                               const point3d& p_end,
                                               const double init_seg_length);

    /**
     * \brief Computes the percentage of intersection of a line of segments with the mesh tets
     *
     * \return A vector of vectors (for each segment) containing pairs <tet, intersection ratio>
     */
    std::vector<intersection_list_t> intersect(const double* points,
                                               int n_points,
                                               int sampling = -1);

    /**
     * \brief Similar to the intersect method but here we deal with independent segments, i.e.
     *        every two points we have a segment not related to previous or following ones.
     *        E.g. seg0 = (points[0], points[1]), seg1 = (points[2], points[3]), etc.
     *
     * \return  A vector of vectors (for each segment) containing pairs <tet, intersection ratio>
     */
    std::vector<intersection_list_t> intersectIndependentSegments(const double* points,
                                                                  int n_points,
                                                                  int sampling = -1);

    /**
     * \brief Gather geometrical entities indices across MPI processes
     *
     * \attention Parallelism: Collective
     *
     * \tparam Entity geometrical entity index strong id type
     * \param entities vector of entities in current process
     * \param datatype MPI data type
     * \return The assembled vector of indices
     */
    template <typename Entity>
    std::vector<Entity> allGatherEntities(const std::vector<Entity>& entities,
                                          MPI_Datatype datatype);

  private:
    std::unordered_map<mesh::diffusion_boundary_name, size_t> diff_bound_name_2_index_;

    void fill_triInfo(const Omega_h::Reals& coords, const Omega_h::Reals& areas);
    void fill_tetInfo(const Omega_h::Reals& coords, const Omega_h::Reals& areas);

    /**
     * \brief get entities owned by this process of the specified class label
     * \tparam Tag strong type tag
     * \param region label of the class in the mesh
     * \param owned owned elements
     * \return entity identifiers
     */
    template <typename Tag>
    osh::LOs getEntitiesImpl(const util::strong_string<Tag>& region, bool owned);

    template <typename Tag, typename Global, typename Local>
    std::vector<Global> getAllEntities(const Tag& tag);

    /**
     *
     * \param verts identifier
     * \return the omega_h class pairs linked to the region.
     */
    std::vector<osh::ClassPair> getClassPairs(const model::vertgroup_id& verts) const;

    /**
     *
     * \param compartment identifier
     * \return the omega_h class pairs linked to the region.
     */
    std::vector<osh::ClassPair> getClassPairs(const model::compartment_id& compartment) const;

    /**
     *
     * \param patch patch identifier
     * \return the omega_h class pairs linked to the region.
     */
    std::vector<osh::ClassPair> getClassPairs(const model::patch_id& patch) const;

    /**
     * \brief Synchronize data across MPI processes
     *
     * \attention Parallelism: Collective
     *
     * \param buff pointer to data
     * \param count number of elements to be synchronized
     * \param datatype MPI data type
     * \param isRoot true if data from this process should be synchronized
     * \param unknownCount true if only the root processes knows the count
     */
    void syncData(void* buff,
                  int count,
                  MPI_Datatype datatype,
                  bool isRoot,
                  bool unknownCount = false) const;

    osh::Mesh mesh_;
    const std::string path_;
    const osh::Real scale_;

    std::unique_ptr<Measure> measure_;
    Measure::element_measure_func measureFunc_;

    /// \brief mask to distinguish owned element from ghost element
    /// owned_elems_mask_[i] != 0 if i is owned
    osh::Bytes owned_elems_mask_;

    /// \brief mask to distinguish owned boundaries from ghost element
    /// owned_bounds_mask_[i] != 0 if i is owned
    osh::Bytes owned_bounds_mask_;

    /// \brief mask to distinguish owned boundaries from ghost element
    /// owned_verts_mask_[i] != 0 if i is owned
    osh::Bytes owned_verts_mask_;

    /** \name Elements
     * \{
     */
    /// Identifiers of the elements owned by this process
    mesh::tetrahedron_ids owned_elems_;
    std::unordered_map<mesh::tetrahedron_global_id_t, mesh::tetrahedron_local_id_t>
        elemGlobal2Local;
    osh::GOs elemLocal2Global;
    /** \} */

    /** \name Boundaries
     * \{
     */
    /// Identifiers of the boundaries owned by this process
    osh::LOs owned_bounds_;
    std::unordered_map<mesh::triangle_global_id_t, mesh::triangle_local_id_t> boundGlobal2Local;
    osh::GOs boundLocal2Global;
    /** \} */

    /** \name Vertices
     * \{
     */
    /// Identifiers of the vertices owned by this process
    osh::LOs owned_verts_;
    std::unordered_map<mesh::vertex_global_id_t, mesh::vertex_local_id_t> vertGlobal2Local;
    osh::GOs vertLocal2Global;
    /** \} */

    // Bounding box of the owned elements
    point3d ownedBBoxMin{};
    point3d ownedBBoxMax{};

    std::vector<TetStruct> tetInfo_;
    std::vector<TriStruct> triInfo_;
    /// provide the number of neighbors of a given element
    osh::Read<osh::LO> neighbors_per_element_;
    /// provide the number of neighbors of an owned element index
    osh::Read<osh::LO> neighbors_per_owned_element_idx_;

    /// get compartment identifier of a given element
    osh::Write<osh::LO> elem2compid_;
    std::unordered_map<mesh::compartment_id, mesh::tetrahedron_ids> compid2elems_;
    /// FIXME TCL: compid2ownedvol should be moved in the `measure` class
    std::vector<osh::Real> compid2ownedvol;
    std::vector<mesh::tetrahedron_ids> comp2owned_elems_;
    std::unordered_map<model::compartment_label, model::compartment_id> compLabelToId;
    std::unordered_map<model::compartment_id, model::compartment_label> compIdtoLabel;
    std::map<model::compartment_id, mesh::compartment_id> apicompid2meshcompid;
    std::vector<model::compartment_id> meshcompid2apicompid;
    std::vector<DistComp*> distcomps;

    std::map<model::patch_id, mesh::patch_id> apipatchid2meshpatchid;
    std::vector<model::patch_id> meshpatchid2apipatchid;
    std::unordered_map<mesh::patch_id, mesh::triangle_ids> patchid2bounds_;
    std::vector<mesh::triangle_ids> patch2owned_bounds_;
    std::vector<DistPatch*> distpatches;

    /// store neighbors distance and area
    util::flat_multimap<osh::Real, 2> tet_neighbors_real_data_;
    /// store neighbors identifier, face identifier, triangle id
    util::flat_multimap<osh::LO, 3> tet_neighbors_int_data_;
    /// neighbor index to be used with \a tet_neighbors_real_data_ and \a
    /// tet_neighbors_int_data_ to only get the neighbors that belong to the
    /// same compartment
    util::flat_multimap<osh::LO, 1> tet_neighbors_in_comp_index_;
    /// all diffusion boundaries
    std::vector<DiffusionBoundary> diffusion_boundaries_;
    /// store the diffusion boundary id
    std::vector<optional_id_t> diffusion_boundary_ids_;
    /// all membranes
    std::map<model::membrane_id, DistMemb*> membranes_;

    /// \brief number of elements in the entire mesh
    const osh::GO total_num_elems_;

    /// \brief number of boundaries in the entire mesh
    const osh::GO total_num_bounds_;

    /// \brief number of vertices in the entire mesh
    const osh::GO total_num_verts_;
};

/*****************************************************************************/
// Implementation:
// "An efficient and robust ray-box intersection algorithm", by
// Amy Williams and Steve Barrus and R. Keith Morley and Peter Shirley
// doi: 10.1145/1198555.1198748
class Ray {
  public:
    using point3d = math::point3d;

    Ray(const point3d& o, const point3d& d)
        : origin(o)
        , direction(d) {
        inv_direction = point3d(1. / d[0], 1. / d[1], 1. / d[2]);
        sign[0] = (inv_direction[0] < 0);
        sign[1] = (inv_direction[1] < 0);
        sign[2] = (inv_direction[2] < 0);
    }

    point3d origin, direction;  // ray origin and direction
    point3d inv_direction;
    std::array<int, 3> sign;
};
// Axis-Aligned Bounding Box in 3D
class AABB3 {
  public:
    using point3d = math::point3d;

    // define a box with ordered corners min and max
    AABB3(const point3d& vmin, const point3d& vmax) {
        bounds[0] = vmin;
        bounds[1] = vmax;
    }

    // Determine whether a ray origin+t*dir and box intersect within the ray's parameterized
    // range (t0,t1)
    bool intersect(const steps::dist::Ray& r, float t0, float t1) const;

    std::array<point3d, 2> bounds;
};
/*****************************************************************************/

}  // namespace steps::dist
