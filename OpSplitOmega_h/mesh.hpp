#pragma once

#include <map>
#include <string>
#include <vector>

#include <Omega_h_array_ops.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>
#include <hadoken/format/format.hpp>
#include <petscdmlabel.h>
#include <petscdmplex.h>


#include "common.hpp"
#include "flat_multimap.hpp"
#include "opsplit/mesh.hpp"
#include "opsplit/vocabulary.hpp"

namespace zee {

struct TetStruct {
    PetscScalar vol;
    osh::Vector<3> centroid;
};

struct TriStruct {
    osh::LO num_neighbors;
    PetscScalar area;
    osh::Vector<2> centroid;
};

/**
 * Provide data and operations dealing with measure of mesh elements i.e
 * surface in 2D and volume in 3D.
 */
class MeasureInfo {
  public:
    /// type of a function taking an element identifier in parameter
    /// and return its measure.
    using element_measure_func = std::function<osh::Real(mesh::element_id)>;

    /**
     *
     * \param num_compartments Number of compartments in the mesh
     * \param t_element_measure_func Functor to get measure of an element
     */
    MeasureInfo(osh::Int num_compartments, const element_measure_func& t_element_measure_func);

    void init(const mesh::element_ids& t_owned_elements, const osh::LOs& t_elem2compid);

    /**
     * \return the sum of all element measures assigned to this MPI rank
     */
    inline osh::Real rank_measure() const noexcept {
        return rank_measure(rank_);
    }

    /**
     * \return the sum of all element measures in the given MPI rank
     */
    inline osh::Real rank_measure(osh::LO rank) const noexcept {
        return rank2measures_[rank * num_measures_per_rank_];
    }

    /**
     * \return the sum of all element measures that belong to the given compartment
     * in the current rank
     */
    inline osh::Real rank_measure(mesh::compartment_id compartment) const noexcept {
        return rank_measure(rank_, compartment);
    }

    /**
     * \return the sum of all element measures part of the given rank and compartment
     */
    inline osh::Real rank_measure(osh::Int rank, mesh::compartment_id compartment) const noexcept {
        return rank2measures_[rank * num_measures_per_rank_ + 1 + compartment.get()];
    }

    /**
     * \return measure of the entire mesh
     */
    inline osh::Real mesh_measure() const noexcept {
        return mesh_measures_[0];
    }

    /**
     * \return measure of the given compartment in the entire mesh
     */
    inline osh::Real mesh_measure(mesh::compartment_id compartment) const noexcept {
        return mesh_measures_[1 + compartment.get()];
    }

    /**
     * \return measure of the given element
     */
    inline osh::Real element_measure(mesh::element_id element) const noexcept {
        return element_measure_func_(element);
    }

    /**
     * Given a total number of molecules to spread within the mesh,
     * this function returns the subset to insert in the given compartment
     * and MPI rank.
     * \param rank MPI rank
     * \param compartment target compartment identifier
     * \param num_molecules number of molecules to spread in the entire mesh
     */
    inline osh::Real molecules_in_rank(osh::Int rank,
                                       mesh::compartment_id compartment,
                                       osh::Real num_molecules) const noexcept {
        return num_molecules * rank_measure(rank, compartment) / mesh_measure(compartment);
    }


    /**
     * Given a total number of molecules to spread in the current rank, this function
     * returns, for a given element and compartment, the number of molecules in regards of its
     * measure compared to the other elements. \param compartment target compartment identifier
     * \param element a mesh element
     * \param num_molecules number of molecules to spread in the current rank
     */
    inline osh::Real molecules_in_element(mesh::compartment_id compartment,
                                          mesh::element_id element,
                                          osh::Real num_molecules) const {
        return num_molecules * element_measure_func_(element) / rank_measure(rank_, compartment);
    }

    /**
     * \return function to get measure of a given element
     */
    inline const element_measure_func& measure_func() const noexcept {
        return element_measure_func_;
    }

  private:
    /// current rank
    const osh::Int rank_;
    /// 1 + number of compartments
    const osh::Int num_measures_per_rank_;
    /// functor to get measure of an element
    const element_measure_func& element_measure_func_;
    /// sum of element measure in the entire mesh
    osh::Reals mesh_measures_;
    /// measures of every compartment for every rank
    osh::Reals rank2measures_;
};

namespace {
/// Functor to compute
template <class Mesh, osh::Int Dim>
struct ElementMeasure {};

template <class Mesh>
struct ElementMeasure<Mesh, 2> {
    PetscScalar operator()(const Mesh& mesh, PetscInt index) {
        return mesh.getTri(index).area;
    }
};

template <class Mesh>
struct ElementMeasure<Mesh, 3> {
    PetscScalar operator()(const Mesh& mesh, PetscInt index) {
        return mesh.getTet(index).vol;
    }
};
}  // namespace

/**
 * \brief "Wrapper class" encapsulates the distribute mesh library.
 *
 * This class encapsulates the distribute mesh library, providing
 * geometry data access of each element in a distributed mesh via
 * the element's index, as well as grouping and ownership data access
 * for compartment, patch and ROI (in the future).
 *
 * Replacement for the Tetmesh class in STEPS.
 */
template <osh::Int Dim>
class OmegaHMesh: public DistMesh {
  public:
    /**
     * \brief Constructor.
     */
    OmegaHMesh(osh::Mesh& mesh,
               const std::string& t_filename,
               PetscScalar t_scale,
               bool t_label_elems);

    void init();

    std::string getBackend() const override;

    /**
     * \brief Get the point index of a tetrahedron from its label index.
     *
     * Point index of an element is the DAG index DMPlex assigns to this element.
     * Label index is the index provided by the Gmsh file (in the case of tetrahedron and vertex),
     * or generated internally during the importing (interpolate layers such as triangle and bar).
     * \param label_idx Label index of the tetrahedron.
     * \return Point index of the tetrahedron.
     */
    // TODO TCL: implement this at some point even if useless in OpSplitOmega_h
    //  PetscInt getTetPointIdx(PetscInt label_idx) const;

    /**
     * \brief Get the point index of a triangle from its label index (Not yet available).
     *
     * Point index of an element is the DAG index DMPlex assigns to this element.
     * Label index is the index provided by the Gmsh file (in the case of tetrahedron and vertex),
     * or generated internally during the importing (interpolate layers such as triangle and bar).
     * \param label_idx Label index of the triangle.
     * \return Point index of the triangle.
     */
    // TODO TCL: implement this at some point even if useless in OpSplitOmega_h
    // PetscInt getTriPointIdx(PetscInt label_idx) const;

    /**
     * \brief Get the point index of a bar from its label index (Not yet available).
     *
     * Point index of an element is the DAG index DMPlex assigns to this element.
     * Label index is the index provided by the Gmsh file (in the case of tetrahedron and vertex),
     * or generated internally during the importing (interpolate layers such as triangle and bar).
     * \param label_idx Label index of the bar.
     * \return Point index of the bar.
     */
    // TODO TCL: implement this at some point even if useless in OpSplitOmega_h
    // PetscInt getBarPointIdx(PetscInt label_idx) const;

    /**
     * \brief Get the point index of a vertex from its label index (Not yet available).
     *
     * Point index of an element is the DAG index DMPlex assigns to this element.
     * Label index is the index provided by the Gmsh file (in the case of tetrahedron and vertex),
     * or generated internally during the importing (interpolate layers such as triangle and bar).
     * \param label_idx Label index of the vertex.
     * \return Point index of the vertex.
     */
    // TODO TCL: implement this at some point even if useless in OpSplitOmega_h
    // PetscInt getVertPointIdx(PetscInt label_idx) const;

    /**
     * \brief Get the label index of a tetrahedron from its point index.
     *
     * Label index is the index provided by the Gmsh file (in the case of tetrahedron and vertex),
     * or generated internally during the importing (interpolate layers such as triangle and bar).
     * Point index of an element is the DAG index DMPlex assigns to this element.
     * \param point_idx Point index of the tetrahedron.
     * \return Label index of the tetrahedron.
     */
    // TODO TCL: implement this at some point  even if useless in OpSplitOmega_h
    // PetscInt getTetLabelIdx(PetscInt point_idx);

    /**
     * \brief Get the label index of a triangle from its point index (Not yet available).
     *
     * Label index is the index provided by the Gmsh file (in the case of tetrahedron and vertex),
     * or generated internally during the importing (interpolate layers such as triangle and bar).
     * Point index of an element is the DAG index DMPlex assigns to this element.
     * \param point_idx Point index of the triangle.
     * \return Label index of the triangle.
     */
    // TODO TCL: implement this at some point  even if useless in OpSplitOmega_h
    //  PetscInt getTriLabelIdx(PetscInt point_idx);

    /**
     * \brief Get the label index of a bar from its point index (Not yet available).
     *
     * Label index is the index provided by the Gmsh file (in the case of tetrahedron and vertex),
     * or generated internally during the importing (interpolate layers such as triangle and bar).
     * Point index of an element is the DAG index DMPlex assigns to this element.
     * \param point_idx Point index of the bar.
     * \return Label index of the bar.
     */
    // TODO TCL: implement this at some point  even if useless in OpSplitOmega_h
    //  PetscInt getBarLabelIdx(PetscInt point_idx);

    /**
     * \brief Get the label index of a vertex from its point index (Not yet available).
     *
     * Label index is the index provided by the Gmsh file (in the case of tetrahedron and vertex),
     * or generated internally during the importing (interpolate layers such as triangle and bar).
     * Point index of an element is the DAG index DMPlex assigns to this element.
     * \param point_idx Point index of the vertex.
     * \return Label index of the vertex.
     */
    // TODO TCL: implement this at some point  even if useless in OpSplitOmega_h
    //  PetscInt getVertLabelIdx(PetscInt point_idx);

    /**
     * \brief Get the geometry data struct of a tetrahedron from its point index.
     *
     * This function returns a TetStruct struct which encapsulates the geometry data
     * of the tetrahedron required in the simulation.
     * \param index Local index of the tetrahedron.
     * \return TetStruct struct of the tetrahedron.
     */
    const TetStruct& getTet(PetscInt index) const;

    /**
     * \brief Get the geometry data struct of a triangle from its point index.
     *
     * This function returns a TriStruct struct which encapsulates the geometry data
     * of the triangle required in the simulation.
     * \param index Local index of the triangle.
     * \return TriStruct struct of the index.
     */
    const TriStruct& getTri(PetscInt index) const;

    /**
     * \brief Check if a point is owned by this process.
     *
     * Giving the point index of a mesh element, this function returns
     * if the element is owned by the process, i.e. not ghost point.
     * \param point_idx Point index of the element.
     * \return True or 't' if the element is owned by the process. False or 'f' if it is not.
     */
    inline bool isOwned(mesh::element_id point_idx) const noexcept {
        assert(point_idx.get() < ownedElemsMask.size());
        return ownedElemsMask[point_idx.get()] != 0;
    }

    /**
     * \brief Return the label of the compartment to which the given tetrahedron belong.
     *
     * \param point_idx Point index of the tetrahedron.
     * \return "Cell Set" label of the compartment.
     */
    // PetscInt getTetCompLabel(PetscInt point_idx) const;

    /**
     * \brief Check if two tetrahedrons are in the same compartment.
     *
     * \param tet0 Point index of the first tetrahedron.
     * \param tet1 Point index of the second tetrahedron.
     * \return True if both tetrahedrons are in the same compartment, false if otherwise.
     */
    bool inSameComp(PetscInt tet0, PetscInt tet1) const;

    /**
     * \brief Return the type of a mesh element.
     *
     * \param point_idx Index of the element point.
     * \return The geometry type of the element.
     */
    ElementType getElementType(PetscInt point_idx) const;

    /**
     * Get measure of an element (area in 2D, volume in 3D)
     * \param index Index of the element point
     * \return measure of the element
     */
    PetscScalar getElementMeasure(PetscInt index) const {
        return ElementMeasure<OmegaHMesh, Dim>()(*this, index);
    }

    PetscScalar getOwnedCompVol(model::compartment_label comp_label) const override;

    mesh::compartment_id getCompID(const model::compartment_id& compartment) noexcept {
        const auto nextid = mesh::compartment_id(
            static_cast<mesh::compartment_id::value_type>(apicompid2meshcompid.size()));
        const auto& status = apicompid2meshcompid.insert({compartment, nextid});
        if (status.second) {
            meshcompid2apicompid.push_back(compartment);
        }
        return status.first->second;
    }

    // TODO TCL: remove this
    //    /**
    //     * \brief Return the DMPlex DM of the distributed mesh.
    //     *
    //     * \return DMPlex DM of the distributed mesh.
    //     */
    //    inline const DM& getDM() const {
    //        return dm;
    //    }

    inline const Omega_h::Mesh& getMesh() const noexcept {
        return mesh;
    }

    inline Omega_h::Mesh& getMesh() noexcept {
        return mesh;
    }

    /**
     * \return number of elements owned by the current rank
     */
    PetscInt getOwnedNumElements() const override;

    std::string createReport() const override;

    const MeasureInfo& getMeasureInfo() const noexcept {
        return *measureInfo;
    }

    inline const mesh::element_ids& getOwnedElems() const noexcept {
        return owned_elems_;
    }

    /**
     * \param patch patch label
     * \return boundaries owned by this process that belong to the given patch
     */
    mesh::boundary_ids getOwnedEntities(const model::patch_id& patch);

    /**
     * \param compartment compartment label
     * \return elements owned by this process that belong to the given compartment
     */
    mesh::element_ids getOwnedEntities(const model::compartment_id& compartment);

    /**
     * \param compartmentId compartment label
     * \return elements in this process (owned or not) that belong to the given compartment
     */
    mesh::element_ids getEntities(const model::compartment_id& compartmentId);

    /**
     * \param compartmentId compartment label
     * \return elements in this process (owned or not) that belong to the given compartment
     */
    mesh::boundary_ids getEntities(const model::patch_id& compartmentId);

    /**
     *
     * @param region patch or compartment identifier
     * @return tuple containing all elements the measure of the region
     */
    std::tuple<osh::LOs, osh::Read<osh::Real>, osh::Real> measure(const model::region_id& region);

    PetscScalar getOwnedPatchArea(model::patch_id patchId) override {
        return std::get<2>(measure(patchId));
    }


    /**
     * \return number of elements owned by this rank
     */
    inline osh::LO numOwnedElems() const noexcept {
        return getOwnedElems().size();
    }

    const mesh::element_ids& getOwnedElems(mesh::compartment_id compartment) const noexcept {
        return comp2owned_elems_[static_cast<size_t>(compartment.get())];
    }

    inline const osh::Read<osh::I8>& getOwnedElemsMask() const noexcept {
        return ownedElemsMask;
    }

    inline const TetStruct& getTetInfo(osh::LO element) const noexcept {
        return tetInfo[static_cast<size_t>(element)];
    }

    /**
     *
     * \return container holding information of triangles
     */
    inline const std::vector<TriStruct>& getTriInfo() const noexcept {
        return triInfo;
    }

    /**
     *
     * \return container holding information of tetrahedrons
     */
    inline const std::vector<TetStruct>& getTetInfo() const noexcept {
        return tetInfo;
    }

    /**
     * \brief Expose floating point data associated to tetrahedrons neighbors.
     * \FIXME(TCL) this should not be exposed directly
     */
    inline const zee::flat_multimap<osh::Real, 2, zee::OSH>& tet_neighbors_real_data() const
        noexcept {
        return tet_neighbors_real_data_;
    }

    /**
     * \brief Expose integral data associated to tetrahedrons neighbors.
     * \FIXME(TCL) this should not be exposed directly
     */
    inline const zee::flat_multimap<osh::LO, 2, zee::OSH>& tet_neighbors_int_data() const noexcept {
        return tet_neighbors_int_data_;
    }

    /**
     * \return Number of neighbors per element identifier
     */
    inline const osh::Read<osh::LO>& neighbors_per_element() const noexcept {
        return neighbors_per_element_;
    }

    /**
     * \return Number of neighbors for a given owned element index in \a owned_elems_ container
     */
    inline const osh::Read<osh::LO> neighbors_per_owned_element_idx() const noexcept {
        return neighbors_per_owned_element_idx_;
    }

    model::compartment_id getCompartment(mesh::element_id element) const noexcept {
        const auto mesh_comp_id(elem2compid[element.get()]);
        return meshcompid2apicompid[static_cast<size_t>(mesh_comp_id)];
    }

  protected:
    void addCompImpl(const model::compartment_id& comp_id,
                     model::compartment_label cell_set_label) override;

  private:
    void importFromFile();
    void fill_triInfo(const Omega_h::Reals& coords, const Omega_h::Reals& areas);
    void fill_tetInfo(const Omega_h::Reals& coords, const Omega_h::Reals& areas);

    inline osh::Int num_compartments() const noexcept {
        return static_cast<osh::Int>(compLabelToId.size());
    }

    /**
     * \brief get entities owned by this process of the specified class label
     * \tparam Tag strong type tag
     * \param region label of the class in the mesh
     * \param owned owned elements
     * \return entity identifiers
     */
    template <typename Tag>
    osh::LOs getEntitiesImpl(const strong_string<Tag>& region, bool owned);

    /**
     *
     * @param label patch or compartment label
     * @return the omega_h class pairs linked to the region.
     */
    template <typename Tag>
    std::vector<osh::ClassPair> getClassPairs(const strong_string<Tag>& label) const;

    PetscErrorCode ierr{};

    Omega_h::Mesh& mesh;

    /// mask to distinguish owned element from ghost element
    /// ownedElemsMask[i] != 0 if i is owned
    osh::Read<osh::I8> ownedElemsMask;
    /// owned elements identifiers
    mesh::element_ids owned_elems_;

    std::vector<mesh::element_ids> comp2owned_elems_;

    std::vector<TetStruct> tetInfo;
    /// store neighbors distance and area
    zee::flat_multimap<osh::Real, 2, zee::OSH> tet_neighbors_real_data_;
    /// store neighbors identifier and face identifier
    zee::flat_multimap<osh::LO, 2, zee::OSH> tet_neighbors_int_data_;
    std::vector<TriStruct> triInfo;
    std::unique_ptr<MeasureInfo> measureInfo;
    MeasureInfo::element_measure_func measureFunc;
    osh::Read<osh::LO> neighbors_per_element_;
    osh::Read<osh::LO> neighbors_per_owned_element_idx_;
    osh::Write<osh::LO> elem2compid;
    std::map<model::compartment_id, mesh::compartment_id> apicompid2meshcompid;
    std::vector<model::compartment_id> meshcompid2apicompid;
    std::vector<PetscScalar> compid2ownedvol;
};

// explicit instantiation declarations
extern template class OmegaHMesh<2>;
extern template std::vector<osh::ClassPair> OmegaHMesh<2>::getClassPairs(
    const model::patch_id& label) const;
extern template std::vector<osh::ClassPair> OmegaHMesh<2>::getClassPairs(
    const model::compartment_id& label) const;

extern template class OmegaHMesh<3>;
extern template std::vector<osh::ClassPair> OmegaHMesh<3>::getClassPairs(
    const model::patch_id& label) const;
extern template std::vector<osh::ClassPair> OmegaHMesh<3>::getClassPairs(
    const model::compartment_id& label) const;


}  // namespace zee
