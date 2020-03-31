#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include <petscdmlabel.h>
#include <petscdmplex.h>

#include "opsplit/vocabulary.hpp"

namespace zee {

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
class DistMesh {
  public:
    /**
     * \brief element types.
     */
    enum class ElementType { Tet, Tri, Bar, Vert, Unknown };

    /**
     * \brief Destructor.
     */
    virtual ~DistMesh();

    /**
     *
     * \return the distributed library name
     */
    virtual std::string getBackend() const = 0;

    /**
     * \brief Import a mesh from a Gmsh file, and distribute its elements.
     *
     * Given the Gmsh filename and the importing scale, this function imports the
     * mesh from the file and distribute its elements across all ranks. Auxiliary
     * data is also generated and stored.
     * \param t_filename The Gmsh file to be imported.
     * \param t_scale The import scale from the Gmsh unit to STEPS unit.
     * e.g. if 1 unit in Gmsh represents 1 um, then scale=1e-6
     * \param t_label_elems Label the element with its rank 0 id,
     * this is very time consuming with large mesh so should be done with caution
     */
    DistMesh(std::string t_filename, PetscScalar t_scale, bool t_label_elems);

    /**
     * \brief Add a compartment to a distributed mesh using label defined in the mesh file.
     *
     * Add a compartment which contains all tetrahedrons with the same "CellSet" label.
     * The value of "CellSet" label of a tetrahedron is the same as its component tag
     * in the Gmsh file.
     * \param compartment Name id of the compartment.
     * \param label Cell Set label of the compartment's elements.
     */
    virtual void addComp(const model::compartment_id& compartment, model::compartment_label label);

    /**
     * \brief Get the "Cell Set" label of a compartment.
     *
     * \param compartment Name id of the compartment.
     * \return "Cell Set" label of the compartment.
     */
    model::compartment_label getCompLabel(const model::compartment_id& compartment) const;

    /**
     * \brief Get the name id of a compartment.
     *
     * \param comp_label "Cell Set" label of the compartment.
     * \return Name id of the compartment.
     */
    const model::compartment_id& getCompID(model::compartment_label comp_label) const;

    /**
     * \brief Get the volume of compartment segment owned by the process.
     *
     * Return the sum of volumes of tetrahedrons that belong to the compartment
     * and owned by the process.
     * \param comp_label "Cell Set" label of the compartment.
     * \return Volume of the compartment segment owned by the process.
     */
    virtual PetscScalar getOwnedCompVol(model::compartment_label comp_label) const = 0;

    /**
     * \brief Get the volume of compartment segment across all processes.
     *
     * Return the sum of volumes of tetrahedrons that belong to the compartment
     * across all processes.
     * \param comp_label "Cell Set" label of the compartment.
     * \return Volume of the compartment.
     */
    PetscScalar getTotalCompVol(model::compartment_label comp_label) const;

    /**
     *
     * \param patchId patch identifier
     * \return Total area of the patch
     */
    PetscScalar getTotalPatchArea(model::patch_id patchId);

    /**
     *
     * \return Area of the owned part of the patch.
     */
    virtual PetscScalar getOwnedPatchArea(model::patch_id /*patchId*/) {
        throw std::logic_error("NOT IMPLEMENTED");
    }

    /**
     * \return number of owned elements in the current rank
     */
    virtual PetscInt getOwnedNumElements() const = 0;

    /**
     * \return number of elements in the entire mesh
     */
    virtual PetscInt getTotalNumElements() const;

    /**
     * \brief Create a report of the local mesh partition.
     *
     * \return A string reporting the content of the local mesh partition.
     */
    virtual std::string createReport() const = 0;

    const std::string filename;
    const PetscScalar scale;
    const bool label_elems;

    /**
     * \brief Get all elements in this partition that belong to a compartment.
     *
     * Return a vector in which all tetrahedral point indices in the mesh partition
     * is stored, if the tetrahedrons belong to the compartment.
     * Note that this includes ghost points.
     * \param comp_label "Cell Set" label of the compartment.
     * \return Vector of the tetrahedral point indices.
     */
    const mesh::element_ids& getCompElements(model::compartment_label comp_label) const;

    /**
     * \brief Get all elements in this partition that belong to a compartment.
     *
     * Return a vector in which all tetrahedral point indices in the mesh partition
     * is stored, if the tetrahedrons belong to the compartment.
     * Note that this includes ghost points.
     * \param comp_id user label of the compartment.
     * \return Vector of the tetrahedral point indices.
     */
    const mesh::element_ids& getCompElements(const model::compartment_id& comp_id) const;

  protected:
    virtual void addCompImpl(const model::compartment_id& comp_id,
                             model::compartment_label cell_set_label) = 0;

    std::unordered_map<model::compartment_label, mesh::element_ids> complabel2elems;
    std::unordered_map<model::compartment_id, model::compartment_label> compIdtoLabel;
    std::unordered_map<model::compartment_label, model::compartment_id> compLabelToId;
};

}  // namespace zee
