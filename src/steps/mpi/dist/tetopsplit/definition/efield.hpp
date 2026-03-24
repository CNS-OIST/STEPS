#pragma once

#include <Omega_h_defines.hpp>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#ifdef USE_PETSC
#include <Omega_h_few.hpp>
#include <petscksp.h>
#endif  // USE_PETSC

#include "complexeventsdef.hpp"
#include "fwd.hpp"
#include "geom/dist/fwd.hpp"
#include "model/fwd.hpp"
#include "model/ghkcurr.hpp"
#include "model/ohmiccurr.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"
#include "util/error.hpp"
#include "util/vocabulary.hpp"


namespace steps::dist {

//------------------------------------------------
/**
 * \brief Base class for the definition of an Ohmic current.
 */
class OhmicCurrdefBase {
  public:
    /**
     * \brief Ohmic current ctor.
     *
     * \param t_conductance conductance that generates the ohmic current, per unit
     * of the channel state, if channel state is specified and per unit surface
     * otherwise.
     * \param t_reversal_potential potential at which current reverses.
     * \param t_channel_state the channel state that enables the ohmic current. If
     * none is provided the ohmic current is enabled.
     */
    OhmicCurrdefBase(const steps::model::OhmicCurrBase& curr)
        : conductance(curr.getG())
        , reversal_potential(curr.getERev()) {}

    OhmicCurrdefBase(const OhmicCurrdefBase&) = delete;
    virtual ~OhmicCurrdefBase() = default;

    osh::Real getReversalPotential(mesh::triangle_id_t triangle) const;
    void setReversalPotential(mesh::triangle_id_t triangle, osh::Real value);
    void reset();

#ifdef USE_PETSC
    /** Boundary condition to get the ohmic current flowing through a triangle an split among the
     * vertexes
     *
     * In case channel_state is not active Avert is used as a quantity of how many "channels" are
     * open
     *
     * check getTriCurrentOnVertex for additional info
     *
     * \param b_id
     * \param mol_state
     * \param Avert
     * \param sim_time
     * \return
     */
    virtual PetscReal getTriBConVertex(const mesh::triangle_id_t& b_id,
                                       const MolState& mol_state,
                                       double Avert,
                                       osh::Real sim_time) const = 0;

    /** get tri current on vertex
     *
     * get the vertex portion of the ohmic current flowing through a triangle. It is 1/3 of the
     * current flowing through the triangle
     *
     * \param potential_on_vertex: voltage on the particular vertex
     * \param patch_id: Id of the patch that contains the triangle
     * \param b_id: triangle id
     * \param mol_state: channel counts
     * \param mesh: required if the channel_state is not specified
     * \param sim_time: required for occupancy
     * \return current
     */
    PetscReal getTriCurrentOnVertex(osh::Real potential_on_vertex,
                                    const mesh::triangle_id_t& b_id,
                                    const MolState& mol_state,
                                    const DistMesh& mesh,
                                    osh::Real sim_time) const;
#else
    static constexpr auto PETSC_ERROR_MSG =
        "STEPS was compiled without PETSc, methods related to membrane potential computations are "
        "not available";
    double getTriBConVertex(const mesh::triangle_id_t& /*b_id*/,
                            const MolState& /*mol_state*/,
                            double /*Avert*/,
                            osh::Real /*sim_time*/) const {
        ArgErrLog(PETSC_ERROR_MSG);
    }
    double getTriCurrentOnVertex(osh::Real /*potential_on_vertex*/,
                                 const mesh::triangle_id_t& /*b_id*/,
                                 const MolState& /*mol_state*/,
                                 const DistMesh& /*mesh*/,
                                 osh::Real /*sim_time*/) const {
        ArgErrLog(PETSC_ERROR_MSG);
    }
#endif  // USE_PETSC

    const osh::Real conductance;

  protected:
    const osh::Real reversal_potential;
    std::unordered_map<mesh::triangle_id_t, osh::Real> reversal_potentials;
};

/**
 * \brief Ohmic current with species channel state.
 */
class OhmicCurrdef: public OhmicCurrdefBase {
  public:
    OhmicCurrdef(const Patchdef& patchdef, const steps::model::OhmicCurr& curr);
    virtual ~OhmicCurrdef() = default;

#ifdef USE_PETSC
    PetscReal getTriBConVertex(const mesh::triangle_id_t& b_id,
                               const MolState& mol_state,
                               double Avert,
                               osh::Real sim_time) const override;
#endif  // USE_PETSC
    friend std::ostream& operator<<(std::ostream& os, OhmicCurrdef const& m);

    container::species_id channel_state;
};

/**
 * \brief Ohmic current with complex channel state.
 */
class ComplexOhmicCurrdef: public OhmicCurrdefBase {
  public:
    ComplexOhmicCurrdef(const Patchdef& patchdef, const steps::model::ComplexOhmicCurr& curr);
    virtual ~ComplexOhmicCurrdef() = default;

#ifdef USE_PETSC
    PetscReal getTriBConVertex(const mesh::triangle_id_t& b_id,
                               const MolState& mol_state,
                               double Avert,
                               osh::Real sim_time) const override;
#endif  // USE_PETSC

    const ComplexFilterDescr channel_state;
    mutable model::complex_filter_occupancy_id occupancy_id;
};

//------------------------------------------------

/**
 * \brief Goldman-Hodgkin-Katz current.
 *
 */
class GHKCurrdefBase {
  public:
    GHKCurrdefBase(const Patchdef& patchdef, const steps::model::GHKcurrBase& curr);

  protected:
    const model::species_name ion_id;
    const osh::I64 valence;
};

class GHKCurrdef: public GHKCurrdefBase {
  public:
    /**
     * \brief GHKCurrent ctor.
     *
     * \param t_ion_channel_state channel state that enables the flow of ions
     * \param t_ion_id the species that flows through the membrane
     * \param t_valence the valence of the species that flows through the membrane
     */
    GHKCurrdef(const Patchdef& patchdef, const steps::model::GHKcurr& curr);

    /// Reset the def-object values to model defaults
    void reset() {}

    const container::species_id channel_state;
};

/**
 * \brief Goldman-Hodgkin-Katz current with a complex channel.
 *
 */
class ComplexGHKCurrdef: public GHKCurrdefBase {
  public:
    /**
     * \brief GHKCurrent ctor.
     *
     * \param t_ion_channel_state channel state that enables the flow of ions
     * \param t_ion_id the species that flows through the membrane
     * \param t_valence the valence of the species that flows through the membrane
     */
    ComplexGHKCurrdef(const Patchdef& patchdef, const steps::model::ComplexGHKcurr& curr);

    /// Reset the def-object values to model defaults
    void reset() {}

    const ComplexFilterDescr channel_state;
};

//------------------------------------------------

/**
 * \brief Membrane definition.
 */
struct Membrane {
    using Stimulus = std::function<osh::Real(osh::Real)>;

    /**
     * \brief Membrane ctor.
     *
     * \param patch patchid of the membrane
     * \param capacitance capacitance of the membrane
     */
    Membrane(Statedef& statedef_, const DistMemb& membrane);

    /// Reset the def-object values to model defaults
    void reset();

    void setStimulus(Stimulus stimulus) noexcept {
        current_ = std::move(stimulus);
    }

    void setConductivity(osh::Real conductivity) noexcept {
        conductivity_ = conductivity;
    }

    void setReversalPotential(osh::Real reversal_potential) noexcept {
        reversal_potential_ = reversal_potential;
    }

    void setCapacitance(osh::Real capacitance) noexcept {
        capacitance_ = capacitance;
    }

    const std::vector<model::patch_id>& getPatchesIds() const noexcept {
        return patches_ids;
    }

    const std::vector<std::reference_wrapper<const Patchdef>>& getPatchdefs() const noexcept {
        return patchdefs;
    }

    osh::Real capacitance() const noexcept {
        return capacitance_;
    }

    const Stimulus& stimulus() const noexcept {
        return current_;
    }

    osh::Real conductivity() const noexcept {
        return conductivity_;
    }

    osh::Real reversal_potential() const noexcept {
        return reversal_potential_;
    }

  private:
    const DistMemb& memb_;
    const Statedef& statedef;
    std::vector<std::reference_wrapper<const Patchdef>> patchdefs;

    std::vector<model::patch_id> patches_ids;
    osh::Real capacitance_;
    osh::Real conductivity_{};
    osh::Real reversal_potential_{};
    Stimulus current_;
};

#ifdef USE_PETSC
/** Small object for preparing values to be inserted in i(), bc(), A0, and diag() in efield
 *
 * It is created with only the capacitance as contribution. Ohmic currents are optional and
 * added after construction
 *
 * Since this struct is related to triangles the numbers 3 and 9 are related to the number of
 * vertexes in a triangle
 */
struct TriMatAndVecs {
    /// local ids of the vertexes of the triangle in PETSc format
    std::array<PetscInt, 3> face_bf2vertsPETSc;
    /// triangle contribution to the stiffness matrix. Size: Nvertexes^2
    std::array<PetscReal, 3 * 3> triStiffnessPETSc;
    /// triangle contribution to the rhs relative to the ohmic currents
    std::array<PetscReal, 3> triBC;
    /// triangle contribution to the rhs relative to the current injections
    std::array<PetscReal, 3> triI;

    /// Ctor
    TriMatAndVecs(const osh::Few<osh::LO, 3>& face_bf2verts,
                  PetscReal tri_capacitance,
                  PetscReal tri_i)
        : face_bf2vertsPETSc{face_bf2verts[0], face_bf2verts[1], face_bf2verts[2]}
        , triStiffnessPETSc{tri_capacitance,
                            0.0,
                            0.0,
                            0.0,
                            tri_capacitance,
                            0.0,
                            0.0,
                            0.0,
                            tri_capacitance}
        , triBC{0.0, 0.0, 0.0}
        , triI{tri_i, tri_i, tri_i} {}


    /// Print: for debugging
    friend std::ostream& operator<<(std::ostream& os, const TriMatAndVecs& obj);
};

#endif  // USE_PETSC

//------------------------------------------------

}  // namespace steps::dist
