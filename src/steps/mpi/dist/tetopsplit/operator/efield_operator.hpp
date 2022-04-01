#pragma once

#include <petscksp.h>

#include "geom/dist/fwd.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"

namespace steps {
namespace dist {

/**
 * Model defined in
 * Iain Hepburn, Robert C. Cannon, Erik De Schutter Efficient calculation of the
 * quasi-static electrical potential on a tetrahedral mesh and its
 * implementation in STEPS Article (PDF Available) in Frontiers in Computational
 * Neuroscience 7:129 · October 2013
 *
 */
class EFieldOperator {
private:
  /// Placeholder for the solution of the linear system [V] (in steps 3: [mV])
  Vec sol_;
  /// Placeholder RHS, boundary conditions, and other currents of the linear system [A] (in steps 3:
  /// [pA])
  Vec rhs_, bc_, i_;
  /// Placeholder for the matrix of the linear system [S] (in steps 3: [nS])
  Mat A_;
  /// Time increment of the solver [s] (in steps 3: [ms])
  osh::Real dt_;
  /// Krylov solver
  KSP ksp_solver_;
  /// Krylov solver preconditioner
  PC ksp_solver_preconditioner_;
  /// Statedef
  const Statedef &state_def;
  /// Mesh
  DistMesh &mesh;
  /// Global index of vertices
  const osh::GOs global_indices_;

  inline Mat &A() noexcept { return A_; }
  inline Vec &rhs() noexcept { return rhs_; }
  inline Vec &bc() noexcept { return bc_; }
  inline Vec& i() noexcept {
      return i_;
  }
  inline Vec &sol() noexcept { return sol_; }

public:
  /**
   * \brief Ctor.
   *
   * \param o_mesh Omega_h mesh
   * \param statedef model definition
   * \param ghk_current_boundaries a mapping from a ghk current reaction index
   * to the mesh triangle id this reaction relates to \param dt maximum
   * e-field step size. Actual step size will depend on other operators as
   * well
   */
  template <typename NumMolecules>
  EFieldOperator(DistMesh& o_mesh,
                 const Statedef& statedef,
                 const std::vector<mesh::triangle_id_t>& ghk_current_boundaries,
                 osh::Real dt,
                 MolState<NumMolecules>& mol_state);

  EFieldOperator(const EFieldOperator &) = delete;

  ~EFieldOperator();

  /**
   *
   * \return get maximum E-Field time step
   */
  inline osh::Real getDt() const noexcept { return dt_; }

  /**
   *
   * \return set maximum E-Field time step
   */
  inline void setDt(const osh::Real dt) noexcept { dt_ = dt; }

    /**
     * \brief Set tolerances (relative and absolute) and on what norm they are used
     *
     * From the manual, the standard tolerances are:
     * - rtol = 1e-5
     * - atol = 1e-50
     *
     *
     */
  void setTolerances(double atol, double rtol, KSPNormType norm_type);

private:
  /**
   * \brief Initialize matrix and vectors
   *
   */
  void setupSystem();

  /** \brief Setup the stiffness matrix
   *
   * @param fixed_voltage_verts: set of indexes of the vertexes for which the voltage must remain
   * constant
   */
  void setupStiffnessMatrix();

  /// Track efield occupancy
  template <typename NumMolecules>
  void setupEfieldOccupancyTracking(MolState<NumMolecules>& mol_state);


  /** Apply ghk currents to i()
   *
   * @param ghk_currents
   */
  void apply_GHKcurrents(const osh::Reals& ghk_currents);

  /** Apply membrane-repated boundary conditions:
   *
   * - capacitance
   * - ohmic currents
   * - current injections
   *
   * @tparam NumMolecules
   * @param A0
   * @param mol_state
   * @param ghk_currents
   * @param sim_time
   * @param dt
   */
  template <typename NumMolecules>
  void apply_membrane_BC(Mat& A0,
                         const MolState<NumMolecules>& mol_state,
                         const osh::Real sim_time,
                         const osh::Real dt,
                         const osh::Write<osh::Real>& potential_on_verts);

  /** Add ohmic currents contributions
   *
   * The real input from the reaction-diffusion part is in the mol_state and the occupancies
   *
   * @tparam NumMolecules
   * @param tri_mat_and_vecs
   * @param membrane
   * @param b_id
   * @param mol_state
   * @param Avert
   * @param sim_time
   */
  template <typename NumMolecules>
  void add_ohmic_currents(TriMatAndVecs& tri_mat_and_vecs,
                          const Membrane& membrane,
                          const mesh::triangle_id_t& b_id,
                          const MolState<NumMolecules>& mol_state,
                          const double Avert,
                          const osh::Real sim_time,
                          const osh::Write<osh::Real>& potential_on_verts) const;

  /// Copy sol -> potential_on_verts
  void get_sol(osh::Write<osh::Real>& potential_on_verts);

  /** Init for the evolve routine
   *
   * Here we put to 0 all the relevant vectors and set A0 to be equal to the coupling matrix A()
   *
   * We do not finalize assembly so we can still add stuff
   *
   * \param A0
   * \param potential_on_verts
   * \param current_on_verts current injection on vertices
   */
  void evolve_init(Mat& A0,
                   const osh::Write<osh::Real>& potential_on_verts,
                   const osh::Read<osh::Real>& current_on_verts);

  /** Finalize assembly of the various vectors and matrices
   *
   * After this we cannot write into a particular element of a vector.
   *
   * @param A0
   */
  void finalize_assembly(Mat& A0);

  /// Combine the various vectors to create the rhs
  void build_rhs();

  /// Zeros appropriate rows and cols so that the voltages for the indexes in fixed_voltage_verts_
  /// remain constant
  void fix_voltages(Mat& A0);

public:
  /**
   * \brief Evolve the E-Field PDE over an interval of time dt
   *
   * \param potential_on_verts potential on vertices
   * \param current_on_verts current injection on vertices
   * \param mol_state molecular state
   * \param ghk_currents mapping from ghk reaction index to ghk current intensity
   * \param sim_time simulation time
   * \param dt time interval
   */
  template <typename NumMolecules>
  void evolve(osh::Write<osh::Real>& potential_on_verts,
              const osh::Read<osh::Real>& current_on_verts,
              const MolState<NumMolecules>& mol_state,
              const osh::Reals& ghk_currents,
              const osh::Real sim_time,
              const osh::Real dt);

  /// pretty printer
  //    template <typename NumMoleculesF>
  friend std::ostream& operator<<(std::ostream& ostr, const EFieldOperator& efo);

private:
  /// owned vertices
  osh::LOs owned_verts_;

  /// triangle to vertices
  osh::LOs tri2verts_;

  /// patch areas
  std::map<model::patch_id, PetscReal> patch_areas_;

  /// patch triangles
  std::map<model::patch_id, mesh::triangle_ids> patch_tris_;

  /// GHK current boundaries
  const std::vector<mesh::triangle_id_t> &ghk_current_boundaries_;

  /// indexes of the fixed voltages. Keep it sorted for log(N) insertions and possible petsc
  /// optimizations
  std::vector<PetscInt> fixed_voltage_verts_;
};

////--------------------------------------------------

// explicit template instantiation declarations

extern template void EFieldOperator::evolve(osh::Write<osh::Real>& potential_on_verts,
                                            const osh::Read<osh::Real>& current_on_verts,
                                            const MolState<osh::I32>& mol_state,
                                            const osh::Reals& ghk_currents,
                                            const osh::Real sim_time,
                                            const osh::Real dt);
extern template void EFieldOperator::evolve(osh::Write<osh::Real>& potential_on_verts,
                                            const osh::Read<osh::Real>& current_on_verts,
                                            const MolState<osh::I64>& mol_state,
                                            const osh::Reals& ghk_currents,
                                            const osh::Real sim_time,
                                            const osh::Real dt);

extern template EFieldOperator::EFieldOperator(
    DistMesh& o_mesh,
    const Statedef& statedef,
    const std::vector<mesh::triangle_id_t>& ghk_current_boundaries,
    osh::Real dt,
    MolState<osh::I32>& mol_state);
extern template EFieldOperator::EFieldOperator(
    DistMesh& o_mesh,
    const Statedef& statedef,
    const std::vector<mesh::triangle_id_t>& ghk_current_boundaries,
    osh::Real dt,
    MolState<osh::I64>& mol_state);

extern template void EFieldOperator::add_ohmic_currents(
    TriMatAndVecs& tri_mat_and_vecs,
    const Membrane& membrane,
    const mesh::triangle_id_t& b_id,
    const MolState<osh::I32>& mol_state,
    const double Avert,
    const osh::Real sim_time,
    const osh::Write<osh::Real>& potential_on_verts) const;

extern template void EFieldOperator::add_ohmic_currents(
    TriMatAndVecs& tri_mat_and_vecs,
    const Membrane& membrane,
    const mesh::triangle_id_t& b_id,
    const MolState<osh::I64>& mol_state,
    const double Avert,
    const osh::Real sim_time,
    const osh::Write<osh::Real>& potential_on_verts) const;


} // namespace dist
} // namespace steps
