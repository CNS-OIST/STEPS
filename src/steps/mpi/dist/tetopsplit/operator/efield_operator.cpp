#include "efield_operator.hpp"

#include <Omega_h_adj.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>
#include <petscmat.h>

#include "../mol_state.hpp"
#include "geom/dist/distmesh.hpp"
#include "util/mesh.hpp"
#include "util/petsc.hpp"
#include "util/profile/profiler_interface.hpp"


#define STRINGIFY(arg) #arg
#define AT             "at " __FILE__ ":" STRINGIFY(__LINE__) ": "

namespace steps::dist {

//----------------------------------------------

EFieldOperator::~EFieldOperator() {
    auto err = VecDestroy(&rhs());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecDestroy(&bc());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecDestroy(&i());
    CHKERRABORT(mesh.comm_impl(), err);
    err = KSPDestroy(&ksp_solver_);
    CHKERRABORT(mesh.comm_impl(), err);
    err = MatDestroy(&A());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecDestroy(&sol());
    CHKERRABORT(mesh.comm_impl(), err);
}

//----------------------------------------------
EFieldOperator::EFieldOperator(DistMesh& o_mesh,
                               const Statedef& statedef,
                               const std::vector<mesh::triangle_id_t>& ghk_current_boundaries,
                               MolState& mol_state)
    : state_def(statedef)
    , mesh(o_mesh)
    , global_indices_(mesh.global_indices(osh::VERT))
    , owned_verts_(osh::collect_marked(mesh.owned_verts_mask()))
    , tri2verts_(mesh.ask_verts_of(osh::FACE))
    , ghk_current_boundaries_(ghk_current_boundaries) {
    setupSystem();
    setupStiffnessMatrix();
    setupEfieldOccupancyTracking(mol_state);
}

//----------------------------------------------

void EFieldOperator::setupSystem() {
    Vec& sol = sol_;
    Vec& rhs = rhs_;
    Mat& A = A_;
    Vec& bc = bc_;
    Vec& i = i_;

    [[maybe_unused]] const auto num_global_vertices = mesh.total_num_verts();
    const auto num_owned_verts = mesh.num_verts();

    // we use PETSC_DETERMINE because global_indices presents only the global indices of owned and
    // ghost vertices, not everything
    auto err = VecCreateMPI(mesh.comm_impl(), num_owned_verts, PETSC_DETERMINE, &sol);
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecDuplicate(sol, &rhs);
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecDuplicate(sol, &bc);
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecDuplicate(sol, &i);
    CHKERRABORT(mesh.comm_impl(), err);
    err = MatCreateAIJ(mesh.comm_impl(),
                       num_owned_verts,
                       num_owned_verts,
                       PETSC_DETERMINE,
                       PETSC_DETERMINE,
                       100,
                       PETSC_NULL,
                       100,
                       PETSC_NULL,
                       &A);
    CHKERRABORT(mesh.comm_impl(), err);

    PetscInt nGlobalVerticesPETSc;
    err = VecGetSize(sol, &nGlobalVerticesPETSc);
    CHKERRABORT(mesh.comm_impl(), err);
    assert(nGlobalVerticesPETSc == num_global_vertices);
    assert(global_indices_.size() == mesh.owned_verts_mask().size());
    {
        ISLocalToGlobalMapping globals_vPETSc;  // SM: IS = Index Set
        err = ISLocalToGlobalMappingCreate(PETSC_COMM_SELF,
                                           1,
                                           global_indices_.size(),
                                           mesh::petsc_pointer(mesh::petsc_cast(global_indices_)),
                                           PETSC_COPY_VALUES,
                                           &globals_vPETSc);
        CHKERRABORT(mesh.comm_impl(), err);
        err = VecSetLocalToGlobalMapping(bc, globals_vPETSc);
        CHKERRABORT(mesh.comm_impl(), err);
        err = VecSetLocalToGlobalMapping(i, globals_vPETSc);
        CHKERRABORT(mesh.comm_impl(), err);
        err = MatSetLocalToGlobalMapping(A, globals_vPETSc, globals_vPETSc);
        CHKERRABORT(mesh.comm_impl(), err);
        err = ISLocalToGlobalMappingDestroy(&globals_vPETSc);
        CHKERRABORT(mesh.comm_impl(), err);
    }
    err = KSPCreate(mesh.comm_impl(), &ksp_solver_);
    CHKERRABORT(mesh.comm_impl(), err);

    // setup default options
    err = KSPSetType(ksp_solver_, KSPPIPECG);
    CHKERRABORT(mesh.comm_impl(), err);
    err = KSPGetPC(ksp_solver_, &ksp_solver_preconditioner_);
    CHKERRABORT(mesh.comm_impl(), err);
    err = PCSetType(ksp_solver_preconditioner_, PCPBJACOBI);
    CHKERRABORT(mesh.comm_impl(), err);

    // to not override CLI options, the following needs
    // to be called after all custom sets
    setPetscOptions();
}


void EFieldOperator::setPetscOptions() {
    auto err = KSPSetFromOptions(ksp_solver_);
    CHKERRABORT(mesh.comm_impl(), err);
    err = PCSetFromOptions(ksp_solver_preconditioner_);
    CHKERRABORT(mesh.comm_impl(), err);
}

std::ostream& operator<<(std::ostream& ostr, const EFieldOperator& efo) {
    return ostr << "A_: " << efo.A_ << "bc_: " << efo.bc_ << "i_: " << efo.i_
                << "rhs_: " << efo.rhs_;
}

void EFieldOperator::add_ohmic_currents(TriMatAndVecs& tri_mat_and_vecs,
                                        const Membrane& membrane,
                                        const mesh::triangle_id_t& b_id,
                                        const MolState& mol_state,
                                        const double Avert,
                                        const osh::Real sim_time,
                                        const osh::Write<osh::Real>& potential_on_verts) const {
    PetscReal avgv = std::accumulate(tri_mat_and_vecs.face_bf2vertsPETSc.begin(),
                                     tri_mat_and_vecs.face_bf2vertsPETSc.end(),
                                     0.0,
                                     [&potential_on_verts](PetscReal curr_sum, PetscInt i) {
                                         return curr_sum + potential_on_verts[i];
                                     });
    avgv /= 3;

    // add ohmic currents
    for (const auto& chan_pair: membrane.channels()) {
        const auto& c = chan_pair.second;
        for (const auto& h: c.ohmic_currents) {
            const PetscReal tri_oc_bc = h.get().getTriBConVertex(b_id, mol_state, Avert, sim_time);
            for (auto ir = 0u; ir < tri_mat_and_vecs.triBC.size(); ++ir) {
                tri_mat_and_vecs.triBC[ir] += tri_oc_bc *
                                              (h.get().getReversalPotential(b_id) - avgv);

                // to switch to an implicit scheme you need to uncomment the following line leave
                // the previous one as it is. Since we could want to do the switch in the future I
                // leave it here, for now, comented Katta
                //
                //  std::for_each(tri_mat_and_vecs.triStiffnessPETSc.begin(),
                //  tri_mat_and_vecs.triStiffnessPETSc.end(), [&tri_oc_bc](PetscReal& i){ i +=
                //  tri_oc_bc/3; });
            }
        }
    }
}

void EFieldOperator::add_leaks(TriMatAndVecs& tri_mat_and_vecs,
                               const double Avert,
                               const double conductivity,
                               const double reversal_potential,
                               const osh::Reals& potential_on_verts) const {
    // Add leaks
    const PetscReal tri_oc_bc = Avert * conductivity;
    // The area used in steps4 (called Avert) is different from the area
    // used in steps 3 (called pSurface).
    // In steps4, Avert is the area of a triangle (touching that vertex) divided by 3.
    // Every vertex touches n triangles and thus it has n areas associated to it.
    // The area has units in m^2.
    // In steps3, pSurface is the sum of the areas of all triangles touching the same vertex,
    // divided by 3. The units are in um^2.

    for (auto ir = 0u; ir < tri_mat_and_vecs.triBC.size(); ++ir) {
        const PetscInt ve_id = tri_mat_and_vecs.face_bf2vertsPETSc[ir];
        tri_mat_and_vecs.triBC[ir] += tri_oc_bc * (reversal_potential - potential_on_verts[ve_id]);
        tri_mat_and_vecs.triStiffnessPETSc[ir * 4] += tri_oc_bc;
    }
}

void EFieldOperator::apply_membrane_BC(
    Mat& A0,
    const MolState& mol_state,
    const osh::Real sim_time,
    const osh::Real dt,
    const osh::Write<osh::Real>& potential_on_verts,
    const osh::Read<osh::Real>& capacitance_on_triangles,
    const osh::Read<osh::Real>& conductivity_on_triangles,
    const osh::Read<osh::Real>& reversal_potential_on_triangles) {
    for (const auto& memb_pair: state_def.membranes()) {
        // IDs
        const auto& membrane = *memb_pair.second;
        const auto& patch_id = membrane.getPatch();
        const auto& patch_tris = patch_tris_[patch_id];
        // useful data required later
        const auto current_density = membrane.stimulus()(sim_time) / patch_areas_[patch_id];

        const auto applyBC = OMEGA_H_LAMBDA(osh::LO triangle_idx) {
            // IDS
            const mesh::triangle_id_t b_id{patch_tris[triangle_idx]};
            const auto& face_bf2verts = osh::gather_verts<3>(tri2verts_, b_id.get());

            // A tri split among the vertexes
            const double Avert = mesh.getTri(b_id).area / 3.0;
            // capacitance
            const PetscReal tri_capacitance = Avert * capacitance_on_triangles[b_id.get()] / dt;
            // current injection
            const auto tri_i = current_density * Avert;
            // create local matrices and vectors

            TriMatAndVecs tri_mat_and_vecs(face_bf2verts, tri_capacitance, tri_i);
            // add ohmic currents
            add_ohmic_currents(
                tri_mat_and_vecs, membrane, b_id, mol_state, Avert, sim_time, potential_on_verts);
            // add leakage
            add_leaks(tri_mat_and_vecs,
                      Avert,
                      conductivity_on_triangles[b_id.get()],
                      reversal_potential_on_triangles[b_id.get()],
                      potential_on_verts);
            // apply
            util::petsc::scalars triBC(tri_mat_and_vecs.triBC);
            auto lerr = VecSetValuesLocal(bc(),
                                          tri_mat_and_vecs.face_bf2vertsPETSc.size(),
                                          tri_mat_and_vecs.face_bf2vertsPETSc.data(),
                                          triBC.data(),
                                          ADD_VALUES);
            CHKERRABORT(mesh.comm_impl(), lerr);
            util::petsc::scalars triI(tri_mat_and_vecs.triI);
            lerr = VecSetValuesLocal(i(),
                                     tri_mat_and_vecs.face_bf2vertsPETSc.size(),
                                     tri_mat_and_vecs.face_bf2vertsPETSc.data(),
                                     triI.data(),
                                     ADD_VALUES);
            CHKERRABORT(mesh.comm_impl(), lerr);
            util::petsc::scalars triStiffnessPETSc(tri_mat_and_vecs.triStiffnessPETSc);
            lerr = MatSetValuesLocal(A0,
                                     tri_mat_and_vecs.face_bf2vertsPETSc.size(),
                                     tri_mat_and_vecs.face_bf2vertsPETSc.data(),
                                     tri_mat_and_vecs.face_bf2vertsPETSc.size(),
                                     tri_mat_and_vecs.face_bf2vertsPETSc.data(),
                                     triStiffnessPETSc.data(),
                                     ADD_VALUES);
            CHKERRABORT(mesh.comm_impl(), lerr);
        };
        osh::parallel_for(patch_tris.size(), applyBC);
    }
}


void EFieldOperator::apply_GHKcurrents(const osh::Reals& ghk_currents) {
    const auto applyGHKCurrents = OMEGA_H_LAMBDA(osh::LO reaction_idx) {
        const mesh::triangle_id_t b_id{ghk_current_boundaries_[static_cast<size_t>(reaction_idx)]};
        const auto& face_bf2verts = osh::gather_verts<3>(tri2verts_, b_id.get());
        const auto tri_i = -ghk_currents[static_cast<size_t>(reaction_idx)] / 3.0;
        std::array<PetscScalar, 3> triI{tri_i, tri_i, tri_i};
        std::array<PetscInt, 3> face_bf2vertsPETSc{static_cast<PetscInt>(face_bf2verts[0]),
                                                   static_cast<PetscInt>(face_bf2verts[1]),
                                                   static_cast<PetscInt>(face_bf2verts[2])};
        auto err =
            VecSetValuesLocal(i(), triI.size(), face_bf2vertsPETSc.data(), triI.data(), ADD_VALUES);
        CHKERRABORT(mesh.comm_impl(), err);
    };
    osh::parallel_for(static_cast<osh::LO>(ghk_current_boundaries_.size()), applyGHKCurrents);
}

void EFieldOperator::get_sol(osh::Write<osh::Real>& potential_on_verts) {
    auto copySolution = OMEGA_H_LAMBDA(osh::LO r) {
        PetscScalar val;
        auto lerr = mesh::petsc_get_values(sol(), global_indices_[owned_verts_[r]], val);
        potential_on_verts[owned_verts_[r]] += util::petsc::to_real(val);
        CHKERRABORT(mesh.comm_impl(), lerr);
    };
    osh::parallel_for(owned_verts_.size(), copySolution);
    const auto& solOmega = mesh.sync_array(osh::VERT, osh::Reals(potential_on_verts), 1);
    std::copy(solOmega.begin(), solOmega.end(), potential_on_verts.begin());
}

void EFieldOperator::evolve_init(Mat& A0,
                                 const osh::Write<osh::Real>& potential_on_verts,
                                 const osh::Read<osh::Real>& current_on_verts) {
    std::vector<PetscInt> idxs(static_cast<size_t>(potential_on_verts.size()));
    std::iota(idxs.begin(), idxs.end(), 0);
    util::petsc::scalars pov(potential_on_verts);
    auto err = VecSetValuesLocal(sol(),
                                 mesh::petsc_int_cast(potential_on_verts.size()),
                                 idxs.data(),
                                 pov.data(),
                                 INSERT_VALUES);
    CHKERRABORT(mesh.comm_impl(), err);

    err = VecZeroEntries(rhs());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecZeroEntries(bc());
    CHKERRABORT(mesh.comm_impl(), err);
    util::petsc::scalars cov(current_on_verts);
    err = VecSetValuesLocal(i(),
                            mesh::petsc_int_cast(potential_on_verts.size()),
                            idxs.data(),
                            cov.data(),
                            INSERT_VALUES);
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyBegin(i());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyEnd(i());
    CHKERRABORT(mesh.comm_impl(), err);

    err = VecSetValuesLocal(sol(),
                            mesh::petsc_int_cast(potential_on_verts.size()),
                            idxs.data(),
                            pov.data(),
                            INSERT_VALUES);
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyBegin(sol());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyEnd(sol());
    CHKERRABORT(mesh.comm_impl(), err);

    err = MatDuplicate(A(), MAT_COPY_VALUES, &A0);
    CHKERRABORT(mesh.comm_impl(), err);
}

void EFieldOperator::finalize_assembly(Mat& A0) {
    auto err = MatAssemblyBegin(A0, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(mesh.comm_impl(), err);
    err = MatAssemblyEnd(A0, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyBegin(bc());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyEnd(bc());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyBegin(i());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyEnd(i());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyBegin(sol());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyEnd(sol());
    CHKERRABORT(mesh.comm_impl(), err);
}

void EFieldOperator::build_rhs() {
    auto err = MatMult(A(), sol(), rhs());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecScale(rhs(), -1.0);
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAXPY(rhs(), 1.0, bc());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAXPY(rhs(), 1.0, i());
    CHKERRABORT(mesh.comm_impl(), err);
}


void EFieldOperator::fix_voltages(Mat& A0) {
    auto err = VecZeroEntries(sol());
    CHKERRABORT(mesh.comm_impl(), err);

    err = VecAssemblyBegin(sol());
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyEnd(sol());
    CHKERRABORT(mesh.comm_impl(), err);

    err = MatZeroRowsColumnsLocal(
        A0, fixed_voltage_verts_.size(), fixed_voltage_verts_.data(), 1.0, sol(), rhs());
    CHKERRABORT(mesh.comm_impl(), err);
}


void EFieldOperator::evolve(osh::Write<osh::Real>& potential_on_verts,
                            const osh::Read<osh::Real>& current_on_verts,
                            const osh::Read<osh::Real>& capacitance_on_triangles,
                            const osh::Read<osh::Real>& conductivity_on_triangles,
                            const osh::Read<osh::Real>& reversal_potential_on_triangles,
                            const MolState& mol_state,
                            const osh::Reals& ghk_currents,
                            const osh::Real sim_time,
                            const osh::Real dt) {
    Instrumentor::phase p("EFieldOperator::evolve()");

    // set A0 = A() (coupling matrix), all the vectors to 0, copy old sol into sol and finalize it
    // to that it can already be used
    Mat A0;
    evolve_init(A0, potential_on_verts, current_on_verts);

    // setup boundary conditions
    apply_membrane_BC(A0,
                      mol_state,
                      sim_time,
                      dt,
                      potential_on_verts,
                      capacitance_on_triangles,
                      conductivity_on_triangles,
                      reversal_potential_on_triangles);
    apply_GHKcurrents(ghk_currents);

    // finalize vector building
    finalize_assembly(A0);

    // build rhs
    build_rhs();

    // override A0 and rhs() so that the voltages for the vertexes in fixed_voltage_vertexes_ are
    // fixed. This must be the last operation on the vectors before solving. The solution vector is
    // zeroed for convenience
    fix_voltages(A0);

    // Update matrix
    auto err = KSPSetOperators(ksp_solver_, A0, A0);
    CHKERRABORT(mesh.comm_impl(), err);

    err = KSPSolve(ksp_solver_, rhs(), sol());
    CHKERRABORT(mesh.comm_impl(), err);
    KSPConvergedReason reason;
    err = KSPGetConvergedReason(ksp_solver_, &reason);
    CHKERRABORT(mesh.comm_impl(), err);

    if (reason <= 0) {
        throw std::logic_error(AT "PETSc Krylov solver not converged.");
    }

    // copy back solution
    get_sol(potential_on_verts);

    // destroy matrix
    err = MatDestroy(&A0);
    CHKERRABORT(mesh.comm_impl(), err);
}

//----------------------------------------------

void EFieldOperator::setupStiffnessMatrix() {
    const auto& coords = mesh.coords();
    osh::Matrix<4, 3> P;
    for (auto ic = 1; ic < 4; ++ic) {
        P[ic - 1] = osh::Vector<4>({-1, 0, 0, 0});
        P[ic - 1][ic] = 1;
    }

    Vec fixed_voltage_verts_petsc;
    auto err = VecDuplicate(rhs(), &fixed_voltage_verts_petsc);
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecSet(fixed_voltage_verts_petsc, 1.0);
    CHKERRABORT(mesh.comm_impl(), err);

    std::unordered_set<model::compartment_id> comp_ids;
    const std::array<PetscScalar, 4> fixed_voltages_mask_stencil = {0.0, 0.0, 0.0, 0.0};
    for (const auto& memb_pair: state_def.membranes()) {
        const auto& membrane = memb_pair.second;
        auto patch_id = membrane->getPatch();
        patch_tris_[patch_id] = mesh.getOwnedEntities(patch_id);
        patch_areas_[patch_id] = mesh.total_measure(patch_id);
        const auto& patchdef = state_def.getPatchdef(patch_id);
        auto inner_comp_id = patchdef.getInnerCompId();
        auto inner_comp_conductivity = state_def.getCompartmentConductivity(inner_comp_id);
        auto owned_tets = mesh.getOwnedEntities(inner_comp_id);
        const auto& tets2verts = mesh.ask_elem_verts();
        auto assembleSystem = OMEGA_H_LAMBDA(osh::LO idx) {
            auto j = owned_tets[idx].get();
            const auto& tet_j2verts = osh::gather_verts<4>(tets2verts, j);
            auto tet_j2x = osh::gather_vectors<4, 3>(coords, tet_j2verts);  // SM: coords of the
            // vertices in tet_j2verts
            auto M = osh::simplex_basis<3, 3>(tet_j2x);
            if (osh::cross(M[0], M[1]) * M[2] <= 0.0) {
                throw std::logic_error(AT "Wrong setup.");
            }
            osh::Matrix<3, 4> N;
            for (auto ic = 1; ic < 4; ++ic) {
                N[ic] = osh::cross(M[ic % 3], M[(ic + 1) % 3]);
            }
            N[0] = -N[1] - N[2] - N[3];  // SM: the divergence theorem on any constant
            // vector field imposes their sum to be 0
            const auto& grad_phi = P * invert(M);  // SM: grad_phi[ic][ir] = dphi_r/dx_c
            std::array<PetscInt, 4> tet_j2vertsPETSc{};
            std::array<PetscScalar, 16> triStiffnessPETSc{};
            triStiffnessPETSc.fill(0);


            for (auto ir = 0u; ir < tet_j2vertsPETSc.size(); ++ir) {
                tet_j2vertsPETSc[ir] = static_cast<PetscInt>(tet_j2verts[static_cast<osh::LO>(ir)]);
                for (auto jr = ir + 1; jr < tet_j2vertsPETSc.size(); ++jr) {
                    auto nDotGrad_phi_times_cond =
                        inner_comp_conductivity * grad_phi *
                        (N[static_cast<osh::LO>(jr)] - N[static_cast<osh::LO>(ir)]) / 24.0;
                    for (auto ic = 0u; ic < tet_j2vertsPETSc.size(); ++ic) {
                        triStiffnessPETSc[ir + 4 * ic] -=
                            nDotGrad_phi_times_cond[static_cast<osh::LO>(ic)];
                        triStiffnessPETSc[jr + 4 * ic] +=
                            nDotGrad_phi_times_cond[static_cast<osh::LO>(ic)];
                    }
                }
            }
            auto err0 = MatSetValuesLocal(A(),
                                          tet_j2vertsPETSc.size(),
                                          tet_j2vertsPETSc.data(),
                                          tet_j2vertsPETSc.size(),
                                          tet_j2vertsPETSc.data(),
                                          triStiffnessPETSc.data(),
                                          ADD_VALUES);
            CHKERRABORT(mesh.comm_impl(), err0);

            // these indexes are touched by an inner compartment. We do not need to mark them as
            // fixed
            err0 = VecSetValuesLocal(fixed_voltage_verts_petsc,
                                     tet_j2vertsPETSc.size(),
                                     tet_j2vertsPETSc.data(),
                                     fixed_voltages_mask_stencil.data(),
                                     INSERT_VALUES);
            CHKERRABORT(mesh.comm_impl(), err0);
        };
        auto comp_ids_insert_result = comp_ids.insert(inner_comp_id);
        if (comp_ids_insert_result.second) {
            osh::parallel_for(owned_tets.size(), assembleSystem);
        }
    }

    err = VecAssemblyBegin(fixed_voltage_verts_petsc);
    CHKERRABORT(mesh.comm_impl(), err);
    err = VecAssemblyEnd(fixed_voltage_verts_petsc);
    CHKERRABORT(mesh.comm_impl(), err);

    // set matrix diagonal to 1s where we have fixed voltages so it does not shrink
    err = MatDiagonalSet(A(), fixed_voltage_verts_petsc, ADD_VALUES);
    CHKERRABORT(mesh.comm_impl(), err);

    // copy results: fixed_voltage_verts_petsc -> fixed_voltage_verts_
    auto copyIntoFixedVoltageVerts = OMEGA_H_LAMBDA(osh::LO idx) {
        PetscScalar val;
        const auto local_idx = owned_verts_[idx];
        const auto global_idx = global_indices_[local_idx];
        auto lerr = mesh::petsc_get_values(fixed_voltage_verts_petsc, global_idx, val);
        CHKERRABORT(mesh.comm_impl(), lerr);
        if (util::petsc::to_real(val) != 0.0) {
            fixed_voltage_verts_.emplace_back(local_idx);
        }
    };
    osh::parallel_for(owned_verts_.size(), copyIntoFixedVoltageVerts);
    std::sort(fixed_voltage_verts_.begin(), fixed_voltage_verts_.end());

    err = MatAssemblyBegin(A(), MAT_FINAL_ASSEMBLY);
    CHKERRABORT(mesh.comm_impl(), err);
    err = MatAssemblyEnd(A(), MAT_FINAL_ASSEMBLY);
    CHKERRABORT(mesh.comm_impl(), err);

    // additional_fixed_voltage_verts_ is just temporary, we destroy it
    err = VecDestroy(&fixed_voltage_verts_petsc);
    CHKERRABORT(mesh.comm_impl(), err);
}

void EFieldOperator::setupEfieldOccupancyTracking(MolState& mol_state) {
    // register channels for ef occupancy tracking
    for (const auto& memb_pair: state_def.membranes()) {
        const auto& membrane = memb_pair.second;
        for (const auto& chan_pair: membrane->channels()) {
            const auto& c = chan_pair.second;
            for (const auto& h: c.ohmic_currents) {
                if (h.get().channel_state) {
                    const auto& patch_id = membrane->getPatch();
                    const auto& patch_tris = mesh.getOwnedEntities(patch_id);
                    for (const auto tri_id: patch_tris) {
                        mol_state.track_occupancy_ef(tri_id, *h.get().channel_state);
                    }
                }
            }
        }
    }
}

}  // namespace steps::dist
