/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */

#pragma once

// STL headers.
#include <algorithm>
#include <memory>
#include <vector>

// STEPS headers.
#include "bdsystem.hpp"
#include "efieldsolver.hpp"
#include "tetmesh.hpp"
#include "vertexconnection.hpp"
#include "vertexelement.hpp"

#include "util/checkpointing.hpp"

namespace steps::solver::efield {

class dVSolverBase: public EFieldSolver {
  public:
    dVSolverBase()
        : pMesh(nullptr)
        , pNVerts(0)
        , pNTris(0) {}

    /** Initialize state with given mesh */
    void initMesh(TetMesh* mesh) override;

    /// checkpoint data
    void checkpoint(std::fstream& cp_file) override;

    /// restore data
    void restore(std::fstream& cp_file) override;

    /** Set membrane conductance and reversal potential (for leak current) */
    void setSurfaceConductance(double g_surface, double v_rev) override;

    /** Get membrane conductance and reversal potential (for leak current) */
    std::pair<double, double> getSurfaceConductance() override;

    /** Set all vertex potentials to v */
    void setPotential(double v) override {
        std::fill(pV.begin(), pV.end(), v);
    }

    /** Retrieve potential at vertex i */
    double getV(vertex_id_t i) const noexcept override {
        return pV[i.get()];
    }

    /** Set potential at vertex i */
    void setV(vertex_id_t i, double v) noexcept override {
        pV[i.get()] = v;
    }

    /** Get voltage clamped status for vertex i */
    bool getClamped(vertex_id_t i) const noexcept override {
        return pVertexClamp[i.get()] != 0;
    }

    /** Set voltage clamped status for vertex i */
    void setClamped(vertex_id_t i, bool clamped) noexcept override {
        pVertexClamp[i.get()] = static_cast<char>(clamped);
    }

    /** Get current through triangle i */
    double getTriI(triangle_local_id i) const noexcept override {
        return -pTriCur[i.get()];
    }

    /** Set current through triangle i to d (pA) */
    void setTriI(triangle_local_id i, double d) noexcept override {
        pTriCur[i.get()] = -d;
    }

    /** Set additional current injection for triangle i to c (pA) */
    void setTriIClamp(triangle_local_id i, double c) noexcept override {
        pTriCurClamp[i.get()] = -c;
    }

    /** Get additional current injection for triangle i (pA) */
    double getTriIClamp(triangle_local_id i) const noexcept override {
        return pTriCurClamp[i.get()];
    }

    /** Set additional current injection for area associated with vertex i to c
     * (pA) */
    void setVertIClamp(vertex_id_t i, double c) noexcept override {
        pVertCurClamp[i.get()] = -c;
    }

    /** Get additional current injection for area associated with vertex i (pA) */
    double getVertIClamp(vertex_id_t i) const noexcept override {
        return pVertCurClamp[i.get()];
    }

  protected:
    /// Generic populate and solve
    template <typename LinSysImpl>
    void _advance(LinSysImpl* L, double dt) {
        // Add up current clamp contributions
        std::copy(pVertCurClamp.begin(), pVertCurClamp.end(), pVertCur.begin());
        for (uint i = 0; i < pNTris; ++i) {
            double c = (pTriCur[i] + pTriCurClamp[i]) / 3.0;

            const auto* triv = pMesh->getTriangle(triangle_local_id(i));
            pVertCur[triv[0].get()] += c;
            pVertCur[triv[1].get()] += c;
            pVertCur[triv[2].get()] += c;
        }

        typename LinSysImpl::matrix_type& A = L->A();
        typename LinSysImpl::vector_type& b = L->b();

        double oodt = 1.0 / dt;

        A.zero();
        for (auto i: vertex_id_t::range(pNVerts)) {
            VertexElement* ve = pMesh->getVertex(i);
            int ind = ve->getIDX();

            if (pVertexClamp[ind]) {
                b.set(ind, 0);
                A.set(ind, ind, 1.0);
            } else {
                double rhs = pVertCur[ind] + pGExt[ind] * (pVExt - pV[ind]);
                double Aii = ve->getCapacitance() * oodt + pGExt[ind];

                for (auto inbr = 0u; inbr < ve->getNCon(); ++inbr) {
                    int k = ve->nbrIdx(inbr);
                    double cc = ve->getCC(inbr);

                    rhs += cc * (pV[k] - pV[ind]);
                    Aii += cc;
                    A.set(ind, k, -cc);
                }
                b.set(ind, rhs);
                A.set(ind, ind, Aii);
            }
        }

        L->solve();

        const typename LinSysImpl::vector_type DV = L->x();
        for (uint i = 0; i < pNVerts; ++i) {
            if (pVertexClamp[i] == false) {
                pV[i] += DV.get(i);
            }
        }

        // reset pTriCur for caller contributions
        std::fill(pTriCur.begin(), pTriCur.end(), 0.0);
    }

    /// Compute matrix half-bw from mesh
    static int meshHalfBW(TetMesh* mesh);

    /// Pointer to the mesh.
    TetMesh* pMesh;

    /// Number of vertices in the mesh, stored locally.
    uint pNVerts;

    /// Number of triangles in the mesh, stored locally.
    uint pNTris;

    /// The local potential over all mesh vertex points.
    std::vector<double> pV;

    /// Leak conductance at vertices
    std::vector<double> pGExt;

    /// Conductance for leak current
    double conductance;

    /// Reversal potential for leak current.
    double pVExt;

    /// Clamped status for each vertex (non-zero => clamped.)
    std::vector<char> pVertexClamp;

    /// Current through each triangle.
    std::vector<double> pTriCur;

    /// Current clamps through each triangle.
    std::vector<double> pTriCurClamp;

    /// Current contribution at each vertex (used only in advance()).
    std::vector<double> pVertCur;

    /// Current clamp through each vertex (adds to any triangle clamps.)
    std::vector<double> pVertCurClamp;
};

class dVSolverBanded: public dVSolverBase {
  public:
    void initMesh(TetMesh* mesh) override {
        dVSolverBase::initMesh(mesh);
        pBDSys.reset(new BDSystem(pNVerts, meshHalfBW(mesh)));
    }

    /// checkpoint data
    void checkpoint(std::fstream& cp_file) override {
        dVSolverBase::checkpoint(cp_file);
        pBDSys->checkpoint(cp_file);
    }

    /// restore data
    void restore(std::fstream& cp_file) override {
        dVSolverBase::restore(cp_file);
        pBDSys->restore(cp_file);
    }

    void advance(double dt) override {
        _advance(pBDSys.get(), dt);
    }

  private:
    std::unique_ptr<BDSystem> pBDSys;
};

}  // namespace steps::solver::efield
