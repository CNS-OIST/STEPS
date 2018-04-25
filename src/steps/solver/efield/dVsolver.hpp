/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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


#ifndef STEPS_SOLVER_EFIELD_DVSOLVER_HPP
#define STEPS_SOLVER_EFIELD_DVSOLVER_HPP 1

// STL headers.
#include <algorithm>
#include <memory>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/efield/bdsystem.hpp"
#include "steps/solver/efield/efieldsolver.hpp"
#include "steps/solver/efield/tetmesh.hpp"
#include "steps/solver/efield/vertexconnection.hpp"
#include "steps/solver/efield/vertexelement.hpp"

namespace steps {
namespace solver {
namespace efield {

class dVSolverBase: public EFieldSolver {
public:
    dVSolverBase(): pMesh(0), pNVerts(0), pNTris(0) {}

    /** Initialize state with given mesh */
    void initMesh(TetMesh *mesh) override;

    /** Set membrane conductance and reversal potential (for leak current) */
    void setSurfaceConductance(double g_surface, double v_rev) override;
    
    /** Set all vertex potentials to v */
    void setPotential(double v) override {
        std::fill(pV.begin(),pV.end(),v);
    }

    /** Retrieve potential at vertex i */
    double getV(int i) const override { return pV[i]; }

    /** Set potential at vertex i */
    void setV(int i, double v) override { pV[i]=v; }

    /** Get voltage clamped status for vertex i */
    bool getClamped(int i) const override { return pVertexClamp[i]; }

    /** Set voltage clamped status for vertex i */
    void setClamped(int i, bool clamped) override { pVertexClamp[i] = clamped; }

    /** Get current through triangle i */
    double getTriI(int i) const override { return -pTriCur[i]; }

    /** Set current through triangle i to d (pA) */
    void setTriI(int i,double d) override { pTriCur[i] = -d; }

    /** Set additional current injection for triangle i to c (pA) */
    void setTriIClamp(int i, double c) override { pTriCurClamp[i] = -c; }

    /** Set additional current injection for area associated with vertex i to c (pA) */
    void setVertIClamp(int i, double c) override { pVertCurClamp[i] = -c; }

protected:
    /// Generic populate and solve
    template <typename LinSysImpl>
    void _advance(LinSysImpl *L, double dt) {
        // Add up current clamp contributions
        std::copy(pVertCurClamp.begin(), pVertCurClamp.end(), pVertCur.begin());
        for (uint i = 0; i < pNTris; ++i) {
            double c = (pTriCur[i] + pTriCurClamp[i]) / 3.0;

            uint *triv = pMesh->getTriangle(i);
            pVertCur[triv[0]] += c;
            pVertCur[triv[1]] += c;
            pVertCur[triv[2]] += c;
        }

        typename LinSysImpl::matrix_type &A=L->A();
        typename LinSysImpl::vector_type &b=L->b();

        double oodt = 1.0/dt;

        A.zero();
        for (uint i = 0; i < pNVerts; ++i) {
            VertexElement * ve = pMesh->getVertex(i);
            int ind = ve->getIDX();

            if (pVertexClamp[ind]) {
                b.set(ind,0);
                A.set(ind,ind,1.0);
            }
            else {
                double rhs = pVertCur[ind] + pGExt[ind] * (pVExt - pV[ind]);
                double Aii = ve->getCapacitance()*oodt + pGExt[ind];

                for (int inbr = 0; inbr < ve->getNCon(); ++inbr) {
                    int k = ve->nbrIdx(inbr);
                    double cc = ve->getCC(inbr);

                    rhs += cc * (pV[k] - pV[ind]);
                    Aii += cc;
                    A.set(ind,k,-cc);
                }
                b.set(ind,rhs);
                A.set(ind,ind,Aii);
            }
        }
        
        L->solve();

        const typename LinSysImpl::vector_type DV=L->x();
        for (uint i = 0; i < pNVerts; ++i)
            if (pVertexClamp[i] == false) pV[i] += DV.get(i);

        // reset pTriCur for caller contributions
        std::fill(pTriCur.begin(), pTriCur.end(), 0.0);
    }

    /// Compute matrix half-bw from mesh
    static int meshHalfBW(TetMesh *mesh);

    /// Pointer to the mesh.
    TetMesh *                   pMesh;

    /// Number of vertices in the mesh, stored locally.
    uint                        pNVerts;

    /// Number of triangles in the mesh, stored locally.
    uint                        pNTris;

    /// The local potential over all mesh vertex points.
    std::vector<double>         pV;

    /// Leak conductance at vertices
    std::vector<double>         pGExt;

    /// Reversal potential for leak current.
    double                      pVExt;

    /// Clamped status for each vertex (non-zero => clamped.)
    std::vector<char>           pVertexClamp;

    /// Current through each triangle.
    std::vector<double>         pTriCur;

    /// Current clamps through each triangle.
    std::vector<double>         pTriCurClamp;

    /// Current contribution at each vertex (used only in advance()).
    std::vector<double>         pVertCur;

    /// Current clamp through each vertex (adds to any triangle clamps.)
    std::vector<double>         pVertCurClamp;
};
    
class dVSolverBanded: public dVSolverBase {
public:
    void initMesh(TetMesh *mesh) override {
        dVSolverBase::initMesh(mesh);
        pBDSys.reset(new BDSystem(pNVerts, meshHalfBW(mesh)));
    }

    void advance(double dt) override {
        _advance(pBDSys.get(), dt);
    }

private:
    std::unique_ptr<BDSystem>  pBDSys;
};


}}} // namespace steps::efield::solver

#endif // ndef STEPS_SIM_EFIELD_DVSOLVER_HPP

// END
