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

#include "mpi/tetvesicle/endocytosis.hpp"

// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include "math/point.hpp"
#include "mpi/tetvesicle/comp_vesraft.hpp"
#include "mpi/tetvesicle/tetvesicle_vesraft.hpp"
#include "mpi/tetvesicle/tri_vesraft.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

Endocytosis::Endocytosis(solver::Endocytosisdef* endodef, std::vector<TriVesRaft*>& tris)
    : pEndocytosisdef(endodef)
    , pCcst(0.0)
    , pKcst(0.0)
    , rExtent(0)
    , pActive(true) {
    AssertLog(pEndocytosisdef != nullptr);

    uint ntris = tris.size();
    AssertLog(ntris > 0);

    pTris.assign(tris.begin(), tris.end());

    double kcst = pEndocytosisdef->kcst();
    AssertLog(kcst >= 0.0);

    pKcst = kcst;

    // It's a type of 1st order reaction
    pCcst = kcst;

    for (uint i = 0; i < ntris; ++i) {
        // Set the position of endocytosis as on the vector of barycentre of the
        // triangle to inner tet barycentre, length radius of vesicle.

        // Now wmvols and tris hold their position, so this should be a piece of
        // piss

        TetVesRaft* tet;
        if (inner()) {
            tet = pTris[i]->iTet();
        } else {
            tet = pTris[i]->oTet();
        }

        math::point3d tet_baryc = tet->position();
        math::point3d tri_baryc = pTris[i]->position();
        math::point3d tri_norm = pTris[i]->norm();

        solver::vesicle_global_id ves_gidx = pEndocytosisdef->rhs_I_ves_uint();

        solver::Vesicledef* vesdef = pEndocytosisdef->statedef()->vesicledef(ves_gidx);

        double ves_diameter = vesdef->diameter();

        // 1) Find tri normal from the tri vertices.
        // 2) Compare dot product of tri normal to baryc-baryc line. If it's
        // negative, flip it. 3) Calculate unit normal 4) Set endo position as
        // radius * unit normal. Original idea was to go along the baryc-baryc line
        // but that usually gives boundary overlap

        // use this function: (something) + math::tri_normal(tri_verts[0],
        // tri_verts[1], tri_verts[2]);

        math::point3d baryc_baryc_vec = tet_baryc - tri_baryc;

        double dotproduct = math::dot(baryc_baryc_vec, tri_norm);

        if (dotproduct < 0.0) {
            tri_norm[0] = -tri_norm[0];
            tri_norm[1] = -tri_norm[1];
            tri_norm[2] = -tri_norm[2];
        }

        // Tri normal is already normalised as returned (mg 1)
        pPos.push_back(tri_baryc + tri_norm * (ves_diameter / 2.0));
    }
}

////////////////////////////////////////////////////////////////////////////////

Endocytosis::~Endocytosis() = default;

////////////////////////////////////////////////////////////////////////////////

void Endocytosis::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pCcst);
    util::checkpoint(cp_file, pKcst);
    util::checkpoint(cp_file, rExtent);
    util::checkpoint(cp_file, pEvents);
    util::checkpoint(cp_file, pActive);
    util::checkpoint(cp_file, pPos);
}

////////////////////////////////////////////////////////////////////////////////

void Endocytosis::restore(std::fstream& cp_file) {
    util::restore(cp_file, pCcst);
    util::restore(cp_file, pKcst);
    util::restore(cp_file, rExtent);
    util::restore(cp_file, pEvents);
    util::restore(cp_file, pActive);
    util::compare(cp_file, pPos);
}

////////////////////////////////////////////////////////////////////////////////

void Endocytosis::reset() {
    resetExtent();
    resetCcst();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void Endocytosis::resetCcst() {
    double kcst = pEndocytosisdef->kcst();
    AssertLog(kcst >= 0.0);

    pKcst = kcst;

    // It's always a 1st order reaction
    pCcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

void Endocytosis::setKcst(double k) {
    AssertLog(k >= 0.0);
    pKcst = k;
    // It's always a 1st order reaction
    pCcst = k;
}
////////////////////////////////////////////////////////////////////////////////

std::vector<solver::EndocytosisEvent> Endocytosis::getEvents() {
    std::vector<solver::EndocytosisEvent> copy(pEvents);
    pEvents.clear();
    return copy;
}

////////////////////////////////////////////////////////////////////////////////


void Endocytosis::setActive(bool active) {
    pActive = active;
}

////////////////////////////////////////////////////////////////////////////////

double Endocytosis::rate() const {
    if (inactive()) {
        return 0.0;
    }

    uint ntris = pTris.size();
    AssertLog(ntris > 0);
    solver::Patchdef* pdef = pTris[0]->patchdef();

    const auto& lhs_s_vec = endodef()->lhs_S();

    for (auto sg: lhs_s_vec.range()) {
        uint lhs = lhs_s_vec[sg];
        if (lhs == 0) {
            continue;
        }

        // We need a spec or specs. Let's sum over the tris
        uint cnt = 0;

        solver::spec_local_id spec_lidx = pdef->specG2L(sg);

        for (uint i = 0; i < ntris; ++i) {
            AssertLog(pTris[i]->patchdef() == pdef);
            cnt += pTris[i]->pools()[spec_lidx];
        }
        // We have summed all available species from each tri. Compare to required
        // lhs
        if (lhs > cnt) {
            //  The required species are not available
            return 0.0;
        }
    }

    // Special kind of reaction that is like a pseudo-first order reaction
    return pCcst;
}

////////////////////////////////////////////////////////////////////////////////

void Endocytosis::apply(TetVesicleVesRaft* solver) {
    AssertLog(solver != nullptr);

    // First see if we can do the endocytosis. If not just exit.
    solver::Vesicledef* ves_def = endodef()->rhs_I_ves();
    solver::vesicle_global_id ves_idx = endodef()->rhs_I_ves_uint();

    uint ntris = pTris.size();
    AssertLog(ntris > 0);
    solver::Patchdef* pdef = pTris[0]->patchdef();
    uint nspecs = pdef->countSpecs();

    // Get the patch and comp pointers but check consistency
    CompVesRaft* comp;
    if (inner()) {
        comp = pTris[0]->iTet()->getCompVesRaft();
    } else {
        comp = pTris[0]->oTet()->getCompVesRaft();
    }

    // Store map of species lidx in patch to population in triangles
    std::map<solver::spec_global_id, int> specs_vesicle;

    bool applied = false;

    // First we need to try the endocytosis:
    solver::vesicle_individual_id ves_unique_index;

    triangle_global_id triApplied;
    for (uint i = 0; i < ntris; ++i) {
        ves_unique_index = comp->addVesicle(ves_def, pPos[i]);

        if (ves_unique_index.valid()) {
            applied = true;
            triApplied = pTris[i]->idx();
            break;
        }
    }

    if (applied) {
        for (uint i = 0; i < ntris; ++i) {
            // If we got here then we have added a vesicle to the compartment. Need to
            // get ALL the species in ALL endocytotic zone triangles and associate
            // them with the vesicle (do that before icomp->addOneVesicle and add to
            // that function??) Then we need to remove all those species from those
            // triangles with SOLVER function setTriCount. pTris[i]->setCount won't do
            // it because the propensites don't get updated.. Could instead find all
            // the dependencies here and adaptively assign pUpdVec?

            AssertLog(pTris[i]->patchdef() == pdef);

            for (auto spec_lidx: solver::spec_local_id::range(nspecs)) {
                uint nmolcs = pTris[i]->pools()[spec_lidx];

                if (nmolcs > 0) {
                    // Think it's better if the vesicles store global species indices
                    solver::spec_global_id spec_gidx = pdef->specL2G(spec_lidx);
                    specs_vesicle[spec_gidx] += nmolcs;

                    // NEED SOMEHOW TO GET SURFACE DIFFUSION INFORMATION IN HERE

                    solver->setTriSpecCount_(pTris[i]->idx(), spec_gidx, 0);
                }
            }
        }

        // Need to add the species in specs_vesicle to the right vesicle in comp

        // When we added the vesicle the index was the total number of vesicles.
        // Ergo the index we need is the new number minus 1

        comp->addVesicleSpecs(ves_idx, ves_unique_index, specs_vesicle);

        addEvent(solver->getTime(), triApplied, ves_unique_index);

        // Try doing a global update for good measure
        // TODO do we need an update, or does solver->_setTriCount take care of all
        // the updates for us??
        // solver->_update();
    }
}

}  // namespace steps::mpi::tetvesicle
