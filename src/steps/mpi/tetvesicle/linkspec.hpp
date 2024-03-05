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

// Standard library & STL headers.
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "math/constants.hpp"
#include "math/point.hpp"
#include "mpi/tetvesicle/linkspecpair.hpp"
#include "mpi/tetvesicle/vesicle.hpp"
#include "solver/linkspecdef.hpp"

namespace steps::mpi::tetvesicle {

class LinkSpec {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    LinkSpec(solver::LinkSpecdef* lspecdef,
             solver::linkspec_individual_id linkspec,
             Vesicle* ves,
             tetrahedron_global_id tet_gidx,
             math::position_rel_to_ves spec_pos_rel);

    LinkSpec(solver::LinkSpecdef* lspecdef, Vesicle* ves, std::fstream& cp_file);

    ~LinkSpec();

    // This is a very important function and must be called before this LinkSpec
    // becomes functional. Can't be part of constructor, though, because LinkSpecs
    // appear in LinkSpecPair consrtuctor, so it's chicken and egg
    void addLinkSpecPair(LinkSpecPair* linkspecpair);

    void removeLinkSpecPair();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore is done by 2nd constructor directly

    ////////////////////////////////////////////////////////////////////////

    inline solver::LinkSpecdef* def() const noexcept {
        return pDef;
    }

    inline solver::linkspec_global_id getGidx() const noexcept {
        return def()->gidx();
    }

    void setDef(solver::LinkSpecdef* linkspec_def) noexcept;

    inline solver::linkspec_individual_id getUniqueID() const noexcept {
        return pUniqueID;
    }

    inline math::position_rel_to_ves const& getPosCartesian_rel() const noexcept {
        return pPosCartesian_rel;
    }

    void setPosCartesian_rel(const math::position_rel_to_ves& pos);

    math::position_abs getPosCartesian_abs() const;

    inline Vesicle* getVesicle() const noexcept {
        return pVesicle;
    }

    LinkSpec* getLinkedSpec() const;

    inline tetrahedron_global_id const& getOverlapTet_gidx() const noexcept {
        return pTetOverlap;
    }

    inline void setOverlapTet_gidx(tetrahedron_global_id tet) {
        pTetOverlap = tet;
    }

    double getLength() const;

    bool movePosAllowed(const math::point3d& move_vector) const;

    void updatePos(double theta, double phi);

    bool withinBounds() const;

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::LinkSpecdef* pDef;

    solver::linkspec_individual_id pUniqueID;
    // position relative to vesicle centre
    math::position_rel_to_ves pPosCartesian_rel;

    math::position_spherical pPosSpherical;

    bool pLinked;
    LinkSpecPair* pPair;

    Vesicle* pVesicle;

    tetrahedron_global_id pTetOverlap;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::mpi::tetvesicle
