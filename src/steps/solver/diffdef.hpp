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
#include <fstream>
#include <string>

// STEPS headers.
#include "statedef.hpp"
#include "util/common.hpp"

#include "solver/fwd.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::solver {

////////////////////////////////////////////////////////////////////////////////
/// Defined diffusion object.
template <typename GlobalId>
class MetaDiffdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the object.
    /// \param d Pointer to the associated Diff object.

    // NOTE: Since one class is used for both volume and surface diffusion,
    // can't think of another way right now but to use unsigned ints and not
    // strong ids for the indices
    MetaDiffdef(Statedef* sd, GlobalId idx, model::Diff* d);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this diffusion rule.
    inline GlobalId gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of this diffusion rule.
    inline const std::string& name() const noexcept {
        return pName;
    }

    /// Return the diffusion constant.
    inline double dcst() const noexcept {
        return pDcst;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: LIGAND
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of the ligand species.
    spec_global_id lig() const;

    int dep(spec_global_id gidx) const;

    bool reqspec(spec_global_id gidx) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup();

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: DIFFUSION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Set the diffusion constant for this diffusion rule.
    ///
    /// \param d Rate constant of the diffusion rule.
    void setDcst(double d);

    /*
    /// Set the ligand of the diffusion rule by its global index
    ///
    /// \param gidx Global index of the diffusion rule.
    void setLig(uint gidx);
    */

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef* pStatedef;

    // The global index of this diffusion rule
    GlobalId pIdx;

    // The string identifier of this diffusion rule
    std::string pName;

    // The diffusion constant
    double pDcst;

    // The chemical species to which this diffusion rule applies,
    // safer to store as a sting, rather than model level spec object pointer
    std::string pLig;

    // The global index of the spec
    spec_global_id ligGIdx;

    bool pSetupdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: LIGAND
    ////////////////////////////////////////////////////////////////////////

    std::vector<int> pSpec_DEP;

    ////////////////////////////////////////////////////////////////////////
};

// explicit template instantiation declarations
extern template class MetaDiffdef<diff_global_id>;
extern template class MetaDiffdef<surfdiff_global_id>;

}  // namespace steps::solver
