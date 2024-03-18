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

#include <iosfwd>
#include <string>

#include "fwd.hpp"
#include "model/fwd.hpp"
#include "model/reac.hpp"

namespace steps::solver {

/// Defined Reaction.
class Reacdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the reaction.
    /// \param r Reference to the associated Reac object.
    Reacdef(Statedef& sd, reac_global_id idx, model::Reac& r);

    Reacdef(const Reacdef&) = delete;
    Reacdef& operator=(const Reacdef&) = delete;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: REACTION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this reaction rule.
    inline reac_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the reaction.
    std::string const name() const noexcept {
        return pName;
    }

    /// Return the order of this reaction.
    uint order() const noexcept {
        return pOrder;
    }

    /// Return the MACROscopic reaction constant.
    double kcst() const noexcept {
        return pKcst;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////
    /// \todo imcompleted.
    uint lhs(spec_global_id gidx) const;
    int dep(spec_global_id gidx) const;
    uint rhs(spec_global_id gidx) const;
    int upd(spec_global_id gidx) const;
    bool reqspec(spec_global_id gidx) const;

    inline solver::spec_global_id_vecCI bgnUpdColl() const noexcept {
        return pSpec_UPD_Coll.begin();
    }
    inline solver::spec_global_id_vecCI endUpdColl() const noexcept {
        return pSpec_UPD_Coll.end();
    }
    inline const spec_global_id_vec& updColl() const noexcept {
        return pSpec_UPD_Coll;
    }

    inline solver::spec_global_id_vec& updColl() noexcept {
        return pSpec_UPD_Coll;
    }

    inline const solver::spec_global_id_vec& UPD_Coll() const noexcept {
        return pSpec_UPD_Coll;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
    void setup(const Statedef& sd);

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    const reac_global_id pIdx;
    const std::string pName;
    const uint pOrder;
    const double pKcst;

    /// The stoichiometry stored as model level Spec objects.
    /// To be used during setup ONLY
    /// \{
    const std::vector<model::Spec*> pLhs;
    const std::vector<model::Spec*> pRhs;
    /// \}

    bool pSetupdone{false};

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    util::strongid_vector<spec_global_id, depT> pSpec_DEP;
    util::strongid_vector<spec_global_id, uint> pSpec_LHS;
    util::strongid_vector<spec_global_id, uint> pSpec_RHS;
    util::strongid_vector<spec_global_id, int> pSpec_UPD;

    spec_global_id_vec pSpec_UPD_Coll;
};

}  // namespace steps::solver
