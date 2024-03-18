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
#include "model/raftsreac.hpp"
#include "solver/fwd.hpp"

namespace steps::solver {

/// Defined Raft Surface Reaction.
class RaftSReacdef {
  public:
    enum orientT  ///< Orientation of the reaction.
    {
        INSIDE = 0,
        OUTSIDE = 1
    };

    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the raft surface reaction.
    /// \param rsr Reference to the RaftSReac object.
    RaftSReacdef(Statedef& sd, raftsreac_global_id idx, model::RaftSReac& rsr);

    RaftSReacdef(const RaftSReacdef&) = delete;
    RaftSReacdef& operator=(const RaftSReacdef&) = delete;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLE SURFACE REACTION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this raft surface reaction rule.
    inline raftsreac_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the raft surface reaction.
    inline std::string const& name() const noexcept {
        return pName;
    }

    /// Return the order of this raft surface reaction.
    inline uint order(void) const noexcept {
        return pOrder;
    }

    /// Return the MACROscopic reaction constant.
    inline double kcst() const noexcept {
        return pKcst;
    }

    inline int immobility() const noexcept {
        return pImmobility;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup(const Statedef& sd);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Returns true if the left hand side of the reaction stoichiometry
    /// involves reactants on the surface and on the inside volume.
    inline bool inside() const noexcept {
        return pOrient == INSIDE;
    }

    /// Returns true if the left hand side of the reaction stoichiometry
    /// involves reactants on the surface and on the outside volume. This
    /// method is mutually exclusive with RaftSReacDef::inside
    ///
    inline bool outside() const noexcept {
        return pOrient == OUTSIDE;
    }

    /// Returns true if any aspect of the surface reaction references
    /// species on the outside volume, regardless of how they are
    /// referenced. Whereas method RaftSReacDef::outside only makes a
    /// statement about the LHS part of the reaction stoichiometry, this
    /// method checks everything, including the right hand side.
    ///
    /// As such, this method will always return true if
    /// RaftSReacDef::outside is true. The converse, however, is not the
    /// the case: RaftSReacDef::outside does not necessarily return true
    /// if this routine returns true.
    ///
    /// It basically polls RaftSReacDef::req_O for each possible species.
    ///
    bool reqOutside() const;

    /// Returns true if any aspect of the surface reaction references
    /// species on the inside volume, regardless of how they are
    /// referenced. Whereas method RaftSReacDef::outside only makes a
    /// statement about the LHS part of the reaction stoichiometry, this
    /// method checks everything, including the right hand side.
    ///
    /// As such, this method will always return true if
    /// RaftSReacDef::outside is true. The converse, however, is not the
    /// the case: RaftSReacDef::outside does not necessarily return true
    /// if this routine returns true.
    ///
    /// It basically polls RaftSReacDef::req_I for each possible species.
    ///
    bool reqInside() const;

    /// Return true if this reaction only involves surface species,
    /// nothing in a volume at all. In that case the reaction constant
    /// should be treated in 2D
    inline bool surf_surf() const noexcept {
        return pSurface_surface;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Returns the number of molecules of species idx required in
    /// the inner volume (_I) , outer volume (_O) , surface patch (_S) or raft
    /// surface (_Rs) to have one occurrence of this surface reaction.
    ///
    uint lhs_S(spec_global_id gidx) const;
    uint lhs_O(spec_global_id gidx) const;
    uint lhs_Rs(spec_global_id gidx) const;
    uint lhs_I(spec_global_id gidx) const;

    uint rsdep(spec_global_id gidx) const;
    uint anti_rsdep(spec_global_id gidx) const;

    /// Returns a description of how an occurrence of this surface reaction
    /// depends on some species, defined by its global index idx, to occur.
    /// See steps/sim/shared/types.hpp for more information on the return
    /// type. This method is distinct from the
    /// SReacDef::req_S and SReacDef::req_O methods.
    ///
    depT dep_S(spec_global_id gidx) const;
    depT dep_O(spec_global_id gidx) const;
    depT dep_Rs(spec_global_id gidx) const;
    depT dep_I(spec_global_id gidx) const;

    /// Returns how many molecules of some species, specified by its
    /// global index, are produced after a single occurrence of this
    /// surface reaction.
    /// '_S' for the surface patch, '_O' for the outer volume,
    /// '_I' for the inner volume and '_Rs' for raft surface.
    ///
    uint rhs_S(spec_global_id gidx) const;
    uint rhs_O(spec_global_id gidx) const;
    uint rhs_Rs(spec_global_id gidx) const;
    uint rhs_I(spec_global_id gidx) const;

    /// Returns how the amount of a species, specified by its global index,
    /// changes as the result of a single occurence of this surface
    /// reaction on the outer volume (_O), raft (_Rs) , inner volume (_I)
    /// surface patch (_S).
    ///
    int upd_S(spec_global_id gidx) const;
    int upd_O(spec_global_id gidx) const;
    int upd_Rs(spec_global_id gidx) const;
    int upd_I(spec_global_id gidx) const;

    /// Returns whether the surface reaction rule references a species,
    /// specified by its global index, on the
    /// outer volume (_O), surface patch (_S), raft (_Rs) or inner volume (_I)
    ///
    bool reqspec_S(spec_global_id gidx) const;
    bool reqspec_O(spec_global_id gidx) const;
    bool reqspec_Rs(spec_global_id gidx) const;
    bool reqspec_I(spec_global_id gidx) const;

    inline const spec_global_id_vec& updColl_S() const noexcept {
        return pSpec_S_UPD_Coll;
    }
    inline const spec_global_id_vec& updColl_O() const noexcept {
        return pSpec_O_UPD_Coll;
    }
    inline const spec_global_id_vec& updColl_Rs() const noexcept {
        return pSpec_Rs_UPD_Coll;
    }
    inline const spec_global_id_vec& updColl_I() const noexcept {
        return pSpec_I_UPD_Coll;
    }

  private:
    inline auto num_specs() const noexcept {
        return pCountSpecs;
    }

    const raftsreac_global_id pIdx;

    const std::string pName;
    const uint pOrder;
    const double pKcst;
    const uint pCountSpecs;

    const int pImmobility;

    /// Does the left-hand side of the stoichiometry involve molecules
    /// on the inside or on the outside?
    const orientT pOrient;

    /// The stoichiometry stored as model level Spec objects.
    /// To be used during setup ONLY
    /// \{
    const std::vector<model::Spec*> pOlhs;
    const std::vector<model::Spec*> pSlhs;
    const std::vector<model::Spec*> pRslhs;
    const std::vector<model::Spec*> pIlhs;

    const std::vector<model::Spec*> pOrhs;
    const std::vector<model::Spec*> pSrhs;
    const std::vector<model::Spec*> pRsrhs;
    const std::vector<model::Spec*> pIrhs;
    /// \}

    std::vector<model::Spec*> pRsdep;
    std::vector<model::Spec*> pAntiRsdep;

    bool pSetupdone{false};

    // Store whether this surface reaction is 2D or not
    bool pSurface_surface{true};

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Array vector describing dependencies for
    /// outer volume (_O_), inner volume (_I_), surface (_S_)
    /// and raft (_Rs_) species. Dependencies can
    /// be stoichiometric or in the rate function (this is not implemented
    /// yet) -- see 'steps/sim/shared/types.hpp'. The vector must be
    /// indexed through global indices, i.e. it runs over all species in
    /// the entire model.
    ///
    util::strongid_vector<spec_global_id, depT> pSpec_S_DEP;
    util::strongid_vector<spec_global_id, depT> pSpec_O_DEP;
    util::strongid_vector<spec_global_id, depT> pSpec_Rs_DEP;
    util::strongid_vector<spec_global_id, depT> pSpec_I_DEP;

    /// Vector describing the left hand (reactant) side of the reaction
    /// stoichiometry, for species in the outer
    /// volume (_O_), inner volume (_I_), surface (_S_) or raft (_Rs_).
    /// The vector must be indexed
    /// through global indices, i.e. it runs over all species in the entire
    /// model.
    ///
    util::strongid_vector<spec_global_id, uint> pSpec_S_LHS;
    util::strongid_vector<spec_global_id, uint> pSpec_O_LHS;
    util::strongid_vector<spec_global_id, uint> pSpec_Rs_LHS;
    util::strongid_vector<spec_global_id, uint> pSpec_I_LHS;

    // Information about the raft surface species dependency
    util::strongid_vector<spec_global_id, uint> pSpec_RsDEP;
    // Information about the raft surface species dependency
    util::strongid_vector<spec_global_id, uint> pSpec_AntiRsDEP;

    /// An array vector describing the right hand (reaction product) side
    /// of the surface reaction stoichiometry, for species in the
    /// outer volume (_O_), inner volume (_I_), surface (_S_) or raft (_Rs_). The
    /// vector must be indexed through global indices, i.e. it runs over
    /// all species in the entire model.
    ///
    util::strongid_vector<spec_global_id, uint> pSpec_S_RHS;
    util::strongid_vector<spec_global_id, uint> pSpec_O_RHS;
    util::strongid_vector<spec_global_id, uint> pSpec_Rs_RHS;
    util::strongid_vector<spec_global_id, uint> pSpec_I_RHS;

    /// An array describing the update vector (i.e. RHS[] - LHS[]) of
    /// the surface reaction, for species in the
    /// outer volume (_O_), inner volume (_I_), patch surface (_S_)
    /// and raft (_Rs_). The vector must be
    /// indexed through global indices, i.e. it runs over all species in
    /// the entire model.
    ///
    util::strongid_vector<spec_global_id, int> pSpec_S_UPD;
    util::strongid_vector<spec_global_id, int> pSpec_O_UPD;
    util::strongid_vector<spec_global_id, int> pSpec_Rs_UPD;
    util::strongid_vector<spec_global_id, int> pSpec_I_UPD;

    /// A vector collecting the global indices of all species that are
    /// updated when this surface reaction rule occurs.
    spec_global_id_vec pSpec_S_UPD_Coll;
    spec_global_id_vec pSpec_O_UPD_Coll;
    spec_global_id_vec pSpec_Rs_UPD_Coll;
    spec_global_id_vec pSpec_I_UPD_Coll;
};

}  // namespace steps::solver
