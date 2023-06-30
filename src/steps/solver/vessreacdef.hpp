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
#include "model/linkspec.hpp"
#include "model/spec.hpp"
#include "model/vessreac.hpp"
#include "solver/api.hpp"
#include "solver/fwd.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"

namespace steps::solver {

/// Defined Vesicle Surface Reaction.
class VesSReacdef {
  public:
    enum orientT  ///< Orientation of the reaction.
    {
        INSIDE = 0,
        OUTSIDE = 1
    };

    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the surface reaction.
    /// \param sr Pointer to the VesSReac object.
    VesSReacdef(Statedef* sd, vessreac_global_id idx, model::VesSReac* sr);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////

    // Need reset function to return to original kcsts
    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLE SURFACE REACTION RULE
    ////////////////////////////////////////////////////////////////////////

    // Added for convenience for solver object access
    inline uint countSpecsGlobal() const noexcept {
        return pStatedef->countSpecs();
    }

    inline uint countLinkSpecsGlobal() const noexcept {
        return pStatedef->countLinkSpecs();
    }

    /// Return the global index of this vesicle surface reaction rule.
    inline vessreac_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the vesicle surface reaction.
    inline std::string const& name() const noexcept {
        return pName;
    }

    /// Return the order of this vesicle surface reaction.
    inline uint order(void) const noexcept {
        return pOrder;
    }

    /// Return the MACROscopic reaction constant.
    inline double kcst() const noexcept {
        return pKcst;
    }

    void setKcst(double k);

    inline int immobility() const noexcept {
        return pImmobility;
    }

    inline double max_distance() const noexcept {
        return pMaxDistance;
    }

    inline uint getExtent() const noexcept {
        return pExtent;
    }

    inline void incExtent() noexcept {
        pExtent++;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Returns true if the left hand side of the reaction stoichiometry
    /// involves reactants on the surface and on the inside volume.
    // inline bool inside() const noexcept
    //{ return (pOrient == INSIDE); }

    /// Returns true if the left hand side of the reaction stoichiometry
    /// involves reactants on the surface and on the outside volume. This
    /// method is mutually exclusive with SReacDef::inside, but not
    /// with SReacDef::insideRef.
    ///
    inline bool outside() const noexcept {
        return (pOrient == OUTSIDE);
    }

    /// Returns true if any aspect of the surface reaction references
    /// species on the outside volume, regardless of how they are
    /// referenced. Whereas method VesSReacDef::outside only makes a
    /// statement about the LHS part of the reaction stoichiometry, this
    /// method checks everything, including the right hand side.
    ///
    /// As such, this method will always return true if
    /// VesSReacDef::outside is true. The converse, however, is not the
    /// the case: VesSReacDef::outside does not necessarily return true
    /// if this routine returns true.
    ///
    /// It basically polls VesSReacDef::req_O for each possible species.
    ///
    bool reqOutside() const;

    /// Returns true if any aspect of the surface reaction references
    /// species on the inside volume, regardless of how they are
    /// referenced. Whereas method VesSReacDef::outside only makes a
    /// statement about the LHS part of the reaction stoichiometry, this
    /// method checks everything, including the right hand side.
    ///
    /// As such, this method will always return true if
    /// VesSReacDef::outside is true. The converse, however, is not the
    /// the case: VesSReacDef::outside does not necessarily return true
    /// if this routine returns true.
    ///
    /// It basically polls VesSReacDef::req_I for each possible species.
    ///
    bool reqInside() const;

    /// True if we need surface
    ///
    bool reqSurface() const;

    /// Return true if this reaction only involves surface species,
    /// nothing in a volume at all. In that case the reaction constant
    /// should be treated in 2D
    inline bool surf_surf(void) const noexcept {
        return pSurface_surface;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Returns the number of molecules of species idx required in
    /// the inner volume (_I) , outer volume (_O) , surface patch (_S) or vesicle
    /// (_V) to have one occurrence of this surface reaction.
    ///
    uint lhs_S(spec_global_id gidx) const;
    uint lhs_O(spec_global_id gidx) const;
    uint lhs_V(spec_global_id gidx) const;

    uint lhs_L(linkspec_global_id gidx) const;

    uint vdep(spec_global_id gidx) const;

    /// Returns a description of how an occurrence of this surface reaction
    /// depends on some species, defined by its global index idx, to occur.
    /// See steps/sim/shared/types.hpp for more information on the return
    /// type. This method is distinct from the
    /// SReacDef::req_S and SReacDef::req_O methods.
    ///
    depT dep_S(spec_global_id gidx);
    depT dep_O(spec_global_id gidx) const;
    depT dep_V(spec_global_id gidx) const;
    depT dep_L(linkspec_global_id gidx) const;

    /// Returns how many molecules of some species, specified by its
    /// global index, are produced after a single occurrence of this
    /// surface reaction.
    /// '_S' for the surface patch, '_O' for the outer volume,
    /// '_I' for the inner volume and '_V' for vesicle.
    ///
    uint rhs_S(spec_global_id gidx) const;
    uint rhs_O(spec_global_id gidx) const;
    uint rhs_V(spec_global_id gidx) const;
    uint rhs_I(spec_global_id gidx) const;
    uint rhs_L(linkspec_global_id gidx) const;

    /// Returns how the amount of a species, specified by its global index,
    /// changes as the result of a single occurence of this surface
    /// reaction on the outer volume (_O), vesicle (_V) , inner volume (_I)
    /// surface patch (_S).
    ///
    int upd_S(spec_global_id gidx) const;
    int upd_O(spec_global_id gidx) const;
    int upd_V(spec_global_id gidx) const;
    uint upd_I(spec_global_id gidx) const;  // I upd always positive
    int upd_L(linkspec_global_id gidx) const;

    /// Returns whether the surface reaction rule references a species,
    /// specified by its global index, on the
    /// outer volume (_O), surface patch (_S), vesicle (_V) or inner volume (_I)
    ///
    bool reqspec_S(spec_global_id gidx) const;
    bool reqspec_O(spec_global_id gidx) const;
    bool reqspec_V(spec_global_id gidx) const;
    bool reqspec_I(spec_global_id gidx) const;
    bool reqspec_L(linkspec_global_id gidx) const;

    inline const spec_global_id_vec& updColl_S() const noexcept {
        return pSpec_S_UPD_Coll;
    }

    inline const spec_global_id_vec& updColl_O() const noexcept {
        return pSpec_O_UPD_Coll;
    }

    inline const spec_global_id_vec& updColl_V() const noexcept {
        return pSpec_V_UPD_Coll;
    }

    inline const spec_global_id_vec& updColl_I() const noexcept {
        return pSpec_I_UPD_Coll;
    }

    inline const linkspec_global_id_vec& updColl_L() const noexcept {
        return pSpec_L_UPD_Coll;
    }

    inline Statedef* statedef() const noexcept {
        return pStatedef;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef* pStatedef;
    vessreac_global_id pIdx;

    uint pExtent;

    // Have to store this to access the original kcsts to allow for reset
    model::VesSReac* pVessreac;

    std::string pName;
    uint pOrder;
    double pKcst;
    int pImmobility;
    double pMaxDistance;

    // The stoichiometry stored as model level Spec objects.
    // To be used during setup ONLY
    model::SpecPVec pOlhs;
    model::SpecPVec pSlhs;
    model::SpecPVec pVlhs;
    model::LinkSpecPVec pLlhs;

    model::SpecPVec pOrhs;
    model::SpecPVec pSrhs;
    model::SpecPVec pVrhs;
    model::SpecPVec pIrhs;
    model::LinkSpecPVec pLrhs;

    model::SpecPVec pVdep;

    bool pSetupdone;

    // Store whether this surface reaction is 2D or not
    bool pSurface_surface;

    /// Does the left-hand side of the stoichiometry involve molecules
    /// on the inside or on the outside?
    orientT pOrient;

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Array vector describing dependencies for
    /// outer volume (_O_), inner volume (_I_), surface (_S_)
    /// and vesicle (_V_) species. Dependencies can
    /// be stoichiometric or in the rate function (this is not implemented
    /// yet) -- see 'steps/sim/shared/types.hpp'. The vector must be
    /// indexed through global indices, i.e. it runs over all species in
    /// the entire model.
    ///

    // Test map instead of sparse array
    // depT                              * pSpec_S_DEP;
    std::map<spec_global_id, depT> pSpec_S_DEP;

    util::strongid_vector<spec_global_id, depT> pSpec_O_DEP;
    util::strongid_vector<spec_global_id, depT> pSpec_V_DEP;
    util::strongid_vector<linkspec_global_id, depT> pSpec_L_DEP;

    /// Vector describing the left hand (reactant) side of the reaction
    /// stoichiometry, for species in the outer
    /// volume (_O_), inner volume (_I_), surface (_S_) or vesicle (_V_).
    /// The vector must be indexed
    /// through global indices, i.e. it runs over all species in the entire
    /// model.
    ///
    util::strongid_vector<spec_global_id, uint> pSpec_S_LHS;
    util::strongid_vector<spec_global_id, uint> pSpec_O_LHS;
    util::strongid_vector<spec_global_id, uint> pSpec_V_LHS;
    util::strongid_vector<linkspec_global_id, uint> pSpec_L_LHS;

    // Information about the vesicle surface species dependency
    util::strongid_vector<spec_global_id, uint> pSpec_VDEP;

    /// An array vector describing the right hand (reaction product) side
    /// of the surface reaction stoichiometry, for species in the
    /// outer volume (_O_), inner volume (_I_), surface (_S_) or vesicle (_V_).
    /// The vector must be indexed through global indices, i.e. it runs over all
    /// species in the entire model.
    ///
    util::strongid_vector<spec_global_id, uint> pSpec_S_RHS;
    util::strongid_vector<spec_global_id, uint> pSpec_O_RHS;
    util::strongid_vector<spec_global_id, uint> pSpec_V_RHS;
    util::strongid_vector<spec_global_id, uint> pSpec_I_RHS;
    util::strongid_vector<linkspec_global_id, uint> pSpec_L_RHS;

    /// An array describing the update vector (i.e. RHS[] - LHS[]) of
    /// the surface reaction, for species in the
    /// outer volume (_O_), inner volume (_I_), patch surface (_S_)
    /// and vesicle (_V_). The vector must be
    /// indexed through global indices, i.e. it runs over all species in
    /// the entire model.
    ///
    util::strongid_vector<spec_global_id, int> pSpec_S_UPD;
    util::strongid_vector<spec_global_id, int> pSpec_O_UPD;
    util::strongid_vector<spec_global_id, int> pSpec_V_UPD;
    util::strongid_vector<spec_global_id, uint> pSpec_I_UPD;
    util::strongid_vector<linkspec_global_id, int> pSpec_L_UPD;

    /// A vector collecting the global indices of all species that are
    /// updated when this surface reaction rule occurs.
    spec_global_id_vec pSpec_S_UPD_Coll;
    spec_global_id_vec pSpec_O_UPD_Coll;
    spec_global_id_vec pSpec_V_UPD_Coll;
    spec_global_id_vec pSpec_I_UPD_Coll;
    linkspec_global_id_vec pSpec_L_UPD_Coll;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::solver
