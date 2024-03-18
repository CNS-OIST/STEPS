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

#include "api.hpp"
#include "fwd.hpp"
#include "model/complexsreac.hpp"
#include "solver/complexeventsdef.hpp"
#include "solver/sreacdef.hpp"

namespace steps::solver {

// Forwards declarations
class Statedef;

/// Defined Surface Reaction.
/// \todo imcompleted.
class ComplexSReacdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the surface reaction.
    /// \param sr Pointer to the ComplexSReac object.
    ComplexSReacdef(Statedef& sd, complexsreac_global_id idx, model::ComplexSReac& sr);

    /// Destructor
    ~ComplexSReacdef() = default;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE REACTION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this surface reaction rule.
    complexsreac_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the surface reaction.
    const std::string& name() const {
        return pName;
    }

    /// Return the order of this surface reaction.
    uint order() const {
        return pOrder;
    }

    /// Return the MACROscopic reaction constant.
    double kcst() const {
        return pKcst;
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
    bool inside() const noexcept {
        return pOrient == SReacdef::INSIDE;
    }

    /// Returns true if any aspect of the surface reaction references
    /// species on the inside volume, regardless of how they are
    /// referenced. Whereas method ComplexSReacDef::inside only makes a
    /// statement about the LHS part of the reaction stoichiometry, this
    /// method checks everything, including the right hand side.
    ///
    /// As such, this method will always return true if
    /// ComplexSReacDef::inside is true. The converse, however, is not the
    /// the case: ComplexSReacDef::inside does not necessarily return true
    /// if this routine returns true.
    ///
    /// It basically polls ComplexSReacDef::req_I for each possible species.

    ///
    bool reqInside() const;

    /// Returns true if the left hand side of the reaction stoichiometry
    /// involves reactants on the surface and on the outside volume. This
    /// method is mutually exclusive with ComplexSReacDef::inside
    ///
    bool outside() const noexcept {
        return pOrient == SReacdef::OUTSIDE;
    }

    /// Returns true if any aspect of the surface reaction references
    /// species on the outside volume, regardless of how they are
    /// referenced. Whereas method ComplexSReacDef::outside only makes a
    /// statement about the LHS part of the reaction stoichiometry, this
    /// method checks everything, including the right hand side.
    ///
    /// As such, this method will always return true if
    /// ComplexSReacDef::outside is true. The converse, however, is not the
    /// the case: ComplexSReacDef::outside does not necessarily return true
    /// if this routine returns true.
    ///
    /// It basically polls ComplexSReacDef::req_O for each possible species.
    ///
    bool reqOutside() const;

    /// Return true if this reaction only involves surface species,
    /// nothing in a volume at all. In that case the reaction constant
    /// should be treated in 2D
    bool surf_surf() const noexcept {
        return pSurface_surface;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Returns the number of molecules of species idx required in
    /// the inner volume (_I), outer volume (_O) or surface patch (_S)
    /// to have one occurence of this surface reaction.
    ///
    uint lhs_I(spec_global_id gidx) const;
    uint lhs_S(spec_global_id gidx) const;
    uint lhs_O(spec_global_id gidx) const;

    /// Returns a description of how an occurence of this surface reaction
    /// depends on some species, defined by its global index idx, to occur.
    /// See steps/sim/shared/types.hpp for more information on the return
    /// type. This method is distinct from the ComplexSReacDef::req_I,
    /// ComplexSReacDef::req_S and ComplexSReacDef::req_O methods.
    ///
    depT dep_I(spec_global_id gidx) const;
    depT dep_S(spec_global_id gidx) const;
    depT dep_O(spec_global_id gidx) const;
    bool complexdep(model::ComplexLocation loc,
                    complex_global_id gidx,
                    complex_substate_id sus) const;

    /// Returns how many molecules of some species, specified by its
    /// global index, are produced after a single occurence of this
    /// surface reaction. '_I' returns this number for the inner volume,
    /// '_S' for the surface patch and '_O' for the outer volume.
    ///
    uint rhs_I(spec_global_id gidx) const;
    uint rhs_S(spec_global_id gidx) const;
    uint rhs_O(spec_global_id gidx) const;

    /// Returns how the amount of a species, specified by its global index,
    /// changes as the result of a single occurence of this surface
    /// reaction on the inside volume (_I), outer volume (_O) or
    /// surface patch (_S).
    ///
    int upd_I(spec_global_id gidx) const;
    int upd_S(spec_global_id gidx) const;
    int upd_O(spec_global_id gidx) const;

    /// Returns whether the surface reaction rule references a species,
    /// specified by its global index, on the inner volume side (_I),
    /// outer volume (_O) or surface patch (_S).
    ///
    bool reqspec_I(spec_global_id gidx) const;
    bool reqspec_S(spec_global_id gidx) const;
    bool reqspec_O(spec_global_id gidx) const;

    const spec_global_id_vec& updColl_I() const noexcept {
        return pSpec_I_UPD_Coll;
    }
    const spec_global_id_vec& updColl_S() const noexcept {
        return pSpec_S_UPD_Coll;
    }
    const spec_global_id_vec& updColl_O() const noexcept {
        return pSpec_O_UPD_Coll;
    }

    const std::map<complex_global_id, std::set<complex_substate_id>>& complexUPDMAP(
        model::ComplexLocation loc) const;

    const std::vector<std::shared_ptr<ComplexUpdateEventdef>>& updEvents(
        model::ComplexLocation loc) const;
    const std::vector<std::shared_ptr<ComplexDeleteEventdef>>& delEvents(
        model::ComplexLocation loc) const;
    const std::vector<std::shared_ptr<ComplexCreateEventdef>>& creEvents(
        model::ComplexLocation loc) const;

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef& pStatedef;
    const complexsreac_global_id pIdx;

    const std::string pName;
    const uint pOrder;
    const double pKcst;

    // The stoichiometry stored as model level Spec objects.
    // To be used during setup ONLY
    const std::vector<model::Spec*> pIlhs;
    const std::vector<model::Spec*> pOlhs;
    const std::vector<model::Spec*> pSlhs;

    const std::vector<model::Spec*> pIrhs;
    const std::vector<model::Spec*> pOrhs;
    const std::vector<model::Spec*> pSrhs;

    bool pSetupdone{false};

    // Store whether this surface reaction is 2D or not
    bool pSurface_surface;

    /// Does the left-hand side of the stoichiometry involve molecules
    /// on the inside or on the outside?
    SReacdef::orientT pOrient;

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Array vector describing dependencies for inner volume (_I_),
    /// outer volume (_O_) and surface (_S_) species. Dependencies can
    /// be stoichiometric or in the rate function (this is not implemented
    /// yet) -- see 'steps/sim/shared/types.hpp'. The vector must be
    /// indexed through global indices, i.e. it runs over all species in
    /// the entire model.
    ///
    util::strongid_vector<spec_global_id, depT> pSpec_I_DEP;
    util::strongid_vector<spec_global_id, depT> pSpec_S_DEP;
    util::strongid_vector<spec_global_id, depT> pSpec_O_DEP;

    /// Vector describing the left hand (reactant) side of the reaction
    /// stoichiometry, for species in the inner volume (_I_), outer
    /// volume (_O_) and surface (_S_) species. The vector must be indexed
    /// through global indices, i.e. it runs over all species in the entire
    /// model.
    ///
    util::strongid_vector<spec_global_id, uint> pSpec_I_LHS;
    util::strongid_vector<spec_global_id, uint> pSpec_S_LHS;
    util::strongid_vector<spec_global_id, uint> pSpec_O_LHS;

    /// An array vector describing the right hand (reaction product) side
    /// of the surface reaction stoichiometry, for species in the inner
    /// volume (_I_), outer volume (_O_) and surface (_S_) species. The
    /// vector must be indexed through global indices, i.e. it runs over
    /// all species in the entire model.
    ///
    util::strongid_vector<spec_global_id, uint> pSpec_I_RHS;
    util::strongid_vector<spec_global_id, uint> pSpec_S_RHS;
    util::strongid_vector<spec_global_id, uint> pSpec_O_RHS;

    /// An array describing the update vector (i.e. RHS[] - LHS[]) of
    /// the surface reaction, for species in the inner volume (_I),
    /// outer volume (_O_) and patch surface (_S_). The vector must be
    /// indexed through global indices, i.e. it runs over all species in
    /// the entire model.
    ///
    util::strongid_vector<spec_global_id, int> pSpec_I_UPD;
    util::strongid_vector<spec_global_id, int> pSpec_S_UPD;
    util::strongid_vector<spec_global_id, int> pSpec_O_UPD;

    /// A vector collecting the global indices of all species that are
    /// updated when this surface reaction rule occurs.
    spec_global_id_vec pSpec_I_UPD_Coll;
    spec_global_id_vec pSpec_S_UPD_Coll;
    spec_global_id_vec pSpec_O_UPD_Coll;

    std::map<model::ComplexLocation, std::vector<std::shared_ptr<ComplexUpdateEventdef>>>
        pComplexUPDEvs;
    std::map<model::ComplexLocation, std::vector<std::shared_ptr<ComplexDeleteEventdef>>>
        pComplexDELEvs;
    std::map<model::ComplexLocation, std::vector<std::shared_ptr<ComplexCreateEventdef>>>
        pComplexCREEvs;

    // location -> {cmplxIdx -> {sub unit states ind}}
    std::map<model::ComplexLocation, std::map<complex_global_id, std::set<complex_substate_id>>>
        pComplex_DEPMAP;
    std::map<model::ComplexLocation, std::map<complex_global_id, std::set<complex_substate_id>>>
        pComplex_UPDMAP;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::solver
