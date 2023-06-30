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
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "fwd.hpp"
#include "util/common.hpp"

namespace steps::model {

// Forward declarations.
class RaftSReac;
class Raftsys;
class Model;
class Spec;

// Auxiliary declarations.
typedef RaftSReac* RaftSReacP;
typedef std::map<std::string, RaftSReacP> RaftSReacPMap;
typedef RaftSReacPMap::iterator RaftSReacPMapI;
typedef RaftSReacPMap::const_iterator RaftSReacPMapCI;
typedef std::vector<RaftSReacP> RaftSReacPVec;
typedef RaftSReacPVec::iterator RaftSReacPVecI;
typedef RaftSReacPVec::const_iterator RaftSReacPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Surface reaction.
///
/// A RaftSReac object describes a reaction which takes place on a raft system.
///
/// \warning Methods start with an underscore are not exposed to Python.
class RaftSReac {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the raft surface reaction.
    /// \param vessurfsys Pointer to the parent raft surface system.
    /// \param ilhs Volume species in the inner compartment
    ///                on the left hand side of the reaction.
    /// \param olhs Volume species in the outer compartment
    ///                on the left hand side of the reaction.
    /// \param slhs Surface species on the left hand side of the reaction.
    /// \param rslhs Species associated with the raft surface, on the left hand
    /// side of the reaction.

    /// \param srhs Surface species on the right hand side of the reaction.
    /// \param orhs Volume species in the outer compartment
    ///             and on the right hand side of the reaction.
    /// \param irhs Volume species in the inner compartment
    ///             and on the right hand side of the reaction.
    /// \param rsrhs Species associated with the raft surface, on the right hand
    /// side of the reaction.

    /// \param rsdeps Dependency on the raft surface
    /// \param anti_rsdeps Anti-Dependency on the raft surface

    /// \param kcst Rate constant of the reaction. \param
    /// \param immobilization A flag used to control the mobility effect of this
    ///			reaction. RaftSReacs can mobilise or immobilize rafts.
    ///

    RaftSReac(std::string const& id,
              Raftsys* raftsys,
              std::vector<Spec*> const& ilhs = {},
              std::vector<Spec*> const& olhs = {},
              std::vector<Spec*> const& slhs = {},
              std::vector<Spec*> const& rslhs = {},
              std::vector<Spec*> const& rsrhs = {},
              std::vector<Spec*> const& srhs = {},
              std::vector<Spec*> const& orhs = {},
              std::vector<Spec*> const& irhs = {},
              std::vector<Spec*> const& rsdeps = {},
              std::vector<Spec*> const& anti_rsdeps = {},
              double kcst = 0.0,
              Immobilization immobilization = NO_EFFECT);

    /// Destructor
    ~RaftSReac();

    ////////////////////////////////////////////////////////////////////////
    // RAFT SURFACE REACTION RULE PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the raft surface reaction rule ID.
    ///
    /// \return ID of the raft surface reaction.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the raft surface reaction rule ID.
    ///
    /// \param id ID of the raft surface reaction.
    void setID(std::string const& id);

    /// Return a pointer to the parent surface system.
    ///
    /// \return Pointer to the surface system.
    inline Raftsys* getRaftsys() const noexcept {
        return pRaftsys;
    }

    /// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent model.
    inline Model* getModel() const noexcept {
        return pModel;
    }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Check if the lhs involves species in the inner compartment.
    /// Will always be false
    /// \return True if ilhs is set.
    ///         False if else.
    inline bool getInner() const noexcept {
        return (!pOuter);
    }

    /// Check if the lhs involves species in the outer compartment,
    /// or there are no volume species on the lhs.
    /// Will always be true
    /// \return True if olhs is set, or neither olhs or ilhs are set.
    ///         False if else.
    inline bool getOuter() const noexcept {
        return pOuter;
    }

    /// Return a list of outer volume species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side outer volume species.
    inline const std::vector<Spec*>& getOLHS() const noexcept {
        return pOLHS;
    }

    /// Set the outer volume species on the left hand side of reaction.
    ///
    /// \param olhs Outer volume species on the left hand side of reaction.
    void setOLHS(std::vector<Spec*> const& olhs);

    /// Return a list of inner volume species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side inner volume species.
    inline const std::vector<Spec*>& getILHS() const noexcept {
        return pILHS;
    }

    /// Set the inner volume species on the left hand side of reaction.
    ///
    /// \param ilhs Inner volume species on the left hand side of reaction.
    void setILHS(std::vector<Spec*> const& ilhs);

    /// Return a list of raft species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side raft species.
    inline const std::vector<Spec*>& getRsLHS() const noexcept {
        return pRsLHS;
    }

    /// Set the raft  species on the left hand side of reaction.
    ///
    /// \param ilhs Rafticle volume species on the left hand side of reaction.
    void setRsLHS(std::vector<Spec*> const& rslhs);

    /// Return a list of surface species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side surface species.
    inline const std::vector<Spec*>& getSLHS() const noexcept {
        return pSLHS;
    }

    /// Set the surface species on the left hand side of reaction.
    ///
    /// \param slhs Surface species on the left hand side of reaction.
    void setSLHS(std::vector<Spec*> const& slhs);

    /// Return a list of raft species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side raft species.
    inline const std::vector<Spec*>& getRsRHS() const noexcept {
        return pRsRHS;
    }

    /// Set the raft species on the right hand side of reaction.
    ///
    /// \param irhs Rafticle species on the right hand side of reaction.
    void setRsRHS(std::vector<Spec*> const& rsrhs);

    /// Return a list of surface species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side surface species.
    inline const std::vector<Spec*>& getSRHS() const noexcept {
        return pSRHS;
    }

    /// Set the surface species on the right hand side of reaction.
    ///
    /// \param srhs Surface species on the right hand side of reaction.
    void setSRHS(std::vector<Spec*> const& srhs);

    /// Return a list of outer volume species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side outer volume species.
    inline const std::vector<Spec*>& getORHS() const noexcept {
        return pORHS;
    }

    /// Set the outer volume species on the right hand side of reaction.
    ///
    /// \param orhs Outer volume species on the right hand side of reaction.
    void setORHS(std::vector<Spec*> const& orhs);

    /// Return a list of inner volume species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side inner volume species.
    inline const std::vector<Spec*>& getIRHS() const noexcept {
        return pIRHS;
    }

    /// Set the inner volume species on the right hand side of reaction.
    ///
    /// \param orhs Inner volume species on the right hand side of reaction.
    void setIRHS(std::vector<Spec*> const& irhs);

    /// Return a list of raft species dependencies for the reaction.
    ///
    /// \return List of pointers of dependencies.
    inline const std::vector<Spec*>& getRsDeps() const noexcept {
        return pRsDeps;
    }

    /// Set the raft species dependencies for the reaction.
    ///
    /// \param rsdeps Raft species dependencies for the reaction.
    void setRsDeps(std::vector<Spec*> const& rsdeps);

    /// Return a list of raft species anti-dependencies for the reaction.
    ///
    /// \return List of pointers of anti-dependencies.
    inline const std::vector<Spec*>& getAntiRsDeps() const noexcept {
        return pAntiRsDeps;
    }

    /// Set the raft species anti-dependencies for the reaction.
    ///
    /// \param anti_rsdeps Raft species anti-dependencies for the reaction.
    void setAntiRsDeps(std::vector<Spec*> const& anti_rsdeps);

    /// Get the order of the raft surface reaction.
    ///
    /// \return Order of the reaction.
    inline uint getOrder() const noexcept {
        return pOrder;
    }

    /// Get the rate constant of the raft surface reaction.
    ///
    /// \return Rate constant of the raft surface reaction.
    inline double getKcst() const noexcept {
        return pKcst;
    }

    /// Set the rate constant of the raft surface reaction.
    ///
    /// \param Rate constant of the raft surface reaction.
    void setKcst(double kcst);

    /// Get the mobility effect of the vesicle surface reaction.
    ///
    /// \return Immobilization of the vesicle surface reaction.
    inline Immobilization getImmobilization() const noexcept {
        return pImmobilization;
    }

    /// Get a list of all species.
    ///
    /// Returns a list of all species involved in this
    /// raft surface reaction, on both the left and righthand side
    /// and does not contain any duplicate members.
    ///
    /// \return List of pointers to the species.
    std::vector<Spec*> getAllSpecs() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////
    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once
    void _handleSelfDelete();

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model* pModel;
    Raftsys* pRaftsys;

    std::vector<Spec*> pILHS;
    std::vector<Spec*> pOLHS;
    std::vector<Spec*> pSLHS;
    std::vector<Spec*> pRsLHS;

    std::vector<Spec*> pRsRHS;
    std::vector<Spec*> pSRHS;
    std::vector<Spec*> pORHS;
    std::vector<Spec*> pIRHS;

    std::vector<Spec*> pRsDeps;
    std::vector<Spec*> pAntiRsDeps;

    Immobilization pImmobilization;

    bool pOuter;
    uint pOrder;
    double pKcst;
};

}  // namespace steps::model
