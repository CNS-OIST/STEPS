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
class VesSReac;
class VesSurfsys;
class Model;
class Spec;
class LinkSpec;

// Auxiliary declarations.
typedef VesSReac* VesSReacP;
typedef std::map<std::string, VesSReacP> VesSReacPMap;
typedef VesSReacPMap::iterator VesSReacPMapI;
typedef VesSReacPMap::const_iterator VesSReacPMapCI;
typedef std::vector<VesSReacP> VesSReacPVec;
typedef VesSReacPVec::iterator VesSReacPVecI;
typedef VesSReacPVec::const_iterator VesSReacPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Surface reaction.
///
/// A VesSReac object describes a reaction which takes place on a vesicle
/// surface system,
///
/// \warning Methods start with an underscore are not exposed to Python.
class VesSReac {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the vesicle surface reaction.
    /// \param vessurfsys Pointer to the parent vesicle surface system.
    /// \param ilhs Volume species in the inner compartment
    ///                on the left hand side of the reaction.
    /// \param olhs Volume species in the outer compartment
    ///                on the left hand side of the reaction.
    /// \param slhs Surface species on the left hand side of the reaction.
    /// \param vlhs Species associated with the vesicle surface, on the left hand
    /// side of the reaction. \param llhs Link species associated with the vesicle
    /// surface, on the left hand side of the reaction.

    /// \param lrhs Link species associated with the vesicle surface, on the right
    /// hand side of the reaction. \param vrhs Species associated with the vesicle
    /// surface, on the right hand side of the reaction.
    /// \param srhs Surface species on the right hand side of the reaction.
    /// \param orhs Volume species in the outer compartment
    ///             and on the right hand side of the reaction.
    /// \param irhs Volume species in the inner compartment
    ///             and on the right hand side of the reaction.
    /// \param vdeps Species dependencies associated with the vesicle surface.
    /// \param kcst Rate constant of the reaction.
    ///
    /// \param immobilization A flag to describe the effect of the vesicle surface
    /// reaction on mobility- can mobilise or immobilize the vesicle
    /// \param max_distance An optional maximum distance at
    /// which the reaction can occur when involving
    ///             reactants on the vesicle surface and a 'patch' surface.

    VesSReac(std::string const& id,
             VesSurfsys* vessurfsys,
             std::vector<Spec*> const& olhs = {},
             std::vector<Spec*> const& slhs = {},
             std::vector<Spec*> const& vlhs = {},
             std::vector<LinkSpec*> const& llhs = {},
             std::vector<LinkSpec*> const& lrhs = {},
             std::vector<Spec*> const& vrhs = {},
             std::vector<Spec*> const& srhs = {},
             std::vector<Spec*> const& orhs = {},
             std::vector<Spec*> const& irhs = {},
             std::vector<Spec*> const& vdeps = {},
             double kcst = 0.0,
             Immobilization immobilization = NO_EFFECT,
             double max_distance = -1.0);

    /// Destructor
    ~VesSReac();

    ////////////////////////////////////////////////////////////////////////
    // VESICLE SURFACE REACTION RULE PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the vesicle surface reaction rule ID.
    ///
    /// \return ID of the vesicle surface reaction.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the vesicle surface reaction rule ID.
    ///
    /// \param id ID of the vesicle surface reaction.
    void setID(std::string const& id);

    /// Return a pointer to the parent surface system.
    ///
    /// \return Pointer to the surface system.
    inline VesSurfsys* getVesSurfsys() const noexcept {
        return pVesSurfsys;
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

    /// Return a list of vesicle species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side vesicle species.
    inline const std::vector<Spec*>& getVLHS() const noexcept {
        return pVLHS;
    }

    /// Set the vesicle  species on the left hand side of reaction.
    ///
    /// \param vlhs Vesicle surface species on the left hand side of reaction.
    void setVLHS(std::vector<Spec*> const& vlhs);

    /// Return a list of link species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side link species.
    inline const std::vector<LinkSpec*>& getLLHS() const noexcept {
        return pLLHS;
    }

    /// Set the link species on the left hand side of reaction.
    ///
    /// \param llhs Vesicle link species on the left hand side of reaction.
    void setLLHS(std::vector<LinkSpec*> const& llhs);

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

    /// Return a list of vesicle species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side vesicle species.
    inline const std::vector<Spec*>& getVRHS() const noexcept {
        return pVRHS;
    }

    /// Set the vesicle species on the right hand side of reaction.
    ///
    /// \param irhs Vesicle species on the right hand side of reaction.
    void setVRHS(std::vector<Spec*> const& vrhs);

    /// Return a list of link species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side link species.
    inline const std::vector<LinkSpec*>& getLRHS() const noexcept {
        return pLRHS;
    }

    /// Set the link species on the right hand side of reaction.
    ///
    /// \param irhs Link species on the right hand side of reaction.
    void setLRHS(std::vector<LinkSpec*> const& lrhs);

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

    /// Return a list of vesicle species dependencies for the reaction.
    ///
    /// \return List of pointers of dependencies.
    inline const std::vector<Spec*>& getVDeps() const noexcept {
        return pVDeps;
    }

    /// Set the vesicle species dependencies for the reaction.
    ///
    /// \param vdeps Vesicle species dependencies for the reaction.
    void setVDeps(std::vector<Spec*> const& vdeps);

    /// Get the order of the vesicle surface reaction.
    ///
    /// \return Order of the reaction.
    inline uint getOrder() const noexcept {
        return pOrder;
    }

    /// Get the rate constant of the vesicle surface reaction.
    ///
    /// \return Rate constant of the vesicle surface reaction.
    inline double getKcst() const noexcept {
        return pKcst;
    }

    /// Set the rate constant of the vesicle surface reaction.
    ///
    /// \param Rate constant of the vesicle surface reaction.
    void setKcst(double kcst);

    /// Get the mobility effect of the vesicle surface reaction.
    ///
    /// \return Immobilization of the vesicle surface reaction.
    inline Immobilization getImmobilization() const noexcept {
        return pImmobilization;
    }

    /// Get the maximim distance of the vesicle surface reaction.
    ///
    /// \return Maximim distance of the vesicle surface reaction.
    inline double getMaxDistance() const noexcept {
        return pMax_distance;
    }

    /// Get a list of all species.
    ///
    /// Returns a list of all species involved in this
    /// vesicle surface reaction, on both the left and righthand side
    /// and does not contain any duplicate members.
    ///
    /// \return List of pointers to the species.
    std::vector<Spec*> getAllSpecs() const;

    /// Get a list of all link species.
    ///
    /// Returns a list of all link species involved in this
    /// vesicle surface reaction, on both the left and righthand side
    /// and does not contain any duplicate members.
    ///
    /// \return List of pointers to the link species.
    std::vector<LinkSpec*> getAllLinkSpecs() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS
    ////////////////////////////////////////////////////////////////////////
    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once
    void _handleSelfDelete();

    /// Set order.
    ///
    void _setOrder();

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model* pModel;
    VesSurfsys* pVesSurfsys;

    std::vector<Spec*> pOLHS;
    std::vector<Spec*> pVLHS;
    std::vector<LinkSpec*> pLLHS;
    std::vector<Spec*> pSLHS;
    std::vector<Spec*> pVRHS;
    std::vector<Spec*> pSRHS;
    std::vector<Spec*> pORHS;
    std::vector<Spec*> pIRHS;
    std::vector<LinkSpec*> pLRHS;

    std::vector<Spec*> pVDeps;

    uint pOrder;
    double pKcst;
    Immobilization pImmobilization;

    double pMax_distance;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::model
