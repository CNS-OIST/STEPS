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
class VesUnbind;
class Volsys;
class Model;
class Spec;
class LinkSpec;
class Vesicle;

// Auxiliary declarations.
typedef VesUnbind* VesUnbindP;
typedef std::map<std::string, VesUnbindP> VesUnbindPMap;
typedef VesUnbindPMap::iterator VesUnbindPMapI;
typedef VesUnbindPMap::const_iterator VesUnbindPMapCI;
typedef std::vector<VesUnbindP> VesUnbindPVec;
typedef VesUnbindPVec::iterator VesUnbindPVecI;
typedef VesUnbindPVec::const_iterator VesUnbindPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Vesicle unbinding reaction in a volume system.
///

class VesUnbind {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the reaction.
    /// \param volsys Pointer to the parent volume system.
    /// \param link The link species that will dissociate, causing unbinding
    /// \param vesicle1 The 1st vesicle for the 1st product species
    /// \param vesicle2 The 2nd vesicle for the 2nd product species
    /// \param product1 The 1st product, located on Vesicle 1
    /// \param product2 The 2nd product, located on Vesicle 2
    /// \param kcst Rate constant for the reaction.
    /// \param immobilization A flag used to control the mobility effect of this
    ///			reaction. VesUnBind can mobilise or immobilize vesicles.

    VesUnbind(std::string const& id,
              Volsys* volsys,
              LinkSpec* link1,
              LinkSpec* link2,
              Vesicle* vesicle1,
              Spec* product1,
              Vesicle* vesicle2,
              Spec* product2,
              double kcst = 0.0,
              Immobilization immobilization = NO_EFFECT);

    /// Destructor
    ~VesUnbind();

    ////////////////////////////////////////////////////////////////////////
    // VESICLE UNBINDING REACTION PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the ID.
    ///
    /// \return ID .
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the ID.
    ///
    /// \param id ID.
    void setID(std::string const& id);

    /// Return a pointer to the parent volume system.
    ///
    /// \return Pointer to the volume system.
    inline Volsys* getVolsys() const noexcept {
        return pVolsys;
    }

    /// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent Model.
    inline Model* getModel() const noexcept {
        return pModel;
    }

    /// Get the mobility effect of the vesicle unbinding reaction.
    ///
    /// \return Immobilization of the vesicle unbinding reaction.
    inline Immobilization getImmobilization() const noexcept {
        return pImmobilization;
    }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Get the reactant pairs
    ///
    /// \return Map

    // NOTE: need to make sure this data is copied to def
    inline const std::pair<Vesicle*, Spec*> getProducts1() const noexcept {
        return pProducts1;
    }

    inline const std::pair<Vesicle*, Spec*> getProducts2() const noexcept {
        return pProducts2;
    }

    // Get the link specs that will trigger the unbinding event
    inline const std::pair<Vesicle*, LinkSpec*> getLinks1() const noexcept {
        return pLinks1;
    }
    inline const std::pair<Vesicle*, LinkSpec*> getLinks2() const noexcept {
        return pLinks2;
    }

    /// Return all species involved in the reaction.
    ///
    /// \return Vector of pointers to the species.
    std::vector<Spec*> getAllSpecs() const;

    /// Return the rate constant of the reaction.
    ///
    /// \return The rate constant of the reaction.
    inline double getKcst() const noexcept {
        return pKcst;
    }

    /// Set or reset the rate constant of the reaction.
    ///
    /// \param kcst The rate constant of the reaction.
    void setKcst(double kcst);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////

    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// \warning Will only be called once
    void _handleSelfDelete();

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model* pModel;
    Volsys* pVolsys;

    std::pair<Vesicle*, Spec*> pProducts1;
    std::pair<Vesicle*, Spec*> pProducts2;

    std::pair<Vesicle*, LinkSpec*> pLinks1;
    std::pair<Vesicle*, LinkSpec*> pLinks2;

    double pKcst;

    Immobilization pImmobilization;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::model
