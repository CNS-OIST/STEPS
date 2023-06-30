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
#include "util/common.hpp"

namespace steps::model {

// Forward declarations.
class Model;
class Spec;
class LinkSpec;
class VesSurfsys;
class VesSReac;
class VesSDiff;
class Exocytosis;

// Auxiliary declarations.
typedef VesSurfsys* VesSurfsysP;
typedef std::map<std::string, VesSurfsysP> VesSurfsysPMap;
typedef VesSurfsysPMap::iterator VesSurfsysPMapI;
typedef VesSurfsysPMap::const_iterator VesSurfsysPMapCI;

////////////////////////////////////////////////////////////////////////////////
/// Vesicle Surface system.
/// Container that collects reactions involving a reactant
/// embedded in a vesicle membrane.
///
/// \warning Methods start with an underscore are not exposed to Python.

class VesSurfsys {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the vesicle surface system.
    /// \param model Pointer to the parent model.
    VesSurfsys(std::string const& id, Model* model);

    /// Destructor
    ~VesSurfsys();

    ////////////////////////////////////////////////////////////////////////
    // VESICLE SURFACE SYSTEM PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the vesicle surface system ID.
    ///
    /// \return ID of the vesicle surface system.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the vesicle surface system ID.
    ///
    /// \param id ID of the vesicle surface system.
    void setID(std::string const& id);

    /// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent model.
    inline Model* getModel() const noexcept {
        return pModel;
    }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return a vesicle surface reaction with name id.
    ///
    /// \param id ID of the surface reaction.
    /// \return Pointer to the surface reaction.
    VesSReac* getVesSReac(std::string const& id) const;

    /// Delete a vesicle surface reaction with name id.
    ///
    /// \param id ID of the vesicle surface reaction.
    void delVesSReac(std::string const& id) const;

    /// Return a list of all vesicle surface reactions.
    ///
    /// \return List of pointers to vesicle surface reactions.
    std::vector<VesSReac*> getAllVesSReacs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): VESICLE SURFACE DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Get a vesicle surface diffusion by its ID.
    ///
    /// \param id ID of the required diffusion.
    /// \return Pointer to the diffusion object.
    VesSDiff* getVesSDiff(std::string const& id) const;

    /// Delete a vesicle surface diffusion by its ID.
    ///
    /// \param id ID of the diffusion to be deleted.
    void delVesSDiff(std::string const& id) const;

    /// Get all vesicle surface diffusions stored in this vesicle surface system.
    ///
    /// \return A vector of pointers to the diffusion objects
    ///         stored in the system.
    std::vector<VesSDiff*> getAllVesSDiffs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): EXOCYTOSIS
    ////////////////////////////////////////////////////////////////////////

    /// Return a exocytotic reaction with name id.
    ///
    /// \param id ID of the exocytotic reaction.
    /// \return Pointer to the exocytotic reaction.
    Exocytosis* getExocytosis(std::string const& id) const;

    /// Delete an exocytotic reaction with name id.
    ///
    /// \param id ID of the exocytotic reaction.
    void delExocytosis(std::string const& id) const;

    /// Return a list of all exocytotic reactions.
    ///
    /// \return List of pointers to exocytotic reactions.
    std::vector<Exocytosis*> getAllExocytosis() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: EXOCYTOSIS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a exocytosis id is occupied.
    ///
    /// \param id ID of the exocytosis.
    void _checkExocytosisID(std::string const& id) const;

    /// Change an exocytosis id from o to n.
    ///
    /// \param o Old id of the exocytosis.
    /// \param n New id of the exocytosis.
    void _handleExocytosisIDChange(std::string const& o, std::string const& n);

    /// Add an exocytosis to the surface system.
    ///
    /// \param Pointer to the exocytosis.
    void _handleExocytosisAdd(Exocytosis* exocyt);

    /// Delete an exocytosis in the surface system.
    ///
    /// \param Pointer to the exocytosis.
    void _handleExocytosisDel(Exocytosis* exocyt);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: VESICLE SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a vesicle surface reaction id is occupied.
    ///
    /// \param id ID of the vesicle surface reaction.
    void _checkVesSReacID(std::string const& id) const;

    /// Change a vesicle surface reaction id from o to n.
    ///
    /// \param o Old id of the vesicle surface reaction.
    /// \param n New id of the vesicle surface reaction.
    void _handleVesSReacIDChange(std::string const& o, std::string const& n);

    /// Add a vesicle surface reaction to the surface system.
    ///
    /// \param vessreac Pointer to the vesicle surface reaction.
    void _handleVesSReacAdd(VesSReac* vessreac);

    /// Delete a vesicle surface reaction in the surface system.
    ///
    /// \param vessreac Pointer to the vesicle surface reaction.
    void _handleVesSReacDel(VesSReac* vessreac);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: VESICLE SURFACE DIFFUSION
    ////////////////////////////////////////////////////////////////////////
    /// Check if a vesicle surface diffusion id is occupied.
    ///
    /// \param id ID of the vesicle surface diffusion.
    void _checkVesSDiffID(std::string const& id) const;

    /// Change the id of a vesicle surface diffusion from o to n.
    ///
    /// \param o Old id of the vesicle surface diffusion.
    /// \param n New id of the vesicle surface diffusion.
    void _handleVesSDiffIDChange(std::string const& o, std::string const& n);

    /// Add a vesicle surface diffusion to the vesicle surface system.
    ///
    /// \param vessdiff Pointer to the vesicle surface diffusion.
    void _handleVesSDiffAdd(VesSDiff* vessdiff);

    /// Delete a vesicle surface diffusion in the vesicle surface system.
    ///
    /// \param vessdiff Pointer to the vesicle surface diffusion.
    void _handleVesSDiffDel(VesSDiff* vessdiff);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return all species in the surface system.
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
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    /// Count the vesicle surface reactions in the surface system.
    ///
    /// \return Number of vesicle surface reactions.
    inline uint _countVesSReacs() const noexcept {
        return pVesSReacs.size();
    }

    /// Get a vesicle surface reaction with index lidx.
    ///
    /// \param lidx Index of the vesicle surface reaction.
    /// \return Pointer to the vesicle surface reaction.
    VesSReac* _getVesSReac(uint lidx) const;

    /// Get all vesicle surface reactions in the vesicle surface system.
    ///
    /// \return Map of vesicle surface reactions.
    inline const std::map<std::string, VesSReac*>& _getAllVesSReacs() const noexcept {
        return pVesSReacs;
    }

    /// Count the vesicle surface diffusion rules in the volume system.
    ///
    /// \return Number of vesicle diffusion rules.
    inline uint _countVesSDiffs() const noexcept {
        return pVesSDiffs.size();
    }

    /// Get a vesicle surface diffusion rule with index lidx.
    ///
    /// \param lidx Index of the vesicle surface diffusion rule.
    /// \return Pointer to the vesicle surface diffusion rule.
    VesSDiff* _getVesSDiff(uint lidx) const;

    /// Get all vesicle surface diffusion rules in the vesicle surface system.
    ///
    /// \return List of pointers to vesicle surface diffusion rules.
    inline const std::map<std::string, VesSDiff*>& _getAllVesSDiffs() const noexcept {
        return pVesSDiffs;
    }

    /// Count the exocytotic reactions in the surface system.
    ///
    /// \return Number of exocytotic reactions.
    inline uint _countExocytosis() const noexcept {
        return pExocytosis.size();
    }

    /// Get a exocytotic reaction with index lidx.
    ///
    /// \param lidx Index of the exocytotic reaction.
    /// \return Pointer to the exocytotic reaction.
    Exocytosis* _getExocytosis(uint lidx) const;

    /// Get all exocytotic reactions in the surface system.
    ///
    /// \return Map of exocytotic reactions.
    inline const std::map<std::string, Exocytosis*>& _getAllExocytosis() const noexcept {
        return pExocytosis;
    }

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
    ////////////////////////////////////////////////////////////////////////
    /// Delete a species in the surface system.
    ///
    /// \param spec Pointer to the species.
    void _handleSpecDelete(Spec* spec);

    /// Delete a link species in the surface system.
    ///
    /// \param spec Pointer to the species.
    void _handleLinkSpecDelete(LinkSpec* spec);

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model* pModel;
    std::map<std::string, VesSReac*> pVesSReacs;
    std::map<std::string, VesSDiff*> pVesSDiffs;
    std::map<std::string, Exocytosis*> pExocytosis;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::model
