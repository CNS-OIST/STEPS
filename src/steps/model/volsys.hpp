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

#include <map>
#include <string>

#include "fwd.hpp"
#include "util/collections.hpp"

namespace steps::model {

/// Volume system
///
/// Container that collects reactions
/// involving a reactant enclosed in a compartment.
///
/// \warning Methods start with an underscore are not exposed to Python.
class Volsys {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the volume system.
    /// \param model Reference to the parent model object.
    Volsys(std::string const& id, Model& model);

    Volsys(const Volsys&) = delete;
    Volsys& operator=(const Volsys&) = delete;

    /// Destructor
    ~Volsys();

    ////////////////////////////////////////////////////////////////////////
    // VOLUME SYSTEM PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the volume system ID.
    ///
    /// \return A string of volume system ID.
    const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the volume system ID.
    ///
    /// \param id The new volume system ID.
    void setID(std::string const& id);

    /// Return a reference to the parent model.
    ///
    /// \return Reference to the parent model.
    Model& getModel() const noexcept {
        return pModel;
    }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Get a reaction by its ID.
    ///
    /// \param id ID of the required reaction.
    /// \return Reference to the reaction object.
    Reac& getReac(std::string const& id) const;

    /// Delete a reaction by its ID.
    ///
    /// \param id ID of the reaction to be deleted.
    void delReac(std::string const& id) const;

    /// Get all reactions stored in this volume system.
    ///
    /// \return A vector of pointers to the reaction objects
    ///         stored in the system.
    util::flat_set<Reac*> getAllReacs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): COMPLEX REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Get all complex reactions stored in this volume system.
    ///
    /// \return A vector of pointers to the complex reaction objects
    ///         stored in the system.
    util::flat_set<ComplexReac*> getAllComplexReacs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Get a difussion by its ID.
    ///
    /// \param id ID of the required difussion.
    /// \return Reference to the diffusion object.
    Diff& getDiff(std::string const& id) const;

    /// Delete a diffusion by its ID.
    ///
    /// \param id ID of the diffusion to be deleted.
    void delDiff(std::string const& id) const;

    /// Get all diffusions stored in this volume system.
    ///
    /// \return A vector of pointers to the diffusion objects
    ///         stored in the system.
    util::flat_set<Diff*> getAllDiffs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): VESICLE BINDINGS
    ////////////////////////////////////////////////////////////////////////

    /// Get a vesicle binding reaction by its ID.
    ///
    /// \param id ID of the required vesicle binding.
    /// \return Reference to the vesicle binding object.
    VesBind& getVesBind(std::string const& id) const;

    /// Delete a vesicle binding reaction by its ID.
    ///
    /// \param id ID of the vesicle binding reaction to be deleted.
    void delVesBind(std::string const& id) const;

    /// Get all vesicle binding reactions stored in this volume system.
    ///
    /// \return A vector of pointers to the vesicle binding reaction objects
    ///         stored in the system.
    util::flat_set<VesBind*> getAllVesBinds() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): VESICLE UNBINDINGS
    ////////////////////////////////////////////////////////////////////////

    /// Get a vesicle unbinding reaction by its ID.
    ///
    /// \param id ID of the required vesicle unbinding.
    /// \return Reference to the vesicle unbinding object.
    VesUnbind& getVesUnbind(std::string const& id) const;

    /// Delete a vesicle unbinding reaction by its ID.
    ///
    /// \param id ID of the vesicle unbinding reaction to be deleted.
    void delVesUnbind(std::string const& id) const;

    /// Get all vesicle unbinding reactions stored in this volume system.
    ///
    /// \return A vector of pointers to the vesicle unbinding reaction objects
    ///         stored in the system.
    util::flat_set<VesUnbind*> getAllVesUnbinds() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Get all species stored in this volume system.
    ///
    /// \return A vector of pointers to the species stored in the system.
    ///
    /// This method returns a list of all species involved in this volume system,
    /// no duplicate member is included.
    util::flat_set<Spec*> getAllSpecs() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////

    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once.
    void _handleSelfDelete();

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a reaction id is occupied.
    ///
    /// \param id ID of the reaction.
    void _checkReacID(std::string const& id) const;

    /// Change the id of a reaction from o to n.
    ///
    /// \param o Old id of the reaction.
    /// \param n New id of the reaction.
    void _handleReacIDChange(std::string const& o, std::string const& n);

    /// Add a reaction to the volume system.
    ///
    /// \param reac Reference to the reaction.
    void _handleReacAdd(Reac& reac);

    /// Delete a reaction in the volume system.
    ///
    /// \param reac Reference to the reaction.
    void _handleReacDel(Reac& reac);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: COMPLEX REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Add a reaction to the volume system.
    ///
    /// \param reac Pointer to the reaction.
    void _handleComplexReacAdd(ComplexReac& reac);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DIFFUSION
    ////////////////////////////////////////////////////////////////////////
    /// Check if a diffusion id is occupied.
    ///
    /// \param id ID of the diffusion.
    void _checkDiffID(std::string const& id) const;

    /// Change the id of a diffusion from o to n.
    ///
    /// \param o Old id of the diffusion.
    /// \param n New id of the diffusion.
    void _handleDiffIDChange(std::string const& o, std::string const& n);

    /// Add a diffusion to the volume system.
    ///
    /// \param diff Reference to the diffusion.
    void _handleDiffAdd(Diff& diff);

    /// Delete a diffusion in the volume system.
    ///
    /// \param diff Reference to the diffusion.
    void _handleDiffDel(Diff& diff);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: VESICLE BINDINGS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a vesicle binding reaction id is occupied.
    ///
    /// \param id ID of the vesicle binding reaction.
    void _checkVesBindID(std::string const& id) const;

    /// Change the id of a vesicle binding reaction from o to n.
    ///
    /// \param o Old id of the vesicle binding reaction.
    /// \param n New id of the vesicle binding reaction.
    void _handleVesBindIDChange(std::string const& o, std::string const& n);

    /// Add a vesicle binding reaction to the volume system.
    ///
    /// \param reac Reference to the vesicle binding reaction.
    void _handleVesBindAdd(VesBind& vesbind);

    /// Delete a vesicle binding reaction in the volume system.
    ///
    /// \param reac Reference to the vesicle binding reaction.
    void _handleVesBindDel(VesBind& vesbind);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: VESICLE UNBINDINGS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a vesicle unbinding reaction id is occupied.
    ///
    /// \param id ID of the vesicle unbinding reaction.
    void _checkVesUnbindID(std::string const& id) const;

    /// Change the id of a vesicle unbinding reaction from o to n.
    ///
    /// \param o Old id of the vesicle unbinding reaction.
    /// \param n New id of the vesicle unbinding reaction.
    void _handleVesUnbindIDChange(std::string const& o, std::string const& n);

    /// Add a vesicle unbinding reaction to the volume system.
    ///
    /// \param vesunbind Reference to the vesicle unbinding reaction.
    void _handleVesUnbindAdd(VesUnbind& vesunbind);

    /// Delete a vesicle unbinding reaction in the volume system.
    ///
    /// \param vesunbind Reference to the vesicle unbinding reaction.
    void _handleVesUnbindDel(VesUnbind& vesunbind);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    /// Count the reactions in the volume system.
    ///
    /// \return Number of reactions.
    inline uint _countReacs() const noexcept {
        return pReacs.size();
    }

    /// Get a reaction with index lidx
    ///
    /// \param lidx Index of the reaction.
    /// \return Reference to the reaction.
    Reac& _getReac(uint lidx) const;

    /// Count the complex reactions in the volume system.
    ///
    /// \return Number of complex reactions.
    uint _countComplexReacs() const noexcept {
        return pComplexReacs.size();
    }

    /// Get a complex reaction with index lidx
    ///
    /// \param lidx Index of the complex reaction.
    /// \return Pointer to the complex reaction.
    ComplexReac& _getComplexReac(uint lidx) const;

    /// Get all species in the volume system.
    ///
    /// \return A list of pointers to the species.
    inline const std::map<std::string, Reac*>& _getAllReacs() const noexcept {
        return pReacs;
    }

    /// Get all complex reacs in the volume system.
    ///
    /// \return A list of pointers to the complex reacs.
    const std::map<std::string, ComplexReac*>& _getAllComplexReacs() const noexcept {
        return pComplexReacs;
    }

    /// Count the diffusion rules in the volume system.
    ///
    /// \return Number of diffusion rules.
    inline uint _countDiffs() const noexcept {
        return pDiffs.size();
    }

    /// Get a diffusion rule with index lidx.
    ///
    /// \param lidx Index of the diffusion rule.
    /// \return Reference to the diffusion rule.
    Diff& _getDiff(uint lidx) const;

    /// Get all diffusion rules in the volume system.
    ///
    /// \return List of pointers to diffusion rules.
    inline const std::map<std::string, Diff*>& _getAllDiffs() const noexcept {
        return pDiffs;
    }

    /// Count the vesicle binding reactions in the volume system.
    ///
    /// \return Number of vesicle binding reactions.
    inline uint _countVesBinds() const noexcept {
        return pVesBinds.size();
    }

    /// Get a vesicle binding reaction with index lidx
    ///
    /// \param lidx Index of the vesicle binding reaction.
    /// \return Reference to the vesicle binding reaction.
    VesBind& _getVesBind(uint lidx) const;

    /// Get all species in the volume system.
    ///
    /// \return A list of pointers to the species.
    inline const std::map<std::string, VesBind*>& _getAllVesBinds() const noexcept {
        return pVesBinds;
    }

    /// Count the vesicle unbinding reactions in the volume system.
    ///
    /// \return Number of vesicle unbinding reactions.
    inline uint _countVesUnbinds() const noexcept {
        return pVesUnbinds.size();
    }

    /// Get a vesicle unbinding reaction with index lidx
    ///
    /// \param lidx Index of the vesicle unbinding reaction.
    /// \return Reference to the vesicle unbinding reaction.
    VesUnbind& _getVesUnbind(uint lidx) const;

    /// Get all species in the volume system.
    ///
    /// \return A list of pointers to the species.
    const std::map<std::string, VesUnbind*>& _getAllVesUnbinds() const noexcept {
        return pVesUnbinds;
    }

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
    ////////////////////////////////////////////////////////////////////////

    /// Delete a species.
    ///
    /// \param Reference to the species.
    void _handleSpecDelete(Spec& spec);

    /// Delete a species.
    ///
    /// \param Reference to the link species.
    void _handleLinkSpecDelete(LinkSpec& lspec);

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model& pModel;
    std::map<std::string, Reac*> pReacs;
    std::map<std::string, ComplexReac*> pComplexReacs;
    std::map<std::string, Diff*> pDiffs;

    std::map<std::string, VesBind*> pVesBinds;
    std::map<std::string, VesUnbind*> pVesUnbinds;
};

}  // namespace steps::model
