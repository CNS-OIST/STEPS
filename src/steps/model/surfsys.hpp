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
#include <vector>

#include "fwd.hpp"
#include "util/collections.hpp"

namespace steps::model {

/// Surface system.
/// Container that collects reactions involving a reactant
/// embedded in a membrane.
///
/// \warning Methods start with an underscore are not exposed to Python.

class Surfsys {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the surface system.
    /// \param model Reference to the parent model.
    Surfsys(std::string const& id, Model& model);

    Surfsys(const Surfsys&) = delete;
    Surfsys& operator=(const Surfsys&) = delete;

    /// Destructor
    ~Surfsys();

    ////////////////////////////////////////////////////////////////////////
    // SURFACE SYSTEM PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the surface system ID.
    ///
    /// \return ID of the surface system.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the surface system ID.
    ///
    /// \param id ID of the surface system.
    void setID(std::string const& id);

    /// Return a reference to the parent model.
    ///
    /// \return Reference to the parent model.
    inline Model& getModel() const noexcept {
        return pModel;
    }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return a surface reaction with name id.
    ///
    /// \param id ID of the surface reaction.
    /// \return Reference to the surface reaction.
    SReac& getSReac(std::string const& id) const;

    /// Delete a surface reaction with name id.
    ///
    /// \param id ID of the surface reaction.
    void delSReac(std::string const& id);

    /// Return a list of all surface reactions.
    ///
    /// \return List of pointers to surface reactions.
    std::vector<SReac*> getAllSReacs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Get a diffusion by its ID.
    ///
    /// \param id ID of the required diffusion.
    /// \return Reference to the diffusion object.
    Diff& getDiff(std::string const& id) const;

    /// Delete a diffusion by its ID.
    ///
    /// \param id ID of the diffusion to be deleted.
    void delDiff(std::string const& id);

    /// Get all diffusions stored in this surface system.
    ///
    /// \return A vector of pointers to the diffusion objects
    ///         stored in the system.
    std::vector<Diff*> getAllDiffs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): VOLTAGE-DEPENDENT REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return a voltage-dependent reaction with name id.
    ///
    /// \param id ID of the voltage-dependent reaction.
    /// \return Reference to the voltage-dependent reaction.
    VDepSReac& getVDepSReac(std::string const& id) const;

    /// Delete a voltage-dependent reaction with name id.
    ///
    /// \param id ID of the voltage-dependent reaction.
    void delVDepSReac(std::string const& id);

    /// Return a list of all voltage-dependent reactions.
    ///
    /// \return List of pointers to voltage-dependent transitions.
    std::vector<VDepSReac*> getAllVDepSReacs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): OHMIC CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return an ohmic current with name id.
    ///
    /// \param id ID of the ohmic current.
    /// \return Reference to the ohmic current.
    OhmicCurr& getOhmicCurr(std::string const& id) const;

    /// Delete an ohmic current with name id.
    ///
    /// \param id ID of the ohmic current.
    void delOhmicCurr(std::string const& id);

    /// Return a list of all ohmic currents.
    ///
    /// \return List of pointers to ohmic currents.
    std::vector<OhmicCurr*> getAllOhmicCurrs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): GHK CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return an ghk current with name id.
    ///
    /// \param id ID of the ghk current.
    /// \return Reference to the ghk current.
    GHKcurr& getGHKcurr(std::string const& id) const;

    /// Delete a ghk current with name id.
    ///
    /// \param id ID of the ghk current.
    void delGHKcurr(std::string const& id);

    /// Return a list of all ghk currents.
    ///
    /// \return List of pointers to ghk currents.
    std::vector<GHKcurr*> getAllGHKcurrs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): RAFT GENESIS
    ////////////////////////////////////////////////////////////////////////

    /// Return a raft genesis with name id.
    ///
    /// \param id ID of the raft genesis.
    /// \return Reference to the raft genesis.
    RaftGen& getRaftGen(std::string const& id) const;

    /// Delete a raft genesis with name id.
    ///
    /// \param id ID of the raft genesis.
    void delRaftGen(std::string const& id);

    /// Return a list of all raft geneses.
    ///
    /// \return List of pointers to raft geneses.
    std::vector<RaftGen*> getAllRaftGens() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): ENDOCYTOSIS
    ////////////////////////////////////////////////////////////////////////

    /// Return a endocytotic reaction with name id.
    ///
    /// \param id ID of the endocytotic reaction.
    /// \return Reference to the endocytotic reaction.
    Endocytosis& getEndocytosis(std::string const& id) const;

    /// Delete an endocytotic reaction with name id.
    ///
    /// \param id ID of the endocytotic reaction.
    void delEndocytosis(std::string const& id);

    /// Return a list of all endocytotic reactions.
    ///
    /// \return List of pointers to endocytotic reactions.
    std::vector<Endocytosis*> getAllEndocytosis() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): COMPLEX SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return a list of all complex reactions.
    ///
    /// \return List of pointers to complex reactions.
    std::vector<ComplexSReac*> getAllComplexSReacs() const;

    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return all species in the surface system.
    ///
    /// \return List of pointers to the species.
    std::vector<Spec*> getAllSpecs() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a surface reaction id is occupied.
    ///
    /// \param id ID of the surface reaction.
    void _checkSReacID(std::string const& id) const;

    /// Change a surface reaction id from o to n.
    ///
    /// \param o Old id of the surface reaction.
    /// \param n New id of the surface reaction.
    void _handleSReacIDChange(std::string const& o, std::string const& n);

    /// Add a surface reaction to the surface system.
    ///
    /// \param Reference to the surface reaction.
    void _handleSReacAdd(SReac& sreac);

    /// Delete a surface reaction in the surface system.
    ///
    /// \param Reference to the surface reaction.
    void _handleSReacDel(SReac& sreac);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: COMPLEX SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Add a complex reaction to the surface system.
    ///
    /// \param reac Pointer to the complex reaction.
    void _handleComplexSReacAdd(ComplexSReac& reac);

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

    /// Add a diffusion to the surface system.
    ///
    /// \param diff Reference to the diffusion.
    void _handleDiffAdd(Diff& diff);

    /// Delete a diffusion in the surface system.
    ///
    /// \param diff Reference to the diffusion.
    void _handleDiffDel(Diff& diff);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: VOLTAGE-DEPENDENT REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a voltage-dependent reaction id is occupied.
    ///
    /// \param id ID of the voltage-dependent reaction.
    void _checkVDepSReacID(std::string const& id) const;

    /// Change a voltage-dependent reaction id from o to n.
    ///
    /// \param o Old id of the voltage-dependent reaction.
    /// \param n New id of the voltage-dependent reaction.
    void _handleVDepSReacIDChange(std::string const& o, std::string const& n);

    /// Add a voltage-dependent reaction to the surface system.
    ///
    /// \param Reference to the voltage-dependent reaction.
    void _handleVDepSReacAdd(VDepSReac& vdepsreac);

    /// Delete a voltage-dependent reaction in the surface system.
    ///
    /// \param Reference to the voltage-dependent reaction.
    void _handleVDepSReacDel(VDepSReac& vdepsreac);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: OHMIC CURRENT
    ////////////////////////////////////////////////////////////////////////

    /// Check if an ohmic current id is occupied.
    ///
    /// \param id ID of the ohmic current.
    void _checkOhmicCurrID(std::string const& id) const;

    /// Change an ohmic current id from o to n.
    ///
    /// \param o Old id of the ohmic current.
    /// \param n New id of the ohmic current.
    void _handleOhmicCurrIDChange(std::string const& o, std::string const& n);

    /// Add a ohmic current to the surface system.
    ///
    /// \param Reference to the ohmic current.
    void _handleOhmicCurrAdd(OhmicCurr& ohmiccurr);

    /// Delete an ohmic current in the surface system.
    ///
    /// \param Reference to the ohmic current.
    void _handleOhmicCurrDel(OhmicCurr& ohmiccurr);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: GHK CURRENT
    ////////////////////////////////////////////////////////////////////////

    /// Check if a GHK current id is occupied.
    ///
    /// \param id ID of the GHK current.
    void _checkGHKcurrID(std::string const& id) const;

    /// Change a GHK current id from o to n.
    ///
    /// \param o Old id of the GHK current.
    /// \param n New id of the GHK current.
    void _handleGHKcurrIDChange(std::string const& o, std::string const& n);

    /// Add a GHK current to the surface system.
    ///
    /// \param Reference to the GHK current.
    void _handleGHKcurrAdd(GHKcurr& ghkcurr);

    /// Delete a GHK current in the surface system.
    ///
    /// \param Reference to the GHK current.
    void _handleGHKcurrDel(GHKcurr& GHKcurr);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: RAFT GENESIS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a raft genesis id is occupied.
    ///
    /// \param id ID of the raft genesis.
    void _checkRaftGenID(std::string const& id) const;

    /// Change a raft genesis id from o to n.
    ///
    /// \param o Old id of the raft genesis.
    /// \param n New id of the raft genesis.
    void _handleRaftGenIDChange(std::string const& o, std::string const& n);

    /// Add a raft genesis to the surface system.
    ///
    /// \param Reference to the raft genesis.
    void _handleRaftGenAdd(RaftGen& raftgen);

    /// Delete a raft genesis in the surface system.
    ///
    /// \param Reference to the raft genesis.
    void _handleRaftGenDel(RaftGen& raftgen);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: ENDOCYTOSIS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a endocytosis id is occupied.
    ///
    /// \param id ID of the endocytosis.
    void _checkEndocytosisID(std::string const& id) const;

    /// Change an endocytosis id from o to n.
    ///
    /// \param o Old id of the endocytosis.
    /// \param n New id of the endocytosis.
    void _handleEndocytosisIDChange(std::string const& o, std::string const& n);

    /// Add an endocytosis to the surface system.
    ///
    /// \param Reference to the endocytosis.
    void _handleEndocytosisAdd(Endocytosis& endocyt);

    /// Delete an endocytosis in the surface system.
    ///
    /// \param Reference to the endocytosis.
    void _handleEndocytosisDel(Endocytosis& endocyt);

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

    /// Count the surface reactions in the surface system.
    ///
    /// \return Number of surface reactions.
    inline uint _countSReacs() const noexcept {
        return pSReacs.size();
    }

    /// Get a surface reaction with index lidx.
    ///
    /// \param lidx Index of the surface reaction.
    /// \return Reference to the surface reaction.
    SReac& _getSReac(uint lidx) const;

    /// Count the complex surface reactions in the surface system.
    ///
    /// \return Number of complex reactions.
    inline uint _countComplexSReacs() const {
        return pComplexSReacs.size();
    }

    /// Get a complex surface reaction with index lidx
    ///
    /// \param lidx Index of the complex reaction.
    /// \return Pointer to the complex reaction.
    ComplexSReac& _getComplexSReac(uint lidx) const;

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

    /// Count the voltage-dependent reactions in the surface system.
    ///
    /// \return Number of voltage-dependent reactions.
    inline uint _countVDepSReacs() const noexcept {
        return pVDepSReacs.size();
    }

    /// Get a voltage-dependent reactions with index lidx.
    ///
    /// \param lidx Index of the voltage-dependent reaction.
    /// \return Reference to the voltage-dependent reaction.
    VDepSReac& _getVDepSReac(uint lidx) const;

    /// Count the ohmic currents in the surface system.
    ///
    /// \return Number of ohmic currents.
    inline uint _countOhmicCurrs() const noexcept {
        return pOhmicCurrs.size();
    }

    /// Get ohmic current object with index lidx.
    ///
    /// \param lidx index of the ohmic current.
    /// \ return Reference to the ohmic current.
    OhmicCurr& _getOhmicCurr(uint lidx) const;

    /// Count the ghk currents in the surface system.
    ///
    /// \return Number of ghk currents.
    inline uint _countGHKcurrs() const noexcept {
        return pGHKcurrs.size();
    }

    /// Get ghk current object with index lidx.
    ///
    /// \param lidx index of the ghk current.
    /// \ return Reference to the ghk current.
    GHKcurr& _getGHKcurr(uint lidx) const;

    /// Count the raft geneses in the surface system.
    ///
    /// \return Number of raft geneses.
    inline uint _countRaftGens() const noexcept {
        return pRaftGens.size();
    }

    /// Get a raft genesis with index lidx.
    ///
    /// \param lidx Index of the raft genesis.
    /// \return Reference to the raft genesis.
    RaftGen& _getRaftGen(uint lidx) const;

    /// Count the endocytotic reactions in the surface system.
    ///
    /// \return Number of endocytotic reactions.
    inline uint _countEndocytosis() const noexcept {
        return pEndocytosis.size();
    }

    /// Get a endocytotic reaction with index lidx.
    ///
    /// \param lidx Index of the endocytotic reaction.
    /// \return Reference to the endocytotic reaction.
    Endocytosis& _getEndocytosis(uint lidx) const;

    /// Get all surface reactions in the surface system.
    ///
    /// \return Map of surface reactions.
    inline const std::map<std::string, SReac*>& _getAllSReacs() const noexcept {
        return pSReacs;
    }

    /// Get all complex surface reactions in the surface system.
    ///
    /// \return Map of complex surface reactions.
    inline const std::map<std::string, ComplexSReac*>& _getAllComplexSReacs() const noexcept {
        return pComplexSReacs;
    }

    /// Get all diffusion rules in the surface system.
    ///
    /// \return List of pointers to diffusion rules.
    inline const std::map<std::string, Diff*>& _getAllDiffs() const noexcept {
        return pDiffs;
    }

    /// Get all voltage-dependent reactions in the system
    ///
    /// \return Map of voltage-dependent reactions
    inline const std::map<std::string, VDepSReac*>& _getAllVDepSReacs() const noexcept {
        return pVDepSReacs;
    }

    /// Get all ohmic currents in the system
    ///
    /// \return Map of ohmic currents
    inline const std::map<std::string, OhmicCurr*>& _getAllOhmicCurrs() const noexcept {
        return pOhmicCurrs;
    }

    /// Get all ghk currents in the system
    ///
    /// \return Map of ghk currents
    inline const std::map<std::string, GHKcurr*>& _getAllGHKcurrs() const noexcept {
        return pGHKcurrs;
    }

    /// Get all raft geneses in the surface system.
    ///
    /// \return Map of raft geneses.
    inline const std::map<std::string, RaftGen*>& _getAllRaftGens() const noexcept {
        return pRaftGens;
    }

    /// Get all endocytotic reactions in the surface system.
    ///
    /// \return Map of endocytotic reactions.
    inline const std::map<std::string, Endocytosis*>& _getAllEndocytosis() const noexcept {
        return pEndocytosis;
    }

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
    ////////////////////////////////////////////////////////////////////////
    /// Delete a species in the surface system.
    ///
    /// \param spec Reference to the species.
    void _handleSpecDelete(Spec& spec);

    /// Delete a channel in the surface system
    ///
    /// \param chan Reference to the channel.
    void _handleChanDelete(Chan& chan);

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model& pModel;
    std::map<std::string, SReac*> pSReacs;
    std::map<std::string, ComplexSReac*> pComplexSReacs;
    std::map<std::string, Diff*> pDiffs;
    std::map<std::string, VDepSReac*> pVDepSReacs;
    std::map<std::string, OhmicCurr*> pOhmicCurrs;
    std::map<std::string, GHKcurr*> pGHKcurrs;
    std::map<std::string, RaftGen*> pRaftGens;
    std::map<std::string, Endocytosis*> pEndocytosis;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::model
