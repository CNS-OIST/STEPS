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

/// Raft system.
/// Container that collects reactions involving a reactant
/// embedded in a raft.
///
/// \warning Methods start with an underscore are not exposed to Python.

class Raftsys {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the raft system.
    /// \param model Reference to the parent model.
    Raftsys(std::string const& id, Model& model);

    Raftsys(const Raftsys&) = delete;
    Raftsys& operator=(const Raftsys&) = delete;

    /// Destructor
    ~Raftsys();

    ////////////////////////////////////////////////////////////////////////
    // RAFT SURFACE SYSTEM PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the surface system ID.
    ///
    /// \return ID of the raft surface system.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the raft surface system ID.
    ///
    /// \param id ID of the raft surface system.
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

    /// Return a raft surface reaction with name id.
    ///
    /// \param id ID of the raft surface reaction.
    /// \return Reference to the raft surface reaction.
    RaftSReac& getRaftSReac(std::string const& id) const;

    /// Delete a raft surface reaction with name id.
    ///
    /// \param id ID of the raft surface reaction.
    void delRaftSReac(std::string const& id) const;

    /// Return a list of all raft surface reactions.
    ///
    /// \return List of pointers to raft surface reactions.
    std::vector<RaftSReac*> getAllRaftSReacs() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: RAFT SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Change a raft surface reaction id from o to n.
    ///
    /// \param o Old id of the raft surface reaction.
    /// \param n New id of the raft surface reaction.
    void _handleRaftSReacIDChange(std::string const& o, std::string const& n);

    /// Add a raft surface reaction to the surface system.
    ///
    /// \param vessreac Reference to the raft surface reaction.
    void _handleRaftSReacAdd(RaftSReac& vessreac);

    /// Delete a raft surface reaction in the surface system.
    ///
    /// \param vessreac Reference to the raft surface reaction.
    void _handleRaftSReacDel(RaftSReac& vessreac);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): RaftEndocytosis
    ////////////////////////////////////////////////////////////////////////

    /// Return a endocytotic reaction with name id.
    ///
    /// \param id ID of the endocytotic reaction.
    /// \return Reference to the endocytotic reaction.
    RaftEndocytosis& getRaftEndocytosis(std::string const& id) const;

    /// Delete an endocytotic reaction with name id.
    ///
    /// \param id ID of the endocytotic reaction.
    void delRaftEndocytosis(std::string const& id) const;

    /// Return a list of all endocytotic reactions.
    ///
    /// \return List of pointers to endocytotic reactions.
    std::vector<RaftEndocytosis*> getAllRaftEndocytosiss() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: RaftEndocytosis
    ////////////////////////////////////////////////////////////////////////

    /// Check if a RaftEndocytosis id is occupied.
    ///
    /// \param id ID of the RaftEndocytosis.
    void _checkRaftEndocytosisID(std::string const& id) const;

    /// Change an RaftEndocytosis id from o to n.
    ///
    /// \param o Old id of the RaftEndocytosis.
    /// \param n New id of the RaftEndocytosis.
    void _handleRaftEndocytosisIDChange(std::string const& o, std::string const& n);

    /// Add an RaftEndocytosis to the surface system.
    ///
    /// \param Reference to the RaftEndocytosis.
    void _handleRaftEndocytosisAdd(RaftEndocytosis& endocyt);

    /// Delete an RaftEndocytosis in the surface system.
    ///
    /// \param Reference to the RaftEndocytosis.
    void _handleRaftEndocytosisDel(RaftEndocytosis& endocyt);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): RAFT DISSOLUTION
    ////////////////////////////////////////////////////////////////////////

    /// Return a raft dissolution with name id.
    ///
    /// \param id ID of the raft dissolution.
    /// \return Reference to the raft dissolution.
    RaftDis& getRaftDis(std::string const& id) const;

    /// Delete a raft dissolution with name id.
    ///
    /// \param id ID of the raft dissolution.
    void delRaftDis(std::string const& id) const;

    /// Return a list of all raft dissolutions.
    ///
    /// \return List of pointers to raft dissolution.
    std::vector<RaftDis*> getAllRaftDiss() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: RAFT DISSOLUTION
    ////////////////////////////////////////////////////////////////////////

    /// Check if a raft dissolution id is occupied.
    ///
    /// \param id ID of the raft dissolution.
    void _checkRaftDisID(std::string const& id) const;

    /// Change a raft dissolution id from o to n.
    ///
    /// \param o Old id of the raft dissolution.
    /// \param n New id of the raft dissolution.
    void _handleRaftDisIDChange(std::string const& o, std::string const& n);

    /// Add a raft dissolution to the raft system.
    ///
    /// \param Reference to the raft dissolution.
    void _handleRaftDisAdd(RaftDis& raftgen);

    /// Delete a raft dissolution in the raft system.
    ///
    /// \param Reference to the raft dissolution.
    void _handleRaftDisDel(RaftDis& raftgen);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return all species in the raft system.
    ///
    /// \return List of pointers to the species.
    util::flat_set<Spec*> getAllSpecs() const;

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

    /// Count the raft surface reactions in the surface system.
    ///
    /// \return Number of raft surface reactions.
    inline uint _countRaftSReacs() const noexcept {
        return pRaftSReacs.size();
    }

    /// Get a raft surface reaction with index lidx.
    ///
    /// \param lidx Index of the raft surface reaction.
    /// \return Reference to the raft surface reaction.
    RaftSReac& _getRaftSReac(uint lidx) const;

    /// Get all raft surface reactions in the raft surface system.
    ///
    /// \return Map of raft surface reactions.
    inline const std::map<std::string, RaftSReac*>& _getAllRaftSReacs() const noexcept {
        return pRaftSReacs;
    }

    /// Count the endocytotic reactions in the raft system.
    ///
    /// \return Number of endocytotic reactions.
    inline uint _countRaftEndocytosis() const noexcept {
        return pRaftEndocytosis.size();
    }

    /// Get a endocytotic reaction with index lidx.
    ///
    /// \param lidx Index of the endocytotic reaction.
    /// \return Reference to the endocytotic reaction.
    RaftEndocytosis& _getRaftEndocytosis(uint lidx) const;

    /// Get all endocytotic reactions in the raft system.
    ///
    /// \return Map of endocytotic reactions.
    inline const std::map<std::string, RaftEndocytosis*>& _getAllRaftEndocytosiss() const noexcept {
        return pRaftEndocytosis;
    }

    /// Count the raft dissolutions in the raft system.
    ///
    /// \return Number of raft dissolutions.
    inline uint _countRaftDiss() const noexcept {
        return pRaftDiss.size();
    }

    /// Get a raft dissolution with index lidx.
    ///
    /// \param lidx Index of the raft dissolution.
    /// \return Reference to the raft dissolution.
    RaftDis& _getRaftDis(uint lidx) const;

    /// Get all raft dissolutions in the raft system.
    ///
    /// \return Map of raft dissolutions.
    inline const std::map<std::string, RaftDis*>& _getAllRaftDiss() const noexcept {
        return pRaftDiss;
    }

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
    ////////////////////////////////////////////////////////////////////////
    /// Delete a species in the surface system.
    ///
    /// \param spec Reference to the species.
    void _handleSpecDelete(Spec& spec);

    ////////////////////////////////////////////////////////////////////////

  private:
    /// Check if a raft surface reaction id is occupied.
    ///
    /// \param id ID of the raft surface reaction.
    void _checkRaftSReacID(std::string const& id) const;

    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model& pModel;
    std::map<std::string, RaftSReac*> pRaftSReacs;
    std::map<std::string, RaftEndocytosis*> pRaftEndocytosis;
    std::map<std::string, RaftDis*> pRaftDiss;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::model
