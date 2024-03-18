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

#include <string>
#include <vector>

#include "fwd.hpp"
#include "util/collections.hpp"

namespace steps::model {

/// Raft dissolution.
///
/// A RaftDis object describes a raft dissolution event which takes place on a
/// surface system, i.e. a patch between two compartments.
///
/// \warning Methods start with an underscore are not exposed to Python.

class RaftDis {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the raft dissolution.
    /// \param raftsys Reference to the parent raft surface system.
    /// \param spec_signature The 'species signature' that- when the given
    //// species counts in the raft are at or below this number- will
    /// trigger the raft dissolution by the specified rate constant
    /// \param kcst The rate constant, first-order

    RaftDis(std::string const& id,
            Raftsys& raftsys,
            std::vector<Spec*> const& spec_signature,
            double kcst = 0.0);

    RaftDis(const RaftDis&) = delete;
    RaftDis& operator=(const RaftDis&) = delete;

    /// Destructor
    ~RaftDis();

    ////////////////////////////////////////////////////////////////////////
    // SURFACE REACTION RULE PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the raft dissolution rule ID.
    ///
    /// \return ID of the raft dissolution.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the raft dissolution rule ID.
    ///
    /// \param id ID of the raft dissolution.
    void setID(std::string const& id);

    /// Return a reference to the parent raft system.
    ///
    /// \return Reference to the raft system.
    inline Raftsys& getRaftsys() const noexcept {
        return pRaftsys;
    }

    /// Return a reference to the parent model.
    ///
    /// \return Reference to the parent model.
    inline Model& getModel() const noexcept {
        return pModel;
    }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return a pointer to the created Raft.
    ///
    /// \return Pointer to the raft.
    inline Raft* getRaft() const noexcept {
        return pRaft;
    }

    /// Return a list, species signature.
    ///
    /// \return List of pointers of species.
    inline const std::vector<Spec*>& getSpecSignature() const noexcept {
        return pSpecSignature;
    }

    /// Get a list of all species, does not contain any duplicate members.
    ///
    /// \return List of pointers to the species.
    util::flat_set<Spec*> getAllSpecs() const;

    /// Get the rate constant of the raft dissolution
    ///
    /// \return Rate constant of the raft dissolution
    inline double getKcst() const noexcept {
        return pKcst;
    }

    /// Set the rate constant of the raft dissolution
    ///
    /// \param Rate constant of the raft dissolution
    void setKcst(double kcst);

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
    Model& pModel;
    Raftsys& pRaftsys;

    std::vector<Spec*> pSpecSignature;
    Raft* pRaft{};

    double pKcst;
};

}  // namespace steps::model
