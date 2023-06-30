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
class RaftGen;
class Surfsys;
class Model;
class Spec;
class Raft;

// Auxiliary declarations.
typedef RaftGen* RaftGenP;
typedef std::map<std::string, RaftGenP> RaftGenPMap;
typedef RaftGenPMap::iterator RaftGenPMapI;
typedef RaftGenPMap::const_iterator RaftGenPMapCI;
typedef std::vector<RaftGenP> RaftGenPVec;
typedef RaftGenPVec::iterator RaftGenPVecI;
typedef RaftGenPVec::const_iterator RaftGenPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Raft genesis.
///
/// A RaftGen object describes a raft genesis event which takes place on a
/// surface system, i.e. a patch between two compartments.
///
/// \warning Methods start with an underscore are not exposed to Python.
class RaftGen {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the raft genesis.
    /// \param surfsys Pointer to the parent surface system.
    /// \param spec_signature The species 'signature' that must be met for a
    /// 			genesis to occur
    /// \param raft Pointer to the created raft
    /// \param kcst Rate constant

    RaftGen(std::string const& id,
            Surfsys* surfsys,
            std::vector<Spec*> const& spec_signature,
            Raft* raft,
            double kcst = 0.0);

    /// Destructor
    ~RaftGen();

    ////////////////////////////////////////////////////////////////////////
    // SURFACE REACTION RULE PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the raft genesis rule ID.
    ///
    /// \return ID of the raft genesis.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the raft genesis rule ID.
    ///
    /// \param id ID of the raft genesis.
    void setID(std::string const& id);

    /// Return a pointer to the parent surface system.
    ///
    /// \return Pointer to the surface system.
    inline Surfsys* getSurfsys() const noexcept {
        return pSurfsys;
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
    std::vector<Spec*> getAllSpecs() const;

    /// Get the rate constant of the raft genesis
    ///
    /// \return Rate constant of the raft genesis
    inline double getKcst() const noexcept {
        return pKcst;
    }

    /// Set the rate constant of the raft genesis
    ///
    /// \param Rate constant of the raft genesis
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
    Model* pModel;
    Surfsys* pSurfsys;

    std::vector<Spec*> pSpecSignature;
    Raft* pRaft;

    double pKcst;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::model
