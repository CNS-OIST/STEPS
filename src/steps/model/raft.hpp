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
#include <set>
#include <string>
#include <vector>

// STEPS headers.
#include "util/common.hpp"

namespace steps::model {

// Forward declarations.
class Model;
class Raft;

// Auxiliary declarations.
typedef Raft* RaftP;
typedef std::map<std::string, RaftP> RaftPMap;
typedef RaftPMap::iterator RaftPMapI;
typedef RaftPMap::const_iterator RaftPMapCI;

typedef std::vector<RaftP> RaftPVec;
typedef RaftPVec::iterator RaftPVecI;
typedef RaftPVec::const_iterator RaftPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Rafts .
/// Component that represents a membrane 'raft' such as a lipid raft or other
/// cluster. Rafts occupy space defined by a diameter projected onto the
/// surface, and can be mobile
///
/// \warning Methods start with an underscore are not exposed to Python.

class Raft {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the rafts.
    /// \param model Pointer to the parent model.
    /// \param diameter Effective diameter of the raft.
    /// \param dcst Diffusion coefficient of the raft.
    Raft(std::string const& id, Model* model, double diameter, double dcst = 0.0);

    /// Destructor
    virtual ~Raft();

    ////////////////////////////////////////////////////////////////////////
    // RAFTIES PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the raft's ID.
    ///
    /// \return ID of the raft.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the raft ID.
    ///
    /// \param id ID of the raft.
    void setID(std::string const& id);

    /// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent model.
    inline Model* getModel() const noexcept {
        return pModel;
    }

    /// Return the raft's diameter.
    ///
    /// \return Diameter of the raft.
    inline double getDiameter() const noexcept {
        return pDiameter;
    }

    /// Get the rate constant of the raft diffusion.
    ///
    /// \return Rate constant of the raft diffusion .
    inline double getDcst() const noexcept {
        return pDcst;
    }

    /// Set the rate constant of the raft diffusion.
    ///
    /// \param Rate constant of the raft diffusion.
    void setDcst(double dcst);

    /// Associate a raft system with name id.
    ///
    /// \param id ID of the raft system.
    void addRaftsys(std::string const& id);

    /// Get a raft surface system.
    ///
    /// \return List of the raft systems associated to the raft.
    inline const std::set<std::string>& getRaftsys() const noexcept {
        return pRaftsys;
    }

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

    // ...

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model* pModel;
    double pDiameter;

    double pDcst;

    std::set<std::string> pRaftsys;
};

}  // namespace steps::model
