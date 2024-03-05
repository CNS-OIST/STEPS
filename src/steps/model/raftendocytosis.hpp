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

////////////////////////////////////////////////////////////////////////////////
/// RaftEndocytosis.
///
/// A specific type of interaction that models endocytosis of a raft to form a
/// vesicle
///
/// \warning Methods start with an underscore are not exposed to Python.
class RaftEndocytosis {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the raft endocytotic reaction.
    /// \param raftsys Reference to the parent raft system.
    /// \param irhs Volume vesicle to appear in the 'inner' compartment
    /// \param orhs Volume vesicle to appear in the 'outer' compartment
    /// \param spec_deps Species dependency for endocytosis on the raft
    /// \param kcst Rate constant
    ///
    RaftEndocytosis(std::string const& id,
                    Raftsys& raftsys,
                    Vesicle* irhs = nullptr,
                    Vesicle* orhs = nullptr,
                    std::vector<Spec*> const& spec_deps = {},
                    double kcst = 0.0);

    RaftEndocytosis(const RaftEndocytosis&) = delete;
    RaftEndocytosis& operator=(const RaftEndocytosis&) = delete;

    /// Destructor
    ~RaftEndocytosis();

    ////////////////////////////////////////////////////////////////////////

    /// Return the raft endocytosis rule ID.
    ///
    /// \return ID of the raft endocytosis.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the raft endocytosis rule ID.
    ///
    /// \param id ID of the raft endocytosis.
    void setID(std::string const& id);

    /// Return a pointer to the parent raft system.
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

    /// Return whether the RaftEndocytosis rule resuts in a
    /// vesicle in the inner compartment (if false then the
    /// resulting vesicle will be in the outer compartment).
    ///
    /// \return Bool, true if 'inner'.
    inline bool getInner() const noexcept {
        return pInner;
    }

    /// Return a list of surface species dependencies.
    ///
    /// \return List of pointers of surface species.
    inline const std::vector<Spec*>& getSpecDeps() const noexcept {
        return pDepSurface;
    }

    /// Return volume vesicle on the right hand side of reaction.
    ///
    /// \return Pointer to right hand side volume species.
    inline Vesicle& getRHS() const noexcept {
        return *pVesicle;
    }

    /// Get the rate constant of the raftendocytotic reaction.
    ///
    /// \return Rate constant of the raftendocytotic reaction.
    inline double getKcst() const noexcept {
        return pKcst;
    }

    /// Set the rate constant of the raftendocytotic reaction.
    ///
    /// \param Rate constant of the raftendocytotic reaction.
    void setKcst(double kcst);

    /// Get all the species involved in the raftendocytotic reaction.
    ///
    /// \return Rate constant of the raftendocytotic reaction.
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

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model& pModel;
    Raftsys& pRaftsys;

    Vesicle* pVesicle;

    std::vector<Spec*> pDepSurface;

    double pKcst;

    bool pInner;
};

}  // namespace steps::model
