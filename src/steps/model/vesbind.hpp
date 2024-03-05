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
#include <utility>
#include <vector>

#include "fwd.hpp"

namespace steps::model {

/// Vesicle binding reaction in a volume system.
class VesBind {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the reaction.
    /// \param volsys Pointer to the parent volume system.
    /// \param vesicle1 The 1st vesicle transporting the 1st reactant
    /// \param vesicle2 The 2nd vesicle transporting the 2nd reactant
    /// \param reactant1 The 1st reactant, located on Vesicle 1
    /// \param reactant2 The 2nd reactant, located on Vesicle 2
    /// \param link The product 'link' species
    /// \param vdeps1 Species dependencies on vesicle 1
    /// \param vdeps2 Species dependencies on vesicle 2
    /// \param ldeps1 Link species dependencies on vesicle 1
    /// \param ldeps2 Link species dependencies on vesicle 2
    /// \param kcst Rate constant for the reaction.
    /// \param immobilization A flag used to control the mobility effect of this
    ///			reaction. VesBind can immobilize vesicles.

    VesBind(std::string const& id,
            Volsys& volsys,
            Vesicle& vesicle1,
            Spec& reactant1,
            Vesicle& vesicle2,
            Spec& reactant2,
            LinkSpec& product1,
            LinkSpec& product2,
            double length_max,
            double length_min,
            std::vector<Spec*> const& vdeps1 = {},
            std::vector<Spec*> const& vdeps2 = {},
            std::vector<LinkSpec*> const& ldeps1 = {},
            std::vector<LinkSpec*> const& ldeps2 = {},
            double kcst = 0.0,
            Immobilization immobilization = NO_EFFECT);

    VesBind(const VesBind&) = delete;
    VesBind& operator=(const VesBind&) = delete;

    /// Destructor
    ~VesBind();

    ////////////////////////////////////////////////////////////////////////
    // REACTION RULE PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the vesbind ID.
    ///
    /// \return ID of the vesbind.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the vesbind  ID.
    ///
    /// \param id ID of the vesbind.
    void setID(std::string const& id);

    /// Return a reference to the parent volume system.
    ///
    /// \return Reference to the volume system.
    inline Volsys& getVolsys() const noexcept {
        return pVolsys;
    }

    /// Return a reference to the parent model.
    ///
    /// \return Reference to the parent Model.
    inline Model& getModel() const noexcept {
        return pModel;
    }

    /// Get the mobility effect of the vesbind.
    ///
    /// \return Immobilization of the vesbind.
    inline Immobilization getImmobilization() const noexcept {
        return pImmobilization;
    }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////
    /// Get the reactant pairs
    ///
    /// \return Map

    //  need to make sure I copy this data to def
    inline const std::pair<Vesicle*, Spec*>& getReactants1() const noexcept {
        return pReactants1;
    }

    inline const std::pair<Vesicle*, Spec*>& getReactants2() const noexcept {
        return pReactants2;
    }

    /// Get the complex product
    ///
    ///    \return Link Species products
    inline const std::pair<Vesicle*, LinkSpec*>& getProducts1() const noexcept {
        return pProducts1;
    }

    inline const std::pair<Vesicle*, LinkSpec*>& getProducts2() const noexcept {
        return pProducts2;
    }

    inline double getLengthMin() const noexcept {
        return pLength_min;
    }

    inline double getLengthMax() const noexcept {
        return pLength_max;
    }

    /// Return all species invloved in the reaction.
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

    /// Return a list of vesicle species dependencies for the reaction.
    ///
    /// \return List of pointers of dependencies.
    const std::vector<Spec*>& getVDeps1() const noexcept {
        return pVDeps1;
    }

    /// Return a list of vesicle species dependencies for the reaction.
    ///
    /// \return List of pointers of dependencies.
    const std::vector<Spec*>& getVDeps2() const noexcept {
        return pVDeps2;
    }

    /// Set the vesicle species dependencies for the reaction.
    ///
    /// \param vdeps Vesicle species dependencies for the reaction.
    void setVDeps1(std::vector<Spec*> const& vdeps);

    /// Set the vesicle species dependencies for the reaction.
    ///
    /// \param vdeps Vesicle species dependencies for the reaction.
    void setVDeps2(std::vector<Spec*> const& vdeps);

    /// Return a list of link species dependencies for the reaction.
    ///
    /// \return List of pointers of dependencies.
    const std::vector<LinkSpec*>& getLDeps1() const noexcept {
        return pLDeps1;
    }

    /// Return a list of link species dependencies for the reaction.
    ///
    /// \return List of pointers of dependencies.
    const std::vector<LinkSpec*>& getLDeps2() const noexcept {
        return pLDeps2;
    }

    /// Set the link species dependencies for the reaction.
    ///
    /// \param vdeps Link species dependencies for the reaction.
    void setLDeps1(std::vector<LinkSpec*> const& ldeps);

    /// Set the link species dependencies for the reaction.
    ///
    /// \param vdeps Link species dependencies for the reaction.
    void setLDeps2(std::vector<LinkSpec*> const& ldeps);

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////

    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// \warning Will only be called once
    void _handleSelfDelete();
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model& pModel;
    Volsys& pVolsys;

    std::pair<Vesicle*, Spec*> pReactants1;
    std::pair<Vesicle*, Spec*> pReactants2;

    std::pair<Vesicle*, LinkSpec*> pProducts1;
    std::pair<Vesicle*, LinkSpec*> pProducts2;

    double pLength_max;
    double pLength_min;

    std::vector<Spec*> pVDeps1;
    std::vector<Spec*> pVDeps2;
    std::vector<LinkSpec*> pLDeps1;
    std::vector<LinkSpec*> pLDeps2;

    double pKcst;

    Immobilization pImmobilization;
};

inline bool operator<(const VesBind& lhs, const VesBind& rhs) {
    return lhs.getID() < rhs.getID();
}

}  // namespace steps::model
