/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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
#include "model/complexevents.hpp"
#include "util/collections.hpp"
#include "util/common.hpp"

namespace steps::model {

////////////////////////////////////////////////////////////////////////////////
/// Complex surface reaction base class.
///
/// A ComplexSReac object describes a reaction which takes place on a surface system,
/// i.e. a patch between two compartments and involves complexes.
///
/// \warning Methods start with an underscore are not exposed to Python.
class ComplexSReacBase {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the surface reaction.
    /// \param surfsys Reference to the parent surface system.
    /// \param olhs Volume species in the outer compartment
    ///                on the left hand side of the reaction.
    /// \param ilhs Volume species in the inner compartment
    ///             and on the left hand side of the reaction.
    /// \param slhs Surface species on the left hand side of the reaction.
    /// \param irhs Volume species in the inner compartment
    ///             and on the right hand side of the reaction.
    /// \param srhs Surface species on the right hand side of the reaction.
    /// \param orhs Volume species in the outer compartment
    ///             and on the right hand side of the reaction.
    /// \param icompEvs Complex events in the inner compartment
    /// \param scompEvs Complex events on the surface
    /// \param ocompEvs Complex events in the outer compartment
    /// \param charge Charge carried by the reaction across an efield membrane.
    ComplexSReacBase(std::string const& id,
                     Surfsys& surfsys,
                     std::vector<Spec*> const& ilhs = {},
                     std::vector<Spec*> const& slhs = {},
                     std::vector<Spec*> const& olhs = {},
                     std::vector<Spec*> const& irhs = {},
                     std::vector<Spec*> const& srhs = {},
                     std::vector<Spec*> const& orhs = {},
                     std::vector<ComplexEvent*> const& icompEvs = {},
                     std::vector<ComplexEvent*> const& scompEvs = {},
                     std::vector<ComplexEvent*> const& ocompEvs = {},
                     int charge = 0);

    ~ComplexSReacBase() = default;

    ////////////////////////////////////////////////////////////////////////
    // REACTION RULE PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    const std::string& getID() const noexcept {
        return pID;
    }

    Surfsys& getSurfsys() const noexcept {
        return pSurfsys;
    }

    Model& getModel() const noexcept {
        return pModel;
    }

    const std::vector<ComplexUpdateEvent*>& getUPDEvents(Location loc) const noexcept;

    const std::vector<ComplexDeleteEvent*>& getDELEvents(Location loc) const noexcept;

    const std::vector<ComplexCreateEvent*>& getCREEvents(Location loc) const noexcept;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Check if the lhs involves species in the inner compartment.
    ///
    /// \return True if ilhs is set.
    ///         False if else.
    bool getInner() const noexcept {
        return !pOuter;
    }

    /// Check if the lhs involves species in the outer compartment,
    /// or there are no volume species on the lhs.
    ///
    /// \return True if olhs is set, or neither olhs or ilhs are set.
    ///         False if else.
    bool getOuter() const noexcept {
        return pOuter;
    }

    /// Check if the lhs involves only species on the patch.
    bool getSurfSurf() const noexcept {
        return pSurfSurf;
    }

    /// Return a list of outer volume species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side outer volume species.
    const std::vector<Spec*>& getOLHS() const noexcept {
        return pOLHS;
    }

    /// Return a list of inner volume species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side inner volume species.
    const std::vector<Spec*>& getILHS() const noexcept {
        return pILHS;
    }

    /// Return a list of surface species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side surface species.
    const std::vector<Spec*>& getSLHS() const noexcept {
        return pSLHS;
    }

    /// Return a list of inner volume species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side inner volume species.
    const std::vector<Spec*>& getIRHS() const noexcept {
        return pIRHS;
    }

    /// Return a list of surface species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side surface species.
    const std::vector<Spec*>& getSRHS() const noexcept {
        return pSRHS;
    }

    /// Return a list of outer volume species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side outer volume species.
    const std::vector<Spec*>& getORHS() const noexcept {
        return pORHS;
    }

    /// Get the charge carried by the reaction across the membrane.
    ///
    /// \return Charge carried by the reaction.
    inline int getCharge() const noexcept {
        return pCharge;
    }

    /// Set the charge carried by the reaction across the membrane.
    ///
    /// \param Charge carried by the reaction.
    void setCharge(int charge);

    /// Get the order of the surface reaction.
    ///
    /// \return Order of the reaction.
    uint getOrder() const noexcept {
        return pOrder;
    }

    /// Get a list of all species.
    ///
    /// Returns a list of all species involved in this
    /// surface reaction, on both the left and righthand side
    /// and does not contain any duplicate members.
    util::flat_set<Spec*> getAllSpecs() const;

  protected:
    void _addEvent(ComplexEvent* ev, Location loc);

    const std::string pID;
    Model& pModel;
    Surfsys& pSurfsys;

    bool pOuter;
    bool pSurfSurf;
    const std::vector<Spec*> pILHS;
    const std::vector<Spec*> pSLHS;
    const std::vector<Spec*> pOLHS;
    const std::vector<Spec*> pIRHS;
    const std::vector<Spec*> pSRHS;
    const std::vector<Spec*> pORHS;
    std::map<Location, std::vector<ComplexUpdateEvent*>> pCompUPD;
    std::map<Location, std::vector<ComplexDeleteEvent*>> pCompDEL;
    std::map<Location, std::vector<ComplexCreateEvent*>> pCompCRE;
    std::map<Location, uint> pLocOrder;
    uint pOrder{};
    int pCharge;
};

inline bool operator<(const ComplexSReacBase& lhs, const ComplexSReacBase& rhs) {
    return lhs.getID() < rhs.getID();
}

////////////////////////////////////////////////////////////////////////////////
/// Complex surface reaction with rate constant
class ComplexSReac: public ComplexSReacBase {
  public:
    ComplexSReac(std::string const& id,
                 Surfsys& surfsys,
                 std::vector<Spec*> const& ilhs = {},
                 std::vector<Spec*> const& slhs = {},
                 std::vector<Spec*> const& olhs = {},
                 std::vector<Spec*> const& irhs = {},
                 std::vector<Spec*> const& srhs = {},
                 std::vector<Spec*> const& orhs = {},
                 std::vector<ComplexEvent*> const& icompEvs = {},
                 std::vector<ComplexEvent*> const& scompEvs = {},
                 std::vector<ComplexEvent*> const& ocompEvs = {},
                 double kcst = 0.0,
                 int charge = 0);

    /// Get the rate constant of the surface reaction.
    ///
    /// \return Rate constant of the surface reaction.
    double getKcst() const noexcept {
        return pKcst;
    }

    /// Set the rate constant of the surface reaction.
    ///
    /// \param Rate constant of the surface reaction.
    void setKcst(double kcst);

  private:
    double pKcst{};
};

////////////////////////////////////////////////////////////////////////////////
/// Complex surface reaction with voltage-dependent rate
class VDepComplexSReac: public ComplexSReacBase {
  public:
    VDepComplexSReac(std::string const& id,
                     Surfsys& surfsys,
                     std::vector<Spec*> const& ilhs = {},
                     std::vector<Spec*> const& slhs = {},
                     std::vector<Spec*> const& olhs = {},
                     std::vector<Spec*> const& irhs = {},
                     std::vector<Spec*> const& srhs = {},
                     std::vector<Spec*> const& orhs = {},
                     std::vector<ComplexEvent*> const& icompEvs = {},
                     std::vector<ComplexEvent*> const& scompEvs = {},
                     std::vector<ComplexEvent*> const& ocompEvs = {},
                     std::vector<double> ktab = {},
                     double vmin = 0.0,
                     double vmax = 0.0,
                     double dv = 0.0,
                     uint tablesize = 0,
                     int charge = 0);

    /// Get the table of transition rates.
    ///
    inline const auto& getK() const noexcept {
        return pK;
    }

    inline const auto& _getK() const noexcept {
        return pK;
    }

    inline double _getVMin() const noexcept {
        return pVMin;
    }

    inline double _getVMax() const noexcept {
        return pVMax;
    }

    inline double _getDV() const noexcept {
        return pDV;
    }

    inline uint _getTablesize() const noexcept {
        return pTablesize;
    }

  private:
    std::vector<double> pK;

    double pVMin;
    double pVMax;
    double pDV;
    uint pTablesize;
};


}  // namespace steps::model
