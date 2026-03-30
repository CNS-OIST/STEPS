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

#include <string>

#include "complexevents.hpp"
#include "fwd.hpp"

namespace steps::model {

/// Ohmic Current base class.
///
/// An OhmicCurr object describes an ohmic current through a channel.
/// Simply current approximation, which is described by channel conductance and
/// reversal potential, with no concentration dependence. Therefore this current
/// does not result in a real transfer of ions between compartments.
///
/// \warning Methods start with an underscore are not exposed to Python.
class OhmicCurrBase {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the ohmic current.
    /// \param surfsys Reference to the parent surface system.
    /// \param erev Reversal potential (in volts).
    /// \param g Single channel conductance (in siemens).
    ///
    OhmicCurrBase(std::string const& id, Surfsys& surfsys, double erev, double g);

    OhmicCurrBase(const OhmicCurrBase&) = delete;
    OhmicCurrBase& operator=(const OhmicCurrBase&) = delete;

    ////////////////////////////////////////////////////////////////////////
    // OHMIC CURRENT PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the ohmic current ID.
    ///
    /// \return ID of the ohmic current.
    const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the ohmic current ID.
    ///
    /// \param id ID of the ohmic current.
    void setID(std::string const& id);

    /// Return a reference to the parent surface system.
    ///
    /// \return Reference to the surface system.
    inline Surfsys& getSurfsys() const noexcept {
        return pSurfsys;
    }

    /// Return a reference to the parent model.
    ///
    /// \return Reference to the parent model.
    inline Model& getModel() const noexcept {
        return pModel;
    }

    /// Return the reversal potential (in volts).
    ///
    /// \return Reversal potential of the ohmic current.
    inline double getERev() const noexcept {
        return pERev;
    }

    /// Change the reversal potential.
    ///
    /// \param erev Reversal potential of the ohmic current.
    void setERev(double erev);

    /// Return the channel conductance (in siemens).
    ///
    /// \return Channel conductance associated with ohmic current.
    inline double getG() const noexcept {
        return pG;
    }

    /// Change the channel conductance.
    ///
    /// \param g Conductance associated with ohmic current.
    void setG(double g);

    ////////////////////////////////////////////////////////////////////////

  protected:
    std::string pID;
    Model& pModel;
    Surfsys& pSurfsys;

    double pERev;
    double pG;
};

////////////////////////////////////////////////////////////////////////
/// Ohmic current with a species-like channel state
class OhmicCurr: public OhmicCurrBase {
  public:
    /// Constructor
    ///
    /// \param id ID of the ohmic current.
    /// \param surfsys Reference to the parent surface system.
    /// \param chanstate The channel state which represents the open state.
    /// \param erev Reversal potential (in volts).
    /// \param g Single channel conductance (in siemens).
    ///
    OhmicCurr(std::string const& id, Surfsys& surfsys, ChanState& chanstate, double erev, double g);

    ~OhmicCurr();

    ////////////////////////////////////////////////////////////////////////
    // OHMIC CURRENT PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return a reference to the associated channel state.
    ///
    /// \return Reference to the channel state.
    inline ChanState& getChanState() const noexcept {
        return *pChanState;
    }

    /// Change the channel state.
    ///
    /// \param chanstate Channel state of the open state.
    void setChanState(ChanState& chanstate);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////
    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once
    void _handleSelfDelete();

  private:
    ChanState* pChanState;
};

////////////////////////////////////////////////////////////////////////
/// Ohmic current with a multi-state complex channel state
class ComplexOhmicCurr: public OhmicCurrBase {
  public:
    /// Constructor
    ///
    /// \param id ID of the ohmic current.
    /// \param surfsys Reference to the parent surface system.
    /// \param cplx The name of the complex type involved in the current
    /// \param filt The complex filter that determines the conducting states
    /// \param erev Reversal potential (in volts).
    /// \param g Single channel conductance (in siemens).
    ///
    ComplexOhmicCurr(std::string const& id,
                     Surfsys& surfsys,
                     std::string const& cplx,
                     const std::vector<std::vector<SubunitStateFilter>>& filt,
                     double erev,
                     double g);

    ////////////////////////////////////////////////////////////////////////
    // OHMIC CURRENT PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return a reference to the associated channel state.
    ///
    /// \return Reference to the channel state.
    inline const ComplexFilterDescr& getChanState() const noexcept {
        return pChanState;
    }

  private:
    const ComplexFilterDescr pChanState;
};

}  // namespace steps::model
