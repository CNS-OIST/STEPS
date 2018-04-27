/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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
 
#ifndef STEPS_MODEL_OHMICCURR_HPP
#define STEPS_MODEL_OHMICCURR_HPP 1

// STL headers.
#include <cassert>
#include <string>
#include <vector>
#include <map>

// STEPS headers.
#include "steps/common.h"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace model {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class OhmicCurr;
class Surfsys;
class Model;
class ChanState;

// Auxiliary declarations.
typedef OhmicCurr *                            OhmicCurrP;
typedef std::map<std::string, OhmicCurrP>   OhmicCurrPMap;
typedef OhmicCurrPMap::iterator             OhmicCurrPMapI;
typedef OhmicCurrPMap::const_iterator       OhmicCurrPMapCI;
typedef std::vector<OhmicCurrP>             OhmicCurrPVec;
typedef OhmicCurrPVec::iterator             OhmicCurrPVecI;
typedef OhmicCurrPVec::const_iterator       OhmicCurrPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Ohmic Current.
///
/// An OhmicCurr object describes an ohmic current through a channel.
/// Simply current approximation, which is described by channel conductance and
/// reversal potential, with no concentration dependence. Therefore this current
/// does not result in a real transfer of ions between compartments.
///
/// \warning Methods start with an underscore are not exposed to Python.

class OhmicCurr
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the ohmic current.
    /// \param surfsys Pointer to the parent surface system.
    /// \param chanstate The channel state which represents the open state.
    /// \param erev Reversal potential (in volts).
    /// \param g Single channel conductance (in siemens).
    ///
    OhmicCurr(std::string const & id, Surfsys * surfsys,
          ChanState * chanstate, double erev, double g);

    /// Destructor
    ~OhmicCurr(void);

    ////////////////////////////////////////////////////////////////////////
    // OHMIC CURRENT PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the ohmic current ID.
    ///
    /// \return ID of the ohmic current.
    std::string getID(void) const
    { return pID; }

    /// Set or change the ohmic current ID.
    ///
    /// \param id ID of the ohmic current.
    void setID(std::string const & id);

    /// Return a pointer to the parent surface system.
    ///
    /// \return Pointer to the surface system.
    Surfsys * getSurfsys(void) const
    { return pSurfsys; }

    /// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent model.
    Model * getModel(void) const
    { return pModel; }

    /// Return a pointer to the associated channel state.
    ///
    /// \return Pointer to the channel state.
    ChanState * getChanState(void) const
    { return pChanState; }

    /// Change the channel state.
    ///
    /// \param chanstate Channel state of the open state.
    void setChanState(ChanState * chanstate);

    /// Return the reversal potential (in volts).
    ///
    /// \return Reversal potential of the ohmic current.
    double getERev(void) const
    { return pERev; }

    /// Change the reversal potential.
    ///
    /// \param erev Reversal potential of the ohmic current.
    void setERev(double erev);

    /// Return the channel conductance (in siemens).
    ///
    /// \return Channel conductance associated with ohmic current.
    double getG(void) const
    { return pG; }

    /// Change the channel conductance.
    ///
    /// \param g Conductance associated with ohmic current.
    void setG(double g);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////
    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once
    void _handleSelfDelete(void);

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    std::string                         pID;
    Model                             * pModel;
    Surfsys                           * pSurfsys;

    ChanState                         * pChanState;
    double                              pERev;
    double                              pG;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_MODEL_OHMICCURR_HPP

// END
