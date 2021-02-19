/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

#ifndef STEPS_MODEL_SURFSYS_HPP
#define STEPS_MODEL_SURFSYS_HPP 1

// STL headers.
#include <cassert>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "steps/common.h"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace model {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Model;
class Spec;
class Surfsys;
class SReac;
class Chan;
class VDepTrans;
class VDepSReac;
class OhmicCurr;
class GHKcurr;
class Diff;

// Auxiliary declarations.
typedef Surfsys *                       SurfsysP;
typedef std::map<std::string, SurfsysP> SurfsysPMap;
typedef SurfsysPMap::iterator           SurfsysPMapI;
typedef SurfsysPMap::const_iterator     SurfsysPMapCI;

////////////////////////////////////////////////////////////////////////////////
/// Surface system.
/// Container that collects reactions involving a reactant
/// embedded in a membrane.
///
/// \warning Methods start with an underscore are not exposed to Python.

class Surfsys
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the surface system.
    /// \param model Pointer to the parent model.
    Surfsys(std::string const & id, Model * model);

    /// Destructor
    ~Surfsys();

    ////////////////////////////////////////////////////////////////////////
    // SURFACE SYSTEM PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the surface system ID.
    ///
    /// \return ID of the surface system.
    std::string getID() const
    { return pID; }
    /// Set or change the surface system ID.
    ///
    /// \param id ID of the surface system.
    void setID(std::string const & id);

    /// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent model.
    Model * getModel() const
    { return pModel; }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return a surface reaction with name id.
    ///
    /// \param id ID of the surface reaction.
    /// \return Pointer to the surface reaction.
    SReac * getSReac(std::string const & id) const;

    /// Delete a surace reaction with name id.
    ///
    /// \param id ID of the surface reaction.
    void delSReac(std::string const & id);

    /// Return a list of all surface reactions.
    ///
    /// \return List of pointers to surface reactions.
    std::vector<SReac *> getAllSReacs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Get a diffusion by its ID.
    ///
    /// \param id ID of the required diffusion.
    /// \return Pointer to the diffusion object.
    Diff * getDiff(std::string const & id) const;

    /// Delete a diffusion by its ID.
    ///
    /// \param id ID of the diffusion to be deleted.
    void delDiff(std::string const & id);

    /// Get all diffusions stored in this surface system.
    ///
    /// \return A vector of pointers to the diffusion objects
    ///         stored in the system.
    std::vector<Diff *> getAllDiffs() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a surface reaction id is occupied.
    ///
    /// \param id ID of the surface reaction.
    void _checkSReacID(std::string const & id) const;

    /// Change a surface reaction id from o to n.
    ///
    /// \param o Old id of the surface reaction.
    /// \param n New id of the surface reaction.
    void _handleSReacIDChange(std::string const & o, std::string const & n);

    /// Add a surface reaction to the surface system.
    ///
    /// \param Pointer to the surface reaction.
    void _handleSReacAdd(SReac * sreac);

    /// Delete a surface reaction in the surface system.
    ///
    /// \param Pointer to the surface reaction.
    void _handleSReacDel(SReac * sreac);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DIFFUSION
    ////////////////////////////////////////////////////////////////////////
    /// Check if a diffusion id is occupied.
    ///
    /// \param id ID of the diffusion.
    void _checkDiffID(std::string const & id) const;

    /// Change the id of a diffusion from o to n.
    ///
    /// \param o Old id of the diffusion.
    /// \param n New id of the diffusion.
    void _handleDiffIDChange(std::string const & o, std::string const & n);

    /// Add a diffusion to the surface system.
    ///
    /// \param diff Pointer to the diffusion.
    void _handleDiffAdd(Diff * diff);

    /// Delete a diffusion in the surface system.
    ///
    /// \param diff Pointer to the diffusion.
    void _handleDiffDel(Diff * diff);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): VOLTAGE-DEPENDENT TRANSITIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return a voltage-dependent transition with name id.
    ///
    /// \param id ID of the voltage-dependent transition.
    /// \return Pointer to the voltage-dependent transition.
    VDepTrans * getVDepTrans(std::string const & id) const;

    /// Delete a voltage-dependent transition with name id.
    ///
    /// \param id ID of the voltage-dependent transition.
    void delVDepTrans(std::string const & id);

    /// Return a list of all voltage-dependent transitions.
    ///
    /// \return List of pointers to voltage-dependent transitions.
    std::vector<VDepTrans *> getAllVDepTrans() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: VOLTAGE-DEPENDENT TRANSITIONS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a voltage-dependent transition id is occupied.
    ///
    /// \param id ID of the voltage-dependent transition.
    void _checkVDepTransID(std::string const & id) const;

    /// Change a voltage-dependent transition id from o to n.
    ///
    /// \param o Old id of the voltage-dependent transition.
    /// \param n New id of the voltage-dependent transition.
    void _handleVDepTransIDChange(std::string const & o, std::string const & n);

    /// Add a voltage-dependent transition to the surface system.
    ///
    /// \param Pointer to the voltage-dependent transition.
    void _handleVDepTransAdd(VDepTrans * vdeptrans);

    /// Delete a voltage-dependent transition in the surface system.
    ///
    /// \param Pointer to the voltage-dependent transition.
    void _handleVDepTransDel(VDepTrans * vdeptrans);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): VOLTAGE-DEPENDENT REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return a voltage-dependent reaction with name id.
    ///
    /// \param id ID of the voltage-dependent reaction.
    /// \return Pointer to the voltage-dependent reaction.
    VDepSReac * getVDepSReac(std::string const & id) const;

    /// Delete a voltage-dependent reaction with name id.
    ///
    /// \param id ID of the voltage-dependent reaction.
    void delVDepSReac(std::string const & id);

    /// Return a list of all voltage-dependent reactions.
    ///
    /// \return List of pointers to voltage-dependent transitions.
    std::vector<VDepSReac *> getAllVDepSReacs() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: VOLTAGE-DEPENDENT REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a voltage-dependent reaction id is occupied.
    ///
    /// \param id ID of the voltage-dependent reaction.
    void _checkVDepSReacID(std::string const & id) const;

    /// Change a voltage-dependent reaction id from o to n.
    ///
    /// \param o Old id of the voltage-dependent reaction.
    /// \param n New id of the voltage-dependent reaction.
    void _handleVDepSReacIDChange(std::string const & o, std::string const & n);

    /// Add a voltage-dependent reaction to the surface system.
    ///
    /// \param Pointer to the voltage-dependent reaction.
    void _handleVDepSReacAdd(VDepSReac * vdepsreac);

    /// Delete a voltage-dependent reaction in the surface system.
    ///
    /// \param Pointer to the voltage-dependent reaction.
    void _handleVDepSReacDel(VDepSReac * vdepsreac);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): OHMIC CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return an ohmic current with name id.
    ///
    /// \param id ID of the ohmic current.
    /// \return Pointer to the ohmic current.
    OhmicCurr * getOhmicCurr(std::string const & id) const;

    /// Delete an ohmic current with name id.
    ///
    /// \param id ID of the ohmic current.
    void delOhmicCurr(std::string const & id);

    /// Return a list of all ohmic currents.
    ///
    /// \return List of pointers to ohmic currents.
    std::vector<OhmicCurr *> getAllOhmicCurrs() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: OHMIC CURRENT
    ////////////////////////////////////////////////////////////////////////

    /// Check if an ohmic current id is occupied.
    ///
    /// \param id ID of the ohmic current.
    void _checkOhmicCurrID(std::string const & id) const;

    /// Change an ohmic current id from o to n.
    ///
    /// \param o Old id of the ohmic current.
    /// \param n New id of the ohmic current.
    void _handleOhmicCurrIDChange(std::string const & o, std::string const & n);

    /// Add a ohmic current to the surface system.
    ///
    /// \param Pointer to the ohmic current.
    void _handleOhmicCurrAdd(OhmicCurr * ohmiccurr);

    /// Delete an ohmic current in the surface system.
    ///
    /// \param Pointer to the ohmic current.
    void _handleOhmicCurrDel(OhmicCurr * ohmiccurr);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): GHK CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return an ghk current with name id.
    ///
    /// \param id ID of the ghk current.
    /// \return Pointer to the ghk current.
    GHKcurr * getGHKcurr(std::string const & id) const;

    /// Delete a ghk current with name id.
    ///
    /// \param id ID of the ghk current.
    void delGHKcurr(std::string const & id);

    /// Return a list of all ghk currents.
    ///
    /// \return List of pointers to ghk currents.
    std::vector<GHKcurr *> getAllGHKcurrs() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: GHK CURRENT
    ////////////////////////////////////////////////////////////////////////

    /// Check if a GHK current id is occupied.
    ///
    /// \param id ID of the GHK current.
    void _checkGHKcurrID(std::string const & id) const;

    /// Change a GHK current id from o to n.
    ///
    /// \param o Old id of the GHK current.
    /// \param n New id of the GHK current.
    void _handleGHKcurrIDChange(std::string const & o, std::string const & n);

    /// Add a GHK current to the surface system.
    ///
    /// \param Pointer to the GHK current.
    void _handleGHKcurrAdd(GHKcurr * ghkcurr);

    /// Delete a GHK current in the surface system.
    ///
    /// \param Pointer to the GHK current.
    void _handleGHKcurrDel(GHKcurr * GHKcurr);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return all species in the surface system.
    ///
    /// \return List of pointers to the species.
    std::vector<Spec *> getAllSpecs() const;

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
    inline uint _countSReacs() const
    { return pSReacs.size(); }

    /// Get a surface reaction with index lidx.
    ///
    /// \param lidx Index of the surface reaction.
    /// \return Pointer to the surface reaction.
    SReac * _getSReac(uint lidx) const;

    /// Count the voltage-dependent transitions in the surface system.
    ///
    /// \return Number of voltage-dependent transitions.
    inline uint _countVDepTrans() const
    { return pVDepTrans.size(); }

    /// Get a voltage-dependent transition with index lidx.
    ///
    /// \param lidx Index of the voltage-dependent transition.
    /// \return Pointer to the voltage-dependent transition.
    VDepTrans * _getVDepTrans(uint lidx) const;

    /// Count the voltage-dependent reactions in the surface system.
    ///
    /// \return Number of voltage-dependent reactions.
    inline uint _countVDepSReacs() const
    { return pVDepSReacs.size(); }

    /// Get a voltage-dependent reactions with index lidx.
    ///
    /// \param lidx Index of the voltage-dependent reaction.
    /// \return Pointer to the voltage-dependent reaction.
    VDepSReac * _getVDepSReac(uint lidx) const;

    /// Count the ohmic currents in the surface system.
    ///
    /// \return Number of ohmic currents.
    inline uint _countOhmicCurrs() const
    { return pOhmicCurrs.size(); }

    /// Get ohmic current object with index lidx.
    ///
    /// \param lidx index of the ohmic current.
    /// \ return Pointer to the ohmic current.
    OhmicCurr * _getOhmicCurr(uint lidx) const;

    /// Count the ghk currents in the surface system.
    ///
    /// \return Number of ghk currents.
    inline uint _countGHKcurrs() const
    { return pGHKcurrs.size(); }

    /// Get ghk current object with index lidx.
    ///
    /// \param lidx index of the ghk current.
    /// \ return Pointer to the ghk current.
    GHKcurr * _getGHKcurr(uint lidx) const;

    /// Get all surface reactions in the surface system.
    ///
    /// \return Map of surface reactions.
    const std::map<std::string, SReac *> & _getAllSReacs() const
    { return pSReacs; }

    /// Get all ohmic currents in the system
    ///
    /// \return Map of ohmic currents
    const std::map<std::string, OhmicCurr *> & _getAllOhmicCurrs() const
    { return pOhmicCurrs; }

    /// Get all ghk currents in the system
    ///
    /// \return Map of ghk currents
    const std::map<std::string, GHKcurr *> & _getAllGHKcurrs() const
    { return pGHKcurrs; }

    /// Get all voltage-dependent transitions in the system
    ///
    /// \return Map of voltage-dependent transitions
    const std::map<std::string, VDepTrans *> & _getAllVDepTrans() const
    { return pVDepTrans; }

    /// Get all voltage-dependent reactions in the system
    ///
    /// \return Map of voltage-dependent reactions
    const std::map<std::string, VDepSReac *> & _getAllVDepSReacs() const
    { return pVDepSReacs; }


    /// Count the diffusion rules in the volume system.
    ///
    /// \return Number of diffusion rules.
    inline uint _countDiffs() const
    { return pDiffs.size(); }

    /// Get a diffusion rule with index lidx.
    ///
    /// \param lidx Index of the diffusion rule.
    /// \return Pointer to the diffusion rule.
    Diff * _getDiff(uint lidx) const;

    /// Get all diffusion rules in the surface system.
    ///
    /// \return List of pointers to diffusion rules.
    const std::map<std::string, Diff *> & _getAllDiffs() const
    { return pDiffs; }

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
    ////////////////////////////////////////////////////////////////////////
    /// Delete a species in the surface system.
    ///
    /// \param spec Pointer to the species.
    void _handleSpecDelete(Spec * spec);

    /// Delete a channel in the surface system
    ///
    /// \param chan Pointer to the channel.
    void _handleChanDelete(Chan * chan);

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    std::string                         pID;
    Model                             * pModel;
    std::map<std::string, SReac *>      pSReacs;
    std::map<std::string, VDepTrans *>  pVDepTrans;
    std::map<std::string, OhmicCurr *>  pOhmicCurrs;
    std::map<std::string, GHKcurr *>    pGHKcurrs;
    std::map<std::string, VDepSReac *>  pVDepSReacs;

    std::map<std::string, Diff *>       pDiffs;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_MODEL_SURFSYS_HPP

// END
