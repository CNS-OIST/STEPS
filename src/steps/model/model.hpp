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

/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

#ifndef STEPS_MODEL_MODEL_HPP
#define STEPS_MODEL_MODEL_HPP 1


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
class Volsys;
class Reac;
class SReac;
class Diff;
class Chan;
class VDepTrans;
class VDepSReac;
class OhmicCurr;
class GHKcurr;

// Auxiliary declarations.

typedef Spec *                           SpecP;
typedef std::vector<SpecP>               SpecPVec;
typedef SpecPVec::iterator               SpecPVecI;
typedef SpecPVec::const_iterator         SpecPVecCI;

typedef Volsys *                         VolsysP;
typedef std::vector<VolsysP>               VolsysPVec;
typedef VolsysPVec::iterator             VolsysPVecI;
typedef VolsysPVec::const_iterator       VolsysPVecCI;

typedef Surfsys *                        SurfsysP;
typedef std::vector<SurfsysP>               SurfsysPVec;
typedef SurfsysPVec::iterator            SurfsysPVecI;
typedef SurfsysPVec::const_iterator      SurfsysPVecCI;

typedef Chan *                           ChanP;
typedef std::vector<ChanP>               ChanPVec;
typedef ChanPVec::iterator               ChanPVecI;
typedef ChanPVec::const_iterator         ChanPVecCI;

typedef std::map<std::string, VolsysP>  VolsysPMap;
typedef VolsysPMap::iterator            VolsysPMapI;
typedef VolsysPMap::const_iterator      VolsysPMapCI;

typedef std::map<std::string, SurfsysP> SurfsysPMap;
typedef SurfsysPMap::iterator           SurfsysPMapI;
typedef SurfsysPMap::const_iterator     SurfsysPMapCI;

////////////////////////////////////////////////////////////////////////////////


/// Top-level container for the objects in a kinetic model.
///
/// A steps::model::Model object is parent to the following objects:
/// <UL>
/// <LI>steps::model::Spec
/// <LI>steps::model::Volsys
/// <LI>steps::model::Surfsys
/// <LI>steps::model::Diff
/// </UL>
/// \sa Spec, Volsys, Surfsys, Diff.
/// \warning Methods start with an underscore are not exposed to Python.
///

class Model
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    Model(void);
    /// Destructor
    ~Model(void);

    // Model * deepcopy(void);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: SPECIES (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////
    /// Return a species with name id.
    ///
    /// \param id ID of the species.
    /// \return Pointer to the species.
    Spec * getSpec(std::string const & id) const;

    /// Delete a species with name id.
    ///
    /// \param id ID of the species.
    void delSpec(std::string const & id);

    /// Return a list of all species in the Model object.
    ///
    /// \return List of pointers to the species in the Model object.
    std::vector<Spec *> getAllSpecs(void) const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: CHANNELS (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////
    /// Return a channel with name id.
    ///
    /// \param id ID of the channel.
    /// \return Pointer to the channel.
    Chan * getChan(std::string const & id) const;

    /* Why overcomplicate things? May get rid of these for all classes actually
    /// Delete a channel with name id.
    ///
    /// \param id ID of the species.
    void delSpec(std::string const & id);
    */

    /// Return a list of all channels in the Model object.
    ///
    /// \return List of pointers to the channels in the Model object.
    std::vector<Chan *> getAllChans(void) const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: VOLSYS (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////

    /// Return a volume system with name id.
    ///
    /// \param id ID of the volume system.
    /// \return Pointer to the volume system.
    Volsys * getVolsys(std::string const & id) const;

    /// Delete a volume system with name id.
    ///
    /// \param id ID of the volume system.
    void delVolsys(std::string const & id);

    /// Return a list of all volume systems in the Model object.
    ///
    /// \return List of pointers to the volume systems in the Model object.
    std::vector<Volsys *> getAllVolsyss(void) const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: SURFSYS (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////

    /// Return a surface system with name id.
    ///
    /// \param id ID of the surface system.
    /// \return Pointer to the surface system.
    Surfsys * getSurfsys(std::string const & id) const;

    /// Delete a surface system with name id.
    ///
    /// \param id ID of the surface system.
    void delSurfsys(std::string const & id);

    /// Return a list of all surface systems in the Model object.
    ///
    /// \return List of pointers to the surface systems in the Model object.
    std::vector<Surfsys *> getAllSurfsyss(void) const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    /// Count the species in the Model object.
    ///
    /// \return Number of species.
    inline uint _countSpecs(void) const
    { return pSpecs.size(); }

    /// Return a species with index gidx.
    ///
    /// \param gidx Index of the species.
    /// \return Pointer to the species.
    Spec * _getSpec(uint gidx) const;

    /// Count the channels in the Model object.
    ///
    /// \return Number of channels.
    inline uint _countChans(void) const
    { return pChans.size(); }

    /// Return a channel with index gidx.
    ///
    /// \param gidx Index of the channel.
    /// \return Pointer to the channel.
    Chan * _getChan(uint gidx) const;

    /// Count the reactions in the Model object.
    ///
    /// \return Number of reactions.
    uint _countReacs(void) const;

    /// Return a reaction with index gidx
    ///
    /// \param gidx Index of the reaction.
    /// \return Pointer to the reaction.
    Reac * _getReac(uint gidx) const;

    /// Count the surface reactions in the Model object.
    ///
    /// \return Number of surface reactions.
    uint _countSReacs(void) const;

    /// Return a surface with index gidx.
    ///
    /// \param gidx Index of the surface reaction.
    /// \return Pointer to the surface reaction.
    SReac * _getSReac(uint gidx) const;

    /// Count the voltage-dependent transitions in the Model object.
    ///
    /// \return Number of voltage-dependent transitions.
    uint _countVDepTrans(void) const;

    /// Return a voltage-dependent reaction with index gidx.
    ///
    /// \param gidx Index of the voltage-dependent reaction.
    /// \return Pointer to the voltage-dependent reaction.
    VDepSReac * _getVDepSReac(uint gidx) const;

    /// Count the voltage-dependent reactions in the Model object.
    ///
    /// \return Number of voltage-dependent reactions.
    uint _countVDepSReacs(void) const;

    /// Return a voltage-dependent transition with index gidx.
    ///
    /// \param gidx Index of the voltage-dependent transition.
    /// \return Pointer to the voltage-dependent transition.
    VDepTrans * _getVDepTrans(uint gidx) const;

    /// Count the ohmic currents in the Model object.
    ///
    /// \return Number of ohmic currents.
    uint _countOhmicCurrs(void) const;

    /// Return an ohmic current with index gidx.
    ///
    /// \param gidx Index of the ohmic current.
    /// \return Pointer to the ohmic current.
    OhmicCurr * _getOhmicCurr(uint gidx) const;

    /// Count the ghk currents in the Model object.
    ///
    /// \return Number of ghk currents.
    uint _countGHKcurrs(void) const;

    /// Return an ghk current with index gidx.
    ///
    /// \param gidx Index of the ghk current.
    /// \return Pointer to the ghk current.
    GHKcurr * _getGHKcurr(uint gidx) const;

    /// Count the volume diffusions in the Model object.
    ///
    /// \return Number of volume diffusions.
    uint _countVDiffs(void) const;

    /// Return a volume diffusion with index gidx.
    ///
    /// \param gidx Index of the volume diffusion.
    /// \return Pointer to the volume diffusion.
    Diff * _getVDiff(uint gidx) const;

    /// Count the surface diffusions in the Model object.
    ///
    /// \return Number of surface diffusions.
    uint _countSDiffs(void) const;

    /// Return a surface diffusion with index gidx.
    ///
    /// \param gidx Index of the surface diffusion.
    /// \return Pointer to the surface diffusion.
    Diff * _getSDiff(uint gidx) const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a species id is occupied.
    ///
    /// \param id ID of the species.
    void _checkSpecID(std::string const & id) const;

    /// Change the id of a species from o to n.
    ///
    /// \param o Old id of the species.
    /// \param n New id of the species.
    void _handleSpecIDChange(std::string const & o, std::string const & n);

    /// Add a species to the Model.
    ///
    /// \param spec Pointer to the species being added.
    void _handleSpecAdd(Spec * spec);

    /// Delete a species in the Model.
    ///
    /// \param spec Pointer to the species being deleted.
    void _handleSpecDel(Spec * spec);


    /// Check if a channel id is occupied.
    ///
    /// \param id ID of the channel.
    void _checkChanID(std::string const & id) const;

    /// Change the id of a channel from o to n.
    ///
    /// \param o Old id of the channel.
    /// \param n New id of the channel.
    void _handleChanIDChange(std::string const & o, std::string const & n);

    /// Add a channel to the Model.
    ///
    /// \param spec Pointer to the channel being added.
    void _handleChanAdd(Chan * chan);

    /// Delete a channel in the Model.
    ///
    /// \param chan Pointer to the channel being deleted.
    void _handleChanDel(Chan * chan);

    /// Check if a volume system id is occupied.
    ///
    /// \param id ID of the volume system.
    void _checkVolsysID(std::string const & id) const;

    /// Change the id of a volume system from o to n.
    ///
    /// \param o Old id of the volume system.
    /// \param n New id of the volume system.
    void _handleVolsysIDChange(std::string const & o, std::string const & n);

    /// Add a volume system to the Model.
    ///
    /// \param volsys Pointer to the volume system being added.
    void _handleVolsysAdd(Volsys * volsys);

    /// Delete a volume system in the Model.
    ///
    /// \param volsys Pointer to the volume system being deleted.
    void _handleVolsysDel(Volsys * volsys);

    /// Check if a surface system id is occupied.
    ///
    /// \param id ID of the surface system.
    void _checkSurfsysID(std::string const & id) const;

    /// Change the id of a surface system from o to n.
    ///
    /// \param o Old id of the surface system.
    /// \param n New id of the surface system.
    void _handleSurfsysIDChange(std::string const & o, std::string const & n);

    /// Add a surface system to the Model.
    ///
    /// \param surfsys Pointer to the surface system being added.
    void _handleSurfsysAdd(Surfsys * surfsys);

    /// Delete a surface system in the Model.
    ///
    /// \param surfsys Pointer to the surface system being deleted.
    void _handleSurfsysDel(Surfsys * surfsys);

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    std::map<std::string, Spec *>       pSpecs;
    std::map<std::string, Chan *>       pChans;
    std::map<std::string, Volsys *>     pVolsys;
    std::map<std::string, Surfsys *>    pSurfsys;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif

// STEPS_MODEL_MODEL_HPP

// END
