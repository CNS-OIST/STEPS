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

#ifndef STEPS_MODEL_VDEPTRANS_HPP
#define STEPS_MODEL_VDEPTRANS_HPP 1

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
class VDepTrans;
class Surfsys;
class Model;
class ChanState;
class Chan;

// Auxiliary declarations.
typedef VDepTrans *                            VDepTransP;
typedef std::map<std::string, VDepTransP>   VDepTransPMap;
typedef VDepTransPMap::iterator             VDepTransPMapI;
typedef VDepTransPMap::const_iterator       VDepTransPMapCI;
typedef std::vector<VDepTransP>             VDepTransPVec;
typedef VDepTransPVec::iterator             VDepTransPVecI;
typedef VDepTransPVec::const_iterator       VDepTransPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Voltage-dependent transition between channel states.
///
///
/// \warning Methods start with an underscore are not exposed to Python.
class VDepTrans
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the voltage-dependant transition.
    /// \param surfsys Pointer to the parent surface system.
    /// \param src The 'source' state, the beginning state of the channel.
    /// \param dst The 'destination' state, the end state of the channel. species in the inner compartment
    /// \param rate A table of the voltage-dependent transition rate.
    ///

    VDepTrans(std::string const & id, Surfsys * surfsys,
             ChanState * src, ChanState * dst,
             std::vector<double> ratetab, double vmin, double vmax,
             double dv, uint tablesize);

    /// Destructor
    ~VDepTrans();

    ////////////////////////////////////////////////////////////////////////
    // VOLTAGE-DEPENDENT TRANSITION PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the voltage-dependent transition ID.
    ///
    /// \return ID of the voltage-dependent transition.
    inline std::string getID() const noexcept
    { return pID; }

    /// Set or change the voltage-dependent transition ID.
    ///
    /// \param id ID of the voltage-dependent transition.
    void setID(std::string const & id);

    /// Return a pointer to the parent surface system.
    ///
    /// \return Pointer to the surface system.
    inline Surfsys * getSurfsys() const noexcept
    { return pSurfsys; }

    /// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent model.
    inline Model * getModel() const noexcept
    { return pModel; }

    /// Return a pointer to the associated channel.
    ///
    /// \return Pointer to associated channel.
    inline Chan * getChan() const noexcept
    { return pChan; }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return a pointer to the 'source' channel state
    ///
    /// \return Pointer of the source channel state
    inline ChanState * getSrc() const noexcept
    { return pSrc; }

    /// Set the 'source' channel state
    ///
    /// \param src Pointer to the channel state
    void setSrc(ChanState * src);

    /// Return a pointer to the 'destination' channel state
    ///
    /// \return Pointer of the destination channel state
    inline ChanState * getDst() const noexcept
    { return pDst; }

    /// Set the 'destination' channel state
    ///
    /// \param src Pointer to the channel state
    void setDst(ChanState * dst);

    /// Get a table of transition rates over the range.
    ///
    std::vector<double> getRate() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER-HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    /// Get the table of transition rates.
    ///
    inline double * _getRate() const noexcept
    { return pRate; }

    inline double _getVMin() const noexcept
    { return pVMin; }

    inline double _getVMax() const noexcept
    { return pVMax; }

    inline double _getDV() const noexcept
    { return pDV; }

    inline uint _getTablesize() const noexcept
    { return pTablesize; }

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////

    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once
    void _handleSelfDelete();

private:

    ////////////////////////////////////////////////////////////////////////

    std::string                         pID;
    Model                             * pModel;
    Surfsys                           * pSurfsys;
    Chan                              * pChan;
    ChanState                         * pSrc;
    ChanState                         * pDst;
    double                            * pRate;

    double                                pVMin;
    double                                 pVMax;
    double                                 pDV;
    uint                                 pTablesize;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_MODEL_VDEPTRANS_HPP

// END
