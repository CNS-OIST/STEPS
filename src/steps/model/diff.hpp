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

#ifndef STEPS_MODEL_DIFF_HPP
#define STEPS_MODEL_DIFF_HPP 1


// STL headers.
#include <cassert>
#include <string>
#include <map>
#include <vector>

// STEPS headers.
#include "steps/common.h"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace model {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Diff;
class Volsys;
class Surfsys;
class Model;
class Spec;

// Auxiliary declarations.
typedef Diff *                            DiffP;
typedef std::map<std::string, DiffP>    DiffPMap;
typedef DiffPMap::iterator              DiffPMapI;
typedef DiffPMap::const_iterator        DiffPMapCI;
typedef std::vector<DiffP>              DiffPVec;
typedef DiffPVec::iterator              DiffPVecI;
typedef DiffPVec::const_iterator        DiffPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Diffusion rule in a volume system.
///
///\warning Methods start with an underscore are not exposed to Python.

class Diff
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the diffusion rule.
    /// \param volsys Volume system which the diffusion rule belongs to.
    /// \param lig Pointers to the species which the diffusion applies to.
    /// \param dcst Diffusion coefficient of the diffusion rule.
    Diff(std::string const & id, Volsys * volsys, Spec * lig, double dcst=0.0);

    /// Constructor
    ///
    /// \param id ID of the diffusion rule.
    /// \param surfsys Surface system which the diffusion rule belongs to.
    /// \param lig Pointers to the species which the diffusion applies to.
    /// \param dcst Diffusion coefficient of the diffusion rule.
    Diff(std::string const & id, Surfsys * surfsys, Spec * lig, double dcst=0.0);

    /// Destructor
    ~Diff();

    ////////////////////////////////////////////////////////////////////////
    // DIFFUSION RULE PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the diffusion rule ID.
    ///
    /// \return ID of the diffusion rule.
    std::string getID() const
    { return pID; }

    /// Set the ID of the diffusion rule.
    ///
    /// \param id ID of the diffusion rule.
    void setID(std::string const & id);

    /// Return a pointer to the parent volume system.
    ///
    /// \return Pointer to the parent volume system.
    Volsys * getVolsys() const
    { return pVolsys; }

    /// Return a pointer to the parent surface system.
    ///
    /// \return Pointer to the parent surface system.
    Surfsys * getSurfsys() const
    { return pSurfsys; }

    /// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent model.
    Model * getModel() const
    { return pModel; }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return a pointer to the species to which this diffusion rule applies
    ///
    /// \return Pointer of the species
    Spec * getLig() const
    { return pLig; }

    /// Set the species which this difusion rule applies to.
    ///
    /// \param lig Pointer to the species
    void setLig(Spec * lig);

    /// Get the rate constant of the diffusion rule.
    ///
    /// \return Rate constant of the diffusion rule.
    double getDcst() const
    { return pDcst; }

    /// Set the rate constant of the diffusion rule.
    ///
    /// \param Rate constant of the diffusion rule.
    void setDcst(double dcst);

    ///  Return a list of all species in this diffusion rule.
    ///
    /// \return List of pointers of species.
    /// \warning Currently will return only one species.
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

private:

    ////////////////////////////////////////////////////////////////////////

    std::string                         pID;
    Model                             * pModel;

    Volsys                            * pVolsys;
    Surfsys                              * pSurfsys;

    Spec                              * pLig;
    double                              pDcst;

    bool                                 pIsvolume;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_MODEL_DIFF_HPP

// END
