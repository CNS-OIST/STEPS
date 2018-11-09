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

#ifndef STEPS_MODEL_SREAC_HPP
#define STEPS_MODEL_SREAC_HPP 1


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
class SReac;
class Surfsys;
class Model;
class Spec;

// Auxiliary declarations.
typedef SReac *                            SReacP;
typedef std::map<std::string, SReacP>   SReacPMap;
typedef SReacPMap::iterator             SReacPMapI;
typedef SReacPMap::const_iterator       SReacPMapCI;
typedef std::vector<SReacP>             SReacPVec;
typedef SReacPVec::iterator             SReacPVecI;
typedef SReacPVec::const_iterator       SReacPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Surface reaction.
///
/// A SReac object describes a reaction which takes place on a surface system,
/// i.e. a patch between two compartments.
///
/// \warning Methods start with an underscore are not exposed to Python.
class SReac
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the surface reaction.
    /// \param surfsys Pointer to the parent surface system.
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
    /// \param kcst Rate constant of the reaction.
    ///
    /// \warning By default, the vlhs are defined in the outer compartment.
    ///          call setInner and SetOuter to change this default setting.
    /// \sa setInner, setOuter.
    SReac(std::string const & id, Surfsys * surfsys,
          std::vector<Spec *> const & olhs = std::vector<Spec *>(),
          std::vector<Spec *> const & ilhs = std::vector<Spec *>(),
          std::vector<Spec *> const & slhs = std::vector<Spec *>(),
          std::vector<Spec *> const & irhs = std::vector<Spec *>(),
          std::vector<Spec *> const & srhs = std::vector<Spec *>(),
          std::vector<Spec *> const & orhs = std::vector<Spec *>(),
          double kcst = 0.0);

    /// Destructor
    ~SReac();

    ////////////////////////////////////////////////////////////////////////
    // SURFACE REACTION RULE PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the surface reaction rule ID.
    ///
    /// \return ID of the surface reaction.
    std::string getID() const
    { return pID; }

    /// Set or change the surface reaction rule ID.
    ///
    /// \param id ID of the surface reaction.
    void setID(std::string const & id);

    /// Return a pointer to the parent surface system.
    ///
    /// \return Pointer to the surface system.
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

    /// Check if the lhs invloves species in the inner compartment.
    ///
    /// \return True if ilhs is set.
    ///         False if else.
    bool getInner() const
    { return (! pOuter); }


    /// Check if the lhs involves species in the outer compartment,
    /// or there are no volume species on the lhs.
    ///
    /// \return True if olhs is set, or neither olhs or ilhs are set.
    ///         False if else.
    bool getOuter() const
    { return pOuter; }


    /// Return a list of outer volume species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side outer volume species.
    const std::vector<Spec *> & getOLHS() const
    { return pOLHS; }

    /// Set the outer volume species on the left hand side of reaction.
    ///
    /// \param olhs Outer volume species on the left hand side of reaction.
    void setOLHS(std::vector<Spec *> const & olhs);

    /// Return a list of inner volume species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side inner volume species.
    const std::vector<Spec *> & getILHS() const
    { return pILHS; }

    /// Set the inner volume species on the left hand side of reaction.
    ///
    /// \param ilhs Inner volume species on the left hand side of reaction.
    void setILHS(std::vector<Spec *> const & ilhs);

    /// Return a list of surface species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side surface species.
    const std::vector<Spec *> & getSLHS() const
    { return pSLHS; }

    /// Set the surface species on the left hand side of reaction.
    ///
    /// \param slhs Surface species on the left hand side of reaction.
    void setSLHS(std::vector<Spec *> const & slhs);

    /// Return a list of inner volume species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side inner volume species.
    const std::vector<Spec *> & getIRHS() const
    { return pIRHS; }

    /// Set the inner volume species on the right hand side of reaction.
    ///
    /// \param irhs Inner volume species on the right hand side of reaction.
    void setIRHS(std::vector<Spec *> const & irhs);

    /// Return a list of surface species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side surface species.
    const std::vector<Spec *> & getSRHS() const
    { return pSRHS; }

    /// Set the surface species on the right hand side of reaction.
    ///
    /// \param srhs Surface species on the right hand side of reaction.
    void setSRHS(std::vector<Spec *> const & srhs);

    /// Return a list of outer volume species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side outer volume species.
    const std::vector<Spec *> & getORHS() const
    { return pORHS; }

    /// Set the outer volume species on the right hand side of reaction.
    ///
    /// \param orhs Outer volume species on the right hand side of reaction.
    void setORHS(std::vector<Spec *> const & orhs);

    /// Get the order of the surface reaction.
    ///
    /// \return Order of the reaction.
    uint getOrder() const
    { return pOrder; }

    /// Get the rate constant of the surface reaction.
    ///
    /// \return Rate constant of the surface reaction.
    double getKcst() const
    { return pKcst; }

    /// Set the rate constant of the surface reaction.
    ///
    /// \param Rate constant of the surface reaction.
    void setKcst(double kcst);

    /// Get a list of all species.
    ///
    /// Returns a list of all species involved in this
    /// surface reaction, on both the left and righthand side
    /// and does not contain any duplicate members.
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

private:

    ////////////////////////////////////////////////////////////////////////

    std::string                         pID;
    Model                             * pModel;
    Surfsys                           * pSurfsys;

    bool                                pOuter;
    std::vector<Spec *>                 pOLHS;
    std::vector<Spec *>                 pILHS;
    std::vector<Spec *>                 pSLHS;
    std::vector<Spec *>                 pIRHS;
    std::vector<Spec *>                 pSRHS;
    std::vector<Spec *>                 pORHS;
    uint                                pOrder;
    double                              pKcst;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_MODEL_SREAC_HPP

// END
