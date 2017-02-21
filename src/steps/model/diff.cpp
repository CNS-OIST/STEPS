/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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


// STL headers.
#include <cassert>
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/model/model.hpp"
#include "steps/model/volsys.hpp"
#include "steps/model/diff.hpp"
#include "steps/model/spec.hpp"

#include "steps/model/surfsys.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

////////////////////////////////////////////////////////////////////////////////

Diff::Diff(string const & id, Volsys * volsys, Spec * lig, double dcst)
: pID(id)
, pModel(0)
, pVolsys(volsys)
, pSurfsys(0)
, pLig(lig)
, pDcst(dcst)
, pIsvolume(true)
{
    if (pVolsys == 0)
    {
        ostringstream os;
        os << "No volsys provided to Diff initializer function.";
        throw steps::ArgErr(os.str());
    }
    if(pDcst < 0.0)
    {
        ostringstream os;
        os << "Diffusion constant can't be negative";
        throw steps::ArgErr(os.str());
    }
    pModel = pVolsys->getModel();
    assert (pModel != 0);

    pVolsys->_handleDiffAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Diff::Diff(string const & id, Surfsys * surfsys, Spec * lig, double dcst)
: pID(id)
, pModel(0)
, pSurfsys(surfsys)
, pVolsys(0)
, pLig(lig)
, pDcst(dcst)
, pIsvolume(false)
{
    if (pSurfsys == 0)
    {
        ostringstream os;
        os << "No surfsys provided to Diff initializer function.";
        throw steps::ArgErr(os.str());
    }
    if(pDcst < 0.0)
    {
        ostringstream os;
        os << "Diffusion constant can't be negative";
        throw steps::ArgErr(os.str());
    }
    pModel = pSurfsys->getModel();
    assert (pModel != 0);

    pSurfsys->_handleDiffAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Diff::~Diff(void)
{
    if (pIsvolume)
    {
        if (pVolsys == 0) return;
    }
    else
    {
        if (pSurfsys == 0) return;
    }
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Diff::_handleSelfDelete(void)
{
    if (pIsvolume)
    {
        pVolsys->_handleDiffDel(this);
        pVolsys = 0;
    }
    else
    {
        pSurfsys->_handleDiffDel(this);
        pSurfsys = 0;
    }
    pDcst = 0.0;
    pLig = 0;
    pModel = 0;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setID(string const & id)
{
    if (pIsvolume)
    {
        assert(pVolsys != 0);
        // The following might raise an exception, e.g. if the new ID is not
        // valid or not unique. If this happens, we don't catch but simply let
        // it pass by into the Python layer.
        pVolsys->_handleDiffIDChange(pID, id);
    }
    else
    {
        assert(pSurfsys != 0);
        // The following might raise an exception, e.g. if the new ID is not
        // valid or not unique. If this happens, we don't catch but simply let
        // it pass by into the Python layer.
        pSurfsys->_handleDiffIDChange(pID, id);
    }
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setDcst(double dcst)
{
    if (pIsvolume)
    {
        assert(pVolsys != 0);
    }
    else
    {
        assert(pSurfsys != 0);
    }
    if(dcst < 0.0)
    {
        ostringstream os;
        os << "Diffusion constant can't be negative";
        throw steps::ArgErr(os.str());
    }
    pDcst = dcst;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setLig(Spec * lig)
{
    assert(lig !=0);
    pLig = lig;
}

////////////////////////////////////////////////////////////////////////////////

vector<Spec *> Diff::getAllSpecs(void) const
{
    SpecPVec specs = SpecPVec();
    specs.push_back(pLig);
    return specs;
}

////////////////////////////////////////////////////////////////////////////////

// END
