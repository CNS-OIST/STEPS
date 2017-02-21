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
#include "steps/model/spec.hpp"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace steps::model;

////////////////////////////////////////////////////////////////////////////////

Spec::Spec(string const & id, Model * model, int valence)
: pID(id)
, pModel(model)
, pValence(valence)
{
    if (pModel == 0)
    {
        ostringstream os;
        os << "No model provided to Spec initializer function";
        throw steps::ArgErr(os.str());
    }
    pModel->_handleSpecAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

Spec::~Spec(void)
{
    if (pModel == 0) return;
    _handleSelfDelete();
}

////////////////////////////////////////////////////////////////////////////////

void Spec::_handleSelfDelete(void)
{
    pModel->_handleSpecDel(this);
    pModel = 0;
}

////////////////////////////////////////////////////////////////////////////////

void Spec::setID(string const & id)
{
    assert(pModel != 0);
    if (id == pID) return;
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pModel->_handleSpecIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

void Spec::setValence(int valence)
{
    assert(pModel != 0);
    pValence = valence;
}

////////////////////////////////////////////////////////////////////////////////

// END
