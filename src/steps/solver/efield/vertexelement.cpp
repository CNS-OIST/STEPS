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


// STL headers.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/efield/vertexconnection.hpp"
#include "steps/solver/efield/vertexelement.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace sefield = steps::solver::efield;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

sefield::VertexElement::VertexElement(uint idx, double * vpos)
: pIDX(idx)
, pXPos(vpos[0])
, pYPos(vpos[1])
, pZPos(vpos[2])
, pSurface(0.0)
, pVolume(0.0)
, pCapacitance(0.0)
, pConnections()
, pNCon(0)
, pNbrs(0)
, pCcs(0)
{
}

////////////////////////////////////////////////////////////////////////////////

sefield::VertexElement::~VertexElement(void)
{
    delete[] pNbrs;
    delete[] pCcs;
}

////////////////////////////////////////////////////////////////////////////////

void sefield::VertexElement::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&pSurface, sizeof(double));
    cp_file.write((char*)&pVolume, sizeof(double));
    cp_file.write((char*)&pCapacitance, sizeof(double));
    cp_file.write((char*)&pNCon, sizeof(uint));
    cp_file.write((char*)pCcs, sizeof(double) * pNCon);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::VertexElement::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&pSurface, sizeof(double));
    cp_file.read((char*)&pVolume, sizeof(double));
    cp_file.read((char*)&pCapacitance, sizeof(double));
    cp_file.read((char*)&pNCon, sizeof(uint));
    cp_file.read((char*)pCcs, sizeof(double) * pNCon);
}

////////////////////////////////////////////////////////////////////////////////

void sefield::VertexElement::fix(void)
{
    pNCon = pConnections.size();
    pNbrs = new VertexElement*[pNCon];
    pCcs = new double[pNCon];

    for (int i = 0; i < pNCon; ++i)
    {
        pNbrs[i] = pConnections[i]->getOther(this);
        pCcs[i] = 0.0;
    }
}

////////////////////////////////////////////////////////////////////////////////

void sefield::VertexElement::applyConductance(double a)
{
    // this has some effect on compilation/optimization: without it,
    // the coupling constants are wrong
    // Iain : what on earth was the following line doing in here?
    // double* uu = new double[pNCon];

    for (int i = 0; i < pNCon; ++i)
    {
        pCcs[i] =  a * pConnections[i]->getGeomCouplingConstant();
    }
}

////////////////////////////////////////////////////////////////////////////////
/*
ostream & operator<< (ostream & os, sefield::VertexElement const & ve)
{
    os << "VertexElement(idx=#" << ve.getIDX() << ", ";
    os << "x=" << ve.getX() << ", ";
    os << "y=" << ve.getY() << ", ";
    os << "z=" << ve.getZ() << ", ";
    os << "ncon=" << ve.getNCon() << ", ";
    os << "con={";
    for (uint i = 0; i < ve.getNCon(); ++i)
    {
        os << "#" << ve.pNbrs[i]->getIDX() << ":";
        os << ve.pCcs[i] << ",";
    }
    os << "}, ";
    os << "cap=" << ve.getCapacitance() << ")";
    return os;
}
*/
////////////////////////////////////////////////////////////////////////////////

// END
