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


#ifndef STEPS_TETEXACT_SREAC_HPP
#define STEPS_TETEXACT_SREAC_HPP 1

////////////////////////////////////////////////////////////////////////////////


// Standard library & STL headers.
#include <map>
#include <string>
#include <vector>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/constants.hpp"
#include "steps/solver/sreacdef.hpp"
#include "steps/tetexact/kproc.hpp"
//#include "tetexact.hpp"


////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetexact {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Tri;
class Tetexact;

////////////////////////////////////////////////////////////////////////////////

class SReac
: public steps::tetexact::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SReac(steps::solver::SReacdef * srdef, steps::tetexact::Tri * tri);
    ~SReac();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    double c() const
    { return pCcst; }
    void resetCcst();

    inline double kcst() const
    { return pKcst; }
    void setKcst(double k);

    double h()
    { return (rate()/pCcst); }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    void setupDeps();
    bool depSpecTet(uint gidx, steps::tetexact::WmVol * tet);
    bool depSpecTri(uint gidx, steps::tetexact::Tri * tri);
    void reset();
    double rate(steps::tetexact::Tetexact * solver = 0);
    std::vector<KProc*> const & apply(steps::rng::RNG * rng, double dt, double simtime);

    uint updVecSize() const
    { return pUpdVec.size(); }

    ////////////////////////////////////////////////////////////////////////

    //inline steps::solver::Reacdef * defr() const
    //{ return pReacdef; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::solver::SReacdef           * pSReacdef;
    steps::tetexact::Tri              * pTri;
    std::vector<KProc*>                 pUpdVec;
    /// Properly scaled reaction constant.
    double                              pCcst;
    // Store the kcst for convenience
    double                              pKcst;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_TETEXACT_SREAC_HPP
