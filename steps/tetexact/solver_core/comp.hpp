////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_TETEXACT_SOLVER_CORE_COMP_HPP
#define STEPS_TETEXACT_SOLVER_CORE_COMP_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
#include <steps/tetexact/solver_core/tet.hpp>

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Comp;

// Auxiliary declarations.
typedef Comp *                          CompP;
typedef std::vector<CompP>              CompPVec;
typedef CompPVec::iterator              CompPVecI;
typedef CompPVec::const_iterator        CompPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Currently, this class is only a lightweight container for grouping 
/// a number of tetrahedrons. There is no internal ordering scheme. Tet
/// objects themselves do not at this point refer back to objects of this 
/// class.
///
/// It also stores the compartment volume as computed by summing the
/// tetrahedra. 
///
/// Later, when we transform the implementation of the Python-to-C-solver
/// interface into a Boost::Python/template based system, this class might
/// become the place in which certain compartment-related concepts are
/// implemented (e.g. getting/setting the concentration in a whole
/// compartment etc). When that happens, this class will probably become
/// related to CompDef through some inheritance mechanism, removing the
/// double storage (and potentially double meaning) of the comparment's 
/// volume.

class Comp
{
    
public:
    
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    
    Comp(steps::sim::CompDef * compdef);
    ~Comp(void);
    
    /// Checks whether the Tet's compdef() corresponds to this object's
    /// CompDef. There is no check whether the Tet object has already
    /// been added to this Comp object before (i.e. no duplicate checking).
    ///
    void addTet(Tet * tet);
    
    /// (Re-)computes the volume of the compartment by adding the 
    /// tetrahedral volumes.
    ///
    void computeVol(void);
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////
    
    inline steps::sim::CompDef * def(void) const
    { return pCompDef; }
    
    inline double vol(void) const
    { return pVol; }
    
    inline uint countTets(void) const
    { return pTets.size(); }
    
    /// If necessary, might be speeded up by pre=computing a rough 
    /// partitioning in Comp::computeVol().
    ///
    Tet * pickTetByVol(double rand01) const;
    
    inline TetPVecCI bgnTet(void) const
    { return pTets.begin(); }
    inline TetPVecCI endTet(void) const
    { return pTets.end(); }
    
    ////////////////////////////////////////////////////////////////////////
    
private:
    
    ////////////////////////////////////////////////////////////////////////
    
    steps::sim::CompDef *       pCompDef;
    double                      pVol;
    TetPVec                     pTets;

    ////////////////////////////////////////////////////////////////////////
    
    // Disable constructors.
    Comp(void);
    Comp(Comp const &);
    
    ////////////////////////////////////////////////////////////////////////
    
};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_TETEXACT_SOLVER_CORE_COMP_HPP

// END
