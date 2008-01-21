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
// $Id:sched.hpp 64 2007-08-20 06:25:41Z stefan $
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_TETEXACT_SOLVER_CORE_SCHED_HPP
#define STEPS_TETEXACT_SOLVER_CORE_SCHED_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <map>
#include <set>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class KProc;
class State;
class Tet;
class Tri;

// Auxiliary declarations.
typedef uint                            SchedIDX;
typedef std::set<SchedIDX>              SchedIDXSet;
typedef SchedIDXSet::iterator           SchedIDXSetI;
typedef SchedIDXSet::const_iterator     SchedIDXSetCI;
typedef std::vector<SchedIDX>           SchedIDXVec;
typedef SchedIDXVec::iterator           SchedIDXVecI;
typedef SchedIDXVec::const_iterator     SchedIDXVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Copies the contents of a set of SchedIDX entries into a vector.
/// The contents of the vector are completely overridden.
///
extern void schedIDXSet_To_Vec(SchedIDXSet const & s, SchedIDXVec & v);

////////////////////////////////////////////////////////////////////////////////

class Sched
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Sched(void);
    
    ~Sched(void);
    
    void addKProc(KProc * kproc);
    
    void build();
    
    ////////////////////////////////////////////////////////////////////////
    // SCHEDULE OPERATIONS
    ////////////////////////////////////////////////////////////////////////
    
    inline double getA0(void) const
    { return pA0; }
    
    KProc * getNext(State * s) const;
    
    /// Asks all kproc's in the schedule to recompute their rates, and
    /// then recomputes the schedule tree.
    ///
    void reset(void);
    
    /// Same as Sched::reset(), but doesn't ask kproc's to recompute
    /// their rates.
    ///
    void recomp(void);
    
    void update(SchedIDXVec const & entries);
    
    /// Update the kproc's of a tet, after a species has been changed.
    /// This also updates kproc's in surrounding triangles.
    ///
    /// Currently doesn't care about the species.
    ///
    void updateSpec(Tet * tet, uint spec_lidx);
    
    /// Update the kproc's of a triangle, after a species has been changed.
    /// This does not need to update the kproc's of any neighbouring
    /// tetrahedrons.
    ///
    /// Currently doesn't care about the species.
    ///
    void updateSpec(Tri * tri, uint spec_lidx);

private:

    ////////////////////////////////////////////////////////////////////////
    // LIST OF KPROCS
    ////////////////////////////////////////////////////////////////////////
    
    std::vector<KProc *>        pKProcs;
    
    ////////////////////////////////////////////////////////////////////////
    // N-ARY TREE
    ////////////////////////////////////////////////////////////////////////
    
    double                      pA0;
    
    std::vector<uint>           pLevelSizes;
    
    std::vector<double*>        pLevels;
    
    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_TETEXACT_SOLVER_CORE_SCHED_HPP

// END
