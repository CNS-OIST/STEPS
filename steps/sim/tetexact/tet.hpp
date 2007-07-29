////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_TETEXACT_COMP_HPP
#define STEPS_SIM_TETEXACT_COMP_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>

////////////////////////////////////////////////////////////////////////////////

#define REAC_LUMP_FACTOR 4

/// Whether the lumping idea positively impacts simulation efficiency
/// should be verified empirically later. To allow this to be changed 
/// easily, we provided accessors for the data fields.
///
class Tet
{

public:

    ////////////////////////////////////////////////////////////////////////

    /// Default constructor.
    Tet(void);
    
    /// Auxiliary pseudo-constructor (objects of type Tet will normally
    /// be build in series).
    void construct
    (
        CompDef * cdef, 
        Tet * n0, Tet * n1, Tet * n2, Tet * n3,
        Tri * t0, Tri * t1, Tri * t2, Tri * t3
    );
    
    /// Destructor.
    ~Tet(void);

    ////////////////////////////////////////////////////////////////////////
    // ACCESS TO SPECIES STUFF
    
    uint poolCount(uint lidx) const
    { return fPoolCount[lidx]; }
    
    uint poolFlags(uint lidx) const
    { return fPoolFlags[lidx]; }
    
    ////////////////////////////////////////////////////////////////////////
    // ACCESS TO REACTION CHANNEL STUFF
    
    double reacKcst(uint lidx) const
    { return fReacStuff[lidx * REAC_LUMP_FACTOR]; }
    
    double reacCcst(uint lidx) const
    { return fReacStuff[(lidx * REAC_LUMP_FACTOR) + 1]; }
    
    double reacHcst(uint lidx) const
    { return fReacStuff[(lidx * REAC_LUMP_FACTOR) + 2]; }
    
    double reacProp(uint lidx) const
    { return fReacStuff[(lidx * REAC_LUMP_FACTOR) + 3]; }
    
    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////
    // SPECIES DATA

    /// Numbers of molecules -- stored as machine word integers.
    ///
    uint *                      pPoolCount;
    
    /// Flags on these pools -- stored as machine word flags.
    ///
    uint *                      pPoolFlags;
    
    ////////////////////////////////////////////////////////////////////////
    // REACTION CHANNELS
    
    KProc **					pReacs;
    
    ////////////////////////////////////////////////////////////////////////
    // DIFFUSION RULES
    
    KProc **
    
};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_TETEXACT_COMP_HPP

// END
