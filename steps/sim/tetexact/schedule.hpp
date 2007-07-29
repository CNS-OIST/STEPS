////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_SCHEDULE_HPP
#define STEPS_SIM_SCHEDULE_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/sim/simenv.hpp>

START_NAMESPACE(steps)
START_NAMESPACE(sim)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.

// See: steps/sim/kproc.hpp
class KProc;

// Auxiliary declarations.

typedef uint                                    ScheduleIDX;

typedef std::vector<ScheduleIDX>                ScheduleIDXVector;
typedef ScheduleIDXVector::iterator             ScheduleIDXVectorIt;
typedef ScheduleIDXVector::const_iterator       ScheduleIDXVectorCtIt;

////////////////////////////////////////////////////////////////////////////////

class Schedule
{

private:

    ////////////////////////////////////////////////////////////////////////
    // LIST OF KPROCS
    ////////////////////////////////////////////////////////////////////////
    
    std::vector<KProc *>                        pKProcs;
    
    ////////////////////////////////////////////////////////////////////////
    // N-ARY TREE
    ////////////////////////////////////////////////////////////////////////
    
    double                                      pA0;
    
    std::vector<uint>                           pLevelSizes;
    
    std::vector<double*>                        pLevels;
    
    ////////////////////////////////////////////////////////////////////////
    // HELPER FUNCTIONS
    ////////////////////////////////////////////////////////////////////////

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Schedule(void);
    
    ~Schedule(void);
    
    void addKProc(KProc * kproc);
    
    void build();
    
    ////////////////////////////////////////////////////////////////////////
    // SCHEDULE OPERATIONS
    ////////////////////////////////////////////////////////////////////////
    
    inline double getA0(void) const
    {
        return pA0;
    }
    
    //KProc * getNext2(double rannum[]) const;
    KProc * getNext(SimEnv & simenv) const;
    
    void reset(void);
    
    void recomp(void);
    
    void update(ScheduleIDXVector const & entries);

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(sim)
END_NAMESPACE(steps)

#endif
// STEPS_SIM_SCHEDULE_HPP

// END
