///////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
///////////////////////////////////////////////////////////////////////////////

%module wmdirect_core

%import "steps/common.h"

%{
//#include <steps/rng/rng.hpp>

class RNG
{
    
    public:
    
	RNG(uint bufsize);
        virtual ~RNG(void);
    
	virtual double getExp(double lambda);
};

#include <steps/sim/wmdirect_core/wmdirect.hpp>
%}

///////////////////////////////////////////////////////////////////////////////
// CREATION & DESTRUCTION
///////////////////////////////////////////////////////////////////////////////

State * create(void);
void destroy(State * s);

void set_rng(State * s, RNG * rng);

uint species_create(State * s, char * name);

uint comp_create(State * s, char * name);
void comp_species_add(State * s, uint cidx, uint sidx); 

void setup(State * s);

///////////////////////////////////////////////////////////////////////////////
// SIMULATION CONTROLS
///////////////////////////////////////////////////////////////////////////////

void reset(State * s);
void step(State * s);
void run(State * s, double endtime);

///////////////////////////////////////////////////////////////////////////////
// SOLVER STATE ACCESS: 
//     GENERAL
///////////////////////////////////////////////////////////////////////////////

double get_time(State * s);
uint get_nsteps(State * s);

///////////////////////////////////////////////////////////////////////////////
// SOLVER STATE ACCESS: 
//     COMPARTMENT
///////////////////////////////////////////////////////////////////////////////

double get_comp_vol(State * s, uint cidx);
void set_comp_vol(State * s, uint cidx, double vol);

uint get_comp_count(State * s, uint cidx, uint sidx);
void set_comp_count(State * s, uint cidx, uint sidx, uint n);

double get_comp_mass(State * s, uint cidx, uint sidx);
void set_comp_mass(State * s, uint cidx, uint sidx, double m);

double get_comp_conc(State * s, uint cidx, uint sidx);
void set_comp_conc(State * s, uint cidx, uint sidx, double c);

bool get_comp_buffer(State * s, uint cidx, uint sidx);
void set_comp_buffer(State * s, uint cidx, uint sidx, bool buf);

///////////////////////////////////////////////////////////////////////////////

// END
