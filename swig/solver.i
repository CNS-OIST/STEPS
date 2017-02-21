/*
 ___license_placeholder___
 */


%module solver_swig

%include "std_string.i"
%include "std_vector.i"

%template(vector_unsigned) std::vector<unsigned int>;


#ifdef WITH_NUMPY
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (unsigned int* IN_ARRAY1, int DIM1) {
    (unsigned int* indices, int input_size)
}
%apply (double* INPLACE_ARRAY1, int DIM1) {
    (double* counts, int output_size)
}

%import "unchecked_stl_seq.i"
UNCHECKED_STL_SEQ_CONVERT(std::vector<unsigned int>,push_back,PyInt_AsUnsignedLongMask)

#endif

////////////////////////////////////////////////////////////////////////////////

%include "error.i"
%import "steps/common.h"
%{
#include "steps/solver/api.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/wmrk4/wmrk4.hpp"
#include "steps/wmdirect/wmdirect.hpp"
#include "steps/tetexact/tetexact.hpp"
#include "steps/tetode/tetode.hpp"
#include "steps/error.hpp"
#include <limits>
%}

#define fwd_api_enum(x) \
%inline %{ enum { x = steps::solver::API::x }; %}

fwd_api_enum(EF_NONE)
fwd_api_enum(EF_DEFAULT)
fwd_api_enum(EF_DV_BDSYS)
fwd_api_enum(EF_DV_SLUSYS)

namespace steps
{
namespace solver
{

class Statedef;

}
}

////////////////////////////////////////////////////////////////////////////////

namespace steps
{
namespace tetode
{

struct structC
{
    unsigned int order;
	unsigned int spec_idx;
};
    
struct structB
{
    std::vector<steps::tetode::structC> info;
};
        
struct structA
{
    double  ccst;
    int upd;
    std::vector<steps::tetode::structB> players;
};

}
}

////////////////////////////////////////////////////////////////////////////////

namespace std
{
    %template(vector_strb) vector<steps::tetode::structB>;
}

////////////////////////////////////////////////////////////////////////////////

%feature("autodoc", "1");

////////////////////////////////////////////////////////////////////////////////

%exception 
{
	try {
		$action
	} catch (steps::ArgErr & ae) {
		PyErr_SetString(PyExc_NameError, ae.getMsg());
		return NULL;
	} catch (steps::NotImplErr & nie) {
		PyErr_SetString(PyExc_NotImplementedError, nie.getMsg());
		return NULL;
	} catch (steps::ProgErr & pe){
        PyErr_SetString(PyExc_RuntimeError, pe.getMsg());
        return NULL;
    }
}

////////////////////////////////////////////////////////////////////////////////

namespace steps
{
namespace solver
{
	
////////////////////////////////////////////////////////////////////////////////

class API
{
public:
	
    %feature("autodoc", "1");
    API(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r = 0);
    %feature("autodoc", "1");
    virtual ~API(void);
	
    %feature("autodoc", 
"
Returns a string of the solver's name.

Syntax::
    
    getSolverName()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverName(void) const = 0;
    
    %feature("autodoc", 
"
Returns a string giving a short description of the solver.

Syntax::
    
    getSolverDesc()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverDesc(void) const = 0;
    
    %feature("autodoc", 
"
Returns a string of the solver authors names.

Syntax::
    
    getSolverAuthors()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverAuthors(void) const = 0;
    
    %feature("autodoc", 
"
Returns a string giving the author's email address.

Syntax::
    
    getSolverEmail()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverEmail(void) const = 0;
    
    %feature("autodoc", 
"
Checkpoint data to a file.
    
Syntax::
    
    checkpoint(file_name)
    
Arguments:
    string file_name
    
Return:
    None
");
    virtual void checkpoint(std::string const & file_name) = 0;
    
    %feature("autodoc", 
"
Restore data from a file.
    
Syntax::
    
    restore(file_name)
    
Arguments:
    string file_name
    
Return:
    None
");
    virtual void restore(std::string const & file_name) = 0;

    %feature("autodoc", 
"
Reset the simulation to the state the solver was initialised to. 
Typically, this resets all concentrations of all chemical species in 
all elements (whether compartments and patches in a well-mixed solver 
or tetrahedrons and triangles in a mesh-based solver) to zero, 
resets the simulation time to zero and resets reaction (and diffusion) 
rates to the default values described in the steps.model objects. 
All reaction (and diffusion) rules are reset to active and all 
compartment volumes and patch areas are reset to default values 
described in steps.geom objects (for well-mixed solvers). 
Usually, this method should be called before starting each simulation iteration.

Syntax::
    
    reset()
    
Arguments:
    None

Return:
    None
");
    virtual void reset(void) = 0;
    
    %feature("autodoc", 
"
Advance the simulation until endtime (given in seconds) is reached. 
The endtime must be larger or equal to the current simulation time.

Syntax::
    
    run(endtime)
    
Arguments:
    float endtime

Return:
    None
");
    virtual void run(double endtime) = 0;
    
    %feature("autodoc", 
"
Returns the current simulation time in seconds.

Syntax::
    
    getTime()
    
Arguments:
    None

Return:
    float
");
    virtual double getTime(void) const = 0;
    
    %feature("autodoc", 
"
Advance the simulation for secs seconds. 

Syntax::
    
    advance(adv)
    
Arguments:
    float adv

Return:
    None
");
	virtual void advance(double adv);
    
    %feature("autodoc", 
"
Advance the simulation for one 'step'. In stochastic solvers this is one 
'realization' of the Gillespie SSA (one reaction 'event'). 
In numerical solvers (currently Wmrk4) this is one time-step, with the 
stepsize defined with the setDT method.

Syntax::
    
    step()
    
Arguments:
    None

Return:
    None
");
    virtual void step(void);
	
    %feature("autodoc", 
"
Set the stepsize for numerical solvers. Must be called before running a 
simulation with these solvers (currently Wmrk4) since there is no default 
stepsize. The deterministic solver Wmrk4 implements a fixed stepsize 
(i.e. not adaptive), although the stepsize can be altered at any point 
during the simulation with this method.

Syntax::
    
    setRk4tDT(dt)
    
Arguments:
    float dt

Return:
    None
");
    virtual void setRk4DT(double dt);
    
    %feature("autodoc", 
"
Returns the stepsize for numerical solvers.

Syntax::
    
    getRk4DT()
    
Arguments:
    None

Return:
    float
");
    virtual double getRk4DT(void) const;

    %feature("autodoc", 
 "
Set the stepsize for numerical solvers. Superceded by setRk4DT, but
Included for backwards compatability.
 
Syntax::
 
    setDT(dt)
 
Arguments:
    float dt
 
Return:
    None
");
    virtual void setDT(double dt);

    %feature("autodoc", 
"
Returns the stepsize for numerical solvers.
Superceded by getRk4DT, but included for backwards compatibility.

Syntax::
    
    getDT()
    
Arguments:
    None

Return:
    float
");
    virtual double getDT(void) const;

    %feature("autodoc", 
"
Set the stepsize for membrane potential solver (default 1us).
This is the time for each voltage calculation step. The SSA will
run until passing this stepsize, so in fact each membrane potential 
time step will vary slightly around the dt so as to be aligned with the SSA.

Syntax::
    
    setEFieldDT(dt)
    
Arguments:
    float dt

Return:
    None
");
    virtual void setEfieldDT(double dt);

    %feature("autodoc", 
"
Set the current simulation time.

Syntax::
    
    setTime(time)
    
Arguments:
    float time

Return:
    None
");
	virtual void setTime(double time);
    

    %feature("autodoc", 
"
Returns the total propensity of the current simulation state 
(the total propensity multiplied by an infinitesimally small 
time dt gives the probability that a reaction will occur in that dt). 
For Tetexact this includes the propensity from the extension of the SSA 
for diffusive flux between tetrahedral elements in the mesh.

Syntax::
    
    getA0()
    
Arguments:
    None

Return:
    float
");		
	virtual double getA0(void) const;

	%feature("autodoc",

"
Set the simulation temperature. Currently, this will only
influence the GHK flux rate, so will only influence simulations
including membrane potential calculation.
	
Syntax::
	
	setTemp(temp)
	
Arguments:
	float temp

Return:
	None
");
    virtual void setTemp(double temp);


	%feature("autodoc",
"
Return the simulation temperature.
			 
Syntax::
			 
	getTemp()
		
Arguments:
	None
	
Return:
	float
");
	virtual double getTemp(void);
	
    %feature("autodoc", 
"
Return the number of 'realizations' of the SSA, the number of reaction 
(and diffusion) events in stochastic solvers.

Syntax::
    
    getNSteps()
    
Arguments:
    None

Return:
    uint
");
	virtual unsigned int getNSteps(void) const;
    
	
    %feature("autodoc", 
"
Set the number of 'realizations' of the SSA, the number of reaction 
(and diffusion) events in stochastic solvers.

Syntax::
    
    setNSteps(nsteps)
    
Arguments:
    uint nsteps

Return:
    None
");
    virtual void setNSteps(unsigned int nsteps);

    %feature("autodoc", 
"
Returns the volume of compartment with identifier string comp (in m^3).

Syntax::
    
    getCompVol(comp)
    
Arguments:
    string comp

Return:
    float
");
    double getCompVol(std::string const & c) const;

    %feature("autodoc", 
"
Set the volume of compartment with identifier string comp (in m^3).

Syntax::
    
    setCompVol(comp, vol)
    
Arguments:
    * string comp
    * float vol

Return:
    None
");
    void setCompVol(std::string const & c, double vol);

    %feature("autodoc", 
"
Returns the number of molecules of a species with identifier string spec 
in compartment with identifier string comp.

In a mesh-based simulation this is the combined count from 
all tetrahedral elements in the compartment.

Syntax::
    
    getCompCount(comp, spec)
    
Arguments:
    * string comp
    * string spec

Return:
    float
");
    double getCompCount(std::string const & c, std::string const & s) const;

    %feature("autodoc", 
"
Set the number of molecules of a species with identifier string spec 
in compartment with identifier string comp.

In a mesh-based simulation this is the combined count from 
all tetrahedral elements in the compartment.

Syntax::
    
    setCompCount(comp, spec, nspec)
    
Arguments:
    * string comp
    * string spec
    * uint nspec

Return:
    None
");
    void setCompCount(std::string const & c, std::string const & s, double n);

    %feature("autodoc", 
"
Returns the amount (in mols) of species with identifier string spec in compartment 
with identifier string comp.

In a mesh-based simulation this is the combined amount from all 
tetrahedral elements in the compartment.

Syntax::
    
    getCompAmount(comp, spec)
    
Arguments:
    * string comp
    * string spec

Return:
    float
");
    double getCompAmount(std::string const & c, std::string const & s) const;

    %feature("autodoc", 
"
Set the amount (in mols) of species with identifier string spec in compartment 
with identifier string comp.

In a mesh-based simulation this is the combined amount from all 
tetrahedral elements in the compartment.

Syntax::
    
    setCompAmount(comp, spec, amount)
    
Arguments:
    * string comp
    * string spec
    * float amount

Return:
    None
");
    void setCompAmount(std::string const & c, std::string const & s, double a);

    %feature("autodoc", 
"
Returns the concentration (in Molar units) of species with identifier string spec 
in compartment with identifier string comp.

Note: in a mesh-based simulation this is calculated from the combined 
number of molecules from all tetrahedral elements in the compartment and the total 
volume of the tetrahedrons.

Syntax::
    
    getCompConc(comp, spec)
    
Arguments:
    * string comp
    * string spec

Return:
    float
");
    double getCompConc(std::string const & c, std::string const & s) const;

    %feature("autodoc", 
"
Sets the concentration (in Molar units) of species with identifier string spec 
in compartment with identifier string comp to conc. In a discrete solver the 
continuous concentration is converted to a discrete number of 
molecules.

Note: in a mesh-based simulation the molecules are divided as 
equally as possible over all tetrahedral elements in the compartment (i.e. a 
uniform distribution).

Syntax::

    setCompConc(comp, spec, conc)
    
Arguments:
    * string comp
    * string spec
    * float conc

Return:
    None
");
    void setCompConc(std::string const & c, std::string const & s, double c);

    %feature("autodoc", 
"
Returns True if species with identifier string spec in compartment with identifier 
string comp is clamped, which means the concentration remains the same 
regardless of reactions that consume or produce molecules of this species. 
Returns False if not.

Note: in a mesh-based simulation it returns True only if the species 
is clamped in all tetrahedral elements of the compartment.

Syntax::
    
    getCompClamped(comp, spec)
    
Arguments:
    * string comp
    * string spec

Return:
    bool
");
    bool getCompClamped(std::string const & c, std::string const & s) const;

    %feature("autodoc", 
"
Sets whether the concentration of species with identifier string spec in compartment 
with identifier string comp is clamped (clamped = True) or not (clamped = False). 
If a species is clamped the concentration stays the same regardless of reactions 
that consume or produce molecules of the species.

Note: in a mesh-based simulation this will set the species to be 
clamped or not in all tetrahedral elements of the compartment.

Syntax::
    
    setCompClamped(comp, spec, clamped)
    
Arguments:
    * string comp
    * string spec
    * bool clamped

Return:
    bool
");
    void setCompClamped(std::string const & c, std::string const & s, bool b);

    %feature("autodoc", 
"
Returns the macroscopic reaction constant of reaction with identifier string reac 
in compartment with identifier string comp. The unit of the reaction constant depends 
on the order of the reaction.

Note: In a mesh-based simulation the value for the compartment is 
returned, although individual tetrahedral elements may have different values 
(set with setTetReacK).

Syntax::
    
    getCompReacK(comp, reac)
    
Arguments:
    * string comp
    * string reac

Return:
    float
");
    double getCompReacK(std::string const & c, std::string const & r) const;

    %feature("autodoc", 
"
Sets the macroscopic reaction constant of reaction with identifier string reac 
in compartment with identifier string comp to kf. The unit of the reaction constant 
depends on the order of the reaction.

Note: In a mesh-based simulation this method sets the reaction 
constant in all tetrahedral elements of the compartment to kf

Note: The default value still comes from the steps.model description, so 
calling reset() will return the reaction constant to that value.

Syntax::
    
    setCompReacK(comp, reac, kf)
    
Arguments:
    * string comp
    * string reac
    * float kf

Return:
    None
");
    void setCompReacK(std::string const & c, std::string const & r, double kf);

    %feature("autodoc", 
"
Returns whether a reaction with identifier string reac in compartment with identifier 
string comp is active (True) or not (False). If it's not active this means that a 
reaction will never occur regardless of whether the reactants are present in 
sufficient numbers or not. 

Note: In a mesh-based simulation this method will return True only 
if the reaction is active in all tetrahedral elements in the compartment. 

Syntax::
    
    getCompReacActive(comp, reac)
    
Arguments:
    * string comp
    * string reac

Return:
    bool
");
    bool getCompReacActive(std::string const & c, std::string const & r) const;

    %feature("autodoc", 
"
Activate (active = True) or deactivate (active = False) a reaction with identifier 
string reac in compartment with identifier string comp. If a reaction is not active 
this means that a reaction will never occur regardless of whether the reactants are 
present in sufficient numbers or not.

Note: In a mesh-based simulation this will activate/deactivate the 
reaction in all tetrahedral elements in the compartment. 

Syntax::
    
    setCompReacActive(comp, reac, active)
    
Arguments:
    * string comp
    * string reac
    * bool active

Return:
    None
");
    void setCompReacActive(std::string const & c, std::string const & r, bool a);
	

    %feature("autodoc", 
"
Returns the diffusion constant of diffusion rule with identifier string diff 
in compartment with identifier string comp. This constant is in units m^2/s.

Note: In a mesh-based solver the value for the compartment is 
returned, although individual or groups of tetrahedral elements may have different 
values (set with setTetDiffD). 

Syntax::
    
    getCompDiffD(comp, diff)
    
Arguments:
    * string comp
    * string diff

Return:
    float
");
    double getCompDiffD(std::string const & c, std::string const & d) const;
	

    %feature("autodoc", 
"
Sets the diffusion constant of diffusion rule with identifier string diff 
in compartment with identifier string comp to dcst (in m^2/s).

Note: This method will set the diffusion constant in all tetrahedral elements 
in the compartment.

Note: The default value still comes from the steps.model description, 
so calling reset() will return the diffusion constants to that value. 

Syntax::
    
    setCompDiffD(comp, diff, dcst)
    
Arguments:
    * string comp
    * string diff
    * float dcst

Return:
    None
");
    void setCompDiffD(std::string const & c, std::string const & d, double dcst);
	
    %feature("autodoc", 
"
Returns whether a diffusion rule with identifier string diff in compartment with 
identifier string comp is active (True) or not (False). If diffusion of a species 
is inactive this means the molecules will remain in place and has the same effect 
as a diffusion constant of zero. 

Syntax::
    
    getCompDiffActive(comp, diff)
    
Arguments:
    * string comp
    * string diff

Return:
    bool
");
    bool getCompDiffActive(std::string const & c, std::string const & d) const;
	
    %feature("autodoc", 
"
Activate (active = True) or deactivate (active = False) a diffusion rule with 
identifier string diff in compartment with identifier string comp. If diffusion 
of a species is inactive this means the molecules will remain in place and is 
effectively the same as setting the diffusion constant to zero

Syntax::
    
    setCompDiffActive(comp, diff, active)
    
Arguments:
    * string comp
    * string diff
    * bool active

Return:
    None
");
    void setCompDiffActive(std::string const & c, std::string const & d, bool act);

    %feature("autodoc", 
"
Returns the 'stochastic reaction constant' (or 'specific probability rate constant') 
of reaction with identifier string reac in compartment with identifier string comp.

The 'stochastic reaction constant' multiplied by infinitesimal time interval dt 
gives the average probability that one reaction channel of this reaction type 
will react accordingly in dt.

Note: in a mesh-based simulation (i.e. Tetexact), the stochastic reaction constant 
is computed as the weighted mean of the stochastic reaction constants in all 
tetrahedral elements of the compartment.

Syntax::
    
    getCompReacC(comp, reac)
    
Arguments:
    * string comp
    * string reac

Return:
    float
");
    double getCompReacC(std::string const & c, std::string const & r) const;

    %feature("autodoc", 
"
Returns h_mu, the distinct number of ways in which reaction with identifier string 
reac can occur in compartment with identifier string comp, by computing the product 
of its reactants. Note: in a mesh-based simulation (i.e. Tetexact), returns the sum 
of the h_mu's over all tetrahedral elements in the compartment. 

Syntax::
    
    getCompReacH(comp, reac)
    
Arguments:
    * string comp
    * string reac

Return:
    float
");
    double getCompReacH(std::string const & c, std::string const & r) const;

    %feature("autodoc", 
"
Returns the propensity of reaction with identifier string reac in compartment 
with identifier string comp. 

The propensity of a reaction is a function of state and is defined as the 
function whose product with infinitesimal time dt gives the probability 
that the reaction will occur in the next dt. It is the 'stochastic reaction 
constant' multiplied by 'h_mu'. 

Note: in a mesh-based simulation (i.e. Tetexact), the propensity of a reaction 
in a compartment is computed as the sum of the propensities in all tetrahedral 
elements of the compartment. 

Syntax::
    
    getCompReacA(comp, reac)
    
Arguments:
    * string comp
    * string reac

Return:
    float
");
    double getCompReacA(std::string const & c, std::string const & r) const;

    %feature("autodoc", 
"
Return the extent of reaction with identifier string reac in compartment with 
identifier string comp, that is the number of times the reaction has occurred up 
to the current simulation time. 

Note: in a mesh-based simulation (i.e. Tetexact), returns the sum of the reaction 
extents in all tetrahedral elements of the compartment.

Syntax::
    
    getCompReacExtent(comp, reac)
    
Arguments:
    * string comp
    * string reac

Return:
    uint
");
    unsigned int getCompReacExtent(std::string const & c, std::string const & r) const;

    %feature("autodoc", 
"
Resets the extent of reaction with identifier string reac in compartment with 
identifier string comp to zero. 

Note: in a mesh-based simulation (i.e. Tetexact), 
resets the extents of the reaction in all tetrahedral elements of the compartment.

Syntax::
    
    resetCompReacExtent(comp, reac)
    
Arguments:
    * string comp
    * string reac

Return:
    None
");
    void resetCompReacExtent(std::string const & c, std::string const & r);

    %feature("autodoc", 
"
Returns the volume (in m^3) of the tetrahedral element with index idx.

Syntax::
    
    getTetVol(idx)
    
Arguments:
    * uint idx

Return:
    float
");
    double getTetVol(unsigned int tidx) const;

    // void setTetVol(unsigned int tidx, double vol);

    %feature("autodoc", 
"
Returns whether species with identifier string spec is defined
in the tetrahedral element with index idx.

Syntax::
    
    getTetSpecDefined(idx, spec)
    
Arguments:
    * uint idx
    * string spec

Return:
    bool
");
    bool getTetSpecDefined(unsigned int tidx, std::string const & s) const;

    %feature("autodoc", 
"
Returns the number of molecules of species with identifier string spec 
in the tetrahedral element with index idx.

Syntax::
    
    getTetCount(idx, spec)
    
Arguments:
    * uint idx
    * string spec

Return:
    uint
");
    double getTetCount(unsigned int tidx, std::string const & s) const;
	
    %feature("autodoc", 
"
Sets the number of molecules of species with identifier string spec in 
tetrahedral element with index idx to n.

Syntax::
    
    setTetCount(idx, spec, n)
    
Arguments:
    * uint idx
    * string spec
    * uint n

Return:
    None
");
    void setTetCount(unsigned int tidx, std::string const & s, double n);
	
    %feature("autodoc", 
"
Returns the amount (in mols) of species with identifier string spec in 
tetrahedral element with index idx.

Syntax::
    
    getTetAmount(idx, spec)
    
Arguments:
    * uint idx
    * string spec

Return:
    float
");
    double getTetAmount(unsigned int tidx, std::string const & s) const;
	
    %feature("autodoc", 
"
Sets the amount (in mols) of species with identifier string spec in tetrahedral 
element with index idx to a. This continuous value must be converted internally 
to a discrete number of molecules by multiplication with Avogadro's 
number. 

Due to the small volumes of tetrahedral elements the difference 
between 'rounding up' and 'rounding down' can be a significant difference in 
concentration.

Syntax::
    
    setTetAmount(idx, spec, a)
    
Arguments:
    * uint idx
    * string spec
    * float a

Return:
    None
");
    void setTetAmount(unsigned int tidx, std::string const & s, double m);
	
    %feature("autodoc", 
"
Returns the concentration (in Molar units) of species with identifier 
string spec in a tetrahedral element with index idx.

Syntax::
    
    getTetConc(idx, spec)
    
Arguments:
    * uint idx
    * string spec

Return:
    float
");
    double getTetConc(unsigned int tidx, std::string const & s) const;

	
    %feature("autodoc", 
"
Sets the concentration (in Molar units) of species with identifier string spec 
in a tetrahedral element with index idx to conc.This continuous value must be 
converted internally to a discrete number of molecules. 

Due to the small volumes of tetrahedral elements the difference between 'rounding 
up' and 'rounding down' can be a large difference in concentration.

Syntax::
    
    setTetConc(idx, spec, conc)
    
Arguments:
    * uint idx
    * string spec
    * conc

Return:
    None
");	
    void setTetConc(unsigned int tidx, std::string const & s, double c);

    %feature("autodoc", 
"
Returns True if concentration of species with identifier string spec in tetrahedral 
element with index idx is clamped, which means the concentration stays the 
same regardless of reactions that consume or produce molecules of this species or 
diffusion of this species into or out of the tetrahedral element. Returns False if 
not.

Syntax::
    
    getTetClamped(idx, spec)
    
Arguments:
    * uint idx
    * string spec

Return:
    bool
");	
    bool getTetClamped(unsigned int tidx, std::string const & s) const;
	
    %feature("autodoc", 
"
Sets whether the concentration of species spec in tetrahedral element with 
index idx is clamped (clamped = True) or not (clamped = False). 
If a species is clamped the concentration stays the same regardless 
of reactions that consume or produce molecules of the species or 
diffusion of the species into or out of the tetrahedral element.

Syntax::
    
    setTetClamped(idx, spec, clamped)
    
Arguments:
    * uint idx
    * string spec
    * bool clamped

Return:
    None
");	
    void setTetClamped(unsigned int tidx, std::string const & s, bool buf);
	
    %feature("autodoc", 
"
Returns the macroscopic reaction constant of reaction with identifier string reac 
in tetrahedral element with index idx. The unit of the reaction constant depends 
on the order of the reaction.

Syntax::
    
    getTetReacK(idx, reac)
    
Arguments:
    * uint idx
    * string reac

Return:
    float
");	
    double getTetReacK(unsigned int tidx, std::string const & r) const;
		
    %feature("autodoc", 
"
Sets the macroscopic reaction constant of reaction with identifier string reac 
in tetrahedral element with index idx to kf. The units of the reaction constant 
depends on the order of the reaction.

Syntax::
    
    setTetReacK(idx, reac, kf)
    
Arguments:
    * uint idx
    * string reac
    * float kf

Return:
    None
");
    void setTetReacK(unsigned int tidx, std::string const & r, double kf);
		
    %feature("autodoc", 
"
Returns whether reaction with identifier string reac in tetrahedral element 
with index idx is active (True) or not (False). If it's not active this means 
that the reaction will never occur regardless of whether reactants are present 
in sufficient numbers or not.

Syntax::
    
    getTetReacActive(idx, reac)
    
Arguments:
    * uint idx
    * string reac

Return:
    bool
");
    bool getTetReacActive(unsigned int tidx, std::string const & r) const;
	
    %feature("autodoc", 
"
Activate (active = True) or deactivate (active = False) a reaction with identifier 
string reac in tetrahedral element with index idx. If it's not active this means 
that the reaction will never occur regardless of whether reactants are present 
in sufficient numbers or not.

Syntax::
    
    setTetReacActive(idx, reac, active)
    
Arguments:
    * uint idx
    * string reac
    * bool active

Return:
    None
");
    void setTetReacActive(unsigned int tidx, std::string const & r, bool act);
	
	
    %feature("autodoc", 
"
Returns the diffusion constant of diffusion rule with identifier string diff 
in tetrahedral element with index idx. This constant is in units m^2/s. If direction_tet is specified, return the diffusion constant towards that direction.

Syntax::
    
    getTetDiffD(idx, diff, direction_tet = UINT_MAX)
    
Arguments:
    * uint idx
    * string diff
    * direction_tet
    
Return:
    float
");
    double getTetDiffD(unsigned int tidx, std::string const & d, uint direction_tet = std::numeric_limits<uint>::max()) const;
	
    %feature("autodoc", 
"
Sets the diffusion constant of diffusion rule with identifier string diff in 
tetrahedral element with index idx to dcst (in m^2/s). Specify direction_tet to set the constant only towards a given tetrahedron direction.
Syntax::
    
    setTetDiffD(idx, diff, dcst, direction_tet = UINT_MAX)
    
Arguments:
    * uint idx
    * string diff
    * dcst
    * direction_tet

Return:
    None
");
    void setTetDiffD(uint tidx, std::string const & d, double dk, uint direction_tet = std::numeric_limits<uint>::max());

    %feature("autodoc", 
"
Returns the diffusion constant of diffusion rule with identifier string diff 
in triangle element with index idx. If direction_tri is specified, return the diffusion constant towards that direction.

Syntax::
    
    getTriDiffD(idx, diff, direction_tri = UINT_MAX)
    
Arguments:
    * uint idx
    * string diff
    * direction_tri
    
Return:
    float
");
    double getTriDiffD(unsigned int tidx, std::string const & d, uint direction_tri = std::numeric_limits<uint>::max()) const;

    %feature("autodoc", 
"
Sets the diffusion constant of diffusion rule with identifier string diff in 
triangle element with index idx to dcst. Specify direction_tri to set the constant only towards a given triangle direction.
Syntax::
    
    setTriDiffD(idx, diff, dcst, direction_tri = UINT_MAX)
    
Arguments:
    * uint idx
    * string diff
    * dcst
    * direction_tri

Return:
    None
");
    void setTriDiffD(uint tidx, std::string const & d, double dk, uint direction_tri = std::numeric_limits<uint>::max());
	
    %feature("autodoc", 
"
Returns whether diffusion with identifier string diff in tetrahedral element 
with index idx is active (True) or not (False). If diffusion of a species 
is inactive this means the molecules will never diffuse out of the tetrahedron 
and has the same effect as a diffusion constant of zero.

Syntax::
    
    getTetDiffActive(idx, diff)
    
Arguments:
    * uint idx
    * string diff

Return:
    bool
");
    bool getTetDiffActive(unsigned int tidx, std::string const & d) const;
	
    %feature("autodoc", 
"
Activate (active = True) or deactivate (active = False) diffusion rule with 
identifier string diff in tetrahedral element with index idx. If diffusion of 
a species is inactive this means the molecules will never diffuse out of the 
tetrahedron and has the same effect as a diffusion constant of zero. 

Syntax::
    
    setTetDiffActive(idx, diff, active)
    
Arguments:
    * uint idx
    * string diff
    * bool active

Return:
    None
");
    void setTetDiffActive(unsigned int tidx, std::string const & d, bool act);

    %feature("autodoc", 
"
Returns the 'stochastic reaction constant' (or 'specific probability rate constant') 
of reaction with identifier string reac in tetrahedral element with index idx.

Syntax::
    
    getTetReacC(idx, reac)
    
Arguments:
    * uint idx
    * string reac

Return:
    float
");
    double getTetReacC(unsigned int tidx, std::string const & r) const;

    %feature("autodoc", 
"
Returns h_mu, the distinct number of ways in which reaction with identifier string 
reac can occur in tetrahedral element with index idx, by computing the product of 
its reactants.

Syntax::
    
    getTetReacH(idx, reac)
    
Arguments:
    * uint idx
    * string reac

Return:
    float
");
    double getTetReacH(unsigned int tidx, std::string const & r) const;

    %feature("autodoc", 
"
Returns the propensity of reaction with identifier string reac in tetrahedral 
element with index idx.

Syntax::
    
    getTetReacA(idx, reac)
    
Arguments:
    * uint idx
    * string reac

Return:
    float
");
    double getTetReacA(unsigned int tidx, std::string const & r) const;
	
    %feature("autodoc", 
"
Returns the propensity of diffusion rule with identifier string diff in 
tetrahedral element with index idx. 

Syntax::
    
    getTetDiffA(idx, reac)
    
Arguments:
    * uint idx
    * string reac

Return:
    float
");
    double getTetDiffA(unsigned int tidx, std::string const & d) const;
	
	%feature("autodoc",
" 
Returns the potential (in volts) of tetrahedral element with index idx, taken at the barycenter.
			
Syntax::
			 
	getTetV(idx)
			 
Arguments:
	uint idx
		
Returns:
	float
");
	double getTetV(unsigned int tidx) const;
	
	%feature("autodoc",
" 
Set the potential (in volts) of tetrahedral element with index idx.
			
Syntax::
			 
	setTetV(idx, v)
			 
Arguments:
	* uint idx
	* float v
			 
Returns:
	None
");	
    void setTetV(unsigned int tidx, double v);

	%feature("autodoc",
" 
Returns true if the potential of tetrahedral element with index idx is clamped
to some voltage.
			
Syntax::
			 
	getTetVClamped(idx)
			 
Arguments:
	uint idx
		
Returns:
	bool
");	
    bool getTetVClamped(unsigned int tidx) const;

	%feature("autodoc",
" 
Sets whether the potential of tetrahedral element with index idx is clamped
(clamped = True) or not (clamped = False).
			
Syntax::
			 
	setTetVClamped(idx, clamped)
			 
Arguments:
	* uint idx
	* bool clamped
			 
Returns:
	None
");		
    void setTetVClamped(unsigned int tidx, bool cl);	

    %feature("autodoc", 
"
Returns the area of patch with identifier string pat (in m^2).

Syntax::
    
    getPatchArea(pat)
    
Arguments:
    * string pat

Return:
    float
");
    double getPatchArea(std::string const & p) const;

    %feature("autodoc", 
"
Sets the area of patch with identifier string pat to area a (in m^2).

Syntax::
    
    setPatchArea(pat, area)
    
Arguments:
    * string pat
    * float area

Return:
    None
");
    void setPatchArea(std::string const & p, double area);

    %feature("autodoc", 
"
Returns the number of molecules of species with identifier string spec in patch 
with identifier string pat.Note: in a mesh-based simulation this 
is the combined count from all triangular elements in the patch. 

Syntax::
    
    getPatchCount(pat, spec)
    
Arguments:
    * string pat
    * string spec

Return:
    float
");
    double getPatchCount(std::string const & p, std::string const & s) const;

    %feature("autodoc", 
"
Sets the number of molecules of species with identifier string spec in patch 
with identifier string pat to n. Note: in a mesh-based simulation the molecules 
are divided as equally as possible over all triangular elements in 
the patch (i.e. a uniform distribution). 

Syntax::
    
    setPatchCount(pat, spec, n)
    
Arguments:
    * string pat
    * string spec
    * uint n

Return:
    float
");
    void setPatchCount(std::string const & p, std::string const & s, double n);

    %feature("autodoc", 
"
Returns the amount (in mols) of species with identifier string spec in patch 
with identifier string pat.

Note: in a mesh-based simulation this is the combined amount 
from all triangular elements in the patch. 

Syntax::
    
    getPatchAmount(pat, spec)
    
Arguments:
    * string pat
    * string spec

Return:
    float
");
    double getPatchAmount(std::string const & p, std::string const & s) const;

    %feature("autodoc", 
"
Sets the amount (in mols) of species with identifier string spec in patch with 
identifier string pat to a. In a discrete solver, such as Wmdirect and Tetexact, 
this continuous value is converted internally into a discrete number of molecules 
by multiplication with Avogadro's number. 

Note: in a mesh-based simulation the molecules are divided as 
equally as possible over all triangular elements in the patch (i.e. a uniform 
distribution).

Syntax::
    
    setPatchAmount(pat, spec, a)
    
Arguments:
    * string pat
    * string spec
    * float a

Return:
    None
");
    void setPatchAmount(std::string const & p, std::string const & s, double a);

    %feature("autodoc", 
"
Sets the amount (in mols) of species with identifier string spec in patch with 
identifier string pat to a. In a discrete solver, such as Wmdirect and Tetexact, 
this continuous value is converted internally into a discrete number of molecules 
by multiplication with Avogadro's number. 

Note: in a mesh-based simulation the molecules are divided as equally 
as possible over all triangular elements in the patch (i.e. a uniform distribution).

Syntax::
    
    getPatchClamped(pat, spec)
    
Arguments:
    * string pat
    * string spec

Return:
    bool
");
    bool getPatchClamped(std::string const & p, std::string const & s) const;

    %feature("autodoc", 
"
Sets whether the species with identifier string spec in patch with identifier 
string pat is clamped (clamped = True) or not (clamped = False). If a species 
is clamped the number of molecules stays the same regardless of surface reactions 
that consume or produce molecules of the species.

Note: in a mesh-based simulation this will set the species to be clamped in all 
triangular elements of the patch.

Syntax::
    
    setPatchClamped(pat, spec, clamped)
    
Arguments:
    * string pat
    * string spec
    * bool clamped

Return:
    None
");
    void setPatchClamped(std::string const & p, std::string const & s, bool buf);

    %feature("autodoc", 
"
Returns the macroscopic reaction constant of surface reaction with identifier 
string sreac in patch with identifier string pat. The unit of the reaction constant 
depends on the order of the reaction.

Note: In a mesh-based solver the value for the patch is returned, 
although individual triangle elements may have different values 
(set with setTriSReacK).

Syntax::
    
    getPatchSReacK(pat, reac)
    
Arguments:
    * string pat
    * string reac

Return:
    float
");
    double getPatchSReacK(std::string const & p, std::string const & r) const;

    %feature("autodoc", 
"
Sets the macroscopic reaction constant of surface reaction with identifier 
string sreac in patch with identifier string pat to kf. The unit of the reaction 
constant depends on the order of the reaction. 

Note: In a mesh-based simulation this method sets the surface 
reaction constant in all triangular elements of the patch to kf.

Note: The default value still comes from the steps.model description, so calling 
reset() will return the surface reaction constant to that value.

Syntax::
    
    setPatchSReacK(pat, reac, kf)
    
Arguments:
    * string pat
    * string reac
    * float kf

Return:
    None
");
    void setPatchSReacK(std::string const & p, std::string const & r, double kf);

    %feature("autodoc", 
"
Returns whether a surface reaction with identifier string sreac in patch with 
identifier string pat is active (True) or not (False). If it's not active this means 
that a surface reaction will never occur regardless of whether the reactants are 
present in sufficient numbers or not. 

Note: In a mesh-based simulation this method will return True only 
if the surface reaction is active in all triangular elements in the patch.

Syntax::
    
    getPatchSReacActive(pat, reac)
    
Arguments:
    * string pat
    * string reac

Return:
    bool
");
    bool getPatchSReacActive(std::string const & p, std::string const & r) const;

    %feature("autodoc", 
"
Activate (active = True) or deactivate (active = False) a surface reaction with 
identifier string sreac in patch with identifier string pat. If a surface reaction 
is not active this means that a reaction will never occur regardless of whether the 
reactants are present in sufficient numbers or not.

Note: In a mesh-based simulation this will activate/ deactivate the 
reaction in all triangular elements in the patch.

Syntax::
    
    setPatchSReacActive(pat, reac, active)
    
Arguments:
    * string pat
    * string reac
    * bool active

Return:
    None
");
    void setPatchSReacActive(std::string const & p, std::string const & r, bool a);

    %feature("autodoc", 
"
Returns the 'stochastic reaction constant' (or 'specific probability rate constant') 
of surface reaction with identifier string sreac in patch with identifier string pat.

Note: in a mesh-based simulation (i.e. Tetexact), the stochastic reaction constant is 
computed as the weighted mean of the stochastic reaction constants in all triangular 
elements of the patch.

Syntax::
    
    getPatchSReacC(pat, reac)
    
Arguments:
    * string pat
    * string reac

Return:
    float
");
    double getPatchSReacC(std::string const & p, std::string const & r) const;

    %feature("autodoc", 
"
Returns h_mu, the distinct number of ways in which surface reaction with identifier 
string sreac can occur in patch with identifier string pat, by computing the product 
of its reactants. Note: in a mesh-based simulation (i.e. Tetexact), returns the sum 
of the h_mu's over all triangular elements in the patch. 

Syntax::
    
    getPatchSReacH(pat, reac)
    
Arguments:
    * string pat
    * string reac

Return:
    float
");
    double getPatchSReacH(std::string const & p, std::string const & r) const;

    %feature("autodoc", 
"
Returns the propensity of surface reaction with identifier string sreac in patch 
with identifier string pat. Note: in a mesh-based simulation (i.e. Tetexact), 
the propensity of a surface reaction in a patch is computed as the sum of the 
propensities in all triangular elements of the patch.

Syntax::
    
    getPatchSReacA(pat, reac)
    
Arguments:
    * string pat
    * string reac

Return:
    float
");
    double getPatchSReacA(std::string const & p, std::string const & r) const;

    %feature("autodoc", 
"
Returns the extent of surface reaction with identifier string sreac in patch 
with identifier string pat, that is the number of times the surface reaction 
has occurred up to the current simulation time. 

Note: in a mesh-based simulation (i.e. Tetexact), returns the sum of the reaction 
extents in all triangular elements of the patch.

Syntax::
    
    getPatchSReacExtent(pat,reac)
    
Arguments:
    * string pat
    * string reac

Return:
    uint
");
    unsigned int getPatchSReacExtent(std::string const & p, std::string const & r) const;

    %feature("autodoc", 
"
Resets the extent of reaction with identifier string sreac in patch with identifier 
string pat to zero. 

Note: in a mesh-based simulation (i.e. Tetexact), resets the extents of the reaction 
in all triangular elements of the patch.

Syntax::
    
    resetPatchSReacExtent(pat, reac)
    
Arguments:
    * string pat
    * string reac

Return:
    None
");
    void resetPatchSReacExtent(std::string const & p, std::string const & r);

    %feature("autodoc", 
"
Returns whether a voltage-dependent surface reaction with identifier string vsreac in patch with 
identifier string pat is active (True) or not (False). If it's not active this means 
that the voltage-dependent surface reaction will never occur regardless of whether the reactants are 
present in sufficient numbers or not. 

Note: In a mesh-based simulation this method will return True only 
if the voltage-dependent surface reaction is active in all triangular elements in the patch.

Syntax::

    getPatchVDepSReacActive(pat, vsreac)

Arguments:
    * string pat
    * string vsreac

Return:
    bool
");
    bool getPatchVDepSReacActive(std::string const & p, std::string const & vsr) const;

    %feature("autodoc", 
"
Activate (active = True) or deactivate (active = False) a voltage-dependent surface reaction with 
identifier string vsreac in patch with identifier string pat. If a voltage-dependent surface reaction 
is not active this means that a reaction will never occur regardless of whether the 
reactants are present in sufficient numbers or not.

Note: In a mesh-based simulation this will activate/ deactivate the 
reaction in all triangular elements in the patch.

Syntax::

    setPatchVDepSReacActive(pat, vsreac, active)

Arguments:
    * string pat
    * string vsreac
    * bool active

Return:
    None
");
    void setPatchVDepSReacActive(std::string const & p, std::string const & vsr, bool a);

    %feature("autodoc",
"
Activates or inactivates diffusion across a diffusion boundary for a species.
             
Syntax::
             
    setDiffBoundaryDiffusionActive(diffb, spec, act)
             
Arguments:
    * string diffb
    * string spec
    * bool act
             
Return:
    None
");
    
    void setDiffBoundaryDiffusionActive(std::string const & db, std::string const & s, bool act);
    
    %feature("autodoc",
"
Returns whether diffusion is active across a diffusion boundary for a species.
             
Syntax::
             
    getDiffBoundaryDiffusionActive(diffb, spec)
             
Arguments:
    * string diffb
    * string spec
             
Return:
    bool
");
    bool getDiffBoundaryDiffusionActive(std::string const & db, std::string const & s) const;

    %feature("autodoc",
"
Set the diffusion constant of tetrahedrons across a diffusion boundary. If direction_comp is provided, only set dcsts of diffusion towards it (Directional dcsts of diffusions in tetrahedrons in the other compartment of the diffusion boundary towards tetrahedons in the direction compartment).
             
Syntax::
             
    setDiffBoundaryDcst(diffb, spec, dcst, direction_comp = '')
             
Arguments:
    * string diffb
    * string spec
    * double dcst
    * direction_comp
             
Return:
    None
");
    void setDiffBoundaryDcst(std::string const & db, std::string const & s, double dcst, std::string const & direction_comp = "");
    
    %feature("autodoc", 
"
Returns the area (in m^2) of the triangular element with index idx.

Syntax::
    
    getTriArea(idx)
    
Arguments:
    * uint idx

Return:
    float
");
    double getTriArea(unsigned int tidx) const;
	
    //void setTriArea(unsigned int tidx, double area);

    %feature("autodoc", 
"
Returns whether species with identifier string spec is defined
in the triangle element with index idx.

Syntax::
    
    getTriSpecDefined(idx, spec)
    
Arguments:
    * uint idx
    * string spec

Return:
    bool
");
    bool getTriSpecDefined(unsigned int tidx, std::string const & s) const;
	
    %feature("autodoc", 
"
Returns the number of molecules of species with identifier string spec 
in the triangular element with index idx.

Syntax::
    
    getTriCount(idx, spec)
    
Arguments:
    * uint idx
    * string spec

Return:
    float
");
    double getTriCount(unsigned int tidx, std::string const & s) const;

    %feature("autodoc", 
"
Sets the number of molecules of species with identifier string spec in 
triangular element with index idx to n. 

Syntax::
    
    setTriCount(idx, spec, n)
    
Arguments:
    * uint idx
    * string spec
    * uint n

Return:
    None
");
    void setTriCount(unsigned int tidx, std::string const & s, double n);

    %feature("autodoc", 
"
Returns the amount (in mols) of species with identifier string spec in triangular 
element with index idx.  

Syntax::
    
    getTriAmount(idx, spec)
    
Arguments:
    * uint idx
    * string spec

Return:
    float
");
    double getTriAmount(unsigned int tidx, std::string const & s) const;

    %feature("autodoc", 
"
Sets the amount (in mols) of species with identifier string spec in triangular 
element with index idx to a. This continuous value must be converted internally 
to a discrete number of molecules by multiplication with Avogadro's number. 

Syntax::
    
    setTriAmount(idx, spec, a)
    
Arguments:
    * uint idx
    * string spec
    * float a

Return:
    None
");
    void setTriAmount(unsigned int tidx, std::string const & s, double m);
    
    %feature("autodoc", 
"
Returns True if the species with identifier string spec in triangular element 
with index idx is clamped, which means the number of molecules stays 
the same regardless of reactions that consume or produce molecules of this species. 
Returns False if not.

Syntax::
    
    getTriClamped(idx, spec)
    
Arguments:
    * uint idx
    * string spec

Return:
    bool
");
    bool getTriClamped(unsigned int tidx, std::string const & s) const;
	
    %feature("autodoc", 
"
Sets whether the concentration of species spec in triangular element with index idx 
is clamped (clamped = True) or not (clamped = False). If a species is clamped the 
concentration stays the same regardless of reactions that consume or produce 
molecules of the species. 

Syntax::
    
    setTriClamped(idx, spec, clamped)
    
Arguments:
    * uint idx
    * string spec
    * bool clamped

Return:
    None
");
    void setTriClamped(unsigned int tidx, std::string const & s, bool buf);

    %feature("autodoc", 
"
Returns the macroscopic reaction constant of surface reaction with identifier 
string sreac in triangular element with index idx. The units of the reaction 
constant depends on the order of the reaction. 

Syntax::
    
    getTriSReacK(idx, reac)
    
Arguments:
    * uint idx
    * string reac

Return:
    float
");
    double getTriSReacK(unsigned int tidx, std::string const & r) const;

    %feature("autodoc", 
"
Sets the macroscopic reaction constant of surface reaction with identifier 
string sreac in triangular element with index idx to kf. The units of the 
reaction constant depends on the order of the reaction.

Syntax::
    
    setTriSReacK(idx, reac, kf)
    
Arguments:
    * uint idx
    * string reac
    * float kf

Return:
    None
");
    void setTriSReacK(unsigned int tidx, std::string const & r, double kf);
	
    %feature("autodoc", 
"
Returns whether surface reaction with identifier string sreac in triangular 
element with index idx is active (True) or not (False). If it's not active 
this means that the surface reaction will never occur regardless of whether 
reactants are present in sufficient numbers or not. 

Syntax::
    
    getTriSReacActive(idx, reac)
    
Arguments:
    * uint idx
    * string reac

Return:
    bool
");
    bool getTriSReacActive(unsigned int tidx, std::string const & r) const;
	
    %feature("autodoc", 
"
Activate (active = True) or deactivate (active = False) a surface reaction 
with identifier string sreac in triangular element with index idx. If it's 
not active this means that the surface reaction will never occur regardless 
of whether reactants are present in sufficient numbers or not.  

Syntax::
    
    setTriSReacActive(idx, reac, active)
    
Arguments:
    * uint idx
    * string reac
    * active

Return:
    None
");
    void setTriSReacActive(unsigned int tidx, std::string const & r, bool act);
	
    %feature("autodoc", 
"
Returns the 'stochastic reaction constant' (or 'specific probability rate constant') 
of surface reaction with identifier string sreac in triangular element with index idx.  

Syntax::
    
    getTriSReacC(idx, reac)
    
Arguments:
    * uint idx
    * string reac

Return:
   float
");
    double getTriSReacC(unsigned int tidx, std::string const & r) const;

    %feature("autodoc", 
"
Returns h_mu, the distinct number of ways in which surface reaction with identifier 
string sreac can occur in triangular element with index idx, by computing the product 
of its reactants. 

Syntax::
    
    getTriSReacH(idx, reac)
    
Arguments:
    * uint idx
    * string reac

Return:
   float
");
    double getTriSReacH(unsigned int tidx, std::string const & r) const;

    %feature("autodoc", 
"
Returns the propensity of surface reaction with identifier string sreac 
in triangular element with index idx. 

Syntax::
    
    getTriSReacA(idx, reac)
    
Arguments:
    * uint idx
    * string reac

Return:
   float
");
    double getTriSReacA(unsigned int tidx, std::string const & r) const;

	
    %feature("autodoc",
" 
Returns the potential (in volts) of triangle element with index idx, taken at the barycenter.
			
Syntax::
			 
	getTriV(idx)
			 
Arguments:
	* uint idx
		
Returns:
	float
");	
	double getTriV(unsigned int tidx) const;
	
	%feature("autodoc",
" 
Set the potential (in volts) of triangle element with index idx.
			
Syntax::
			 
	setTriV(idx, v)
			 
Arguments:
	* uint idx
	* float v
			 
Returns:
	None
");		
    void setTriV(unsigned int tidx, double v);

	%feature("autodoc",
" 
Returns true if the potential of triangle element with index idx is clamped
to some voltage.
			
Syntax::
			 
	getTriVClamped(idx)
			 
Arguments:
	* uint idx
		
Returns:
	bool
");		
    bool getTriVClamped(unsigned int tidx) const;

	%feature("autodoc",
" 
Sets whether the potential of triangle element with index idx is clamped
(clamped = True) or not (clamped = False).
			
Syntax::
			 
	setTriVClamped(idx, clamped)
			 
Arguments:
	* uint idx
	* bool clamped
			 
Returns:
	None
");			
    void setTriVClamped(unsigned int tidx, bool cl);

%feature("autodoc",
" 
Returns the ohmic current of triangle element with index idx, in amps.
                       
Syntax::
                        
       getTriOhmicI(idx)
                        
Arguments:
       * uint idx
               
Returns:
       float
");            
    double getTriOhmicI(unsigned int tidx);

%feature("autodoc",
" 
Returns the ohmic current with identifier string oc of triangle element with index idx, in amps.
             
Syntax::
                        
    getTriOhmicI(idx, oc)
                        
Arguments:
    * uint idx
    * string oc

Returns:
    float
");            
    double getTriOhmicI(unsigned int tidx, std::string const & oc);

%feature("autodoc",
"
Returns the GHK current of triangle element with index idx, in amps.
             
Syntax::
			 
    getTriGHKI(idx)
			 
Arguments:
    * uint idx
             
Returns:
    float
");		
    double getTriGHKI(unsigned int tidx);
    
%feature("autodoc",
" 
    Returns the ghk current with identifier string ghk of triangle element with index idx, in amps.
             
Syntax::
			 
    getTriGHKI(idx, ghk)
			 
Arguments:
    * uint idx
    * string ghk
             
Returns:
    float
");		
    double getTriGHKI(unsigned int tidx, std::string const & ghk);
    
%feature("autodoc",
" 
Returns the current of triangle element with index idx, in amps, 
at the last EField calculation step.
			 
Syntax::
			 
    getTriI(idx)
			 
Arguments:
    * uint idx
			 
Returns:
    float
");		
    double getTriI(unsigned int tidx) const;

%feature("autodoc",
" 
Set current clamp to triangle element with index idx to current i (amps).
NOTE: Convention is maintained that a positive current clamp is depolarizing, a negative current clamp is hyperpolarizing.
			
Syntax::
			 
	setTriIClamp(idx, i)
			 
Arguments:
	* uint idx
	* float i
			 
Returns:
	None
");			
    void setTriIClamp(unsigned int tidx, double i);

    %feature("autodoc", 
"
Returns whether voltage-dependent surface reaction with identifier string vsreac in triangular 
element with index idx is active (True) or not (False). If it's not active 
this means that the voltage-dependent surface reaction will never occur regardless of whether 
reactants are present in sufficient numbers or not. 

Syntax::

    getTriVDepSReacActive(idx, reac)

Arguments:
    * uint idx
    * string vsreac

Return:
    bool
");
    bool getTriVDepSReacActive(unsigned int tidx, std::string const & vsr) const;

    %feature("autodoc", 
"
Activate (active = True) or deactivate (active = False) a voltage-dependent surface reaction 
with identifier string vsreac in triangular element with index idx. If it's 
not active this means that the voltage-dependent surface reaction will never occur regardless 
of whether reactants are present in sufficient numbers or not.  

Syntax::

    setTriVDepSReacActive(idx, vsreac, active)

Arguments:
    * uint idx
    * string vsreac
    * active

Return:
    None
");
    void setTriVDepSReacActive(unsigned int tidx, std::string const & vsr, bool act);


    %feature("autodoc",
" 
Returns the potential (in volts) of vertex element with index idx.
			
Syntax::
			 
	getVertV(idx)
			 
Arguments:
	* uint idx
		
Returns:
	float
");	
    double getVertV(unsigned int vidx) const;

	%feature("autodoc",
" 
Set the potential (in volts) of vertex element with index idx.
			
Syntax::
			 
	setVertV(idx, v)
			 
Arguments:
	* uint idx
	* float v
			 
Returns:
	None
");		
    void setVertV(unsigned int vidx, double v);

	%feature("autodoc",
" 
Returns true if the potential of vertex element with index idx is clamped
to some voltage.
			
Syntax::
			 
	getVertVClamped(idx)
			 
Arguments:
	* uint idx
		
Returns:
	bool
");		
    bool getVertVClamped(unsigned int vidx) const;

	%feature("autodoc",
" 
Sets whether the potential of vertex element with index idx is clamped
(clamped = True) or not (clamped = False).
			
Syntax::
			 
	setVertVClamped(idx, clamped)
			 
Arguments:
	* uint idx
	* bool clamped
			 
Returns:
	None
");		
    void setVertVClamped(unsigned int vidx, bool cl);

    %feature("autodoc",
" 
Set current clamp to vertex element with index idx to current i (Amps). 
NOTE: Convention is maintained that a positive current clamp is depolarizing, a negative current clamp is hyperpolarizing.

Syntax::

    setVertIClamp(idx, i)
			 
Arguments:
    * uint idx
    * float i
			 
Returns:
    None
");			
    void setVertIClamp(unsigned int vidx, double i);
    
	%feature("autodoc", 
"
Sets the potential (in volts) of membrane with string identifier memb.
NOTE: This method will set the potential of all nodes in the volume conductor
to the same value.

Syntax::
			
	setMembPotential(memb, v)
		
Arguments:
	* string memb
	* double v

Returns:
	None
"
);
    void setMembPotential(std::string memb, double v);
	
	%feature("autodoc", 
"
Sets the specific membrane capacitance (in farad / m^2) of membrane with string identifier memb.
			 
Syntax::
			 
    setMembCapac(memb, cm)
			 
Arguments:
    * string memb
    * double cm
			 
Returns:
    None
"
);
    void setMembCapac(std::string memb, double cm);		 

	%feature("autodoc", 
"
Sets the bulk electrical resistivity (in ohm.m) of the volume conductor 
assocaited with membrane with string identifier memb.
			 
Syntax::
			 
    setMembVolRes(memb, ro)
			 
Arguments:
    * string memb
    * double ro
			 
Returns:
    None
"
);
    void setMembVolRes(std::string memb, double ro);
	
    
	%feature("autodoc", 
"
Sets the  electrical resistivity (in ohm/m^2) of the membrane with string identifier memb. Reversal potential vrev is required in Volts.
			 
Syntax::
			 
    setMembRes(memb, ro, vrev)
			 
Arguments:
    * string memb
    * double ro
    * double vrev
			 
Returns:
    None
"
);
    void setMembRes(std::string memb, double ro, double vrev);                     

%feature("autodoc",
"
Return the number of compartments in the solver.

Syntax::

    getNComps()

Arguments:
    None

Returns:
    integer
"
);
    uint getNComps(void) const;

%feature("autodoc",
"
Return the number of patches in the solver.

Syntax::

    getNPatches()

Arguments:
    None

Returns:
    integer
"
);
    uint getNPatches(void) const;
    
	%feature("autodoc",
"
Return the name of a compartment by its index in the solver.
    
Syntax::
    
    getCompName(c_idx)
    
Arguments:
    * integer c_idx
    
Returns:
    string
"
);
    std::string getCompName(uint c_idx) const;
    
	%feature("autodoc",
"
Return the name of a patch by its index in the solver.

Syntax::

    getPatchName(p_idx)

Arguments:
    * integer p_idx

Returns:
    string
"
);
    std::string getPatchName(uint p_idx) const;
    
	%feature("autodoc",
"
Get number of species in a compartment.

Syntax::

    getNCompSpecs(c_idx)

Arguments:
    * integer c_idx

Returns:
    integer
"
);
    
    uint getNCompSpecs(uint c_idx) const;
    
	%feature("autodoc",
"
Get number of species in a patch.

Syntax::

    getNPatchSpecs(p_idx)

Arguments:
    * integer p_idx

Returns:
    integer
"
);
    uint getNPatchSpecs(uint p_idx) const;
    
	%feature("autodoc",
"
Get the name of a species in a compartment.

Syntax::

    getCompSpecName(c_idx, s_idx)

Arguments:
    * integer c_idx
    * integer s_idx

Returns:
    string
"
);
    std::string getCompSpecName(uint c_idx, uint s_idx) const;
    
	%feature("autodoc",
"
Get the name of a species in a patch.

Syntax::

    getPatchSpecName(p_idx, s_idx)

Arguments:
    * integer p_idx
    * integer s_idx

Returns:
    string
"
);
    std::string getPatchSpecName(uint p_idx, uint s_idx) const;
    
             ////////////////////////////////////////////////////////////////////////
             // Batch Data Access
             ////////////////////////////////////////////////////////////////////////
             
    %feature("autodoc",
"
Get the counts of a species s in a list of tetrahedrons.

Syntax::

    getBatchTetCounts(tets, s)

Arguments:
    * list<uint> tets
    * string s

Return:
    list<double>
"
);
    
    virtual std::vector<double> getBatchTetCounts(std::vector<uint> const & tets, std::string const & s) const;
    
    %feature("autodoc",
"
Get the counts of a species s in a list of triangles.

Syntax::

    getBatchTriCounts(tris, s)

Arguments:
    * list<uint> tris
    * string s

Return:
    list<double>
"
);
    
    virtual std::vector<double> getBatchTriCounts(std::vector<uint> const & tris, std::string const & s) const;
	
    
    %feature("autodoc",
"
Get the counts of a species s in a list of tetrahedrons.

Syntax::
    getBatchTetCountsNP(indices, s, counts)

Arguments:
    * numpy.array<uint> indices
    * string s
    * numpy.array<double, lengh = len(indices)>

Return:
    None
"
);
    virtual void getBatchTetCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const;
             
    %feature("autodoc",
"
Get the counts of a species s in a list of triangles.

Syntax::
    getBatchTriCountsNP(indices, s, counts)

Arguments:
    * numpy.array<uint> indices
    * string s
    * numpy.array<double, lengh = len(indices)>

Return:
    None
"
);
    virtual void getBatchTriCountsNP(unsigned int* indices, int input_size, std::string const & s, double* counts, int output_size) const;
    
             ////////////////////////////////////////////////////////////////////////
             // ROI Data Access
             ////////////////////////////////////////////////////////////////////////

    %feature("autodoc",
"
Get the counts of a species s in tetrehedrons of a ROI.

Syntax::

    getROITetCounts(ROI_id, s)

Arguments:
    * string ROI_id
    * string s

Return:
    list<double>
"
);
    virtual std::vector<double> getROITetCounts(std::string ROI_id, std::string const & s) const;

             
    %feature("autodoc",
"
Get the counts of a species s in triangles of a ROI.

Syntax::

    getROITriCounts(ROI_id, s)

Arguments:
    * string ROI_id
    * string s

Return:
    list<double>
"
);

    virtual std::vector<double> getROITriCounts(std::string ROI_id, std::string const & s) const;
	
    
    %feature("autodoc",
"
Get the counts of a species s in tetrehedrons of a ROI.

Syntax::

    getROITetCountsNP(ROI_id, s, counts)

Arguments:
    * string ROI_id
    * string s
    * numpy.array<double, lengh = len(indices)>

Return:
    None
"
);
    virtual void getROITetCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const;
    %feature("autodoc",
"
Get the counts of a species s in triangles of a ROI.

Syntax::

    getROITriCountsNP(ROI_id, s, counts)

Arguments:
    * string ROI_id
    * string s
    * numpy.array<double, lengh = len(indices)>

Return:
    None
"
);
    virtual void getROITriCountsNP(std::string ROI_id, std::string const & s, double* counts, int output_size) const;

    %feature("autodoc",
"
Get the volume of a ROI.

Syntax::

    getROIVol(ROI_id)

Arguments:
    * string ROI_id

Return:
    double
"
);
 	virtual double getROIVol(std::string ROI_id) const;

    %feature("autodoc",
"
Get the area of a ROI.

Syntax::

    getROIArea(ROI_id)

Arguments:
    * string ROI_id

Return:
    double
"
);
	virtual double getROIArea(std::string ROI_id) const;
    
    %feature("autodoc",
"
Get the count of a species in a ROI.

Syntax::

    getROICount(ROI_id, s)

Arguments:
    * string ROI_id
    * string s

Return:
    double
"
);
 	virtual double getROICount(std::string ROI_id, std::string const & s) const;
    
    %feature("autodoc",
"
Set the count of a species in a ROI.

Syntax::

    setROICount(ROI_id, s, count)

Arguments:
    * string ROI_id
    * string s
    * double count

Return:
    None
"
);
    virtual void setROICount(std::string ROI_id, std::string const & s, double count);
    
    %feature("autodoc",
"
Get the amount of a species in a ROI.

Syntax::

    getROIAmount(ROI_id, s, count)

Arguments:
    * string ROI_id
    * string s

Return:
    double
"
);
 	virtual double getROIAmount(std::string ROI_id, std::string const & s) const;
    
    %feature("autodoc",
"
Get the concentration of a species in a ROI.

Syntax::

    getROIConc(ROI_id, s, count)

Arguments:
    * string ROI_id
    * string s

Return:
    double
"
);
	virtual double getROIConc(std::string ROI_id, std::string const & s) const;
    
    %feature("autodoc",
"
Set the concentration of a species in a ROI.

Syntax::

    setROIConc(ROI_id, s, conc)

Arguments:
    * string ROI_id
    * string s
    * double conc

Return:
    None
"
);
    virtual void setROIConc(std::string ROI_id, std::string const & s, double conc);
    
    %feature("autodoc",
"
Set a species in a ROI to be clamped or not. The count of species s in the ROI is clamped if
b is True, not clamped if b is False.

Syntax::

    setROIClamped(ROI_id, s, b)

Arguments:
    * string ROI_id
    * string s
    * boolean b

Return:
    None
"
);
	virtual void setROIClamped(std::string ROI_id, std::string const & s, bool b);
    
    %feature("autodoc",
"
Sets the macroscopic reaction constant of reaction with identifier string r
in a ROI with identifier string ROI_id to kf. The unit of the reaction constant
depends on the order of the reaction.

Note: The default value still comes from the steps.model description, so 
calling reset() will return the reaction constant to that value.

Syntax::

    setROIReacK(ROI_id, r, kf)

Arguments:
    * string ROI_id
    * string r
    * double kf

Return:
    None
"
);
	virtual void setROIReacK(std::string ROI_id, std::string const & r, double kf);
    
    %feature("autodoc",
"
Sets the macroscopic reaction constant of surface reaction with identifier string sr
in a ROI with identifier string ROI_id to kf. The unit of the reaction constant
depends on the order of the reaction.

Note: The default value still comes from the steps.model description, so 
calling reset() will return the reaction constant to that value.

Syntax::

    setROISReacK(ROI_id, sr, kf)

Arguments:
    * string ROI_id
    * string sr
    * double kf

Return:
    None
"
);
  	virtual void setROISReacK(std::string ROI_id, std::string const & sr, double kf);
    
    %feature("autodoc",
"
Sets the macroscopic diffusion constant of diffusion with identifier string d
in a ROI with identifier string ROI_id to dk.

Note: The default value still comes from the steps.model description, so 
calling reset() will return the diffusion constant to that value.

Syntax::

    setROIDiffD(ROI_id, d, dk)

Arguments:
    * string ROI_id
    * string d
    * double dk

Return:
    None
"
);
	virtual void setROIDiffD(std::string ROI_id, std::string const & d, double dk);
    
    %feature("autodoc",
"
Set reaction r in a ROI to be active or not.

Syntax::

    setROIReacActive(ROI_id, r, a)

Arguments:
    * string ROI_id
    * string r
    * bool a

Return:
    None
"
);
	virtual void setROIReacActive(std::string ROI_id, std::string const & r, bool a);
    
    %feature("autodoc",
"
Set surface reaction sr in a ROI to be active or not.

Syntax::

    setROISReacActive(ROI_id, sr, a)

Arguments:
    * string ROI_id
    * string sr
    * bool a

Return:
    None
"
);
 	virtual void setROISReacActive(std::string ROI_id, std::string const & sr, bool a);

    %feature("autodoc",
"
Set diffusion d in a ROI to be active or not.

Syntax::

    setROIDiffActive(ROI_id, sr, a)

Arguments:
    * string ROI_id
    * string sr
    * bool a

Return:
    None
"
);
	virtual void setROIDiffActive(std::string ROI_id, std::string const & d, bool a);
    
    %feature("autodoc",
"
Set voltage dependent surface reaction vsr in a ROI to be active or not.

Syntax::

    setROIVDepSReacActive(ROI_id, vsr, a)

Arguments:
    * string ROI_id
    * string vsr
    * bool a

Return:
    None
"
);
    virtual void setROIVDepSReacActive(std::string ROI_id, std::string const & vsr, bool a);

    %feature("autodoc",
"
Return the extent of reaction with identifier string r in ROI with
identifier string ROI_id, that is the number of times the reaction has occurred up
to the current simulation time.

Syntax::

    getROIReacExtent(ROI_id, r)

Arguments:
    * string ROI_id
    * string r

Return:
    uint
"
);

    virtual uint getROIReacExtent(std::string ROI_id, std::string const & r) const;
    
    %feature("autodoc",
"
Reset the extent of reaction with identifier string r in ROI with
identifier string ROI_id, that is the number of times the reaction has occurred up
to the current simulation time, to 0.

Syntax::

    resetROIReacExtent(ROI_id, r)

Arguments:
    * string ROI_id
    * string r

Return:
    None
"
);
    virtual void resetROIReacExtent(std::string ROI_id, std::string const & r);
    
    %feature("autodoc",
"
Return the extent of surface reaction with identifier string sr in ROI with
identifier string ROI_id, that is the number of times the reaction has occurred up
to the current simulation time.

Syntax::

    getROISReacExtent(ROI_id, sr)

Arguments:
    * string ROI_id
    * string sr

Return:
    uint
"
);
    virtual uint getROISReacExtent(std::string ROI_id, std::string const & sr) const;
    
    %feature("autodoc",
"
Reset the extent of surface reaction with identifier string r in ROI with
identifier string ROI_id, that is the number of times the reaction has occurred up
to the current simulation time, to 0.

Syntax::

    resetROISReacExtent(ROI_id, r)

Arguments:
    * string ROI_id
    * string sr

Return:
    None
"
);
    virtual void resetROISReacExtent(std::string ROI_id, std::string const & sr);
    
    %feature("autodoc",
"
Return the extent of diffusion with identifier string d in ROI with
identifier string ROI_id, that is the number of times the diffusion has occurred up
to the current simulation time.

Syntax::

    getROIDiffExtent(ROI_id, d)

Arguments:
    * string ROI_id
    * string d

Return:
    uint
"
);
    virtual uint getROIDiffExtent(std::string ROI_id, std::string const & d) const;
    
    %feature("autodoc",
"
Reset the extent of diffusion with identifier string d in ROI with
identifier string ROI_id, that is the number of times the diffusion has occurred up
to the current simulation time, to 0.

Syntax::

    resetROIDiffExtent(ROI_id, d)

Arguments:
    * string ROI_id
    * string d

Return:
    None
"
);
    virtual void resetROIDiffExtent(std::string ROI_id, std::string const & d);

};  

////////////////////////////////////////////////////////////////////////////////
	
} // end namespace solver
} // end namespace steps

////////////////////////////////////////////////////////////////////////////////

%feature("autodoc", "1");

namespace steps
{
namespace wmrk4
{
	
class Wmrk4 : public steps::solver::API
{
	
public:
	%feature("autodoc", "1");
	Wmrk4(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);
    %feature("autodoc", "1");
	~Wmrk4(void);
    
    %feature("autodoc", 
"
Returns a string of the solver's name.

Syntax::
    
    getSolverName()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverName(void) const;
    
    %feature("autodoc", 
"
Returns a string giving a short description of the solver.

Syntax::
    
    getSolverDesc()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverDesc(void) const;
    
    %feature("autodoc", 
"
Returns a string of the solver authors names.

Syntax::
    
    getSolverAuthors()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverAuthors(void) const;
    
    %feature("autodoc", 
"
Returns a string giving the author's email address.

Syntax::
    
    getSolverEmail()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverEmail(void) const;

    %feature("autodoc", 
"
Checkpoint data to a file.
    
Syntax::
    
    checkpoint(file_name)
    
Arguments:
    string file_name
    
Return:
    None
");
    virtual void checkpoint(std::string const & file_name);
    
    %feature("autodoc", 
"
Restore data from a file.
    
Syntax::
    
    restore(file_name)
    
Arguments:
    string file_name
    
Return:
    None
");
    virtual void restore(std::string const & file_name);
    
    %feature("autodoc", 
"
Reset the simulation to the state the solver was initialised to. 
Typically, this resets all concentrations of all chemical species in 
all elements (whether compartments and patches in a well-mixed solver 
or tetrahedrons and triangles in a mesh-based solver) to zero, 
resets the simulation time to zero and resets reaction (and diffusion) 
rates to the default values described in the steps.model objects. 
All reaction (and diffusion) rules are reset to active and all 
compartment volumes and patch areas are reset to default values 
described in steps.geom objects (for well-mixed solvers). 
Usually, this method should be called before starting each simulation iteration.

Syntax::
    
    reset()
    
Arguments:
    None

Return:
    None
");
    virtual void reset(void);
    
    %feature("autodoc", 
"
Advance the simulation until endtime (given in seconds) is reached. 
The endtime must be larger or equal to the current simulation time.

Syntax::
    
    run(endtime)
    
Arguments:
    float endtime

Return:
    None
");
    virtual void run(double endtime);
    
    %feature("autodoc", 
"
Returns the current simulation time in seconds.

Syntax::
    
    getTime()
    
Arguments:
    None

Return:
    float
");
    virtual double getTime(void) const;


	////////////////////////////////////////////////////////////////////////			
		
};

////////////////////////////////////////////////////////////////////////////////
	
} // end namespace wmrk4
} // end namespace steps

////////////////////////////////////////////////////////////////////////////////

namespace steps
{
namespace wmdirect
{

class Wmdirect : public steps::solver::API
{

public:
    %feature("autodoc", "1");
	Wmdirect(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);
	%feature("autodoc", "1");
    ~Wmdirect(void);
     %feature("autodoc", 
"
Returns a string of the solver's name.

Syntax::
    
    getSolverName()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverName(void) const;
    
    %feature("autodoc", 
"
Returns a string giving a short description of the solver.

Syntax::
    
    getSolverDesc()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverDesc(void) const;
    
    %feature("autodoc", 
"
Returns a string of the solver authors names.

Syntax::
    
    getSolverAuthors()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverAuthors(void) const;
    
    %feature("autodoc", 
"
Returns a string giving the author's email address.

Syntax::
    
    getSolverEmail()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverEmail(void) const;

    %feature("autodoc", 
    "
    Checkpoint data to a file.
    
    Syntax::
    
    checkpoint(file_name)
    
    Arguments:
    string file_name
    
    Return:
    None
    ");
    virtual void checkpoint(std::string const & file_name);
    
    %feature("autodoc", 
"
Restore data from a file.
    
Syntax::
    
    restore(file_name)
    
Arguments:
    string file_name
    
Return:
    None
");
    virtual void restore(std::string const & file_name);
    
    %feature("autodoc", 
"
Reset the simulation to the state the solver was initialised to. 
Typically, this resets all concentrations of all chemical species in 
all elements (whether compartments and patches in a well-mixed solver 
or tetrahedrons and triangles in a mesh-based solver) to zero, 
resets the simulation time to zero and resets reaction (and diffusion) 
rates to the default values described in the steps.model objects. 
All reaction (and diffusion) rules are reset to active and all 
compartment volumes and patch areas are reset to default values 
described in steps.geom objects (for well-mixed solvers). 
Usually, this method should be called before starting each simulation iteration.

Syntax::
    
    reset()
    
Arguments:
    None

Return:
    None
");
    virtual void reset(void);
    
    %feature("autodoc", 
"
Advance the simulation until endtime (given in seconds) is reached. 
The endtime must be larger or equal to the current simulation time.

Syntax::
    
    run(endtime)
    
Arguments:
    float endtime

Return:
    None
");
    virtual void run(double endtime);
    
    %feature("autodoc", 
"
Returns the current simulation time in seconds.

Syntax::
    
    getTime()
    
Arguments:
    None

Return:
    float
");
    virtual double getTime(void) const;


	
	////////////////////////////////////////////////////////////////////////			
	
};
		
////////////////////////////////////////////////////////////////////////////////

} // end namespace wmdirect
} // end namespace steps

////////////////////////////////////////////////////////////////////////////////

namespace steps
{
namespace tetexact
{

class Tetexact : public steps::solver::API
{	

public:
    %feature("autodoc", "1");
    Tetexact(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r, int calcMembPot = EF_NONE);
    %feature("autodoc", "1");
    ~Tetexact(void);
    %feature("autodoc", 
"
Returns a string of the solver's name.

Syntax::
    
    getSolverName()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverName(void) const;
    
    %feature("autodoc", 
"
Returns a string giving a short description of the solver.

Syntax::
    
    getSolverDesc()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverDesc(void) const;
    
    %feature("autodoc", 
"
Returns a string of the solver authors names.

Syntax::
    
    getSolverAuthors()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverAuthors(void) const;
    
    %feature("autodoc", 
"
Returns a string giving the author's email address.

Syntax::
    
    getSolverEmail()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverEmail(void) const;

    %feature("autodoc", 
"
Checkpoint data to a file.
    
Syntax::
    
    checkpoint(file_name)
    
Arguments:
    string file_name
    
Return:
    None
");
    virtual void checkpoint(std::string const & file_name);
    
    %feature("autodoc", 
"
Restore data from a file.
    
Syntax::
    
    restore(file_name)
    
Arguments:
    string file_name
    
Return:
    None
");
    virtual void restore(std::string const & file_name);
    
    %feature("autodoc", 
"
Reset the simulation to the state the solver was initialised to. 
Typically, this resets all concentrations of all chemical species in 
all elements (whether compartments and patches in a well-mixed solver 
or tetrahedrons and triangles in a mesh-based solver) to zero, 
resets the simulation time to zero and resets reaction (and diffusion) 
rates to the default values described in the steps.model objects. 
All reaction (and diffusion) rules are reset to active and all 
compartment volumes and patch areas are reset to default values 
described in steps.geom objects (for well-mixed solvers). 
Usually, this method should be called before starting each simulation iteration.

Syntax::
    
    reset()
    
Arguments:
    None

Return:
    None
");
    virtual void reset(void);
    
    %feature("autodoc", 
"
Advance the simulation until endtime (given in seconds) is reached. 
The endtime must be larger or equal to the current simulation time.

Syntax::
    
    run(endtime)
    
Arguments:
    float endtime

Return:
    None
");
    virtual void run(double endtime);
    
    %feature("autodoc", 
"
Returns the current simulation time in seconds.

Syntax::
    
    getTime()
    
Arguments:
    None

Return:
    float
");
    virtual double getTime(void) const;

%feature("autodoc", 
"
Saves the vertex optimization in the Efield structure.
             
Syntax::
             
    saveMembOpt()
             
Arguments:
    string filename
             
Return:
    None
"
);
    void saveMembOpt(std::string const & opt_file_name);

};
	
////////////////////////////////////////////////////////////////////////////////
	
} // end namespace tetexact
} // end namespace steps


////////////////////////////////////////////////////////////////////////////////

namespace steps
{
namespace tetode
{

class TetODE : public steps::solver::API
{	

public:

    TetODE(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r, int calcMembPot = EF_NONE);
    ~TetODE(void);

    %feature("autodoc", 
"
Returns a string of the solver's name.

Syntax::

    getSolverName()

Arguments:
    None

Return:
    string
");
    virtual std::string getSolverName(void) const;

    %feature("autodoc", 
"
Returns a string giving a short description of the solver.

Syntax::

    getSolverDesc()

Arguments:
    None

Return:
    string
");
    virtual std::string getSolverDesc(void) const;

%feature("autodoc", 
"
Returns a string of the solver authors names.

Syntax::

    getSolverAuthors()

Arguments:
    None

Return:
    string
");
    virtual std::string getSolverAuthors(void) const;

    %feature("autodoc", 
"
Returns a string giving the author's email address.

Syntax::

    getSolverEmail()

Arguments:
    None

Return:
    string
");
    virtual std::string getSolverEmail(void) const;

    %feature("autodoc", 
"
Checkpoint data to a file.

Syntax::

    checkpoint(file_name)

Arguments:
    string file_name

Return:
    None
");
    virtual void checkpoint(std::string const & file_name);

    %feature("autodoc", 
"
Restore data from a file.

Syntax::

    restore(file_name)

Arguments:
    string file_name

Return:
    None
");
    virtual void restore(std::string const & file_name);

%feature("autodoc", 
"
Reset the simulation to the state the solver was initialised to. 
Typically, this resets all concentrations of all chemical species in 
all elements (whether compartments and patches in a well-mixed solver 
or tetrahedrons and triangles in a mesh-based solver) to zero, 
resets the simulation time to zero and resets reaction (and diffusion) 
rates to the default values described in the steps.model objects. 
All reaction (and diffusion) rules are reset to active and all 
compartment volumes and patch areas are reset to default values 
described in steps.geom objects (for well-mixed solvers). 
Usually, this method should be called before starting each simulation iteration.

Syntax::

    reset()

Arguments:
    None

Return:
    None
");

    %feature("autodoc", 
"
Reset the simulation to the state the solver was initialised to. 
Typically, this resets all concentrations of all chemical species in 
all elements (whether compartments and patches in a well-mixed solver 
or tetrahedrons and triangles in a mesh-based solver) to zero, 
resets the simulation time to zero and resets reaction (and diffusion) 
rates to the default values described in the steps.model objects. 
All reaction (and diffusion) rules are reset to active and all 
compartment volumes and patch areas are reset to default values 
described in steps.geom objects (for well-mixed solvers). 
Usually, this method should be called before starting each simulation iteration.

Syntax::

    reset()

Arguments:
    None

Return:
    None
");
    virtual void reset(void);

    %feature("autodoc", 
"
Advance the simulation until endtime (given in seconds) is reached. 
The endtime must be larger or equal to the current simulation time.

Syntax::

    run(endtime)

Arguments:
    float endtime

Return:
    None
");
    virtual void run(double endtime);

    %feature("autodoc", 
"
Returns the current simulation time in seconds.

Syntax::

    getTime()

Arguments:
    None

Return:
    float
");
    virtual double getTime(void) const;
    
%feature("autodoc", 
"
Sets the maximum number of steps in CVODE per call to run().
Default is 10000 if this function is not called.
             
Syntax::
             
    setMaxNumSteps()
             
Arguments:
    unsigned integer maxn
             
Return:
    None
");
    void setMaxNumSteps(unsigned int maxn);

%feature("autodoc", 
"
Set the absolute tolerance and the relative tolerance for CVODE.
             
Syntax::
             
    setTolerance(atol, rtol)
             
Arguments:
    float atol
    float rtol
             
Return:
    None
");
    void setTolerances(double atol, double rtol);


////////////////////////////////////////////////////////////////////////			

};

////////////////////////////////////////////////////////////////////////////////

} // end namespace tetode
} // end namespace steps
