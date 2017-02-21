/*
 ___license_placeholder___
 */


%module model_swig

%include "python/std_map.i"
%include "python/std_string.i"
%include "python/std_vector.i"
%include "error.i"
%import "steps/common.h"

%{
#include "steps/model/model.hpp"
#include "steps/model/diff.hpp"
#include "steps/model/reac.hpp"
#include "steps/model/spec.hpp"
#include "steps/model/sreac.hpp"
#include "steps/model/surfsys.hpp"
#include "steps/model/volsys.hpp"
#include "steps/model/vdeptrans.hpp"
#include "steps/model/vdepsreac.hpp"
#include "steps/model/ohmiccurr.hpp"
#include "steps/model/ghkcurr.hpp"
#include "steps/model/chan.hpp"
#include "steps/model/chanstate.hpp"
%}

////////////////////////////////////////////////////////////////////////////////

namespace steps 
{ 
namespace model 
{
	
class Model;
class Diff;
class Reac;
class Spec;
class SReac;
class Surfsys;
class Volsys;
class Chan;
class ChanState;	
class VDepTrans;
class VDepSReac;
class OhmicCurr;
class GHKcurr;
}
}

////////////////////////////////////////////////////////////////////////////////

namespace std
{
	

%template(vector_dbl) vector<double>;
%template(vector_chn) vector<steps::model::Chan *>;
%template(vector_cst) vector<steps::model::ChanState *>;
%template(vector_spc) vector<steps::model::Spec *>;
%template(vector_rec) vector<steps::model::Reac *>;
%template(vector_src) vector<steps::model::SReac *>;
%template(vector_dif) vector<steps::model::Diff *>;
%template(vector_vdt) vector<steps::model::VDepTrans *>;
%template(vector_vdsr) vector<steps::model::VDepSReac *>;
%template(vector_ohc) vector<steps::model::OhmicCurr *>;
%template(map_str_dbl) map<std::string, double>;
%template(vector_vsys) vector<steps::model::Volsys *>;
%template(vector_ssys) vector<steps::model::Surfsys *>;
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
%feature("autodoc", "1");


///////////////////////////////////////////////////////////////////////////////	

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

///////////////////////////////////////////////////////////////////////////////

namespace steps 
{ 
namespace model 
{

////////////////////////////////////////////////////////////////////////////////

//bool isValidID(std::string const & id);
//void checkID(std::string const & id);

////////////////////////////////////////////////////////////////////////////////


class Model
{

public:

	Model(void);
	~Model(void);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the steps.model.Spec species object with 
identifier string spec_id (if defined).

Syntax::

    getSpec(spec_id)

Arguments:
    string spec_id
             
Return:
    steps.model.Spec
");
	Spec * getSpec(std::string const & id) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Remove the steps.model.Spec species object with identifier 
string spec_id (if defined) from the model.

Syntax::

    delSpec(spec_id)

Arguments:
    string spec_id

Return:
    None
");
	void delSpec(std::string const & id);
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of steps.model.Spec object references of all species in the model.

Syntax::

    getAllSpecs()
    
Arguments:
    None
             
Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getAllSpecs(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the steps.model.Chan channel object with
identifier string chan_id (if defined).

Syntax::

    getSpec(chan_id)

Arguments:
    string chan_id
             
Return:
    steps.model.Chan
");
	Chan * getChan(std::string const & id) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of steps.model.Chan object references of all channels in the model.

Syntax::

    getAllChans()
    
Arguments:
    None
             
Return:
    list<steps.model.Chan>
");
	std::vector<Chan *> getAllChans(void) const;

    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the steps.model.Volsys volume system object with 
identifier string vsys_id (if defined).

Syntax::

    getVolsys(vsys_id)

Arguments:
    string vsys_id

Return:
    steps.model.Volsys
");
	Volsys * getVolsys(std::string const & id) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Remove the steps.model.Volsys volume system object with identifier string 
vsys_id (if defined) from the model.

Syntax::

    delVolsys(vsys_id)

Arguments:
    string vsys_id

Return:
    None
");
	void delVolsys(std::string const & id);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of steps.model.Volsys object references of all volume systems in the model.

Syntax::

    getAllVolsyss()
    
Arguments:
    None
             
Return:
    list<steps.model.Volsys>
");
	std::vector<Volsys *> getAllVolsyss(void) const;
	    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the steps.model.Surfsys surface system object with 
identifier string ssys_id (if defined).

Syntax::

    getSurfsys(ssys_id)

Arguments:
    string ssys_id

Return:
    steps.model.Surfsys
");
	Surfsys * getSurfsys(std::string const & id) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Remove the steps.model.Surfsys surface system object with identifier string 
ssys_id (if defined) from the model.

Syntax::

    delSurfsys(ssys_id)

Arguments:
    string ssys_id

Return:
    None
");
	void delSurfsys(std::string const & id);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of steps.model.Surfsys object references of all surface systems in the model.

Syntax::

    getAllSurfsyss()
    
Arguments:
    None
             
Return:
    list<steps.model.Surfsys>
");
	std::vector<Surfsys *> getAllSurfsyss(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

class Spec
{
	
public:
	
	Spec(std::string const & id, Model * model, int valence = 0);
	~Spec(void);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the identifier string of the species.

Syntax::

    getID()

Arguments:
    None

Return:
    string
");
	std::string getID(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set the identifier string of the species.

Syntax::

    setID(name)

Arguments:
    string name

Return:
    None
");
	void setID(std::string const & id);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Syntax::

    getModel()
    
Arguments:
    None

Return:
    steps.model.Model
    
Attribute:
    model
");	
	Model * getModel(void) const;
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc",

"
Set the valence of the species.
			 
Syntax::
			 
    setValence(valence)

Arguments:
	int valence

Return:
	None
");
	void setValence(int valence);	
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc",
"
Returns the valence of the species.
	
Syntax::
	
	getValence()
			 
Arguments:
	None
		
Return:
	int
");
	int getValence(void) const;	
};

	
////////////////////////////////////////////////////////////////////////////////

class Chan
{
	
public:
	
	Chan(std::string const & id, Model * model);
	~Chan(void);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the identifier string of the channel.

Syntax::

    getID()

Arguments:
	None
			 
Return:
    string
");		
	std::string getID(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set the identifier string of the channel.

Syntax::

    setID(name)

Arguments:
    string name
			 
Return:
    None
");	
	void setID(std::string const & id);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Syntax::
			 
    getModel()
			 
Arguments:
    None
			 
Return:
    steps.model.Model
			 
Attribute:
    model
");		
	Model * getModel(void) const;

////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc",
"
Returns a reference to channel state of the channel with string identifier id.
	
Syntax::
			 
	getChanState(id)
	
Arguments:
	string id
		
Return:
	steps.model.Chanstate
");
	ChanState * getChanState(std::string const & id) const;	
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc",
"
Returns a list of steps.model.Chanstate object references of all channel states in the channel.
			 
Syntax::
	
	getAllChanStates()

Arguments:
	None
		
Return:
	list<steps.model.ChanState>
");
	std::vector<ChanState *> getAllChanStates(void) const;
	
};	

////////////////////////////////////////////////////////////////////////////////

class ChanState : public steps::model::Spec
{
	
public:
	
	ChanState(std::string const & id, Model * model, Chan * chan);
	~ChanState(void);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set the identifier string of the channel state.
			 
Syntax::
			 
    setID(name)
			 
Arguments:
    string name
			 
Return:
    None
");		
	void setID(std::string const & id);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Chan container object.
			 
Syntax::
			 
    getChan()
			 
Arguments:
    None
			 
Return:
    steps.model.Chan
");		
	Chan * getChan(void) const;
	
};	

////////////////////////////////////////////////////////////////////////////////

class Surfsys
{
	
public:
	
	Surfsys(std::string const & id, Model * model);
	~Surfsys(void);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the identifier string of the surface system.

Syntax::
    
    getID()

Arguments:
    None

Return:
    string
");		
	std::string getID(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set the identifier string of the surface system.

Syntax::

    setID(name)
    
Arguments:
    string name

Return:
    None
");	
	void setID(std::string const & id);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Syntax::

    getModel()

Arguments:
    None

Return:
    steps.model.Model
");		
	Model * getModel(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the steps.model.SReac surface-reaction object 
with identifier sreac_id (if defined in the surface system).

Syntax::

    getSReac(id)

Arguments:
    string id

Return:
    steps.model.SReac
");		
	SReac * getSReac(std::string const & id) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Remove the steps.model.SReac surface-reaction object with identifier 
id from the surface system.

Syntax::

    delSReac(id)
    
Arguments:
    string id

Return:
    None
");	
	void delSReac(std::string const & id);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of references to all steps.model.SReac surface-reaction 
objects defined in the surface system.

Syntax::

    getAllSReacs()

Arguments:
    None

Return:
    list<steps.model.SReac>
");	
	std::vector<SReac *> getAllSReacs(void) const;

    
////////////////////////////////////////////////////////////////////////////////
%feature("autodoc", 
"
Returns a reference to the steps.model.Diff diffusion-rule object with 
identifier diff_id (if defined in the surface system).

Syntax::

    getDiff(diff_id)

Arguments:
    string diff_id

Return:
    steps.model.Diff
");
    Diff * getDiff(std::string const & id) const;

////////////////////////////////////////////////////////////////////////////////
%feature("autodoc", 
"
Remove the steps.model.Diff diffusion-rule object with identifier diff_id 
from the surface system.

Syntax::

    delDiff(diff_id)

Arguments:
    string diff_id

Return:
    None
");
    void delDiff(std::string const & id);

////////////////////////////////////////////////////////////////////////////////
%feature("autodoc", 
"
Returns a list of references to all steps.model.Diff diffusion-rule objects 
defined in the surface system.

Syntax::

    getAllDiffs()

Arguments:
    None

Return:
    list<steps.model.Diff>
");
    std::vector<Diff *> getAllDiffs(void) const;
    
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Returns a reference to the steps.model.VDepTrans voltage-dependent transition object 
with identifier id (if defined in the surface system).
	
Syntax::
	
    getVDepTrans(id)
	
Arguments:
    string id
	
Return:
    steps.model.VDepTrans
");		
	VDepTrans * getVDepTrans(std::string const & id) const;
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Remove the steps.model.VDepTrans voltage-dependent transition object with identifier 
id from the surface system.
	
Syntax::
	
    delVDepTrans(id)
    
Arguments:
    string id
	
Return:
    None
");	
	void delVDepTrans(std::string const & id);
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Returns a list of references to all steps.model.VDepTrans voltage-dependent transition 
objects defined in the surface system.
			 
Syntax::
			 
    getAllVDepTrans()
			 
Arguments:
    None
			 
Return:
    list<steps.model.VDepTrans>
");	
	std::vector<VDepTrans *> getAllVDepTrans(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 	
"
Returns a reference to the steps.model.OhmicCurr ohmic current object 
with identifier id (if defined in the surface system).
	
Syntax::
	
    getOhmicCurr(id)
	
Arguments:
    string id
	
Return:
    steps.model.OhmicCurr
");		
	OhmicCurr * getOhmicCurr(std::string const & id) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 	
"
Remove the steps.model.OhmicCurr ohmic current object with identifier 
id from the surface system.
	
Syntax::
	
    delOhmicCurr(id)
    
Arguments:
    string id
	
Return:
    None
");	
	void delOhmicCurr(std::string const & id);
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Returns a list of references to all steps.model.OhmicCurr ohmic current 
objects defined in the surface system.
			 
Syntax::
			 
    getAllOhmicCurrs()
			 
Arguments:
    None
			 
Return:
    list<steps.model.OhmicCurr>
");	
	std::vector<OhmicCurr *> getAllOhmicCurrs(void) const;
	
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Returns a reference to the steps.model.GHKcurr ghk current object 
with identifier id (if defined in the surface system).
			 
Syntax::
			 
    getGHKcurr(id)
			 
Arguments:
    string id
			 
Return:
    steps.model.GHKcurr
");		
	GHKcurr * getGHKcurr(std::string const & id) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 	
"
Remove the steps.model.GHKcurr ghk current object with identifier 
id from the surface system.
			 
Syntax::
			 
    delGHKcurr(id)
			 
Arguments:
    string id
			 
Return:
    None
");		
	void delGHKcurr(std::string const & id);
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Returns a list of references to all steps.model.GHKcurr ghk current 
objects defined in the surface system.
			 
Syntax::
			 
    getAllGHKcurrs()
			 
Arguments:
    None
			 
Return:
    list<steps.model.GHKcurrs>
");		
	std::vector<GHKcurr *> getAllGHKcurrs(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Spec species objects included 
in the surface system; that is all reactants and products in the surface 
reactions belonging to this surface system. No duplicate member is included.

Syntax::

    getAllSpecs()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");		
	std::vector<Spec *> getAllSpecs(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

class Volsys
{
	
public:
	
	Volsys(std::string const & id, Model * model);
	~Volsys(void);
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the identifier string of the volume system.

Syntax::

    getID()
    
Arguments:
    None

Return:
    string
");		
	std::string getID(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set the identifier string of the volume system.

Syntax::

    setID(name)
    
Arguments:
    string name

Return:
    None
");	
	void setID(std::string const & id);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Syntax::

    getModel()
    
Arguments:
    None

Return:
    steps.model.Model
");
	Model * getModel(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the steps.model.Reac reaction-rule object with 
identifier string reac_id (if defined in the volume system).

Syntax::

    getReac(reac_id)
    
Arguments:
    string reac_id

Return:
    steps.model.Reac
");
	Reac * getReac(std::string const & id) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Remove the steps.model.Reac reaction-rule object with identifier reac_id 
(if defined) from the volume system.

Syntax::

    delReac(reac_id)

Arguments:
    string reac_id

Return:
    None
");
	void delReac(std::string const & id);
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Reac objects in this volume 
system; that is all reaction rules belonging to this volume system. No duplicate 
member is included.

Syntax::

    getAllReacs()
    
Arguments:
    None

Return:
    list<steps.model.Reac>
");
	std::vector<Reac *> getAllReacs(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the steps.model.Diff diffusion-rule object with 
identifier diff_id (if defined in the volume system).

Syntax::

    getDiff(diff_id)
    
Arguments:
    string diff_id

Return:
    steps.model.Diff
");
	Diff * getDiff(std::string const & id) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Remove the steps.model.Diff diffusion-rule object with identifier diff_id 
from the volume system.

Syntax::

    delDiff(diff_id)
    
Arguments:
    string diff_id

Return:
    None
");
	void delDiff(std::string const & id);
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Diff diffusion-rule objects 
defined in the volume system.

Syntax::

    getAllDiffs()
    
Arguments:
    None

Return:
    list<steps.model.Diff>
");
	std::vector<Diff *> getAllDiffs(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Spec objects in this volume system; 
that is all reactants, products or diffusing species in the reaction and diffusion 
rules belonging to this volume system. No duplicate member is included.

Syntax::

    getAllSpecs()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getAllSpecs(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//%feature("kwargs") Diff::Diff;

class Diff
{
public:
	
	Diff(std::string const & id, Volsys * volsys, Spec * lig, double dcst=0.0);
	Diff(std::string const & id, Surfsys * surfsys, Spec * lig, double dcst=0.0);
    ~Diff(void);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the identifier string of the diffusion rule.

Syntax::

    getID()
    
Arguments:
    None

Return:
    string
");
	std::string getID(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set the identifier string of the diffusion rule.

Syntax::

    setID(name)
    
Arguments:
    string name

Return:
    None
");
	void setID(std::string const & id);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Volsys volume system object if a volume diffusion object.

Syntax::

    getVolsys()
    
Arguments:
    None

Return:
    steps.model.Volsys
");
	Volsys * getVolsys(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Surfsys surface system object if surface diffusion object.

Syntax::

    getSurfys()

Arguments:
    None

Return:
    steps.model.Surfsys
");
    Surfsys * getSurfsys(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Syntax::

    getModel()
    
Arguments:
    None

Return:
    steps.model.Model
");
	Model * getModel(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
get a reference to the steps.model.Spec species object to which this 
diffusion rule is applied.

Syntax::

    getLig()
    
Arguments:
    None

Return:
    steps.model.Spec
");
	Spec * getLig(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a reference to the steps.model.Spec species object to which this 
diffusion rule is applied.

Syntax::

    setLig(lig)
    
Arguments:
    steps.model.Spec lig

Return:
    None
");
	void setLig(Spec * lig);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the diffusion constant for the diffusion rule, in s.i. units.

Syntax::

    getDcst()
    
Arguments:
    None

Return:
    float
");
	double getDcst(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set the diffusion constant for the diffusion rule, in s.i. units.

Syntax::

    setDcst(dcst)
    
Arguments:
    float dcst

Return:
    None
");
	void setDcst(double dcst);
	
};

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
%feature("kwargs") Reac::Reac;

class Reac
{
public:
	
	Reac(std::string const & id, Volsys * volsys, 
		 std::vector<Spec *> const lhs = std::vector<Spec *>(),
		 std::vector<Spec *> const rhs = std::vector<Spec *>(), 
		 double kcst = 0.0);
	~Reac(void);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the identifier string of the reaction rule.

Syntax::

    getID()
    
Arguments:
    None

Return:
    string
");
	std::string getID(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set the identifier string of the reaction rule.

Syntax::

    setID(name)
    
Arguments:
    string name

Return:
    None
");
	void setID(std::string const & id);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Volsys volume system object.

Syntax::

    getVolsys()
    
Arguments:
    None

Return:
    steps.model.Volsys
");
	Volsys * getVolsys(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Syntax::

    getModel()
    
Arguments:
    None

Return:
    steps.model.Model
");
	Model * getModel(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects on the 
left hand side of the reaction: the reactants.

Syntax::

    getLHS()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getLHS(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects on the 
left hand side of the reaction: the reactants.

Syntax::

    setLHS(lhs)

Arguments:
    list<steps.model.Spec> lhs

Return:
    None
");
	void setLHS(std::vector<Spec *> const & lhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects on the 
right hand side of the reaction: the reactants.

Syntax::

    getRHS()

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getRHS(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects on the 
right hand side of the reaction: the reactants.

Syntax::

    setRHS(rhs)
    
Arguments:
    list<steps.model.Spec> rhs

Return:
    None
");
	void setRHS(std::vector<Spec *> const & rhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns the order of this reaction.

Syntax::

    getOrder()
    
Arguments:
    None

Return:
    int
");
	int getOrder(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the kinetic reaction rate constant, in s.i. units, 
where the actual units depend on the order of the reaction.

Syntax::
    
    getKcst()

Arguments:
    None

Return:
    float
");
	double getKcst(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set the kinetic reaction rate constant, in s.i. units, 
where the actual units depend on the order of the reaction.

Syntax::
    
    setKcst(kcst)

Arguments:
    float kcst

Return:
    None
");
	void setKcst(double kcst);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Spec species objects in 
the reaction; that is all reactants and products. No duplicate member 
is included.

Syntax::

    getAllSpecs()

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getAllSpecs(void) const;
};


////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
%feature("kwargs") SReac::SReac;

class SReac
{
public:
	
	SReac(std::string const & id, Surfsys * surfsys, 
		  std::vector<Spec *> const & olhs = std::vector<Spec *>(),
		  std::vector<Spec *> const & ilhs = std::vector<Spec *>(),
		  std::vector<Spec *> const & slhs = std::vector<Spec *>(), 
		  std::vector<Spec *> const & irhs = std::vector<Spec *>(),
		  std::vector<Spec *> const & srhs = std::vector<Spec *>(), 
		  std::vector<Spec *> const & orhs = std::vector<Spec *>(),
		  double kcst = 0.0);
	~SReac(void);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the identifier string of the surface reaction rule.

Syntax::

    getID()
    
Arguments:
    None

Return:
    string
");
	std::string getID(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set the identifier string of the surface reaction rule.

Syntax::

    setID(name)
    
Arguments:
    string name

Return:
    None
");
	void setID(std::string const & id);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Surfsys surface system object.

Syntax::

    getSurfsys()
    
Arguments:
    None

Return:
    steps.model.Surfsys
");
	Surfsys * getSurfsys(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object of parent surface system object.

Syntax::

    getModel()
    
Arguments:
    None

Return:
    steps.model.Model
");
	Model * getModel(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", "Obsolete");
	bool getInner(void) const;
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", "Obsolete");
	bool getOuter(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the left hand side outer volume reactants.

Syntax::

    getOLHS()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getOLHS(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the left hand side outer volume reactants.

Syntax::

    setOLHS(olhs)
    
Arguments:
    list<steps.model.Spec) olhs

Return:
    None
");
	void setOLHS(std::vector<Spec *> const & olhs);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the left hand side inner volume reactants.

Syntax::
    
    getILHS()

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getILHS(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the left hand side inner volume reactants.

Syntax::

    setILHS(ilhs)

Arguments:
    list<steps.model.Spec> ilhs

Return:
    None
");
	void setILHS(std::vector<Spec *> const & ilhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the left hand side surface reactants.

Syntax::

    getSLHS()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getSLHS(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the left hand side surface reactants.

Syntax::

    setSLHS(slhs)
    
Arguments:
    list<steps.model.Spec> slhs

Return:
    None
");
	void setSLHS(std::vector<Spec *> const & slhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side inner volume reactants.

Syntax::

    getIRHS()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getIRHS(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the right hand side inner volume reactants.

Syntax::

    setIRHS(irhs)
    
Arguments:
    list<steps.model.Spec> irhs

Return:
    None
");
	void setIRHS(std::vector<Spec *> const & irhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side surface reactants.

Syntax::
    
    getSRHS()

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getSRHS(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the right hand side surface reactants.

Syntax::

    setSRHS(srhs)
    
Arguments:
    list<steps.model.Spec> srhs

Return:
    None
");
	void setSRHS(std::vector<Spec *> const & srhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side outer volume reactants.

Syntax::

    getORHS()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getORHS(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side outer volume reactants.

Syntax::

    setORHS(orhs)
    
Arguments:
    list<steps.model.Spec> orhs

Return:
    None
");
	void setORHS(std::vector<Spec *> const & orhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns the order of this surface reaction.

Syntax::

    getOrder()
    
Arguments:
    None

Return:
    int
");
	int getOrder(void) const;
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the kinetic reaction rate constant, in s.i. units, 
where the actual units depend on the order of the surface reaction.

Syntax::

    getKcst()
    
Arguments:
    None

Return:
    float
");
	double getKcst(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set the kinetic reaction rate constant, in s.i. units, 
where the actual units depend on the order of the surface reaction.

Syntax::

    setKcst(kcst)
    
Arguments:
    float kcst

Return:
    None
");
	void setKcst(double kcst);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Spec species objects in 
the surface reaction; that is all reactants and products. No duplicate member 
is included.

Syntax::

    getAllSpecs()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getAllSpecs(void) const;
	
};

	
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
%feature("kwargs") VDepTrans::VDepTrans;

class VDepTrans
{
	
public:
	
	VDepTrans(std::string const & id, Surfsys * surfsys,
			  ChanState * src, ChanState * dst,
			  std::vector<double> ratetab, double vmin, double vmax, 
			  double dv, uint tablesize);
	~VDepTrans(void);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the identifier string of the voltage-dependent transition.
			 
Syntax::
			 
   getID()
			 
Arguments:
    None
			 
Return:
    string
");	
	std::string getID(void) const;
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Set the identifier string of the voltage-dependent transition.
			 
Syntax::
			 
    setID(name)
			 
Arguments:
    string name
			 
Return:
    None
");
	void setID(std::string const & id);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Surfsys surface system object.
			 
Syntax::
			 
    getSurfsys()
			 
Arguments:
    None
			 
Return:
    steps.model.Surfsys
");
	Surfsys * getSurfsys(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object of parent surface system object.

Syntax::

    getModel()
    
Arguments:
    None

Return:
    steps.model.Model
");	
	Model * getModel(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the steps.model.Chan container object of source and destination channel states.

Syntax::

    getChan()
    
Arguments:
    None

Return:
    steps.model.Chan
");
	Chan * getChan(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the 'source' (left-hand side) steps.model.ChanState object.

Syntax::

    getSrc()
    
Arguments:
    None

Return:
    steps.model.ChanState
");
	ChanState * getSrc(void) const;

////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Set the 'source' (left-hand side) channel state.
			 
Syntax::
			 
    setSrc(src)
			 
Arguments:
    steps.model.ChanState src
			 
Return:
    None
");
	void setSrc(ChanState * src);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the 'destination' (right-hand side) steps.model.ChanState object.

Syntax::

    getDst()
    
Arguments:
    None

Return:
    steps.model.ChanState
");
	ChanState * getDst(void) const;
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Set the 'destination' (right-hand side) channel state.
			 
Syntax::
			 
    setDst(dst)
			 
Arguments:
    steps.model.ChanState dst
			 
Return:
    None
");	
	void setDst(ChanState * dst);

////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc",
"
Return a list of transition rates in the default voltage range.
			 
Syntax::
    
    getRate()
			 
Arguments:
    None
	
Return:
	list<float>
");
	std::vector<double> getRate(void) const;
	
};	

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
%feature("kwargs") VDepSReac::VDepSReac;

class VDepSReac
{
	
public:
	
	VDepSReac(std::string const & id, Surfsys * surfsys,
			  std::vector<Spec *> const & olhs = std::vector<Spec *>(),
			  std::vector<Spec *> const & ilhs = std::vector<Spec *>(),
			  std::vector<Spec *> const & slhs = std::vector<Spec *>(),
			  std::vector<Spec *> const & irhs = std::vector<Spec *>(),
			  std::vector<Spec *> const & srhs = std::vector<Spec *>(),
			  std::vector<Spec *> const & orhs = std::vector<Spec *>(),
			  std::vector<double> ktab = std::vector<double>(),
			  double vmin = 0.0, double vmax = 0.0,
			  double dv = 0.0, uint tablesize = 0);
	~VDepSReac(void);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the identifier string of the voltage-dependent reaction.
			 
Syntax::
			 
   getID()
			 
Arguments:
    None
			 
Return:
    string
");	
	std::string getID(void) const;
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Set the identifier string of the voltage-dependent reaction.
			 
Syntax::
			 
    setID(name)
			 
Arguments:
    string name
			 
Return:
    None
");
	void setID(std::string const & id);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Surfsys surface system object.
			 
Syntax::
			 
    getSurfsys()
			 
Arguments:
    None
			 
Return:
    steps.model.Surfsys
");
	Surfsys * getSurfsys(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object of parent surface system object.

Syntax::

    getModel()
    
Arguments:
    None

Return:
    steps.model.Model
");	
	Model * getModel(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the left hand side outer volume reactants.

Syntax::

    getOLHS()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getOLHS(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the left hand side outer volume reactants.

Syntax::

    setOLHS(olhs)
    
Arguments:
    list<steps.model.Spec) olhs

Return:
    None
");
	void setOLHS(std::vector<Spec *> const & olhs);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the left hand side inner volume reactants.

Syntax::
    
    getILHS()

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getILHS(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the left hand side inner volume reactants.

Syntax::

    setILHS(ilhs)

Arguments:
    list<steps.model.Spec> ilhs

Return:
    None
");
	void setILHS(std::vector<Spec *> const & ilhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the left hand side surface reactants.

Syntax::

    getSLHS()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getSLHS(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the left hand side surface reactants.

Syntax::

    setSLHS(slhs)
    
Arguments:
    list<steps.model.Spec> slhs

Return:
    None
");
	void setSLHS(std::vector<Spec *> const & slhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side inner volume reactants.

Syntax::

    getIRHS()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getIRHS(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the right hand side inner volume reactants.

Syntax::

    setIRHS(irhs)
    
Arguments:
    list<steps.model.Spec> irhs

Return:
    None
");
	void setIRHS(std::vector<Spec *> const & irhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side surface reactants.

Syntax::
    
    getSRHS()

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getSRHS(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the right hand side surface reactants.

Syntax::

    setSRHS(srhs)
    
Arguments:
    list<steps.model.Spec> srhs

Return:
    None
");
	void setSRHS(std::vector<Spec *> const & srhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side outer volume reactants.

Syntax::

    getORHS()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getORHS(void) const;
    
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side outer volume reactants.

Syntax::

    setORHS(orhs)
    
Arguments:
    list<steps.model.Spec> orhs

Return:
    None
");
	void setORHS(std::vector<Spec *> const & orhs);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns the order of this voltage-dependent reaction.

Syntax::

    getOrder()
    
Arguments:
    None

Return:
    int
");
	int getOrder(void) const;

////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc",
"
Return a list of reaction 'constants' in the default voltage range.
			 
Syntax::
    
    getK()
			 
Arguments:
    None
	
Return:
	list<float>
");
	std::vector<double> getK(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Spec species objects in 
the voltage-dependent reaction; that is all reactants and products. No duplicate member 
is included.

Syntax::

    getAllSpecs()
    
Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getAllSpecs(void) const;
	
};	

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
%feature("kwargs") OhmicCurr::OhmicCurr;

class OhmicCurr
{
	
public:
	OhmicCurr(std::string const & id, Surfsys * surfsys, 
			  ChanState * chanstate, double erev, double g);
	~OhmicCurr(void);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the identifier string of the ohmic current.
			 
Syntax::
			 
   getID()
			 
Arguments:
    None
			 
Return:
    string
");	
	std::string getID(void) const;
			 
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Set the identifier string of the ohmic current.
			 
Syntax::
			 
    setID(name)
			 
Arguments:
    string name
			 
Return:
    None
");			 
	void setID(std::string const & id);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Surfsys surface system object.
			 
Syntax::
			 
    getSurfsys()
			 
Arguments:
    None
			 
Return:
    steps.model.Surfsys
");
	Surfsys * getSurfsys(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object of parent surface system object.

Syntax::

    getModel()
    
Arguments:
    None

Return:
    steps.model.Model
");				 
	Model * getModel(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the steps.model.ChanState channel state object.

Syntax::

    getChanState()
    
Arguments:
    None

Return:
    steps.model.ChanState
");			 
	ChanState * getChanState(void) const;
			 
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Set the channel state for this ohmic current.
			 
Syntax::
			 
    setChanState(chanstate)
			 
Arguments:
    steps.model.ChanState chanstate
			 
Return:
    None
");				 
	void setChanState(ChanState * chanstate);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns the reversal potential of this ohmic current in volts.

Syntax::

    getERev()
    
Arguments:
    None

Return:
    float
");				 
	double getERev(void) const;
			 
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Set the reveral potential for this ohmic current in volts.
			 
Syntax::
			 
    setERev(erev)
			 
Arguments:
    flaot erev
			 
Return:
    None
");				 
	void setERev(double erev);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns the single-channel conductance for this ohmic current in Siemens.

Syntax::

    getG()
    
Arguments:
    None

Return:
    float
");	
	double getG(void) const;
			 
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Set the single-channel conductance for this ohmic current in Siemens.
			 
Syntax::
			 
    setG(g)
			 
Arguments:
    flaot g
			 
Return:
    None
");	
	void setG(double g);
	
};

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
%feature("kwargs") GHKcurr::GHKcurr;
////////////////////////////////////////////////////////////////////////////////
%feature("kwargs") GHKcurr::setPInfo;	

class GHKcurr
{
	
public:
	GHKcurr(std::string const & id, Surfsys * surfsys,
			ChanState * chanstate, Spec * ion, bool computeflux = false,
            double virtual_oconc = -1, double vshift = 0.0);
	~GHKcurr(void);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Get the identifier string of the ghk current.
			 
Syntax::
			 
   getID()
			 
Arguments:
    None
			 
Return:
    string
");				 
	std::string getID(void) const;
			 
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Set the identifier string of the ghk current.
			 
Syntax::
			 
    setID(name)
			 
Arguments:
    string name
			 
Return:
    None
");	
	void setID(std::string const & id);

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Surfsys surface system object.
			 
Syntax::
			 
    getSurfsys()
			 
Arguments:
    None
			 
Return:
    steps.model.Surfsys
");			 
	Surfsys * getSurfsys(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object of parent surface system object.

Syntax::

    getModel()
    
Arguments:
    None

Return:
    steps.model.Model
");				 			 
	Model * getModel(void) const;

////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc", 
"
Returns a reference to the steps.model.ChanState channel state object.

Syntax::

    getChanState()
    
Arguments:
    None

Return:
    steps.model.ChanState
");				 
	ChanState * getChanState(void) const;
			 
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Set the channel state for this ghk current.
			 
Syntax::
			 
    setChanState(chanstate)
			 
Arguments:
    steps.model.ChanState chanstate
			 
Return:
    None
");		
	void setChanState(ChanState * chanstate);
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Returns a reference to steps.model.Spec object- the ion of this ghk current.

Syntax::

    getIon()
    
Arguments:
    None

Return:
    steps.model.Spec
");				 
	Spec * getIon(void) const;
	
////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc", 
"
Set the ion for this ghk current.
			 
Syntax::
			 
    setIon(ion)
			 
Arguments:
    steps.model.Spec ion
			 
Return:
    None
");			
	void setIon(Spec * ion);

	// something getPInfo(void) const;

////////////////////////////////////////////////////////////////////////////////
	%feature("autodoc",
"
Supply information from a channel measurement in order to find the permeability. 
A measured single-channel coonductance (in Siemens) should be supplied, along with the potential (in volts),
temperature (in Kelvins), the 'outer' concentration and 'inner' concentration of ion (in molar units).

Syntax::
	
	setPInfo(g, V, T, oconc, iconc)
	
Arguments:
	* float g (the conductance)
	* float V (the voltage)
	* float T (the temperature)
	* float oconc (the 'outer' concentration)
	* float iconc (the 'inner' concentration)

Return:
	None
");
	void setPInfo(double g, double V, double T, double oconc, double iconc);
	
////////////////////////////////////////////////////////////////////////////////
    %feature("autodoc",
"
Set the single-channel permeability directly (units: cubic meters / second). 

Syntax::

    setP(p)

Arguments:
    * float p (the single-channel permeability)

Return:
    None
");
	void setP(double p);
    
};

////////////////////////////////////////////////////////////////////////////////

} // namespace model
} // namespace steps

// END
