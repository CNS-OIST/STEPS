////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009 Okinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006 University of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

%module model_swig

%include "python/std_map.i"
%include "python/std_string.i"
%include "python/std_vector.i"
%include "error.i"
%import "cpp/steps/common.h"

%{
#include "../cpp/steps/model/model.hpp"
#include "../cpp/steps/model/diff.hpp"
#include "../cpp/steps/model/reac.hpp"
#include "../cpp/steps/model/spec.hpp"
#include "../cpp/steps/model/sreac.hpp"
#include "../cpp/steps/model/surfsys.hpp"
#include "../cpp/steps/model/volsys.hpp"
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

}
}

////////////////////////////////////////////////////////////////////////////////

namespace std
{
	
%template(vector_str) vector<std::string>;
%template(vector_spc) vector<steps::model::Spec *>;
%template(vector_rec) vector<steps::model::Reac *>;
%template(vector_src) vector<steps::model::SReac *>;
%template(vector_dif) vector<steps::model::Diff *>;

}

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
	}
}

///////////////////////////////////////////////////////////////////////////////

namespace steps 
{ 
namespace model 
{

////////////////////////////////////////////////////////////////////////////////

bool isValidID(std::string const & id);
void checkID(std::string const & id);

////////////////////////////////////////////////////////////////////////////////


class Model
{
	
public:

	Model(void);
	~Model(void);
	
    %feature("autodoc", 
"
Returns a reference to the steps.model.Spec species object with 
identifier string spec_id (if defined).

Arguments:
    string spec_id
             
Return:
    steps.model.Spec
");
	Spec * getSpec(std::string const & id) const;
    
    %feature("autodoc", 
"
Remove the steps.model.Spec species object with identifier 
string spec_id (if defined) from the model.

Arguments:
    string spec_id

Return:
    None
");
	void delSpec(std::string const & id);
    
    %feature("autodoc", 
"
Returns a list of steps.model.Spec object references of all species in the model.

Arguments:
    None
             
Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getAllSpecs(void) const;
    
    %feature("autodoc", 
"
Returns a reference to the steps.model.Volsys volume system object with 
identifier string vsys_id (if defined).

Arguments:
    string vsys_id

Return:
    steps.model.Volsys
");
	Volsys * getVolsys(std::string const & id) const;
    
    %feature("autodoc", 
"
Remove the steps.model.Volsys volume system object with identifier string 
vsys_id (if defined) from the model.

Arguments:
    string vsys_id

Return:
    None
");
	void delVolsys(std::string const & id);
    
    %feature("autodoc", 
"
Returns a reference to the steps.model.Surfsys surface system object with 
identifier string ssys_id (if defined).

Arguments:
    string ssys_id

Return:
    steps.model.Surfsys
");
	Surfsys * getSurfsys(std::string const & id) const;
    
    %feature("autodoc", 
"
Remove the steps.model.Surfsys surface system object with identifier string 
ssys_id (if defined) from the model.

Arguments:
    string ssys_id

Return:
    None
");
	void delSurfsys(std::string const & id);
	
};

////////////////////////////////////////////////////////////////////////////////

class Spec
{
	
public:
	
	Spec(std::string const & id, Model * model);
	~Spec(void);

    %feature("autodoc", 
"
Get the identifier string of the species.

Arguments:
    None

Return:
    string
");
	std::string getID(void) const;
    
    %feature("autodoc", 
"
Set the identifier string of the species.

Arguments:
    string name

Return:
    None
");
	void setID(std::string const & id);

    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Arguments:
    None

Return:
    steps.model.Model
    
Attribute:
    model
");	
	Model * getModel(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

class Surfsys
{
	
public:
	
	Surfsys(std::string const & id, Model * model);
	~Surfsys(void);

    %feature("autodoc", 
"
Get the identifier string of the surface system.

Arguments:
    None

Return:
    string
");		
	std::string getID(void) const;
    
    %feature("autodoc", 
"
Set the identifier string of the surface system.

Arguments:
    string name

Return:
    None
");	
	void setID(std::string const & id);

    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Arguments:
    None

Return:
    steps.model.Model
");		
	Model * getModel(void) const;

    %feature("autodoc", 
"
Returns a reference to the steps.model.SReac surface-reaction object 
with identifier sreac_id (if defined in the surface system.

Arguments:
    string sreac_id

Return:
    steps.model.SReac
");		
	SReac * getSReac(std::string const & id) const;
    
    %feature("autodoc", 
"
Remove the steps.model.SReac surface-reaction object with identifier 
sreac_id from the surface system.

Arguments:
    string sreac_id

Return:
    None
");	
	void delSReac(std::string const & id);
    
    %feature("autodoc", 
"
Returns a list of references to all steps.model.SReac surface-reaction 
objects defined in the surface system.

Arguments:
    None

Return:
    list<steps.model.SReac>
");	
	std::vector<SReac *> getAllSReacs(void) const;

    %feature("autodoc", 
"
Returns a list of references to all steps.model.Spec species objects included 
in the surface system; that is all reactants and products in the surface 
reactions belonging to this surface system. No duplicate member is included.

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
    
    %feature("autodoc", 
"
Get the identifier string of the volume system.

Arguments:
    None

Return:
    string
");		
	std::string getID(void) const;
    
    %feature("autodoc", 
"
Set the identifier string of the volume system.

Arguments:
    string name

Return:
    None
");	
	void setID(std::string const & id);

    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Arguments:
    None

Return:
    steps.model.Model
");
	Model * getModel(void) const;
	
    %feature("autodoc", 
"
Returns a reference to the steps.model.Reac reaction-rule object with 
identifier string reac_id (if defined in the volume system).

Arguments:
    string reac_id

Return:
    steps.model.Reac
");
	Reac * getReac(std::string const & id) const;
    
    %feature("autodoc", 
"
Remove the steps.model.Reac reaction-rule object with identifier reac_id 
(if defined) from the volume system.

Arguments:
    string reac_id

Return:
    None
");
	void delReac(std::string const & id);
    
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Spec objects in this volume 
system; that is all reactants, products or diffusing species in the reaction 
and diffusion rules belonging to this volume system. No duplicate member is 
included.

Arguments:
    None

Return:
    list<steps.model.Reac>
");
	std::vector<Reac *> getAllReacs(void) const;
	
    %feature("autodoc", 
"
Returns a reference to the steps.model.Diff diffusion-rule object with 
identifier diff_id (if defined in the volume system).

Arguments:
    string diff_id

Return:
    steps.model.Diff
");
	Diff * getDiff(std::string const & id) const;
    
    %feature("autodoc", 
"
Remove the steps.model.Diff diffusion-rule object with identifier diff_id 
from the volume system.

Arguments:
    string diff_id

Return:
    None
");
	void delDiff(std::string const & id);
    
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Diff diffusion-rule objects 
defined in the volume system.

Arguments:
    None

Return:
    list<steps.model.Diff>
");
	std::vector<Diff *> getAllDiffs(void) const;
	
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Spec objects in this volume system; 
that is all reactants, products or diffusing species in the reaction and diffusion 
rules belonging to this volume system. No duplicate member is included.

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getAllSpecs(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

%feature("kwargs") Diff::Diff;

class Diff
{
public:
	
	Diff(std::string const & id, Volsys * volsys, Spec * lig, double dcst=0.0);
	~Diff(void);
	
    %feature("autodoc", 
"
Get the identifier string of the diffusion rule.

Arguments:
    None

Return:
    string
");
	std::string getID(void) const;
    
    %feature("autodoc", 
"
Set the identifier string of the diffusion rule.

Arguments:
    string name

Return:
    None
");
	void setID(std::string const & id);
	
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Volsys volume system object.

Arguments:
    None

Return:
    steps.model.Volsys
");
	Volsys * getVolsys(void) const;
	
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Arguments:
    None

Return:
    steps.model.Model
");
	Model * getModel(void) const;
	
    %feature("autodoc", 
"
get a reference to the steps.model.Spec species object to which this 
diffusion rule is applied.

Arguments:
    None

Return:
    steps.model.Spec
");
	Spec * getLig(void) const;
    
    %feature("autodoc", 
"
Set a reference to the steps.model.Spec species object to which this 
diffusion rule is applied.

Arguments:
    steps.model.Spec lig

Return:
    None
");
	void setLig(Spec * lig);
	
    %feature("autodoc", 
"
Get the diffusion constant for the diffusion rule, in s.i. units.

Arguments:
    None

Return:
    float
");
	double getDcst(void) const;
    
    %feature("autodoc", 
"
Set the diffusion constant for the diffusion rule, in s.i. units.

Arguments:
    float dcst

Return:
    None
");
	void setDcst(double dcst);
	
};

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
	
    %feature("autodoc", 
"
Get the identifier string of the reaction rule.

Arguments:
    None

Return:
    string
");
	std::string getID(void) const;

    %feature("autodoc", 
"
Set the identifier string of the reaction rule.

Arguments:
    string name

Return:
    None
");
	void setID(std::string const & id);
	
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Volsys volume system object.

Arguments:
    None

Return:
    steps.model.Volsys
");
	Volsys * getVolsys(void) const;
	
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Arguments:
    None

Return:
    steps.model.Model
");
	Model * getModel(void) const;
	
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects on the 
left hand side of the reaction: the reactants.

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getLHS(void) const;
    
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects on the 
left hand side of the reaction: the reactants.

Arguments:
    list<steps.model.Spec> lhs

Return:
    None
");
	void setLHS(std::vector<Spec *> const & lhs);
	
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects on the 
right hand side of the reaction: the reactants.

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getRHS(void) const;

    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects on the 
right hand side of the reaction: the reactants.

Arguments:
    list<steps.model.Spec> rhs

Return:
    None
");
	void setRHS(std::vector<Spec *> const & rhs);
	
    %feature("autodoc", 
"
Returns the order of this reaction.

Arguments:
    None

Return:
    int
");
	int getOrder(void) const;
	
    %feature("autodoc", 
"
Get the kinetic reaction rate constant, in s.i. units, 
where the actual units depend on the order of the reaction.

Arguments:
    None

Return:
    float
");
	double getKcst(void) const;
    
    %feature("autodoc", 
"
Set the kinetic reaction rate constant, in s.i. units, 
where the actual units depend on the order of the reaction.

Arguments:
    float kcst

Return:
    None
");
	void setKcst(double kcst);
	
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Spec species objects in 
the reaction; that is all reactants and products. No duplicate member 
is included.

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getAllSpecs(void) const;
};


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
	
    %feature("autodoc", 
"
Get the identifier string of the surface reaction rule.

Arguments:
    None

Return:
    string
");
	std::string getID(void) const;
    
    %feature("autodoc", 
"
Set the identifier string of the surface reaction rule.

Arguments:
    string name

Return:
    None
");
	void setID(std::string const & id);
	
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Surfsys surface system object.

Arguments:
    None

Return:
    steps.model.Surfsys
");
	Surfsys * getSurfsys(void) const;
	
    %feature("autodoc", 
"
Returns a reference to the parent steps.model.Model container object.

Arguments:
    None

Return:
    steps.model.Model
");
	Model * getModel(void) const;
	
    %feature("autodoc", "Obsolete");
	bool getInner(void) const;
    %feature("autodoc", "Obsolete");
	bool getOuter(void) const;
	
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the left hand side outer volume reactants.

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getOLHS(void) const;

    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the left hand side outer volume reactants.

Arguments:
    list<steps.model.Spec) olhs

Return:
    None
");
	void setOLHS(std::vector<Spec *> const & olhs);

    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the left hand side inner volume reactants.

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getILHS(void) const;
    
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the left hand side inner volume reactants.

Arguments:
    list<steps.model.Spec> ilhs

Return:
    None
");
	void setILHS(std::vector<Spec *> const & ilhs);
	
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the left hand side surface reactants.

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getSLHS(void) const;

    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the left hand side surface reactants.

Arguments:
    list<steps.model.Spec> slhs

Return:
    None
");
	void setSLHS(std::vector<Spec *> const & slhs);
	
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side inner volume reactants.

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getIRHS(void) const;
    
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the right hand side inner volume reactants.

Arguments:
    list<steps.model.Spec> irhs

Return:
    None
");
	void setIRHS(std::vector<Spec *> const & irhs);
	
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side surface reactants.

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getSRHS(void) const;
    
    %feature("autodoc", 
"
Set a list of references to steps.model.Spec species objects; 
the right hand side surface reactants.

Arguments:
    list<steps.model.Spec> srhs

Return:
    None
");
	void setSRHS(std::vector<Spec *> const & srhs);
	
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side outer volume reactants.

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getORHS(void) const;
    
    %feature("autodoc", 
"
Get a list of references to steps.model.Spec species objects; 
the right hand side outer volume reactants.

Arguments:
    list<steps.model.Spec> orhs

Return:
    None
");
	void setORHS(std::vector<Spec *> const & orhs);
	
    %feature("autodoc", 
"
Returns the order of this surface reaction.

Arguments:
    None

Return:
    int
");
	int getOrder(void) const;
	
    %feature("autodoc", 
"
Get the kinetic reaction rate constant, in s.i. units, 
where the actual units depend on the order of the surface reaction.

Arguments:
    None

Return:
    float
");
	double getKcst(void) const;
    
    %feature("autodoc", 
"
Set the kinetic reaction rate constant, in s.i. units, 
where the actual units depend on the order of the surface reaction.

Arguments:
    float kcst

Return:
    None
");
	void setKcst(double kcst);
	
    %feature("autodoc", 
"
Returns a list of references to all steps.model.Spec species objects in 
the surface reaction; that is all reactants and products. No duplicate member 
is included.

Arguments:
    None

Return:
    list<steps.model.Spec>
");
	std::vector<Spec *> getAllSpecs(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

} // namespace model
} // namespace steps

// END
