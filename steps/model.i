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

%module model_swig

%include "python/std_map.i"
%include "python/std_string.i"
%include "python/std_vector.i"
%include "error.i"
%import "cpp/steps/common.h"

%{
// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif
#include <steps/model/model.hpp>
#include <steps/model/diff.hpp>
#include <steps/model/reac.hpp>
#include <steps/model/spec.hpp>
#include <steps/model/sreac.hpp>
#include <steps/model/surfsys.hpp>
#include <steps/model/volsys.hpp>
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
	
	Spec * getSpec(std::string const & id) const;
	void delSpec(std::string const & id);
	std::vector<Spec *> getAllSpecs(void) const;
	
	Volsys * getVolsys(std::string const & id) const;
	void delVolsys(std::string const & id);
	
	Surfsys * getSurfsys(std::string const & id) const;
	void delSurfsys(std::string const & id);
	
};

////////////////////////////////////////////////////////////////////////////////

class Spec
{
	
public:
	
	Spec(std::string const & id, Model * model);
	~Spec(void);
	
	std::string getID(void) const;
	void setID(std::string const & id);
	
	Model * getModel(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

class Surfsys
{
	
public:
	
	Surfsys(std::string const & id, Model * model);
	~Surfsys(void);
	
	std::string getID(void) const;
	void setID(std::string const & id);
	
	Model * getModel(void) const;
	
	SReac * getSReac(std::string const & id) const;
	void delSReac(std::string const & id);
	std::vector<SReac *> getAllSReacs(void) const;
	
	std::vector<Spec *> getAllSpecs(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

class Volsys
{
	
public:
	
	Volsys(std::string const & id, Model * model);
	~Volsys(void);
	
	std::string getID(void) const;
	void setID(std::string const & id);
	
	Model * getModel(void) const;
	
	Reac * getReac(std::string const & id) const;
	void delReac(std::string const & id);
	std::vector<Reac *> getAllReacs(void) const;
	
	Diff * getDiff(std::string const & id) const;
	void delDiff(std::string const & id);
	std::vector<Diff *> getAllDiffs(void) const;
	
	std::vector<Spec *> getAllSpecs(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

%feature("kwargs") Diff::Diff;

class Diff
{
public:
	
	Diff(std::string const & id, Volsys * volsys, Spec * lig, double dcst=0.0);
	~Diff(void);
	
	std::string getID(void) const;
	void setID(std::string const & id);
	
	Volsys * getVolsys(void) const;
	
	Model * getModel(void) const;
	
	Spec * getLig(void) const;
	void setLig(Spec * lig);
	
	double getDcst(void) const;
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
	
	std::string getID(void) const;
	void setID(std::string const & id);
	
	Volsys * getVolsys(void) const;
	
	Model * getModel(void) const;
	
	std::vector<Spec *> getLHS(void) const;
	void setLHS(std::vector<Spec *> const & lhs);
	
	std::vector<Spec *> getRHS(void) const;
	void setRHS(std::vector<Spec *> const & rhs);
	
	int getOrder(void) const;
	
	double getKcst(void) const;
	void setKcst(double kcst);
	
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
	
	std::string getID(void) const;
	void setID(std::string const & id);
	
	Surfsys * getSurfsys(void) const;
	
	Model * getModel(void) const;
	
	bool getInner(void) const;
	bool getOuter(void) const;
	
	std::vector<Spec *> getOLHS(void) const;
	void setOLHS(std::vector<Spec *> const & olhs);

	std::vector<Spec *> getILHS(void) const;
	void setILHS(std::vector<Spec *> const & ilhs);
	
	std::vector<Spec *> getSLHS(void) const;
	void setSLHS(std::vector<Spec *> const & slhs);
	
	std::vector<Spec *> getIRHS(void) const;
	void setIRHS(std::vector<Spec *> const & irhs);
	
	std::vector<Spec *> getSRHS(void) const;
	void setSRHS(std::vector<Spec *> const & srhs);
	
	std::vector<Spec *> getORHS(void) const;
	void setORHS(std::vector<Spec *> const & orhs);
	
	int getOrder(void) const;
	
	double getKcst(void) const;
	void setKcst(double kcst);
	
	std::vector<Spec *> getAllSpecs(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

} // namespace model
} // namespace steps

// END
