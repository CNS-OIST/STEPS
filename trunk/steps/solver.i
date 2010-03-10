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

%module solver_swig

%include "python/std_string.i"
%include "python/std_vector.i"
%include "error.i"
%import "cpp/steps/common.h"

%{
#include "../cpp/steps/solver/api.hpp"
#include "../cpp/steps/solver/statedef.hpp"
#include "../cpp/steps/wmrk4/wmrk4.hpp"
#include "../cpp/steps/wmdirect/wmdirect.hpp"
#include "../cpp/steps/tetexact/tetexact.hpp"
#include "../cpp/steps/error.hpp"
%}

////////////////////////////////////////////////////////////////////////////////

namespace steps
{
namespace solver
{

class Statedef;

}
}

////////////////////////////////////////////////////////////////////////////////

namespace std
{
///////// any vector definitions here /////
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
	
    API(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);
    virtual ~API(void);
	
    virtual std::string getSolverName(void) const = 0;
    virtual std::string getSolverDesc(void) const = 0;
    virtual std::string getSolverAuthors(void) const = 0;
    virtual std::string getSolverEmail(void) const = 0;

	
    virtual void reset(void) = 0;
    virtual void run(double endtime) = 0;
	virtual void advance(double adv);
    virtual void step(void);
	
    virtual void setDT(double dt);

    virtual double getTime(void) const = 0;
	virtual void setTime(double time);
    virtual double getDT(void) const;
	
    virtual double getA0(void) const;
	
	virtual unsigned int getNSteps(void) const;
    virtual void setNSteps(unsigned int nsteps);

    double getCompVol(std::string const & c) const;

    void setCompVol(std::string const & c, double vol);

    double getCompCount(std::string const & c, std::string const & s) const;

    void setCompCount(std::string const & c, std::string const & s, double n);

    double getCompAmount(std::string const & c, std::string const & s) const;

    void setCompAmount(std::string const & c, std::string const & s, double a);

    double getCompConc(std::string const & c, std::string const & s) const;

    void setCompConc(std::string const & c, std::string const & s, double c);

    bool getCompClamped(std::string const & c, std::string const & s) const;

    void setCompClamped(std::string const & c, std::string const & s, bool b);

    double getCompReacK(std::string const & c, std::string const & r) const;

    void setCompReacK(std::string const & c, std::string const & r, double kf);

    bool getCompReacActive(std::string const & c, std::string const & r) const;

    void setCompReacActive(std::string const & c, std::string const & r, bool a);
	
    double getCompDiffD(std::string const & c, std::string const & d) const;
	
    void setCompDiffD(std::string const & c, std::string const & d, double dcst);
	
    bool getCompDiffActive(std::string const & c, std::string const & d) const;
	
    void setCompDiffActive(std::string const & c, std::string const & d, bool act);

    double getCompReacC(std::string const & c, std::string const & r) const;

    double getCompReacH(std::string const & c, std::string const & r) const;

    double getCompReacA(std::string const & c, std::string const & r) const;

    unsigned int getCompReacExtent(std::string const & c, std::string const & r) const;

    void resetCompReacExtent(std::string const & c, std::string const & r);

    double getTetVol(unsigned int tidx) const;

    // void setTetVol(unsigned int tidx, double vol);

    double getTetCount(unsigned int tidx, std::string const & s) const;
	
    void setTetCount(unsigned int tidx, std::string const & s, double n);
	
    double getTetAmount(unsigned int tidx, std::string const & s) const;
	
    void setTetAmount(unsigned int tidx, std::string const & s, double m);
	
    double getTetConc(unsigned int tidx, std::string const & s) const;
	
    void setTetConc(unsigned int tidx, std::string const & s, double c);

    bool getTetClamped(unsigned int tidx, std::string const & s) const;
	
    void setTetClamped(unsigned int tidx, std::string const & s, bool buf);
	
    double getTetReacK(unsigned int tidx, std::string const & r) const;
	
    void setTetReacK(unsigned int tidx, std::string const & r, double kf);
	
    bool getTetReacActive(unsigned int tidx, std::string const & r) const;
	
    void setTetReacActive(unsigned int tidx, std::string const & r, bool act);
	
    double getTetDiffD(unsigned int tidx, std::string const & d) const;
	
    void setTetDiffD(unsigned int tidx, std::string const & d, double dk);
	
    bool getTetDiffActive(unsigned int tidx, std::string const & d) const;
	
    void setTetDiffActive(unsigned int tidx, std::string const & d, bool act);

    double getTetReacC(unsigned int tidx, std::string const & r) const;

    double getTetReacH(unsigned int tidx, std::string const & r) const;

    double getTetReacA(unsigned int tidx, std::string const & r) const;
	
    double getTetDiffA(unsigned int tidx, std::string const & d) const;

    double getPatchArea(std::string const & p) const;

    void setPatchArea(std::string const & p, double area);

    double getPatchCount(std::string const & p, std::string const & s) const;

    void setPatchCount(std::string const & p, std::string const & s, double n);

    double getPatchAmount(std::string const & p, std::string const & s) const;

    void setPatchAmount(std::string const & p, std::string const & s, double a);

    bool getPatchClamped(std::string const & p, std::string const & s) const;

    void setPatchClamped(std::string const & p, std::string const & s, bool buf);

    double getPatchSReacK(std::string const & p, std::string const & r) const;

    void setPatchSReacK(std::string const & p, std::string const & r, double kf);

    bool getPatchSReacActive(std::string const & p, std::string const & r) const;

    void setPatchSReacActive(std::string const & p, std::string const & r, bool a);

    double getPatchSReacC(std::string const & p, std::string const & r) const;

    double getPatchSReacH(std::string const & p, std::string const & r) const;

    double getPatchSReacA(std::string const & p, std::string const & r) const;

    unsigned int getPatchSReacExtent(std::string const & p, std::string const & r) const;

    void resetPatchSReacExtent(std::string const & p, std::string const & r);

    double getTriArea(unsigned int tidx) const;
	
    // void setTriArea(unsigned int tidx, double area);
	
    double getTriCount(unsigned int tidx, std::string const & s) const;

    void setTriCount(unsigned int tidx, std::string const & s, double n);

    double getTriAmount(unsigned int tidx, std::string const & s) const;

    void setTriAmount(unsigned int tidx, std::string const & s, double m);
    
    bool getTriClamped(unsigned int tidx, std::string const & s) const;
	
    void setTriClamped(unsigned int tidx, std::string const & s, bool buf);

    double getTriSReacK(unsigned int tidx, std::string const & r) const;

    void setTriSReacK(unsigned int tidx, std::string const & r, double kf);
	
    bool getTriSReacActive(unsigned int tidx, std::string const & r) const;
	
    void setTriSReacActive(unsigned int tidx, std::string const & r, bool act);
	
    double getTriSReacC(unsigned int tidx, std::string const & r) const;

    double getTriSReacH(unsigned int tidx, std::string const & r) const;

    double getTriSReacA(unsigned int tidx, std::string const & r) const;

};  

////////////////////////////////////////////////////////////////////////////////
	
} // end namespace solver
} // end namespace steps

////////////////////////////////////////////////////////////////////////////////

namespace steps
{
namespace wmrk4
{
	
class Wmrk4 : public steps::solver::API
{
	
public:
	
	Wmrk4(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);
	~Wmrk4(void);

	std::string getSolverName(void) const;
	std::string getSolverDesc(void) const;
	std::string getSolverAuthors(void) const;
	std::string getSolverEmail(void) const;


	void reset(void);
	void run(double endtime);
	void step(void);
		
	void setDT(double dt);
		
	double getTime(void) const;
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
	Wmdirect(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);
	~Wmdirect(void);
	
	std::string getSolverName(void) const;
	std::string getSolverDesc(void) const;
	std::string getSolverAuthors(void) const;
	std::string getSolverEmail(void) const;
	
	
	void reset(void);
	void run(double endtime);
    void advance(double adv);
	void step(void);

	double getTime(void) const;
	
    double getA0(void) const;
	
    unsigned int getNSteps(void) const;
    void setTime(double time);
    void setNSteps(unsigned int nsteps);	
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
    Tetexact(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);
    ~Tetexact(void);

	std::string getSolverName(void) const;
	std::string getSolverDesc(void) const;
	std::string getSolverAuthors(void) const;
	std::string getSolverEmail(void) const;
	
	
	void reset(void);
	void run(double endtime);
	void step(void);
    void advance(double adv);
	
	double getTime(void) const;
	
    double getA0(void) const;
	
    unsigned int getNSteps(void) const;
    void setTime(double time);
    void setNSteps(unsigned int nsteps);
	
	////////////////////////////////////////////////////////////////////////			
	
};
	
////////////////////////////////////////////////////////////////////////////////
	
} // end namespace tetexact
} // end namespace steps
    