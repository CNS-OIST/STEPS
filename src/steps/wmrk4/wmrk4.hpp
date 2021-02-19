/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */


#ifndef STEPS_SOLVER_WMRK4_HPP
#define STEPS_SOLVER_WMRK4_HPP 1


// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/api.hpp"
#include "steps/solver/statedef.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace wmrk4  {
using namespace steps::solver;

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.

// Auxiliary declarations.
typedef std::vector<double>             dVec;
typedef dVec::iterator                    dVecI;

typedef std::vector<uint>                uiVec;
typedef uiVec::iterator                    uiVecI;

////////////////////////////////////////////////////////////////////////////////

struct Reactant{
    uint globalIndex;
    uint order;
    Reactant(uint gidx, uint o) : globalIndex(gidx), order(o) {}
};

struct SpecieInReaction{
    uint globalIndex;
    int populationChange;
    SpecieInReaction(uint gidx, int v) : globalIndex(gidx), populationChange(v) {}
};

struct Reaction{
    std::vector<Reactant> reactants;
    std::vector<SpecieInReaction> affectedSpecies;
    bool isActivated;
    double c; // constant relating populations, not concentrations !
    void addSpecies(uint gidx, uint order, int populationChange)
    {
        if (order > 0)
            reactants.push_back(Reactant(gidx, order));
        if (populationChange != 0)
            affectedSpecies.push_back(SpecieInReaction(gidx, populationChange));
    }
};

class Wmrk4: public API
{

public:

    Wmrk4(steps::model::Model *m, steps::wm::Geom *g, const rng::RNGptr &r);
    ~Wmrk4() override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER INFORMATION
    ////////////////////////////////////////////////////////////////////////

    std::string getSolverName() const override;
    std::string getSolverDesc() const override;
    std::string getSolverAuthors() const override;
    std::string getSolverEmail() const override;


    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS
    ////////////////////////////////////////////////////////////////////////

    void reset() override;
    void run(double endtime) override;
    void advance(double adv) override;
    void step() override;

    inline void setDT(double dt) override
    { setRk4DT(dt); }

    void setRk4DT(double dt) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      GENERAL
    ////////////////////////////////////////////////////////////////////////

    double getTime() const override;

    void checkpoint(std::string const & file_name) override;
    void restore(std::string const & file_name) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

     double _getCompVol(uint cidx) const override;
    void _setCompVol(uint cidx, double vol) override;

     double _getCompCount(uint cidx, uint sidx) const override;
     void _setCompCount(uint cidx, uint sidx, double n) override;

     double _getCompAmount(uint cidx, uint sidx) const override;
    void _setCompAmount(uint cidx, uint sidx, double a) override;

    double _getCompConc(uint cidx, uint sidx) const override;
     void _setCompConc(uint cidx, uint sidx, double c) override;

    bool _getCompClamped(uint cidx, uint sidx) const override;
    void _setCompClamped(uint cidx, uint sidx, bool b) override;

    double _getCompReacK(uint cidx, uint ridx) const override;
    void _setCompReacK(uint cidx, uint ridx, double kf) override;

     bool _getCompReacActive(uint cidx, uint ridx) const override;
    void _setCompReacActive(uint cidx, uint ridx, bool a) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

    double _getPatchArea(uint pidx) const override;
    void _setPatchArea(uint pidx, double area) override;

     double _getPatchCount(uint pidx, uint sidx) const override;
    void _setPatchCount(uint pidx, uint sidx, double n) override;

    double _getPatchAmount(uint pidx, uint sidx) const override;
     void _setPatchAmount(uint pidx, uint sidx, double a) override;

    bool _getPatchClamped(uint pidx, uint sidx) const override;
    void _setPatchClamped(uint pidx, uint sidx, bool buf) override;

    double _getPatchSReacK(uint pidx, uint ridx) const override;
      void _setPatchSReacK(uint pidx, uint ridx, double kf) override;

     bool _getPatchSReacActive(uint pidx, uint ridx) const override;
     void _setPatchSReacActive(uint pidx, uint ridx, bool a) override;

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////
    // WMRK4 SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

    /// initialises vectors and builds the reaction matrix.
    ///
    void _setup();

    /// this function refills the values and flags vectors
    /// called if a flag or a count changed
    ///
    void _refill();

    /// refill the Ccst vectors
    /// called if reaction constants are changed during simulation
    ///
    void _refillCcst();

    /// returns properly scaled reaction constant
    ///
    double _ccst(double kcst, double vol, uint order);

    /// returns properly scaled reaction constant for surface-surface case
    ///
    double _ccst2D(double kcst, double area, uint order);

    /// the Runge-Kutta algorithm
    ///
    void _rk4(double pdt);

    /// the simple stepper
    ///
    void _rksteps(double t1, double t2);

    /// the derivatives calculator
    ///
    void _setderivs(dVec& vals, dVec& dydx);

    /// update local values vector,
    /// then update state with computed counts
    ///
    void _update();

    ////////////////////////////////////////////////////////////////////////
    // WMRK4 SOLVER MEMBERS
    ////////////////////////////////////////////////////////////////////////

    /// number of species total: all species in all comps and patches
    uint                                pSpecs_tot;

    /// number of reactions total: all reactions and surface reactions
    /// in each comp and patch
    uint                                pReacs_tot;

    /// vector holding current molecular counts (as doubles)
    dVec                                pVals;

    /// vector holding flags on species
    uiVec                                pSFlags;

    /// vector holding new, calculated counts
    dVec                                pNewVals;

    /// vector of present derivatives
    dVec                                pDyDx;

    /// the time step
    double                                pDT;

    /// objects to contain temporary values important in algorithm
    dVec                                yt;
    dVec                                dyt;
    dVec                                dym;

    std::vector<Reaction> 				reactions;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_SOLVER_WMRK4_HPP

// END
