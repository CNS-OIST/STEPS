/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#pragma once

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include "solver/api.hpp"
#include "solver/statedef.hpp"

#include "util/checkpointing.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace steps::wmrk4 {

// Auxiliary declarations.
typedef std::vector<double> dVec;
typedef dVec::iterator dVecI;

typedef std::vector<uint> uiVec;
typedef uiVec::iterator uiVecI;

struct Reactant {
    uint globalIndex;
    uint order;
    Reactant(uint gidx, uint o)
        : globalIndex(gidx)
        , order(o) {}
};

struct SpecieInReaction {
    uint globalIndex;
    int populationChange;
    SpecieInReaction(uint gidx, int v)
        : globalIndex(gidx)
        , populationChange(v) {}
};

struct Reaction {
    std::vector<Reactant> reactants;
    std::vector<SpecieInReaction> affectedSpecies;
    bool isActivated;
    double c;  // constant relating populations, not concentrations !
    void addSpecies(uint gidx, uint order, int populationChange) {
        if (order > 0) {
            reactants.push_back(Reactant(gidx, order));
        }
        if (populationChange != 0) {
            affectedSpecies.push_back(SpecieInReaction(gidx, populationChange));
        }
    }

    void checkpoint(std::fstream& cp_file) {
        util::checkpoint(cp_file, isActivated);
        util::checkpoint(cp_file, c);
    }

    void restore(std::fstream& cp_file) {
        util::restore(cp_file, isActivated);
        util::restore(cp_file, c);
    }
};

class Wmrk4: public solver::API {
  public:
    Wmrk4(model::Model* m, wm::Geom* g, const rng::RNGptr& r);
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

    inline void setDT(double dt) override {
        setRk4DT(dt);
    }

    void setRk4DT(double dt) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      GENERAL
    ////////////////////////////////////////////////////////////////////////

    double getTime() const override;

    void checkpoint(std::string const& file_name) override;
    void restore(std::string const& file_name) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

    double _getCompVol(solver::comp_global_id cidx) const override;
    void _setCompVol(solver::comp_global_id cidx, double vol) override;

    double _getCompSpecCount(solver::comp_global_id cidx,
                             solver::spec_global_id sidx) const override;
    void _setCompSpecCount(solver::comp_global_id cidx,
                           solver::spec_global_id sidx,
                           double n) override;

    double _getCompSpecAmount(solver::comp_global_id cidx,
                              solver::spec_global_id sidx) const override;
    void _setCompSpecAmount(solver::comp_global_id cidx,
                            solver::spec_global_id sidx,
                            double a) override;

    double _getCompSpecConc(solver::comp_global_id cidx,
                            solver::spec_global_id sidx) const override;
    void _setCompSpecConc(solver::comp_global_id cidx,
                          solver::spec_global_id sidx,
                          double c) override;

    bool _getCompSpecClamped(solver::comp_global_id cidx,
                             solver::spec_global_id sidx) const override;
    void _setCompSpecClamped(solver::comp_global_id cidx,
                             solver::spec_global_id sidx,
                             bool b) override;

    double _getCompReacK(solver::comp_global_id cidx, solver::reac_global_id ridx) const override;
    void _setCompReacK(solver::comp_global_id cidx,
                       solver::reac_global_id ridx,
                       double kf) override;

    bool _getCompReacActive(solver::comp_global_id cidx,
                            solver::reac_global_id ridx) const override;
    void _setCompReacActive(solver::comp_global_id cidx,
                            solver::reac_global_id ridx,
                            bool a) override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

    double _getPatchArea(solver::patch_global_id pidx) const override;
    void _setPatchArea(solver::patch_global_id pidx, double area) override;

    double _getPatchSpecCount(solver::patch_global_id pidx,
                              solver::spec_global_id sidx) const override;
    void _setPatchSpecCount(solver::patch_global_id pidx,
                            solver::spec_global_id sidx,
                            double n) override;

    double _getPatchSpecAmount(solver::patch_global_id pidx,
                               solver::spec_global_id sidx) const override;
    void _setPatchSpecAmount(solver::patch_global_id pidx,
                             solver::spec_global_id sidx,
                             double a) override;

    bool _getPatchSpecClamped(solver::patch_global_id pidx,
                              solver::spec_global_id sidx) const override;
    void _setPatchSpecClamped(solver::patch_global_id pidx,
                              solver::spec_global_id sidx,
                              bool buf) override;

    double _getPatchSReacK(solver::patch_global_id pidx,
                           solver::sreac_global_id ridx) const override;
    void _setPatchSReacK(solver::patch_global_id pidx,
                         solver::sreac_global_id ridx,
                         double kf) override;

    bool _getPatchSReacActive(solver::patch_global_id pidx,
                              solver::sreac_global_id ridx) const override;
    void _setPatchSReacActive(solver::patch_global_id pidx,
                              solver::sreac_global_id ridx,
                              bool a) override;

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
    uint pSpecs_tot;

    /// number of reactions total: all reactions and surface reactions
    /// in each comp and patch
    uint pReacs_tot;

    /// vector holding current molecular counts (as doubles)
    dVec pVals;

    /// vector holding flags on species
    uiVec pSFlags;

    /// vector holding new, calculated counts
    dVec pNewVals;

    /// vector of present derivatives
    dVec pDyDx;

    /// the time step
    double pDT;

    /// objects to contain temporary values important in algorithm
    dVec yt;
    dVec dyt;
    dVec dym;

    std::vector<Reaction> reactions;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::wmrk4
