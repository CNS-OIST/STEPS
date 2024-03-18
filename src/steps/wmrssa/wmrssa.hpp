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
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

// STEPS headers.
#include "comp.hpp"
#include "kproc.hpp"
#include "patch.hpp"
#include "solver/api.hpp"
#include "solver/compdef.hpp"
#include "solver/fwd.hpp"
#include "solver/patchdef.hpp"
#include "solver/statedef.hpp"

namespace steps::wmrssa {

// Auxiliary declarations.
typedef solver::kproc_global_id SchedIDX;
typedef std::set<SchedIDX> SchedIDXSet;
typedef SchedIDXSet::iterator SchedIDXSetI;
typedef SchedIDXSet::const_iterator SchedIDXSetCI;
typedef std::vector<SchedIDX> SchedIDXVec;
typedef SchedIDXVec::iterator SchedIDXVecI;
typedef SchedIDXVec::const_iterator SchedIDXVecCI;

/// Copies the contents of a set of SchedIDX entries into a vector.
/// The contents of the vector are completely overridden.
///
extern void schedIDXSet_To_Vec(SchedIDXSet const& s, SchedIDXVec& v);

class Wmrssa: public solver::API {
  public:
    Wmrssa(model::Model* m, wm::Geom* g, const rng::RNGptr& r);
    ~Wmrssa() override;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::string const& file_name) override;

    /// restore data
    void restore(std::string const& file_name) override;

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

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      GENERAL
    ////////////////////////////////////////////////////////////////////////

    double getTime() const override;

    uint getNSteps() const override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      ADVANCE
    //      Developer only
    ////////////////////////////////////////////////////////////////////////

    void setTime(double time) override;
    void setNSteps(uint nsteps) override;

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

    double _getCompReacK(solver::comp_global_id cidx, solver::reac_global_id sridx) const override;
    void _setCompReacK(solver::comp_global_id cidx,
                       solver::reac_global_id sridx,
                       double kf) override;

    bool _getCompReacActive(solver::comp_global_id cidx,
                            solver::reac_global_id sridx) const override;
    void _setCompReacActive(solver::comp_global_id cidx,
                            solver::reac_global_id sridx,
                            bool a) override;

    double _getCompReacC(solver::comp_global_id cidx, solver::reac_global_id sridx) const override;

    double _getCompReacH(solver::comp_global_id cidx, solver::reac_global_id sridx) const override;

    long double _getCompReacA(solver::comp_global_id cidx,
                              solver::reac_global_id sridx) const override;

    unsigned long long _getCompReacExtent(solver::comp_global_id cidx,
                                          solver::reac_global_id sridx) const override;
    void _resetCompReacExtent(solver::comp_global_id cidx, solver::reac_global_id sridx) override;

    ////////////////////////////////////////////////////////////////////////
    /*
        double _getCompDiffD(solver::comp_global_id cidx, uint didx);
        void _setCompDiffD(solver::comp_global_id cidx, uint didx);

        bool _getCompDiffActive(solver::comp_global_id cidx, uint didx);
        void _setCompDiffActive(solver::comp_global_id cidx, uint didx, bool act);
    */

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
                           solver::sreac_global_id sridx) const override;
    void _setPatchSReacK(solver::patch_global_id pidx,
                         solver::sreac_global_id sridx,
                         double kf) override;

    bool _getPatchSReacActive(solver::patch_global_id pidx,
                              solver::sreac_global_id sridx) const override;
    void _setPatchSReacActive(solver::patch_global_id pidx,
                              solver::sreac_global_id sridx,
                              bool a) override;

    double _getPatchSReacC(solver::patch_global_id pidx,
                           solver::sreac_global_id sridx) const override;

    double _getPatchSReacH(solver::patch_global_id pidx,
                           solver::sreac_global_id sridx) const override;

    double _getPatchSReacA(solver::patch_global_id pidx,
                           solver::sreac_global_id sridx) const override;

    unsigned long long _getPatchSReacExtent(solver::patch_global_id pidx,
                                            solver::sreac_global_id sridx) const override;
    void _resetPatchSReacExtent(solver::patch_global_id pidx,
                                solver::sreac_global_id sridx) override;

    ////////////////////////////////////////////////////////////////////////

    // Called from local Comp or Patch objects. Ass KProc to this object
    void addKProc(wmrssa::KProc* kp);

    inline uint countKProcs() const noexcept {
        return pKProcs.size();
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////
    // WMRSSA SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

    uint _addComp(solver::Compdef* cdef);

    uint _addPatch(solver::Patchdef* pdef);

    // called when local comp, patch, reac, sreac objects have been created
    // by constructor
    void _setup();

    void _build();

    uint _getNext();

    void _reset();

    void _update(SchedIDXVec const& entries);

    void _executeStep(wmrssa::KProc* kp, double dt);

    ////////////////////////////////////////////////////////////////////////
    // LIST OF WMRSSA SOLVER OBJECTS
    ////////////////////////////////////////////////////////////////////////

    std::vector<wmrssa::KProc*> pKProcs;

    std::vector<wmrssa::Comp*> pComps;
    std::map<solver::Compdef*, Comp*> pCompMap;

    std::vector<wmrssa::Patch*> pPatches;

    /// \brief sum of propensities
    double pA0{0.0};

    ////////////////////////////////////////////////////////////////////////
    // N-ARY TREE
    ////////////////////////////////////////////////////////////////////////

    std::vector<std::vector<double>> pLevels;

    ////////////////////////////////////////////////////////////////////////

    // Keeps track of whether _build() has been called
    bool pBuilt{false};

    ////////////////////////////////////////////////////////////////////////

    // Tables to hold update vector indices and random numbers respectively,
    // to be re-used each step.
    std::vector<uint> pIndices;
    uint pMaxUpSize{0};
    std::vector<double> pRannum;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::wmrssa
