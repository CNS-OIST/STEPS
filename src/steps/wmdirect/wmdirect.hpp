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
#include <unordered_set>
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

namespace steps::wmdirect {


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

class Wmdirect: public solver::API {
  public:
    Wmdirect(model::Model* m, wm::Geom* g, const rng::RNGptr& r);
    ~Wmdirect() override;

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

    inline double getA0() const noexcept override {
        return pA0;
    }

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
    //      GENERAL
    ////////////////////////////////////////////////////////////////////////

    inline solver::complex_individual_id getNextComplexInd() {
        return pNextComplexStateInd++;
    }

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

    double _getCompComplexCount(
        solver::comp_global_id cidx,
        solver::complex_global_id sidx,
        const std::vector<util::strongid_vector<solver::complex_substate_id,
                                                model::SubunitStateFilter>>& f) const override;
    void _setCompComplexCount(solver::comp_global_id cidx,
                              solver::complex_global_id sidx,
                              const util::strongid_vector<solver::complex_substate_id, uint>& f,
                              double n) override;

    double _getCompComplexAmount(
        solver::comp_global_id cidx,
        solver::complex_global_id sidx,
        const std::vector<util::strongid_vector<solver::complex_substate_id,
                                                model::SubunitStateFilter>>& f) const override;
    void _setCompComplexAmount(solver::comp_global_id cidx,
                               solver::complex_global_id sidx,
                               const util::strongid_vector<solver::complex_substate_id, uint>& f,
                               double a) override;

    double _getCompComplexConc(
        solver::comp_global_id cidx,
        solver::complex_global_id sidx,
        const std::vector<util::strongid_vector<solver::complex_substate_id,
                                                model::SubunitStateFilter>>& f) const override;
    void _setCompComplexConc(solver::comp_global_id cidx,
                             solver::complex_global_id sidx,
                             const util::strongid_vector<solver::complex_substate_id, uint>& f,
                             double c) override;

    double _getCompComplexSUSCount(
        solver::comp_global_id cidx,
        solver::complex_global_id sidx,
        const std::vector<
            util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f,
        solver::complex_substate_id m) const override;

    double _getCompComplexSUSConc(
        solver::comp_global_id cidx,
        solver::complex_global_id sidx,
        const std::vector<
            util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f,
        solver::complex_substate_id m) const override;

    double _getCompComplexSUSAmount(
        solver::comp_global_id cidx,
        solver::complex_global_id sidx,
        const std::vector<
            util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f,
        solver::complex_substate_id m) const override;

    bool _getCompReacActive(solver::comp_global_id cidx,
                            solver::reac_global_id ridx) const override;
    void _setCompReacActive(solver::comp_global_id cidx,
                            solver::reac_global_id ridx,
                            bool a) override;

    double _getCompReacC(solver::comp_global_id cidx, solver::reac_global_id ridx) const override;

    double _getCompReacH(solver::comp_global_id cidx, solver::reac_global_id ridx) const override;

    long double _getCompReacA(solver::comp_global_id cidx,
                              solver::reac_global_id ridx) const override;

    unsigned long long _getCompReacExtent(solver::comp_global_id cidx,
                                          solver::reac_global_id ridx) const override;
    void _resetCompReacExtent(solver::comp_global_id cidx, solver::reac_global_id ridx) override;

    unsigned long long _getCompComplexReacExtent(solver::comp_global_id cidx,
                                                 solver::complexreac_global_id ridx) const override;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

    double _getPatchArea(solver::patch_global_id cidx) const override;
    void _setPatchArea(solver::patch_global_id cidx, double area) override;

    double _getPatchSpecCount(solver::patch_global_id cidx,
                              solver::spec_global_id sidx) const override;
    void _setPatchSpecCount(solver::patch_global_id cidx,
                            solver::spec_global_id sidx,
                            double n) override;

    double _getPatchSpecAmount(solver::patch_global_id cidx,
                               solver::spec_global_id sidx) const override;
    void _setPatchSpecAmount(solver::patch_global_id cidx,
                             solver::spec_global_id sidx,
                             double a) override;

    bool _getPatchSpecClamped(solver::patch_global_id cidx,
                              solver::spec_global_id sidx) const override;
    void _setPatchSpecClamped(solver::patch_global_id cidx,
                              solver::spec_global_id sidx,
                              bool buf) override;

    double _getPatchComplexCount(
        solver::patch_global_id cidx,
        solver::complex_global_id sidx,
        const std::vector<util::strongid_vector<solver::complex_substate_id,
                                                model::SubunitStateFilter>>& f) const override;
    void _setPatchComplexCount(solver::patch_global_id cidx,
                               solver::complex_global_id sidx,
                               const util::strongid_vector<solver::complex_substate_id, uint>& f,
                               double n) override;

    double _getPatchComplexAmount(
        solver::patch_global_id cidx,
        solver::complex_global_id sidx,
        const std::vector<util::strongid_vector<solver::complex_substate_id,
                                                model::SubunitStateFilter>>& f) const override;
    void _setPatchComplexAmount(solver::patch_global_id cidx,
                                solver::complex_global_id sidx,
                                const util::strongid_vector<solver::complex_substate_id, uint>& f,
                                double a) override;

    double _getPatchComplexSUSCount(
        solver::patch_global_id cidx,
        solver::complex_global_id sidx,
        const std::vector<
            util::strongid_vector<solver::complex_substate_id, model::SubunitStateFilter>>& f,
        solver::complex_substate_id m) const override;

    double _getPatchSReacK(solver::patch_global_id cidx,
                           solver::sreac_global_id ridx) const override;
    void _setPatchSReacK(solver::patch_global_id cidx,
                         solver::sreac_global_id ridx,
                         double kf) override;

    bool _getPatchSReacActive(solver::patch_global_id cidx,
                              solver::sreac_global_id ridx) const override;
    void _setPatchSReacActive(solver::patch_global_id cidx,
                              solver::sreac_global_id ridx,
                              bool a) override;

    double _getPatchSReacC(solver::patch_global_id cidx,
                           solver::sreac_global_id ridx) const override;

    double _getPatchSReacH(solver::patch_global_id cidx,
                           solver::sreac_global_id ridx) const override;

    double _getPatchSReacA(solver::patch_global_id cidx,
                           solver::sreac_global_id ridx) const override;

    unsigned long long _getPatchSReacExtent(solver::patch_global_id cidx,
                                            solver::sreac_global_id ridx) const override;
    void _resetPatchSReacExtent(solver::patch_global_id cidx,
                                solver::sreac_global_id ridx) override;

    unsigned long long _getPatchComplexSReacExtent(
        solver::patch_global_id pidx,
        solver::complexsreac_global_id ridx) const override;

    ////////////////////////////////////////////////////////////////////////

    // Called from local Comp or Patch objects. Ass KProc to this object
    void addKProc(steps::wmdirect::KProc* kp);

    inline uint countKProcs() const noexcept {
        return pKProcs.size();
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////
    // WMDIRECT SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

    uint _addComp(solver::Compdef* cdef);

    uint _addPatch(solver::Patchdef* pdef);

    // called when local comp, patch, reac, sreac objects have been created
    // by constructor
    void _setup();

    void _build();

    steps::wmdirect::KProc* _getNext();

    void _reset();

    void _update(SchedIDXVec const& entries);

    void _clearComplexFilterUpdates();

    void _executeStep(steps::wmdirect::KProc* kp, double dt);

    ////////////////////////////////////////////////////////////////////////
    // LIST OF WMDIRECT SOLVER OBJECTS
    ////////////////////////////////////////////////////////////////////////

    std::vector<steps::wmdirect::KProc*> pKProcs;

    std::vector<steps::wmdirect::Comp*> pComps;
    std::map<solver::Compdef*, Comp*> pCompMap;

    std::vector<steps::wmdirect::Patch*> pPatches;

    /// \brief sum of propensities
    double pA0{0.0};

    solver::complex_individual_id pNextComplexStateInd{};

    ////////////////////////////////////////////////////////////////////////
    // N-ARY TREE
    ////////////////////////////////////////////////////////////////////////

    std::vector<uint> pLevelSizes;

    std::vector<double*> pLevels;

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

}  // namespace steps::wmdirect
