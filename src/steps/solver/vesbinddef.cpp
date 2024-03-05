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

#include "solver/vesbinddef.hpp"

#include "solver/statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

VesBinddef::VesBinddef(Statedef& sd, vesbind_global_id idx, model::VesBind& vb)
    : pIdx(idx)
    , pName(vb.getID())
    , pKcst(vb.getKcst())
    , pMaxDistance(vb.getLengthMax())
    , pMinDistance(vb.getLengthMin())
    , pReactants1(vb.getReactants1())
    , pReactants2(vb.getReactants2())
    , pProducts1(vb.getProducts1())
    , pProducts2(vb.getProducts2())
    , pVdep1(vb.getVDeps1())
    , pVdep2(vb.getVDeps2())
    , pLdep1(vb.getLDeps1())
    , pLdep2(vb.getLDeps2())
    , pImmobility(vb.getImmobilization()) {
    uint nspecs = sd.countSpecs();

    pSpec_VDEP1.container().resize(nspecs);
    pSpec_VDEP2.container().resize(nspecs);

    uint nlspecs = sd.countLinkSpecs();

    pSpec_LDEP1.container().resize(nlspecs);
    pSpec_LDEP2.container().resize(nlspecs);
}

////////////////////////////////////////////////////////////////////////////////

void VesBinddef::checkpoint(std::fstream& /*cp_file*/) const {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void VesBinddef::restore(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void VesBinddef::setup(const Statedef& sd) {
    AssertLog(pSetupdone == false);
    const auto& reactants1 = pReactants1;
    const auto& reactants2 = pReactants2;
    const auto& products1 = pProducts1;
    const auto& products2 = pProducts2;

    // The vesicles '1' and '2' should match
    AssertLog(reactants1.first == products1.first);
    AssertLog(reactants2.first == products2.first);

    pVesicle_1_idx = sd.getVesicleIdx(*reactants1.first);
    pSpec_1_gidx = sd.getSpecIdx(*reactants1.second);

    pVesicle_2_idx = sd.getVesicleIdx(*reactants2.first);
    pSpec_2_gidx = sd.getSpecIdx(*reactants2.second);

    pProduct1_gidx = sd.getLinkSpecIdx(*products1.second);
    pProduct2_gidx = sd.getLinkSpecIdx(*products2.second);

    pProduct1def = &sd.linkspecdef(pProduct1_gidx);
    pProduct2def = &sd.linkspecdef(pProduct2_gidx);

    for (auto const& vdep1: pVdep1) {
        spec_global_id sidx = sd.getSpecIdx(*vdep1);
        pSpec_VDEP1[sidx] += 1;
    }

    for (auto const& vdep2: pVdep2) {
        spec_global_id sidx = sd.getSpecIdx(*vdep2);
        pSpec_VDEP2[sidx] += 1;
    }

    for (auto const& ldep1: pLdep1) {
        linkspec_global_id lsidx = sd.getLinkSpecIdx(*ldep1);
        pSpec_LDEP1[lsidx] += 1;
    }

    for (auto const& ldep2: pLdep2) {
        linkspec_global_id lsidx = sd.getLinkSpecIdx(*ldep2);
        pSpec_LDEP2[lsidx] += 1;
    }
    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

vesicle_global_id VesBinddef::getVes1idx() const {
    AssertLog(pSetupdone == true);
    return pVesicle_1_idx;
}

////////////////////////////////////////////////////////////////////////////////

vesicle_global_id VesBinddef::getVes2idx() const {
    AssertLog(pSetupdone == true);
    return pVesicle_2_idx;
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id VesBinddef::getSpec1gidx() const {
    AssertLog(pSetupdone == true);
    return pSpec_1_gidx;
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id VesBinddef::getSpec2gidx() const {
    AssertLog(pSetupdone == true);
    return pSpec_2_gidx;
}

////////////////////////////////////////////////////////////////////////////////

linkspec_global_id VesBinddef::getLinkSpec1gidx() const {
    AssertLog(pSetupdone == true);
    return pProduct1_gidx;
}

////////////////////////////////////////////////////////////////////////////////

linkspec_global_id VesBinddef::getLinkSpec2gidx() const {
    AssertLog(pSetupdone == true);
    return pProduct2_gidx;
}

////////////////////////////////////////////////////////////////////////////////

LinkSpecdef& VesBinddef::getLinkSpec1def() const {
    AssertLog(pSetupdone == true);
    return *pProduct1def;
}

////////////////////////////////////////////////////////////////////////////////

LinkSpecdef& VesBinddef::getLinkSpec2def() const {
    AssertLog(pSetupdone == true);
    return *pProduct2def;
}

////////////////////////////////////////////////////////////////////////////////

uint VesBinddef::vdep1(spec_global_id gidx) const {
    return pSpec_VDEP1.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint VesBinddef::vdep2(spec_global_id gidx) const {
    return pSpec_VDEP2.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint VesBinddef::ldep1(linkspec_global_id gidx) const {
    return pSpec_LDEP1.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint VesBinddef::ldep2(linkspec_global_id gidx) const {
    return pSpec_LDEP2.at(gidx);
}

}  // namespace steps::solver
