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

#include "mpi/tetvesicle/raftproxy.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "math/point.hpp"
#include "mpi/mpi_common.hpp"
#include "mpi/tetvesicle/patch_rdef.hpp"
#include "mpi/tetvesicle/raftsreac.hpp"
#include "solver/raftdef.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

RaftProxy::RaftProxy(solver::Raftdef* raftdef,
                     TriRDEF* central_tri,
                     solver::raft_individual_id unique_index)
    : pDef(raftdef)
    , pIndex(unique_index)
    , pTri(central_tri)
    , pImmobilityUpdate(0) {
    AssertLog(pDef != nullptr);
    AssertLog(pTri != nullptr);

    pRaftIndex = pDef->gidx();

    pPoolCount.resize(def()->countSpecs_global());
}

////////////////////////////////////////////////////////////////////////////////

void RaftProxy::checkpoint(std::fstream& /*cp_file*/) {
    // Reserve. Nothing to do here because only created when rafts are created
}

////////////////////////////////////////////////////////////////////////////////

void RaftProxy::restore(std::fstream& /*cp_file*/) {
    // Reserve.
}

////////////////////////////////////////////////////////////////////////////////

void RaftProxy::setSpecCountByLidx(solver::spec_local_id slidx, uint count) {
    AssertLog(slidx < def()->countSpecs());
    solver::spec_global_id spec_gidx = def()->specL2G(slidx);
    pPoolCount[spec_gidx.get()] = count;
}

////////////////////////////////////////////////////////////////////////////////

uint RaftProxy::getSpecCountByLidx(solver::spec_local_id slidx) {
    AssertLog(slidx < def()->countSpecs());
    solver::spec_global_id spec_gidx = def()->specL2G(slidx);
    return pPoolCount[spec_gidx.get()];
}

////////////////////////////////////////////////////////////////////////////////

void RaftProxy::setSpecCountByGidx(solver::spec_global_id sgidx, uint count) {
    AssertLog(sgidx < def()->countSpecs_global());
    pPoolCount[sgidx.get()] = count;
}

////////////////////////////////////////////////////////////////////////////////

uint RaftProxy::getSpecCountByGidx(solver::spec_global_id sgidx) {
    AssertLog(sgidx < def()->countSpecs_global());
    return pPoolCount[sgidx.get()];
}

////////////////////////////////////////////////////////////////////////////////

std::map<steps::index_t, uint> RaftProxy::getSpecs() {
    std::map<steps::index_t, uint> specs;
    for (auto spec_gidx: solver::spec_global_id::range(def()->countSpecs_global())) {
        uint count = pPoolCount[spec_gidx.get()];
        if (count > 0) {
            specs[spec_gidx.get()] = count;
        }
    }
    return specs;
}

////////////////////////////////////////////////////////////////////////////////

void RaftProxy::updImmobility(int mob_upd) {
    pImmobilityUpdate += mob_upd;
}

////////////////////////////////////////////////////////////////////////////////

bool RaftProxy::getRaftSReacActive(solver::raftsreac_global_id rsreacidx) const {
    return pRaftSReac_inactive.find(rsreacidx) == pRaftSReac_inactive.end();
}

////////////////////////////////////////////////////////////////////////////////

void RaftProxy::setRaftSReacInActive(solver::raftsreac_global_id rsreacidx) {
    pRaftSReac_inactive.insert(rsreacidx);
}

////////////////////////////////////////////////////////////////////////////////

}  // namespace steps::mpi::tetvesicle
