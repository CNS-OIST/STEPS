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

#include "mpi/mpi_common.hpp"
#include "util/strong_id.hpp"

#include <mpi.h>

#include "solver/endocytosisdef.hpp"
#include "solver/exocytosisdef.hpp"
#include "solver/raftendocytosisdef.hpp"

namespace steps::mpi::tetvesicle {

// Auxiliary declarations.
typedef solver::kproc_global_id SchedIDX;
typedef std::set<SchedIDX> SchedIDXSet;
typedef SchedIDXSet::iterator SchedIDXSetI;
typedef SchedIDXSet::const_iterator SchedIDXSetCI;
typedef std::vector<SchedIDX> SchedIDXVec;
typedef SchedIDXVec::iterator SchedIDXVecI;
typedef SchedIDXVec::const_iterator SchedIDXVecCI;

enum SubVolType { SUB_TET, SUB_TRI };
enum SyncDirection { VESRAFT_TO_RDEF, RDEF_TO_VESRAFT };

// 'V2R' reflects VesRaft -> RDEF communication,
// 'R2V' reflects RDEF -> VesRaft communication
// 'Sync' reflects bidirectional communication

////////////// Regular species tri and tet pools, sync /////////////////////////

struct PoolCountSync {
    index_t container_global_index{};
    solver::spec_global_id spec_global_index{};
    uint count{};
};

////////////// Tet overlap ////////////////////////////////////////////

struct TetV2R {
    tetrahedron_global_id tetrahedron_global_index{};
    double overlap{};
};

/////////// Vesicles and vesicle species and link species, V2R /////////////////

struct VesProxyV2R {
    tetrahedron_global_id tetrahedron_global_index{};
    solver::vesicle_global_id vesicle_global_index{};
    solver::vesicle_individual_id vesicle_individual_index{};
    bool contains_link{};
    std::array<double, 3> vesicle_central_position{};
};

struct VesSurfSpecV2R {
    tetrahedron_global_id tetrahedron_global_index{};
    solver::vesicle_individual_id vesicle_individual_index{};
    solver::spec_global_id surface_spec_global_index{};
    solver::pointspec_individual_id surface_spec_individual_index{};
    std::array<double, 3> surface_spec_position_abs{};
};

/* Not necessary because inner species are not available for SSA interactions
struct VesInnerSpecV2R {
  index_t tetrahedron_global_index;
  index_t vesicle_individual_index;
  index_t inner_spec_global_index;
  uint count;
};
*/

struct VesLinkSpecV2R {
    tetrahedron_global_id tetrahedron_global_index{};
    solver::vesicle_individual_id vesicle_individual_index{};
    solver::linkspec_global_id linkspec_global_index{};
    solver::linkspec_individual_id linkspec_individual_index{};
    std::array<double, 3> linkspec_position_abs{};
};

/////////// Vesicles and vesicle species and link species, R2V /////////////////

struct VesProxyR2V {
    tetrahedron_global_id tetrahedron_global_index{};
    solver::vesicle_individual_id vesicle_individual_index{};
    solver::exocytosis_global_id exo_applied_global_index{};
    int immobility_update{};
};

struct VesSurfSpecR2V {
    tetrahedron_global_id tetrahedron_global_index{};
    solver::vesicle_individual_id vesicle_individual_index{};
    solver::spec_global_id surface_spec_global_index{};
    solver::pointspec_individual_id surface_spec_individual_index{};
    std::array<double, 3> surface_spec_position_abs{};
};

struct VesInnerSpecR2V {
    tetrahedron_global_id tetrahedron_global_index{};
    solver::vesicle_individual_id vesicle_individual_index{};
    solver::spec_global_id inner_spec_global_index{};
    uint count{};
};

// new link specs and corresponding pairs can be created by RDEF during
// vesbinding events. Important that this comes first in update so that VesRaft
// knows about these new LinkSpecs
struct VesLinkSpecPairR2V {
    tetrahedron_global_id tetrahedron_global_index{};
    solver::vesicle_individual_id vesicle1_individual_index{};
    solver::vesicle_individual_id vesicle2_individual_index{};
    solver::linkspec_global_id linkspec1_global_index{};
    solver::linkspec_global_id linkspec2_global_index{};
    solver::linkspec_individual_id linkspec1_individual_index{};
    solver::linkspec_individual_id linkspec2_individual_index{};
    std::array<double, 3> linkspec1_position_abs{};
    std::array<double, 3> linkspec2_position_abs{};
    double min_length{};
    double max_length{};
};

struct VesLinkSpecR2V {
    tetrahedron_global_id tetrahedron_global_index{};
    solver::linkspec_global_id linkspec_global_index{};
    solver::linkspec_individual_id linkspec_individual_index{};
    solver::vesicle_individual_id vesicle_individual_index{};
    // double linkspec_position_abs[3];
};

/////////// Raft and raft species, V2R /////////////////

struct RaftProxyV2R {
    triangle_global_id triangle_global_index;
    solver::raft_global_id raft_global_index;
    solver::raft_individual_id raft_individual_index;
};

struct RaftSpecV2R {
    triangle_global_id triangle_global_index;
    solver::raft_individual_id raft_individual_index;
    solver::spec_global_id spec_global_index;
    double count;
};

struct RaftSReacInactiveV2R {
    triangle_global_id triangle_global_index;
    solver::raft_individual_id raft_individual_index;
    solver::raftsreac_global_id raftsreac_global_index;
};

/////////// Raft and raft species, V2R /////////////////
struct RaftProxyR2V {
    triangle_global_id triangle_global_index;
    solver::raft_individual_id raft_individual_index;
    int immobility_update;
};

struct RaftSpecR2V {
    triangle_global_id triangle_global_index;
    solver::raft_individual_id raft_individual_index;
    solver::spec_global_id spec_global_index;
    uint count;
};

////////////// Raft generation, R2V ////////////////////////////////////////////

struct RaftGenCountR2V {
    triangle_global_id triangle_global_index{};
    solver::raftgen_global_id raftgen_global_index{};
    uint count{};
};


struct MPIDataTypeUtil {
    MPI_Datatype MPI_PoolCountSync;
    MPI_Datatype MPI_TetV2R;
    MPI_Datatype MPI_VesProxyV2R;
    MPI_Datatype MPI_VesSurfSpecV2R;
    MPI_Datatype MPI_VesLinkSpecV2R;
    MPI_Datatype MPI_VesProxyR2V;
    MPI_Datatype MPI_VesSurfSpecR2V;
    MPI_Datatype MPI_VesInnerSpecR2V;
    MPI_Datatype MPI_VesLinkSpecPairR2V;
    MPI_Datatype MPI_VesLinkSpecR2V;
    MPI_Datatype MPI_RaftProxyV2R;
    MPI_Datatype MPI_RaftSpecV2R;
    MPI_Datatype MPI_RaftSReacInactiveV2R;
    MPI_Datatype MPI_RaftProxyR2V;
    MPI_Datatype MPI_RaftSpecR2V;
    MPI_Datatype MPI_RaftGenCountR2V;
    MPI_Datatype MPI_ExocytosisEventSync;
    MPI_Datatype MPI_RaftEndocytosisEventSync;
    MPI_Datatype MPI_EndocytosisEventSync;


    void commitAllDataTypes() {
        commitPoolCountSync(MPI_PoolCountSync);
        commitTetV2R(MPI_TetV2R);
        commitVesProxyV2R(MPI_VesProxyV2R);
        commitVesSurfSpecV2R(MPI_VesSurfSpecV2R);
        commitVesLinkSpecV2R(MPI_VesLinkSpecV2R);
        commitVesProxyR2V(MPI_VesProxyR2V);
        commitVesSurfSpecR2V(MPI_VesSurfSpecR2V);
        commitVesInnerSpecR2V(MPI_VesInnerSpecR2V);
        commitVesLinkSpecPairR2V(MPI_VesLinkSpecPairR2V);
        commitVesLinkSpecR2V(MPI_VesLinkSpecR2V);
        commitRaftProxyV2R(MPI_RaftProxyV2R);
        commitRaftSpecV2R(MPI_RaftSpecV2R);
        commitRaftSReacInactiveV2R(MPI_RaftSReacInactiveV2R);
        commitRaftProxyR2V(MPI_RaftProxyR2V);
        commitRaftSpecR2V(MPI_RaftSpecR2V);
        commitRaftGenCountR2V(MPI_RaftGenCountR2V);
        commitExocytosisEventSync(MPI_ExocytosisEventSync);
        commitRaftEndocytosisEventSync(MPI_RaftEndocytosisEventSync);
        commitEndocytosisEventSync(MPI_EndocytosisEventSync);
    }

    void freeAllDataTypes() {
        MPI_Type_free(&MPI_PoolCountSync);
        MPI_Type_free(&MPI_TetV2R);
        MPI_Type_free(&MPI_VesProxyV2R);
        MPI_Type_free(&MPI_VesSurfSpecV2R);
        MPI_Type_free(&MPI_VesLinkSpecV2R);
        MPI_Type_free(&MPI_VesProxyR2V);
        MPI_Type_free(&MPI_VesSurfSpecR2V);
        MPI_Type_free(&MPI_VesInnerSpecR2V);
        MPI_Type_free(&MPI_VesLinkSpecPairR2V);
        MPI_Type_free(&MPI_VesLinkSpecR2V);
        MPI_Type_free(&MPI_RaftProxyV2R);
        MPI_Type_free(&MPI_RaftSpecV2R);
        MPI_Type_free(&MPI_RaftSReacInactiveV2R);
        MPI_Type_free(&MPI_RaftProxyR2V);
        MPI_Type_free(&MPI_RaftSpecR2V);
        MPI_Type_free(&MPI_RaftGenCountR2V);
        MPI_Type_free(&MPI_ExocytosisEventSync);
        MPI_Type_free(&MPI_RaftEndocytosisEventSync);
        MPI_Type_free(&MPI_EndocytosisEventSync);
    }

  private:
    void commitPoolCountSync(MPI_Datatype& new_type) {
        PoolCountSync temp{};
        constexpr size_t num_blocks = 3;
        int blocklens[num_blocks] = {1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_UNSIGNED};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.container_global_index, addr);
        MPI_Get_address(&temp.spec_global_index, addr + 1);
        MPI_Get_address(&temp.count, addr + 2);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitTetV2R(MPI_Datatype& new_type) {
        TetV2R temp{};
        constexpr size_t num_blocks = 2;
        int blocklens[num_blocks] = {1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX, MPI_DOUBLE};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.tetrahedron_global_index, addr);
        MPI_Get_address(&temp.overlap, addr + 1);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitVesProxyV2R(MPI_Datatype& new_type) {
        VesProxyV2R temp{};
        constexpr size_t num_blocks = 5;
        int blocklens[num_blocks] = {1, 1, 1, 1, 3};
        MPI_Datatype old_types[num_blocks] = {
            MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_C_BOOL, MPI_DOUBLE};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.tetrahedron_global_index, addr);
        MPI_Get_address(&temp.vesicle_global_index, addr + 1);
        MPI_Get_address(&temp.vesicle_individual_index, addr + 2);
        MPI_Get_address(&temp.contains_link, addr + 3);
        MPI_Get_address(temp.vesicle_central_position.data(), addr + 4);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitVesSurfSpecV2R(MPI_Datatype& new_type) {
        VesSurfSpecV2R temp{};
        constexpr size_t num_blocks = 5;
        int blocklens[num_blocks] = {1, 1, 1, 1, 3};
        MPI_Datatype old_types[num_blocks] = {
            MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_DOUBLE};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.tetrahedron_global_index, addr);
        MPI_Get_address(&temp.vesicle_individual_index, addr + 1);
        MPI_Get_address(&temp.surface_spec_global_index, addr + 2);
        MPI_Get_address(&temp.surface_spec_individual_index, addr + 3);
        MPI_Get_address(temp.surface_spec_position_abs.data(), addr + 4);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitVesLinkSpecV2R(MPI_Datatype& new_type) {
        VesLinkSpecV2R temp{};
        constexpr size_t num_blocks = 5;
        int blocklens[num_blocks] = {1, 1, 1, 1, 3};
        MPI_Datatype old_types[num_blocks] = {
            MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_DOUBLE};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.tetrahedron_global_index, addr);
        MPI_Get_address(&temp.vesicle_individual_index, addr + 1);
        MPI_Get_address(&temp.linkspec_global_index, addr + 2);
        MPI_Get_address(&temp.linkspec_individual_index, addr + 3);
        MPI_Get_address(temp.linkspec_position_abs.data(), addr + 4);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitVesProxyR2V(MPI_Datatype& new_type) {
        VesProxyR2V temp{};
        constexpr size_t num_blocks = 4;
        int blocklens[num_blocks] = {1, 1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_INT};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.tetrahedron_global_index, addr);
        MPI_Get_address(&temp.vesicle_individual_index, addr + 1);
        MPI_Get_address(&temp.exo_applied_global_index, addr + 2);
        MPI_Get_address(&temp.immobility_update, addr + 3);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitVesSurfSpecR2V(MPI_Datatype& new_type) {
        VesSurfSpecR2V temp{};
        constexpr size_t num_blocks = 5;
        int blocklens[num_blocks] = {1, 1, 1, 1, 3};
        MPI_Datatype old_types[num_blocks] = {
            MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_DOUBLE};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.tetrahedron_global_index, addr);
        MPI_Get_address(&temp.vesicle_individual_index, addr + 1);
        MPI_Get_address(&temp.surface_spec_global_index, addr + 2);
        MPI_Get_address(&temp.surface_spec_individual_index, addr + 3);
        MPI_Get_address(temp.surface_spec_position_abs.data(), addr + 4);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitVesInnerSpecR2V(MPI_Datatype& new_type) {
        VesInnerSpecR2V temp{};
        constexpr size_t num_blocks = 4;
        int blocklens[num_blocks] = {1, 1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_UNSIGNED};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.tetrahedron_global_index, addr);
        MPI_Get_address(&temp.vesicle_individual_index, addr + 1);
        MPI_Get_address(&temp.inner_spec_global_index, addr + 2);
        MPI_Get_address(&temp.count, addr + 3);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitVesLinkSpecPairR2V(MPI_Datatype& new_type) {
        VesLinkSpecPairR2V temp{};
        constexpr size_t num_blocks = 11;
        int blocklens[num_blocks] = {1, 1, 1, 1, 1, 1, 1, 3, 3, 1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_DOUBLE,
                                              MPI_DOUBLE,
                                              MPI_DOUBLE,
                                              MPI_DOUBLE};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.tetrahedron_global_index, addr);
        MPI_Get_address(&temp.vesicle1_individual_index, addr + 1);
        MPI_Get_address(&temp.vesicle2_individual_index, addr + 2);
        MPI_Get_address(&temp.linkspec1_global_index, addr + 3);
        MPI_Get_address(&temp.linkspec2_global_index, addr + 4);
        MPI_Get_address(&temp.linkspec1_individual_index, addr + 5);
        MPI_Get_address(&temp.linkspec2_individual_index, addr + 6);
        MPI_Get_address(temp.linkspec1_position_abs.data(), addr + 7);
        MPI_Get_address(temp.linkspec2_position_abs.data(), addr + 8);
        MPI_Get_address(&temp.min_length, addr + 9);
        MPI_Get_address(&temp.max_length, addr + 10);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitVesLinkSpecR2V(MPI_Datatype& new_type) {
        VesLinkSpecR2V temp{};
        constexpr size_t num_blocks = 4;
        int blocklens[num_blocks] = {1, 1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.tetrahedron_global_index, addr);
        MPI_Get_address(&temp.linkspec_global_index, addr + 1);
        MPI_Get_address(&temp.linkspec_individual_index, addr + 2);
        MPI_Get_address(&temp.vesicle_individual_index, addr + 3);
        // MPI_Get_address(temp.linkspec_position_abs, addr + 4);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitRaftProxyV2R(MPI_Datatype& new_type) {
        RaftProxyV2R temp{};
        constexpr size_t num_blocks = 3;
        int blocklens[num_blocks] = {1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_STEPS_INDEX};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.triangle_global_index, addr);
        MPI_Get_address(&temp.raft_global_index, addr + 1);
        MPI_Get_address(&temp.raft_individual_index, addr + 2);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitRaftSpecV2R(MPI_Datatype& new_type) {
        RaftSpecV2R temp{};
        constexpr size_t num_blocks = 4;
        int blocklens[num_blocks] = {1, 1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_DOUBLE};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.triangle_global_index, addr);
        MPI_Get_address(&temp.raft_individual_index, addr + 1);
        MPI_Get_address(&temp.spec_global_index, addr + 2);
        MPI_Get_address(&temp.count, addr + 3);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitRaftSReacInactiveV2R(MPI_Datatype& new_type) {
        RaftSReacInactiveV2R temp{};
        constexpr size_t num_blocks = 3;
        int blocklens[num_blocks] = {1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_STEPS_INDEX};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.triangle_global_index, addr);
        MPI_Get_address(&temp.raft_individual_index, addr + 1);
        MPI_Get_address(&temp.raftsreac_global_index, addr + 2);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitRaftProxyR2V(MPI_Datatype& new_type) {
        RaftProxyR2V temp{};
        constexpr size_t num_blocks = 3;
        int blocklens[num_blocks] = {1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_INT};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.triangle_global_index, addr);
        MPI_Get_address(&temp.raft_individual_index, addr + 1);
        MPI_Get_address(&temp.immobility_update, addr + 2);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitRaftSpecR2V(MPI_Datatype& new_type) {
        RaftSpecR2V temp{};
        const size_t num_blocks = 4;
        int blocklens[num_blocks] = {1, 1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_STEPS_INDEX,
                                              MPI_UNSIGNED};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.triangle_global_index, addr);
        MPI_Get_address(&temp.raft_individual_index, addr + 1);
        MPI_Get_address(&temp.spec_global_index, addr + 2);
        MPI_Get_address(&temp.count, addr + 3);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitRaftGenCountR2V(MPI_Datatype& new_type) {
        RaftGenCountR2V temp{};
        const size_t num_blocks = 3;
        int blocklens[num_blocks] = {1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {MPI_STEPS_INDEX, MPI_STEPS_INDEX, MPI_UNSIGNED};
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.triangle_global_index, addr);
        MPI_Get_address(&temp.raftgen_global_index, addr + 1);
        MPI_Get_address(&temp.count, addr + 2);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitExocytosisEventSync(MPI_Datatype& new_type) {
        steps::solver::ExocytosisEvent temp{};
        constexpr size_t num_blocks = 4;
        int blocklens[num_blocks] = {1, 1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {
            MPI_DOUBLE,
            MPI_STEPS_INDEX,
            MPI_STEPS_INDEX,
            MPI_STEPS_INDEX,
        };
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.time, addr);
        MPI_Get_address(&temp.vidx, addr + 1);
        MPI_Get_address(&temp.tidx, addr + 2);
        MPI_Get_address(&temp.ridx, addr + 3);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitRaftEndocytosisEventSync(MPI_Datatype& new_type) {
        steps::solver::RaftEndocytosisEvent temp{};
        constexpr size_t num_blocks = 4;
        int blocklens[num_blocks] = {1, 1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {
            MPI_DOUBLE,
            MPI_STEPS_INDEX,
            MPI_STEPS_INDEX,
            MPI_STEPS_INDEX,
        };
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.time, addr);
        MPI_Get_address(&temp.ridx, addr + 1);
        MPI_Get_address(&temp.tidx, addr + 2);
        MPI_Get_address(&temp.vidx, addr + 3);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void commitEndocytosisEventSync(MPI_Datatype& new_type) {
        steps::solver::EndocytosisEvent temp{};
        constexpr size_t num_blocks = 3;
        int blocklens[num_blocks] = {1, 1, 1};
        MPI_Datatype old_types[num_blocks] = {
            MPI_DOUBLE,
            MPI_STEPS_INDEX,
            MPI_STEPS_INDEX,
        };
        MPI_Aint addr[num_blocks];
        MPI_Aint disp1[num_blocks];
        MPI_Get_address(&temp.time, addr);
        MPI_Get_address(&temp.tidx, addr + 1);
        MPI_Get_address(&temp.vidx, addr + 2);
        _process(new_type, blocklens, old_types, addr, disp1, num_blocks);
    }

    void _process(MPI_Datatype& new_type,
                  int blocklens[],
                  MPI_Datatype old_types[],
                  MPI_Aint addr[],
                  MPI_Aint disp1[],
                  size_t num_blocks) {
        MPI_Aint lb, extent;
        MPI_Datatype tmp_type;
        disp1[0] = 0;
        for (size_t i = 1; i < num_blocks; i++) {
            disp1[i] = addr[i] - addr[0];
        }
        MPI_Type_create_struct(num_blocks, blocklens, disp1, old_types, &tmp_type);
        MPI_Type_get_extent(tmp_type, &lb, &extent);
        MPI_Type_create_resized(tmp_type, lb, extent, &new_type);
        MPI_Type_commit(&new_type);
    }
};

template <typename I, typename T>
void strong_map_to_vecs(const std::map<I, T>& p_i,
                        std::vector<steps::index_t>& idx_o,
                        std::vector<T>& v_o) {
    idx_o.clear();
    idx_o.reserve(p_i.size());
    v_o.clear();
    v_o.reserve(p_i.size());
    for (const auto& [key, value]: p_i) {
        idx_o.push_back(key.get());
        v_o.push_back(value);
    }
}

template <typename I, typename T>
void vecs_to_strong_map(const std::vector<steps::index_t>& idx_i,
                        const std::vector<T>& v_i,
                        std::map<I, T>& p_o) {
    const std::size_t nentries = idx_i.size();
    p_o.clear();
    for (std::size_t e = 0; e < nentries; e++) {
        p_o.insert_or_assign(I(idx_i[e]), v_i[e]);
    }
}

template <typename T>
void flatten_vecvec(const std::vector<std::vector<T>>& vv,
                    std::vector<T>& fv,
                    std::vector<std::size_t>& sizes) {
    sizes.reserve(sizes.size() + vv.size());
    for (const auto& v: vv) {
        std::size_t vec_entries = v.size();
        sizes.push_back(vec_entries);
        fv.insert(fv.end(), v.begin(), v.end());
    }
}

template <typename T>
void restruct_vecvec(const std::vector<T>& fv,
                     const std::vector<std::size_t>& sizes,
                     std::vector<std::vector<T>>& vv) {
    auto itr = fv.begin();
    std::size_t nsubvecs = sizes.size();
    vv.resize(nsubvecs);
    for (std::size_t v = 0; v < nsubvecs; v++) {
        vv[v].assign(itr, itr + sizes[v]);
        itr += sizes[v];
    }
}

template <typename T>
inline void MPI_BcastVec(std::vector<T>& vec,
                         MPI_Datatype datatype,
                         int source_rank,
                         int my_rank,
                         MPI_Comm communicator,
                         bool fixed_size = false) {
    int n_entries = vec.size();
    if (!fixed_size) {
        MPI_Bcast(&n_entries, 1, MPI_INT, source_rank, communicator);
        if (my_rank != source_rank) {
            vec.resize(n_entries);
        }
    }
    if (n_entries > 0) {
        MPI_Bcast(vec.data(), n_entries, datatype, source_rank, communicator);
    }
}

template <typename T>
inline void MPI_GatherVec(std::vector<T>& vec,
                          MPI_Datatype datatype,
                          int recv_rank,
                          int my_rank,
                          int n_hosts,
                          MPI_Comm communicator) {
    if (my_rank == recv_rank) {
        std::vector<int> n_entries(n_hosts, 0);
        MPI_Gather(MPI_IN_PLACE, 1, MPI_INT, n_entries.data(), 1, MPI_INT, recv_rank, communicator);
        std::vector<int> entries_per_host(n_hosts);
        std::vector<int> dspl(n_hosts);
        uint total_entries = 0;
        for (int p = 0; p < n_hosts; p++) {
            entries_per_host[p] = n_entries[p];
            dspl[p] = total_entries;
            total_entries += entries_per_host[p];
        }
        vec.resize(total_entries);
        MPI_Gatherv(MPI_IN_PLACE,
                    0,
                    datatype,
                    vec.data(),
                    entries_per_host.data(),
                    dspl.data(),
                    datatype,
                    recv_rank,
                    communicator);
    } else {
        int n_entries = vec.size();
        MPI_Gather(&n_entries, 1, MPI_INT, nullptr, 1, MPI_INT, recv_rank, communicator);
        MPI_Gatherv(vec.data(),
                    n_entries,
                    datatype,
                    nullptr,
                    nullptr,
                    nullptr,
                    datatype,
                    recv_rank,
                    communicator);
    }
}

template <typename T>
inline void MPI_ConditionalReduce(T* input,
                                  T* output,
                                  int count,
                                  MPI_Datatype datatype,
                                  MPI_Op op,
                                  bool sync_output = false,
                                  int output_rank = 0,
                                  MPI_Comm communicator = MPI_COMM_WORLD) {
    if (sync_output) {
        MPI_Allreduce(input, output, count, datatype, op, communicator);
    } else {
        MPI_Reduce(input, output, count, datatype, op, output_rank, communicator);
    }
}

template <typename T>
inline void MPI_ConditionalReduce(std::vector<T>& input_vec,
                                  std::vector<T>& output_vec,
                                  MPI_Datatype datatype,
                                  MPI_Op op,
                                  bool sync_output = false,
                                  int output_rank = 0,
                                  MPI_Comm communicator = MPI_COMM_WORLD) {
    output_vec.resize(input_vec.size());
    MPI_ConditionalReduce<T>(input_vec.data(),
                             output_vec.data(),
                             input_vec.size(),
                             datatype,
                             op,
                             sync_output,
                             output_rank,
                             communicator);
}

template <typename T>
inline T MPI_ConditionalReduce(T input,
                               MPI_Datatype datatype,
                               MPI_Op op,
                               bool sync_output = false,
                               int output_rank = 0,
                               MPI_Comm communicator = MPI_COMM_WORLD) {
    T output;
    MPI_ConditionalReduce<T>(
        &input, &output, 1, datatype, op, sync_output, output_rank, communicator);
    return output;
}

template <typename T, typename Tag>
inline util::strong_id<T, Tag> MPI_ConditionalReduce(util::strong_id<T, Tag> input,
                                                     MPI_Datatype datatype,
                                                     MPI_Op op,
                                                     bool sync_output = false,
                                                     int output_rank = 0,
                                                     MPI_Comm communicator = MPI_COMM_WORLD) {
    return util::strong_id<T, Tag>(
        MPI_ConditionalReduce(input.get(), datatype, op, sync_output, output_rank, communicator));
}

template <typename T>
inline void MPI_ConditionalBcast(T* data,
                                 size_t count,
                                 MPI_Datatype datatype,
                                 int source_host,
                                 int my_rank,
                                 bool sync_output = false,
                                 int output_rank = 0,
                                 MPI_Comm communicator = MPI_COMM_WORLD) {
    if (sync_output) {
        MPI_Bcast(data, count, datatype, source_host, communicator);
    } else if (source_host != output_rank) {
        if (my_rank == source_host) {
            MPI_Send(data, count, datatype, output_rank, MPI_CONDITIONAL_BCAST, communicator);
        } else if (my_rank == output_rank) {
            MPI_Recv(data,
                     count,
                     datatype,
                     source_host,
                     MPI_CONDITIONAL_BCAST,
                     communicator,
                     MPI_STATUS_IGNORE);
        }
    }
}

template <typename T>
inline void MPI_ConditionalBcast(std::vector<T>& data_vec,
                                 size_t count,
                                 MPI_Datatype datatype,
                                 int source_host,
                                 int my_rank,
                                 bool sync_output = false,
                                 int output_rank = 0,
                                 MPI_Comm communicator = MPI_COMM_WORLD) {
    if (my_rank != source_host) {
        data_vec.resize(count);
    }
    MPI_ConditionalBcast<T>(data_vec.data(),
                            count,
                            datatype,
                            source_host,
                            my_rank,
                            sync_output,
                            output_rank,
                            communicator);
}

template <typename T>
inline T MPI_ConditionalBcast(T data,
                              MPI_Datatype datatype,
                              int source_host,
                              int my_rank,
                              bool sync_output = false,
                              int output_rank = 0,
                              MPI_Comm communicator = MPI_COMM_WORLD) {
    MPI_ConditionalBcast<T>(
        &data, 1, datatype, source_host, my_rank, sync_output, output_rank, communicator);
    return data;
}

template <typename T, typename Tag>
inline util::strong_id<T, Tag> MPI_ConditionalBcast(util::strong_id<T, Tag> id,
                                                    MPI_Datatype datatype,
                                                    int source_host,
                                                    int my_rank,
                                                    bool sync_output = false,
                                                    int output_rank = 0,
                                                    MPI_Comm communicator = MPI_COMM_WORLD) {
    MPI_ConditionalBcast<T>(
        &id.reference(), 1, datatype, source_host, my_rank, sync_output, output_rank, communicator);
    return id;
}

template <typename T>
inline void MPI_ConditionalBcast(std::vector<T>& data_vec,
                                 MPI_Datatype datatype,
                                 int source_host,
                                 int my_rank,
                                 bool sync_output = false,
                                 int output_rank = 0,
                                 MPI_Comm communicator = MPI_COMM_WORLD) {
    auto count = MPI_ConditionalBcast(data_vec.size(),
                                      MPI_STD_SIZE_T,
                                      source_host,
                                      my_rank,
                                      sync_output,
                                      output_rank,
                                      communicator);
    if (my_rank != source_host) {
        data_vec.resize(count);
    }
    MPI_ConditionalBcast<T>(data_vec.data(),
                            count,
                            datatype,
                            source_host,
                            my_rank,
                            sync_output,
                            output_rank,
                            communicator);
}

template <typename I, typename T>
inline void MPI_ConditionalBcast(std::map<I, T>& data_map,
                                 MPI_Datatype datatype,
                                 int source_host,
                                 int my_rank,
                                 bool sync_output = false,
                                 int output_rank = 0,
                                 MPI_Comm communicator = MPI_COMM_WORLD) {
    std::vector<index_t> indices;
    std::vector<T> values;
    std::size_t nentries = MPI_ConditionalBcast<std::size_t>(data_map.size(),
                                                             MPI_STD_SIZE_T,
                                                             source_host,
                                                             my_rank,
                                                             sync_output,
                                                             output_rank,
                                                             communicator);


    if (my_rank == source_host) {
        strong_map_to_vecs<I, T>(data_map, indices, values);
    } else {
        indices.resize(nentries);
        values.resize(nentries);
    }
    MPI_ConditionalBcast<steps::index_t>(indices.data(),
                                         nentries,
                                         MPI_STEPS_INDEX,
                                         source_host,
                                         my_rank,
                                         sync_output,
                                         output_rank,
                                         communicator);
    MPI_ConditionalBcast<T>(values.data(),
                            nentries,
                            datatype,
                            source_host,
                            my_rank,
                            sync_output,
                            output_rank,
                            communicator);
    if (my_rank != source_host) {
        vecs_to_strong_map<I, T>(indices, values, data_map);
    }
}

template <typename T = steps::index_t, typename P>
inline void MPI_ConditionalBcast(std::vector<util::strong_id<T, P>>& strong_data_vec,
                                 int source_host,
                                 int my_rank,
                                 bool sync_output = false,
                                 int output_rank = 0,
                                 MPI_Comm communicator = MPI_COMM_WORLD) {
    if (my_rank == source_host) {
        std::vector<steps::index_t> sync_vec = strong_type_to_value_type(strong_data_vec);
        auto nentries = MPI_ConditionalBcast<std::size_t>(sync_vec.size(),
                                                          MPI_STD_SIZE_T,
                                                          source_host,
                                                          my_rank,
                                                          sync_output,
                                                          output_rank,
                                                          communicator);
        MPI_ConditionalBcast<steps::index_t>(sync_vec.data(),
                                             nentries,
                                             MPI_STEPS_INDEX,
                                             source_host,
                                             my_rank,
                                             sync_output,
                                             output_rank,
                                             communicator);
    } else {
        auto nentries = MPI_ConditionalBcast<std::size_t>(
            0u, MPI_STD_SIZE_T, source_host, my_rank, sync_output, output_rank, communicator);
        std::vector<steps::index_t> sync_vec(nentries);
        MPI_ConditionalBcast<steps::index_t>(sync_vec.data(),
                                             nentries,
                                             MPI_STEPS_INDEX,
                                             source_host,
                                             my_rank,
                                             sync_output,
                                             output_rank,
                                             communicator);
        strong_data_vec.clear();
        strong_data_vec.reserve(sync_vec.size());
        for (auto e: sync_vec) {
            strong_data_vec.emplace_back(e);
        }
    }
}

template <typename T>
inline void MPI_ConditionalBcast(std::vector<std::vector<T>>& mv,
                                 MPI_Datatype datatype,
                                 int source_host,
                                 int my_rank,
                                 bool sync_output = false,
                                 int output_rank = 0,
                                 MPI_Comm communicator = MPI_COMM_WORLD) {
    std::vector<T> fv;
    std::vector<std::size_t> sizes;
    if (my_rank == source_host) {
        flatten_vecvec(mv, fv, sizes);
        std::array<std::size_t, 2> nentries{fv.size(), sizes.size()};
        MPI_ConditionalBcast<std::size_t>(nentries.data(),
                                          nentries.size(),
                                          MPI_STD_SIZE_T,
                                          source_host,
                                          my_rank,
                                          sync_output,
                                          output_rank,
                                          communicator);
        MPI_ConditionalBcast<T>(fv.data(),
                                nentries[0],
                                datatype,
                                source_host,
                                my_rank,
                                sync_output,
                                output_rank,
                                communicator);
        MPI_ConditionalBcast<std::size_t>(sizes.data(),
                                          nentries[1],
                                          MPI_STD_SIZE_T,
                                          source_host,
                                          my_rank,
                                          sync_output,
                                          output_rank,
                                          communicator);
    } else {
        std::array<std::size_t, 2> nentries{0, 0};
        MPI_ConditionalBcast<std::size_t>(nentries.data(),
                                          nentries.size(),
                                          MPI_STD_SIZE_T,
                                          source_host,
                                          my_rank,
                                          sync_output,
                                          output_rank,
                                          communicator);
        fv.resize(nentries[0]);
        sizes.resize(nentries[1]);
        MPI_ConditionalBcast<T>(fv.data(),
                                nentries[0],
                                datatype,
                                source_host,
                                my_rank,
                                sync_output,
                                output_rank,
                                communicator);
        MPI_ConditionalBcast<std::size_t>(sizes.data(),
                                          nentries[1],
                                          MPI_STD_SIZE_T,
                                          source_host,
                                          my_rank,
                                          sync_output,
                                          output_rank,
                                          communicator);
        restruct_vecvec(fv, sizes, mv);
    }
}

}  // namespace steps::mpi::tetvesicle
