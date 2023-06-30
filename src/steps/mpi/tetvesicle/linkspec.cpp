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

#include "mpi/tetvesicle/linkspec.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "math/point.hpp"
#include "mpi/mpi_common.hpp"
#include "mpi/tetvesicle/tetvesicle_vesraft.hpp"
#include "solver/linkspecdef.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::mpi::tetvesicle {

////////////////////////////////////////////////////////////////////////////////

// Just chuck this here for now
void cross_prod(double a1,
                double a2,
                double a3,
                double b1,
                double b2,
                double b3,
                double& p1,
                double& p2,
                double& p3) {
    p1 = (a2 * b3 - a3 * b2);
    p2 = (a3 * b1 - a1 * b3);
    p3 = (a1 * b2 - a2 * b1);
}

double dot_prod(double a1, double a2, double a3, double b1, double b2, double b3) {
    return a1 * b1 + a2 * b2 + a3 * b3;
}

////////////////////////////////////////////////////////////////////////////////

LinkSpec::LinkSpec(solver::LinkSpecdef* lspecdef,
                   solver::linkspec_individual_id linkspec,
                   Vesicle* ves,
                   tetrahedron_global_id tet_gidx,
                   math::position_rel_to_ves spec_pos_rel)
    : pDef(lspecdef)
    , pUniqueID(linkspec)
    , pPosCartesian_rel(spec_pos_rel)
    , pLinked(false)
    , pPair(nullptr)
    , pVesicle(ves)
    , pTetOverlap(tet_gidx) {
    AssertLog(pDef != nullptr);
    AssertLog(pVesicle != nullptr);

    pPosSpherical.setRadius(ves->getDiam() / 2.0);
    pPosSpherical.setTheta(atan2(pPosCartesian_rel[1], pPosCartesian_rel[0]));
    pPosSpherical.setPhi(acos(pPosCartesian_rel[2] / pPosSpherical.getRadius()));

    pVesicle->comp()->solverVesRaft()->recordLinkSpec_(pUniqueID, this);
}

////////////////////////////////////////////////////////////////////////////////

LinkSpec::LinkSpec(solver::LinkSpecdef* lspecdef, Vesicle* ves, std::fstream& cp_file)
    : pDef(lspecdef)
    , pLinked(false)
    , pPair(nullptr)
    , pVesicle(ves) {
    AssertLog(pDef != nullptr);
    AssertLog(pVesicle != nullptr);

    util::restore(cp_file, pUniqueID);
    util::restore(cp_file, pPosCartesian_rel);
    util::restore(cp_file, pTetOverlap);

    pPosSpherical.setRadius(ves->getDiam() / 2.0);
    pPosSpherical.setTheta(atan2(pPosCartesian_rel[1], pPosCartesian_rel[0]));
    pPosSpherical.setPhi(acos(pPosCartesian_rel[2] / pPosSpherical.getRadius()));

    pVesicle->comp()->solverVesRaft()->recordLinkSpec_(pUniqueID, this);
}

////////////////////////////////////////////////////////////////////////////////

LinkSpec::~LinkSpec() {
    pVesicle->comp()->solverVesRaft()->removeLinkSpec_(pUniqueID, this);
    if (pLinked) {
        auto ls = getLinkedSpec();
        // Isolate the other end of the link so that the pair is only removed once
        ls->removeLinkSpecPair();  // sets pLinked = false on other LinkSpec
        // Remove the link spec pair from the solver
        pVesicle->comp()->solverVesRaft()->removeLinkSpecPair_(pPair);
        // Remove the link species from the vesicle
        ls->getVesicle()->remLinkSpec(ls->getUniqueID(), ls);
    }
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpec::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pUniqueID);
    util::checkpoint(cp_file, pPosCartesian_rel);
    util::checkpoint(cp_file, pTetOverlap);
    // pLinked and pPair are taken care of by LinkSpecPair checkpoint/restore
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpec::setDef(solver::LinkSpecdef* linkspec_def) noexcept {
    pDef = linkspec_def;
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpec::addLinkSpecPair(LinkSpecPair* linkspecpair) {
    if (pLinked) {
        ProgErrLog("Trying to add partner to LinkSpecies that is already paired.");
    }
    pPair = linkspecpair;
    pLinked = true;
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpec::removeLinkSpecPair() {
    if (!pLinked) {
        ProgErrLog("Trying to remove partner from LinkSpecies that is not paired.");
    }
    pPair = nullptr;
    pLinked = false;
}

////////////////////////////////////////////////////////////////////////////////

math::position_abs LinkSpec::getPosCartesian_abs() const {
    math::position_abs vespos = pVesicle->getPosition();
    math::position_abs pos_abs{pPosCartesian_rel[0] + vespos[0],
                               pPosCartesian_rel[1] + vespos[1],
                               pPosCartesian_rel[2] + vespos[2]};
    return pos_abs;
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpec::setPosCartesian_rel(const math::position_rel_to_ves& pos_rel) {
    double radius = pPosSpherical.getRadius();
    pPosSpherical.setTheta(atan2(pos_rel[1], pos_rel[0]));
    pPosSpherical.setPhi(acos(pos_rel[2] / radius));

    pPosCartesian_rel = pos_rel;

    return;
}

////////////////////////////////////////////////////////////////////////////////

double LinkSpec::getLength() const {
    auto posA_abs = getPosCartesian_abs();
    auto posB_abs = getLinkedSpec()->getPosCartesian_abs();

    return math::distance(posA_abs, posB_abs);
}

////////////////////////////////////////////////////////////////////////////////

bool LinkSpec::movePosAllowed(const math::point3d& move_vector) const {
    if (!pLinked) {
        ProgErrLog("LinkSpecies has not been paired.");
    }

    auto posA_abs = getPosCartesian_abs();

    math::position_abs posA_abs_new{posA_abs[0] + move_vector[0],
                                    posA_abs[1] + move_vector[1],
                                    posA_abs[2] + move_vector[2]};

    auto posB_abs = getLinkedSpec()->getPosCartesian_abs();

    double newlength = math::distance(posA_abs_new, posB_abs);

    return (newlength <= pPair->max_length() && newlength >= pPair->min_length());
}

////////////////////////////////////////////////////////////////////////////////

LinkSpec* LinkSpec::getLinkedSpec() const {
    if (!pLinked) {
        ProgErrLog("LinkSpecies has not been paired.");
    }
    return pPair->getPairedLinkSpec(this);
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpec::updatePos(double theta, double phi) {
    double radius = pPosSpherical.getRadius();
    // double theta_old = pPosSpherical[1];
    double phi_old = pPosSpherical.getPhi();

    double x_start = pPosCartesian_rel[0];
    double y_start = pPosCartesian_rel[1];
    double z_start = pPosCartesian_rel[2];

    // find the cartesians coordinates of the new positions relative to north pole
    double x_rot = radius * sin(phi) * cos(theta);
    double y_rot = radius * sin(phi) * sin(theta);
    double z_rot = radius * cos(phi);

    // Find the unit rotation axis, which is the cross product of the 'starting'
    // positions and the z-axis (if phi is relative to z-axis)
    double rotationaxis_x, rotationaxis_y,
        rotationaxis_z;  // z is overkill because it will always be zero

    cross_prod(x_start, y_start, z_start, 0, 0, 1, rotationaxis_x, rotationaxis_y, rotationaxis_z);

    double mod_rot_axis = sqrt((rotationaxis_x * rotationaxis_x) +
                               (rotationaxis_y * rotationaxis_y) +
                               (rotationaxis_z * rotationaxis_z));

    if (mod_rot_axis != 0.0) {
        rotationaxis_x = rotationaxis_x / mod_rot_axis;
        rotationaxis_y = rotationaxis_y / mod_rot_axis;
        rotationaxis_z = rotationaxis_z / mod_rot_axis;
    }

    // The _rot values are a starting vector that must be rotated around the
    // rotation axis to get the final position. T do that we need to find the
    // cross product and dot product:

    double cross_x, cross_y, cross_z;
    cross_prod(x_rot,
               y_rot,
               z_rot,
               rotationaxis_x,
               rotationaxis_y,
               rotationaxis_z,
               cross_x,
               cross_y,
               cross_z);
    double dot_p = dot_prod(rotationaxis_x, rotationaxis_y, rotationaxis_z, x_rot, y_rot, z_rot);

    // Apply Rodrigues' formula
    double v_x = x_rot * cos(phi_old) + cross_x * (sin(phi_old)) +
                 rotationaxis_x * dot_p * (1 - cos(phi_old));
    double v_y = y_rot * cos(phi_old) + cross_y * (sin(phi_old)) +
                 rotationaxis_y * dot_p * (1 - cos(phi_old));
    double v_z = z_rot * cos(phi_old) + cross_z * (sin(phi_old));

    pPosCartesian_rel[0] = v_x;
    pPosCartesian_rel[1] = v_y;
    pPosCartesian_rel[2] = v_z;

    pPosSpherical.setTheta(atan2(v_y, v_x));
    pPosSpherical.setPhi(acos(v_z / radius));
}

////////////////////////////////////////////////////////////////////////////////

bool LinkSpec::withinBounds() const {
    if (!pLinked) {
        ProgErrLog("LinkSpecies has not been paired.");
    }

    double newlength = getLength();
    return (newlength <= pPair->max_length() && newlength >= pPair->min_length());
}

}  // namespace steps::mpi::tetvesicle
