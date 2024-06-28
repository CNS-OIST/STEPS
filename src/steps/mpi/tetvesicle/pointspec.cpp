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

#include "mpi/tetvesicle/pointspec.hpp"

// logging

// STEPS headers.
#include "math/constants.hpp"
#include "math/point.hpp"
#include "math/tools.hpp"
#include "util/checkpointing.hpp"
#include "util/error.hpp"

namespace steps::mpi::tetvesicle {

////////////////////////////////////////////////////////////////////////////////

PointSpec::PointSpec(solver::spec_global_id gidx,
                     double radius,
                     solver::pointspec_individual_id idx)
    : pIdx(idx)
    , pSpec_gidx(gidx) {
    pPosSpherical.setRadius(radius);
}

////////////////////////////////////////////////////////////////////////////////

PointSpec::PointSpec(solver::spec_global_id gidx, double radius, std::fstream& cp_file)
    : pSpec_gidx(gidx) {
    util::restore(cp_file, pIdx);
    util::restore(cp_file, pSpec_gidx);
    util::restore(cp_file, pPosSpherical);
    util::restore(cp_file, pPosCartesian);
    util::restore(cp_file, pTetOverlap);
}

////////////////////////////////////////////////////////////////////////////////

PointSpec::~PointSpec() = default;

////////////////////////////////////////////////////////////////////////////////

void PointSpec::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pIdx);
    util::checkpoint(cp_file, pSpec_gidx);
    util::checkpoint(cp_file, pPosSpherical);
    util::checkpoint(cp_file, pPosCartesian);
    util::checkpoint(cp_file, pTetOverlap);
}

////////////////////////////////////////////////////////////////////////////////

void PointSpec::setPosCartesian(math::position_rel_to_ves const& pos) noexcept {
    double radius = pPosSpherical.getRadius();
    pPosSpherical.setTheta(atan2(pos[1], pos[0]));
    pPosSpherical.setPhi(acos(pos[2] / radius));

    pPosCartesian = pos;
}

////////////////////////////////////////////////////////////////////////////////

void PointSpec::setPosSpherical(double theta, double phi) {
    if (theta < -math::PI or theta > math::PI) {
        ProgErrLog("Theta is outside range -pi <-> pi.\n");
    }
    if (phi < 0.0 or phi > math::PI) {
        ProgErrLog("Phi is outside range 0<->pi.\n");
    }

    pPosSpherical.setTheta(theta);
    pPosSpherical.setPhi(phi);

    pPosCartesian[0] = pPosSpherical.getRadius() * sin(phi) * cos(theta);
    pPosCartesian[1] = pPosSpherical.getRadius() * sin(phi) * sin(theta);
    pPosCartesian[2] = pPosSpherical.getRadius() * cos(phi);
}

////////////////////////////////////////////////////////////////////////////////

void PointSpec::updatePos(double theta, double phi) {
    // if (theta < -math::PI || theta > math::PI) {
    if (theta < 0.0 || theta > 2 * math::PI) {
        ProgErrLog("Theta is outside range -0 <-> 2pi.\n");
    }
    if (phi < 0.0 || phi > math::PI) {
        ProgErrLog("Phi is outside range 0 <-> pi.\n");
    }

    double radius = pPosSpherical.getRadius();
    // double theta_old = pPosSpherical[1];
    double phi_old = pPosSpherical.getPhi();

    double x_start = pPosCartesian[0];
    double y_start = pPosCartesian[1];
    double z_start = pPosCartesian[2];

    // find the cartesians coordinates of the new positions relative to north pole
    double x_rot = radius * sin(phi) * cos(theta);
    double y_rot = radius * sin(phi) * sin(theta);
    double z_rot = radius * cos(phi);

    // Find the unit rotation axis, which is the cross product of the 'starting'
    // positions and the z-axis (if phi is relative to z-axis)
    double rotationaxis_x, rotationaxis_y,
        rotationaxis_z;  // z is overkill because it will always be zero

    math::cross_product(
        x_start, y_start, z_start, 0, 0, 1, rotationaxis_x, rotationaxis_y, rotationaxis_z);

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
    math::cross_product(x_rot,
                        y_rot,
                        z_rot,
                        rotationaxis_x,
                        rotationaxis_y,
                        rotationaxis_z,
                        cross_x,
                        cross_y,
                        cross_z);
    double dot_p =
        math::dot_product(rotationaxis_x, rotationaxis_y, rotationaxis_z, x_rot, y_rot, z_rot);

    // Apply Rodrigues' formula
    double v_x = x_rot * cos(phi_old) + cross_x * (sin(phi_old)) +
                 rotationaxis_x * dot_p * (1 - cos(phi_old));
    double v_y = y_rot * cos(phi_old) + cross_y * (sin(phi_old)) +
                 rotationaxis_y * dot_p * (1 - cos(phi_old));
    double v_z = z_rot * cos(phi_old) + cross_z * (sin(phi_old));

    pPosCartesian[0] = v_x;
    pPosCartesian[1] = v_y;
    pPosCartesian[2] = v_z;

    pPosSpherical.setTheta(atan2(v_y, v_x));

    pPosSpherical.setPhi(acos(v_z / radius));
}

}  // namespace steps::mpi::tetvesicle
