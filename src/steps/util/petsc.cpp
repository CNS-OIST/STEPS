/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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

#include "petsc.hpp"
#include <cassert>
#include <memory>
#include <petscsys.h>

namespace steps::util::petsc {

static std::unique_ptr<PetscFixture> fixture{nullptr};

PetscFixture::PetscFixture(int* argc, char*** argv, const char file[], const char help[]) {
    PetscErrorCode ierr = PetscInitialize(argc, argv, file, help);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

PetscFixture::~PetscFixture() {
    int flag;
    MPI_Finalized(&flag);
    if (flag == 0) {
        PetscErrorCode ierr = PetscFinalize();
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
    }
}

void Init(int* argc, char*** argv, const char file[], const char help[]) {
    if (not fixture) {
        fixture = std::make_unique<PetscFixture>(argc, argv, file, help);
    }
}

void Finalize() {
    fixture = nullptr;
}

}  // namespace steps::util::petsc
