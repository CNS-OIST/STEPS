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

#ifndef STEPS_EXTERN
#define STEPS_EXTERN extern
#endif

#include <cstdint>

#ifdef __cplusplus
STEPS_EXTERN "C" {
#endif
/* __cplusplus */

/* ************************************************************************** */

/* General contact. */
#define STEPS_CONTACT_GENERAL "steps.dev@gmail.com"

/* Primary author. */
#define STEPS_CONTACT_FIRSTAUTHOR "steps.dev@gmail.com"

/* Use this address for composing (error) messages that advise a user to
 * contact the development ...uh... 'team'. */
#define STEPS_CONTACT_BUGREPORT "steps.dev@gmail.com"

/* Official website URL. */
#define STEPS_CONTACT_HOMEPAGE "steps.sourceforge.net"

    /* ************************************************************************** */

    /* Abbrevations for unsigned versions of integral types. */
    typedef unsigned char uchar;
    typedef unsigned short int ushort;
    typedef unsigned int uint;
    typedef unsigned long int ulong;

    /* ************************************************************************** */

#ifdef __cplusplus
}
#endif
/* __cplusplus */

#ifdef __cplusplus

#ifdef STEPS_USE_DIST_MESH
#include <Omega_h_defines.hpp>

#define INITIAL_COMPARTMENT_ID -1

namespace steps {
namespace osh = Omega_h;
}  // namespace steps

#endif  // !STEPS_USE_DIST_MESH

#if defined(__clang__)
#define STEPS_FALLTHROUGH [[clang::fallthrough]]
/* Test for GCC >= 7.0.0 */
#elif defined(__GNUC__) && (__GNUC__ > 7 || __GNUC__ == 7)
#define STEPS_FALLTHROUGH [[gnu::fallthrough]]
#else
#define STEPS_FALLTHROUGH ((void) 0)
#endif

namespace steps {

#ifdef STEPS_USE_64BITS_INDICES
using index_t = std::uint64_t;
#else
using index_t = std::uint32_t;
#endif

inline const double SPERM_WHALE_BRAIN_VOLUME = 1e-2;  // cubic meters
inline const double DEFAULT_MEMB_POT = -65e-3;        // Default membrane potential in V

}  // namespace steps

#endif  //__cplusplus
