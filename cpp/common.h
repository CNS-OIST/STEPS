/*
////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
*/

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

#ifndef STEPS_COMMON_H
#define STEPS_COMMON_H 1

/* ************************************************************************** */

#ifndef STEPS_EXTERN
#define STEPS_EXTERN extern
#endif

#ifdef __cplusplus
STEPS_EXTERN "C" {
#endif
/* __cplusplus */

/* ************************************************************************** */

/* General contact. */
#define STEPS_CONTACT_GENERAL       "steps.dev@gmail.com"

/* Primary author. */
#define STEPS_CONTACT_FIRSTAUTHOR   "steps.dev@gmail.com"

/* Use this address for composing (error) messages that advise a user to
 * contact the development ...uh... 'team'. */
#define STEPS_CONTACT_BUGREPORT     "steps.dev@gmail.com"

/* Official website URL. */
#define STEPS_CONTACT_HOMEPAGE      "steps.sourceforge.net"

/* ************************************************************************** */

#ifdef __cplusplus
#define START_NAMESPACE(X) namespace X {
#define END_NAMESPACE(X) }
#define NAMESPACE_ALIAS(X,Y) namespace Y = X
#define USING(X,Y) using X::Y
#define USING_NAMESPACE(X) using namespace X
#endif
/* __cplusplus */

/* ************************************************************************** */

/* Abbrevations for unsigned versions of integral types. */
typedef unsigned char                           uchar;
typedef unsigned short int                      ushort;
typedef unsigned int                            uint;
typedef unsigned long int                       ulong;

/* ************************************************************************** */

#ifdef __cplusplus
}
#endif
/* __cplusplus */

#endif
/* STEPS_COMMON_H */

/* END */
