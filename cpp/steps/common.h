/*
 * STEPS - STochastic Engine for Pathway Simulation
 * Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
 *
 * This file is part of STEPS.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
 * 
 * $Id$
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
#define STEPS_CONTACT_GENERAL       "wils@oist.jp"

/* Primary author. */
#define STEPS_CONTACT_FIRSTAUTHOR   "wils@oist.jp"

/* Use this address for composing (error) messages that advise a user to
 * contact the development ...uh... 'team'. */
#define STEPS_CONTACT_BUGREPORT     "wils@oist.jp"

/* Official website URL. */
#define STEPS_CONTACT_HOMEPAGE      "http://www.tnb.ua.ac.be/steps"

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
