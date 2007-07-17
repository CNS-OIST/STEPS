/*
 * STEPS - STochastic Engine for Pathway Simulation
 * Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
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
#define STEPS_CONTACT_GENERAL       "stefan@tnb.ua.ac.be"

/* Primary author. */
#define STEPS_CONTACT_FIRSTAUTHOR   "stefan@tnb.ua.ac.be"

/* Use this address for composing (error) messages that advise a user to
 * contact the development ...uh... 'team'. */
#define STEPS_CONTACT_BUGREPORT     "stefan@tnb.ua.ac.be"

/* Official website URL. */
#define STEPS_CONTACT_HOMEPAGE      "http://www.tnb.ua.ac.be/steps"

/* ************************************************************************** */

#ifdef __cplusplus
#define START_NAMESPACE(X) namespace X {
#define END_NAMESPACE(X) }
#define NAMESPACE_ALIAS(X,Y) namespace Y = X
#define USING(X,Y) using X::Y
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
