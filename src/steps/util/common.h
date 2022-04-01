/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

#ifdef __cplusplus

#include <istream>
#include <ostream>
#include <vector>

#ifdef STEPS_USE_DIST_MESH
#include <Omega_h_defines.hpp>

#define INITIAL_COMPARTMENT_ID -1

namespace steps {
namespace osh = Omega_h;
} // namespace steps

#endif // !STEPS_USE_DIST_MESH

#if defined(__clang__)
#define STEPS_FALLTHROUGH [[clang::fallthrough]]
/* Test for GCC >= 7.0.0 */
#elif defined(__GNUC__) && (__GNUC__ > 7 || __GNUC__ == 7)
#define STEPS_FALLTHROUGH [[gnu::fallthrough]]
#else
#define STEPS_FALLTHROUGH ((void)0)
#endif

namespace steps {

template<typename T>
void checkpoint(std::ostream &istr, std::vector<T> &v, bool with_size = true) {
  if (with_size) {
    auto size = v.size();
    istr.write(reinterpret_cast<char *>(size), sizeof(decltype(size)));
  }
  istr.write(reinterpret_cast<char *>(v.data()), sizeof(T) * v.size());
}

template<typename T>
void restore(std::istream &istr, uint nelems, std::vector<T> &v) {
  v.resize(nelems);
  istr.read(reinterpret_cast<char *>(v.data()), sizeof(T) * nelems);
}

template<typename T>
void restore(std::istream &istr, std::vector<T> &v) {
  size_t nelems;
  istr.read(reinterpret_cast<char *>(&nelems), sizeof(size_t));
  restore(istr, nelems, v);
}

} // namespace steps

#endif //__cplusplus

#endif
/* STEPS_COMMON_H */

/* END */
