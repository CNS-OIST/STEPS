/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_UTIL_CHECKID_HPP
#define STEPS_UTIL_CHECKID_HPP 1

#include <string>

namespace steps {
namespace util {

/** Test string for validity as steps object identifier.
 *
 * \param s  ID string.
 * \return   True if valid.
 *
 * Valid ID strings consist of an alphabetic character [A-Za-z_]
 * followed by a possibly empty sequence of alphanumeric characters
 * [A-Za-z_0-9].
 */
bool isValidID(const char *s);

/** Test string for validity as steps object identifier.
 *
 * \param s  ID string.
 * \return   True if valid.
 *
 * See isValidID(const char *)
 */
bool isValidID(const std::string &s);

/** Throw exception if string is an invalid identifier.
 *
 * \param s  ID string.
 *
 * Throws steps::ArgErr if isValidID(s) is false.
 */
void checkID(const char *s);

/** Throw exception if string is an invalid identifier.
 *
 * \param s  ID string.
 *
 * Throws steps::ArgErr if isValidID(s) is false.
 */
void checkID(const std::string &s);

}} // namespace steps::util


#endif // ndef STEPS_UTIL_CHECKID_HPP
