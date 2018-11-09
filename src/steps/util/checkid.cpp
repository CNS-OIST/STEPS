/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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

#include <string>

#include "steps/error.hpp"
#include "steps/util/checkid.hpp"
// logging
#include "easylogging++.h"
namespace steps {
namespace util {

static inline bool ascii_is_alpha(char c) {
    return c >= 'A' && c <= 'Z' || c >= 'a' && c <= 'z' || c=='_';
}

static inline bool ascii_is_alphanum(char c) {
    return ascii_is_alpha(c) || c >='0' && c <='9';
}

bool isValidID(const char *s) {
    if (!ascii_is_alpha(*s)) { return false;
}
    while (*++s != 0) { if (!ascii_is_alphanum(*s)) { return false;
}
}

    return true;
}

bool isValidID(const std::string &s) {
    return isValidID(s.c_str());
}

void checkID(const char *s) {
    if (!isValidID(s)) {ArgErrLog("'"+std::string(s)+"' is not a valid id.");}
}

void checkID(const std::string &s) {
    if (!isValidID(s)) {ArgErrLog("'"+s+"' is not a valid id.");}
}

}  // namespace util
}  // namespace steps
