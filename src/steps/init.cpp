/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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

#include "third_party/easylogging++.h"
#include "steps/init.hpp"

_INITIALIZE_EASYLOGGINGPP

void steps::init(void) {

    // easylogging initialization
    el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);
    el::Loggers::addFlag(el::LoggingFlag::CreateLoggerAutomatically);
    
    // STEPS DEBUG Logger
    el::Configurations debug_conf;
    debug_conf.set(el::Level::Debug, el::ConfigurationType::Format, "[%datetime][%logger][%loc]: %msg");
    debug_conf.set(el::Level::Debug,
             el::ConfigurationType::ToStandardOutput, "false");
    
    // DEBUG level logging only applies if NDEBUG not defined, so
    // create log file only if NDEBUG not defined (or for Windows, _DEBUG is defined)

#if !defined(NDEBUG) || defined(_DEBUG)
    debug_conf.set(el::Level::Debug,
             el::ConfigurationType::ToFile, "true");
    debug_conf.set(el::Level::Debug,
                   el::ConfigurationType::Filename, ".logs/debug_log");
#endif
    
    el::Loggers::getLogger("steps_debug");
    el::Loggers::reconfigureLogger("steps_debug", debug_conf);
}

