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

#include "init.hpp"
#include <easylogging++.h>

#ifdef OMEGA_H_USE_GMSH
#include <gmsh.h>
#endif  // OMEGA_H_USE_GMSH

#include "util/profile/profiler_interface.hpp"

/*
HOWTO: Logging, Assertion and Exception

1. To enable assertion logging, make sure ENABLE_ASSERTLOG is defined in
src/CMakeLists.txt
2. Make sure "easylogging++.h" and "steps/error.hpp" are included as header
files
3. Use AssertLog(assert_condition) to add assertion, for example
            int a = -1;
            AssertLog(a > 0);
4. Use error logs like ErrLog(msg), NotImplErrLog(msg), ArgErrLog(msg),
   ProgErrLog(msg), SysErrLog(msg) and IOErrLog(msg), or condition error logs like
   ErrLogIf(error_condition, msg), NotImplErrLogIf(error_condition, msg),
   ArgErrLogIf(error_condition, msg), ProgErrLogIf(error_condition, msg),
   SysErrLogIf(error_condition, msg) and IOErrLogIf(error_condition, msg)
   to log and throw associated exception, for example
            NotImplErrLog("This function is not implemented in this solver.");
            int i = 4;
            ArgErrLogIf(i != 2, "i should be 2");
5. Use CLOG(WARNING, "general_log") to log warning message, for example
            CLOG(WARNING, "general_log") << "This is a warning message.";
6. Use CLOG(INFO, "general_log") to log information message, for example
            CLOG(INFO, "general_log") << "This is a general info message.";

NOTE: WARNING and INFO mesages do not terminate the simulation.
      INFO messages are not displayed in terminal under parallel mode.
*/

INITIALIZE_EASYLOGGINGPP

void steps::init() {
    // easylogging initialization
    el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);
    el::Loggers::addFlag(el::LoggingFlag::CreateLoggerAutomatically);
    el::Loggers::addFlag(el::LoggingFlag::AutoSpacing);
    el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);

    // This is the default log for serial simulation
    // This configuration will be rewritten if steps.mpi module is imported in
    // parallel simulation
    el::Configurations serial_conf;

    serial_conf.set(el::Level::Global,
                    el::ConfigurationType::Format,
                    "[%datetime][%level][%loc][%func]: %msg");
    serial_conf.set(el::Level::Global, el::ConfigurationType::ToStandardOutput, "false");
    serial_conf.set(el::Level::Global, el::ConfigurationType::ToFile, "true");
    serial_conf.set(el::Level::Global, el::ConfigurationType::Filename, ".logs/general_log_0.txt");
    serial_conf.set(el::Level::Global, el::ConfigurationType::MaxLogFileSize, "2097152");

    serial_conf.set(el::Level::Fatal, el::ConfigurationType::ToStandardOutput, "true");
    serial_conf.set(el::Level::Error, el::ConfigurationType::ToStandardOutput, "true");
    serial_conf.set(el::Level::Warning, el::ConfigurationType::ToStandardOutput, "true");
    serial_conf.set(el::Level::Info, el::ConfigurationType::ToStandardOutput, "true");

    el::Loggers::getLogger("general_log");
    el::Loggers::reconfigureLogger("general_log", serial_conf);

#ifdef OMEGA_H_USE_GMSH
    ::gmsh::initialize();
    // change the verbosity level
    // (0: silent except for fatal errors, 1: +errors, 2: +warnings, 3: +direct, 4: +information, 5:
    // +status, 99: +debug)
    ::gmsh::option::setNumber("General.Verbosity", 2);
#endif  // OMEGA_H_USE_GMSH

    steps::Instrumentor::init_profile();
}
