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

#ifndef STEPS_ERROR_HPP
#define STEPS_ERROR_HPP 1

#include <easylogging++.h>

// Standard library & STL headers.
#include <exception>
#include <sstream>
#include <string>

// STEPS headers.
#include "common.h"

#ifdef ENABLE_ASSERTLOG
inline std::string compose_assert_msg(const char *file, int line,
                                      const char *condition_str) {
  std::stringstream ss;
  ss << "Assertion Fail [" << file << ":" << line << "]: " << condition_str;
  return ss.str();
}

#define AssertLog(assert_condition)                                            \
  if (!(assert_condition))                                                     \
  CLOG(ERROR, "general_log")                                                   \
      << std::string("Assertion Fail: ") + (#assert_condition),                \
      throw steps::AssertErr(                                                  \
          compose_assert_msg(__FILE__, __LINE__, (#assert_condition)))
#else
#define AssertLog(assert_condition)                                            \
  {}
#endif

#define ErrLog(msg)                                                            \
  CLOG(ERROR, "general_log") << std::string("GeneralErr: ") + (msg),           \
      throw steps::Err(std::string("GeneralErr: ") + (msg))

#define NotImplErrLog(msg)                                                     \
  CLOG(ERROR, "general_log") << std::string("NotImplErr: ") + (msg),           \
      throw steps::NotImplErr(std::string("NotImplErr: ") + (msg))

#define ArgErrLog(msg)                                                         \
  CLOG(ERROR, "general_log") << std::string("ArgErr: ") + (msg),               \
      throw steps::ArgErr(std::string("ArgErr: ") + (msg))

#define ProgErrLog(msg)                                                        \
  CLOG(ERROR, "general_log") << std::string("ProgErr: ") + (msg),              \
      throw steps::ProgErr(std::string("ProgErr: ") + (msg))

#define SysErrLog(msg)                                                         \
  CLOG(ERROR, "general_log") << std::string("SysErr: ") + (msg),               \
      throw steps::SysErr(std::string("SysErr: ") + (msg))

#define IOErrLog(msg)                                                          \
  CLOG(ERROR, "general_log") << std::string("IOErr: ") + (msg),                \
      throw steps::IOErr(std::string("IOErr: ") + (msg))

inline std::string compose_err_msg(const char *error_type,
                                   const std::string &msg,
                                   const char *condition_str) {
  std::stringstream ss;
  ss << error_type << ": " << msg << "\n[Error Condition] " << condition_str;
  return ss.str();
}

#define ErrLogIf(error_condition, msg)                                         \
  if (error_condition)                                                         \
  CLOG(ERROR, "general_log")                                                   \
      << compose_err_msg("GeneralErr", msg, (#error_condition)),               \
      throw steps::Err(compose_err_msg("GeneralErr", msg, (#error_condition)))

#define NotImplErrLogIf(error_condition, msg)                                  \
  if (error_condition)                                                         \
  CLOG(ERROR, "general_log")                                                   \
      << compose_err_msg("NotImplErr", msg, (#error_condition)),               \
      throw steps::NotImplErr(                                                 \
          compose_err_msg("NotImplErr", msg, (#error_condition)))

#define ArgErrLogIf(error_condition, msg)                                      \
  if (error_condition)                                                         \
  CLOG(ERROR, "general_log")                                                   \
      << compose_err_msg("ArgErr", msg, (#error_condition)),                   \
      throw steps::ArgErr(compose_err_msg("ArgErr", msg, (#error_condition)))

#define ProgErrLogIf(error_condition, msg)                                     \
  if (error_condition)                                                         \
  CLOG(ERROR, "general_log")                                                   \
      << compose_err_msg("ProgErr", msg, (#error_condition)),                  \
      throw steps::ProgErr(                                                    \
          compose_err_msg("ProgErr", msg, (#error_condition)))

#define SysErrLogIf(error_condition, msg)                                      \
  if (error_condition)                                                         \
  CLOG(ERROR, "general_log")                                                   \
      << compose_err_msg("SysErr", msg, (#error_condition)),                   \
      throw steps::SysErr(compose_err_msg("SysErr", msg, (#error_condition)))

#define IOErrLogIf(error_condition, msg)                                       \
  if (error_condition)                                                         \
  CLOG(ERROR, "general_log")                                                   \
      << compose_err_msg("IOErr", msg, (#error_condition)),                    \
      throw steps::IOErr(compose_err_msg("IOErr", msg, (#error_condition)))

namespace steps {

////////////////////////////////////////////////////////////////////////////////

/// Base STEPS exception class. All 'real' exceptions are derived from this.
///
struct Err : public std::exception {
  Err(std::string const &msg = "") : pMessage(msg) {}

  const char *getMsg() const noexcept;
  const char *what() const noexcept override { return getMsg(); }

private:
  std::string pMessage;
};

////////////////////////////////////////////////////////////////////////////////

struct AssertErr : public Err {
  AssertErr(std::string const &msg = "") : Err(msg) {}
};

////////////////////////////////////////////////////////////////////////////////

/// This exception gets thrown when some STEPS interface function is not
/// implemented, usually in a solver module. This can be because the
/// function does not make sense for that particular solver, or because
/// the programmer took a shortcut...
///
struct NotImplErr : public Err {
  NotImplErr(std::string const &msg = "") : Err(msg) {}
};

////////////////////////////////////////////////////////////////////////////////

/// This exception usually gets thrown when some caller-provided argument
/// doesn't make sense immediately. This 'caller' should be restricted
/// to the actual user as much as possible.
///
struct ArgErr : public Err {
  ArgErr(std::string const &msg = "") : Err(msg) {}
};

////////////////////////////////////////////////////////////////////////////////

/// This exception gets thrown whenever the program knows it has failed,
/// in other words not because of some mistake by the user, but e.g.
/// because of some hard coded limit (running out of index numbers,...)
///
/// In principle, it should be only be thrown in absurdly rare cases.
///
struct ProgErr : public Err {
  ProgErr(std::string const &msg = "") : Err(msg) {}
};

////////////////////////////////////////////////////////////////////////////////

/// Generic base class for any system-induced error. These exceptions signal
/// truly unexpected situations, such as out-of-memory or I/O problems.
///
struct SysErr : public Err {
  SysErr(std::string const &msg = "") : Err(msg) {}
};

////////////////////////////////////////////////////////////////////////////////

/// Gets thrown whenever there is a system failure involving I/O, such
/// as dealing with files.
///
struct IOErr : public SysErr {
  IOErr(std::string const &msg = "") : SysErr(msg) {}
};

////////////////////////////////////////////////////////////////////////////////
} // namespace steps

#endif
// STEPS_ERROR_HPP

// END
