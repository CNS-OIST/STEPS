
dnl $Id$

dnl ----------------------------------------------------------------------------
dnl     BOILERPLATE
dnl
dnl     Standard boilerplate code.
dnl ----------------------------------------------------------------------------

dnl Initialize the autotools system.
AC_INIT([steps], [0.5.0], [wils@oist.jp])

AC_PREREQ(2.57)

dnl Copyright notice.
AC_COPYRIGHT([[STEPS -- STochastic Engine for Pathway Simulation
Copyright (C) 2005-2008 by Stefan Wils. All rights reserved.
See the LICENSE file in the root directory for more information.
]])

dnl Specify the location for the configuration header.
AM_CONFIG_HEADER(cpp/steps/config.h)


dnl Surround the config.h file with a header include guard.
AH_TOP([
#ifndef STEPS_CONFIG_H
#define STEPS_CONFIG_H 1
])
dnl Header include guard (cont'd) + type definition shortcuts.
AH_BOTTOM([
/* Find int16 and uint16. */
#if SIZEOF_SHORT_INT == 2
typedef short int int16;
typedef unsigned short int uint16;
#elif SIZEOF_INT == 2
typedef int int16;
typedef unsigned int uint16;
#elif SIZEOF_LONG_INT == 2
typedef long int int16;
typedef unsigned long int uint16;
#else
/* 16 BIT INTEGER NOT FOUND! */
#endif

/* Find int32 and uint32. */
#if SIZEOF_INT == 4
typedef int int32;
typedef unsigned int uint32;
#elif SIZEOF_SHORT_INT == 4
typedef short int int32;
typedef unsigned short int uint32;
#elif SIZEOF_LONG_INT == 4
typedef long int int32;
typedef unsigned long int uint32;
#else
/* 32 BIT INTEGER NOT FOUND! */
#endif

/* Find int64 and uint64. */
#if SIZEOF_INT == 8
typedef int int64;
typedef unsigned int uint64;
#elif SIZEOF_LONG_INT == 8
typedef long int int64;
typedef unsigned long int uint64;
#elif SIZEOF_SHORT_INT == 8
typedef short int int64;
typedef unsigned short int uint64;
#else
/* 64 BIT INTEGER NOT FOUND! */
#endif

#endif
])

dnl Initialize automake.
AM_INIT_AUTOMAKE([AC_PACKAGE_NAME], [AC_PACKAGE_VERSION])

dnl Add maintainer mode option.
dnl AM_MAINTAINER_MODE

dnl Uncomment to not build static (or shared) libraries by default.
dnl AC_DISABLE_SHARED
AC_DISABLE_STATIC

dnl ----------------------------------------------------------------------------
dnl Dynamic loading
dnl ----------------------------------------------------------------------------

dnl Build dlls on Windows. For this to work all exported library symbols  
dnl must be properly marked with __declspec(dllexport). 
AC_LIBTOOL_WIN32_DLL
dnl Use libtool to build the system (we use AC_PROG_LIBTOOL instead of 
dnl AM_PROG_LIBTOOL, since we're using libtool/automake >= 1.4).
AC_PROG_LIBTOOL

dnl Perform a build-time sanity check.
AM_SANITY_CHECK

dnl ----------------------------------------------------------------------------
dnl     OPTIONS
dnl
dnl     Macros which add command-line options to configure, such as
dnl     AC_ARG_ENABLE. It is typical to put support code for the option in this
dnl     section as well, if it is short enough.
dnl ----------------------------------------------------------------------------

dnl ----------------------------------------------------------------------------
dnl     PROGRAMS
dnl
dnl     Checks for programs that are either needed by the configure process,
dnl     build process, or by one of the programs being built. This usually
dnl     involves calls to macros like AC_CHECK_PROG and AC_PATH_TOOL.
dnl ----------------------------------------------------------------------------

dnl ----------------------------------------------------------------------------
dnl Compiler: C
dnl ----------------------------------------------------------------------------

dnl Check for a C compiler.
AC_PROG_CC

dnl ----------------------------------------------------------------------------
dnl Compiler: C++
dnl ----------------------------------------------------------------------------

dnl Check for a C++ compiler.
AC_PROG_CXX

dnl Does it support namespaces? Gives an error if not found.
dnl AC_CXX_NAMESPACES

dnl Check for a normal C/C++ preprocessor.
AC_PROG_CPP

dnl Check for C and C++ warning options.
dnl AC_COMPILE_WARNINGS

dnl ----------------------------------------------------------------------------
dnl Python
dnl ----------------------------------------------------------------------------

AM_PATH_PYTHON(2.5)
AC_PYTHON_DEVEL(2.5)

dnl ----------------------------------------------------------------------------
dnl SWIG
dnl ----------------------------------------------------------------------------

AC_PROG_SWIG(1.3.20)
SWIG_ENABLE_CXX
dnl SWIG_MULTI_MODULE_SUPPORT
SWIG_PYTHON

dnl ----------------------------------------------------------------------------
dnl Other programs.
dnl ----------------------------------------------------------------------------

dnl Check for common build/install tools.
AC_PROG_INSTALL

dnl If `ln -s' works on the current file system (the operating system and
dnl file system support symbolic links), set the output variable LN_S to 
dnl `ln -s'; otherwise, if `ln' works, set LN_S to `ln', and otherwise set
dnl it to `cp -p'.
AC_PROG_LN_S

dnl If grep -E is available, 
AC_PROG_EGREP

dnl Define an output variable SET_MAKE that is empty if the version of make 
dnl installed on the host automatically defines the variable $MAKE. If this
dnl make does not define $MAKE, SET_MAKE equals "MAKE=make". Hence, @SET_MAKE@ 
dnl should be added to each Makefile.am that calls make for subdirectories
dnl (or in every Makefile.am).
AC_PROG_MAKE_SET

dnl ----------------------------------------------------------------------------
dnl     LIBRARIES
dnl
dnl     Checks for libraries come before checks for other objects visible to
dnl     C/C++. This is necessary because some other checks work by trying to
dnl     link or run a program; by checking for libraries first you ensure that
dnl     the resulting programs can be linked.
dnl ----------------------------------------------------------------------------

dnl ----------------------------------------------------------------------------
dnl NumPy stuff.
dnl ----------------------------------------------------------------------------

AC_PYTHON_MODULE([numpy])
AX_NUMPY_INCLUDE

dnl ----------------------------------------------------------------------------
dnl     HEADERS
dnl
dnl     Next come checks for existence of headers.
dnl ----------------------------------------------------------------------------

dnl Check for basic standard c library.
AC_HEADER_STDC

dnl ----------------------------------------------------------------------------
dnl     TYPEDEFS AND STRUCTURES
dnl
dnl     Checks for typedefs come after checking for headers for the simple
dnl     reason that typedefs appear in headers, and we need to know which 
dnl     headers we can use before we look inside them.
dnl ----------------------------------------------------------------------------


dnl ----------------------------------------------------------------------------
dnl     FUNCTIONS
dnl
dnl     Finally we check for functions. These come last because functions have
dnl     dependencies on the preceding items: when searching for functions, 
dnl     libraries are needed in order to correctly link, headers are needed
dnl     in order to find prototypes (this is especially important for C++,
dnl     which has stricter prototyping rules than C), and typedefs are needed
dnl     for those functions which use or return types which are not built in.
dnl ----------------------------------------------------------------------------


dnl ----------------------------------------------------------------------------
dnl     OTHER SYSTEM/COMPILER CHARACTERISTICS
dnl
dnl     Check other system/compiler characteristics.
dnl ----------------------------------------------------------------------------

dnl See which path style we're dealing with here.
dnl AC_ARG_WITH_PATH_STYLE

dnl Check endianess even when crosscompiling (>2.53).
AC_C_BIGENDIAN(
[AC_DEFINE([IS_BIG_ENDIAN],, [Big endian system])], 
[AC_DEFINE([IS_LITTLE_ENDIAN],, [Little endian system])])

dnl Pre-check size of some integer types.
AC_CHECK_SIZEOF([char])
AC_CHECK_SIZEOF([short int])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long int])

dnl Pre-check size of generic pointer.
AC_CHECK_SIZEOF([void*])

dnl Pre-check size of some real types.
AC_CHECK_SIZEOF([float])
AC_CHECK_SIZEOF([double])
AC_C_LONG_DOUBLE

dnl ----------------------------------------------------------------------------
dnl     OUTPUT
dnl
dnl     Generate the autotools system for this package.
dnl ----------------------------------------------------------------------------

AC_OUTPUT([
cpp/steps/compat/Makefile
cpp/steps/console/Makefile
cpp/steps/geom/Makefile
cpp/steps/math/Makefile
cpp/steps/model/Makefile
cpp/steps/rng/Makefile
cpp/steps/solver/Makefile
cpp/steps/tetexact/Makefile
cpp/steps/tools/Makefile
cpp/steps/wmdirect/Makefile
cpp/steps/wmrk4/Makefile
cpp/steps/Makefile
cpp/Makefile
steps/Makefile
Makefile
])
dnl PUT THIS IN THE ABOVE LIST cpp/steps/tetmesh/Makefile

