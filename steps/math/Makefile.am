## ./steps/math/Makefile.am
##
## Process this file with automake to produce Makefile.in

## $Id$

## Set $MAKE if necessary.
@SET_MAKE@

## VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

## Headers that should be included in the distribution and installed.
pkginclude_HEADERS = \
    aabbox3.hpp \
    accumulate.hpp \
    color.hpp \
    constants.hpp \
    geometry.hpp \
    iota.hpp \
    linalg.hpp \
    matrix4.hpp \
    normal3.hpp \
    point3.hpp \
    tools.hpp \
    vector3.hpp

## Compiler flags.
AM_CXXFLAGS = \
	-I$(top_srcdir)

## Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak
    
## END
