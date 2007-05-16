## ./steps/Makefile.am
##
## Process this file with automake to produce Makefile.in

## $Id$

## Set $MAKE if necessary.
@SET_MAKE@

## VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

## Process subdirectories.
SUBDIRS	= \
    examples \
    geom \
    math \
    model \
    rng \
    sim \
    tests \
    .

## Python package.
stepsdir = $(pythondir)/steps
steps_PYTHON = \
    __init__.py \
    error.py \
    tools.py

## Headers that should be included in the distribution and installed.
pkginclude_HEADERS = \
    common.h

## Compiler flags.
AM_CXXFLAGS = \
	-I$(top_srcdir)

## Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak

## END
