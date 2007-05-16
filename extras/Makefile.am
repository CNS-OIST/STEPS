## ./extras/Makefile.am
##
## Process this file with automake to produce Makefile.in

## $Id$

## Set $MAKE if necessary.
@SET_MAKE@

## VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

## Things that should be added to the distribution.
EXTRA_DIST = \
    metrics.pl \
    untab.sh

## Installation directives.
extrasdir = $(pkgdatadir)/extras
nobase_extras_DATA = \
    metrics.pl \
    untab.sh

## Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak

## END
