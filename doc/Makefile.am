## ./doc/Makefile.am
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
    devnotes.doc

## Installation directives.
docdir = $(pkgdatadir)/doc
nobase_doc_DATA = \
    devnotes.doc

## Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak

## END
