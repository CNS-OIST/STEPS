## ./steps/tests/Makefile.am
##
## Process this file with automake to produce Makefile.in

## $Id$

## Set $MAKE if necessary.
@SET_MAKE@

## VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

## Python modules in local packages.
steps_testsdir = $(pythondir)/steps/tests
steps_tests_PYTHON = \
    __init__.py \
    test_misc.py \
    test_model.py \
    tools.py

## Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak
    
## END
