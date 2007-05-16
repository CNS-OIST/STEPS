## ./steps/examples/Makefile.am
##
## Process this file with automake to produce Makefile.in

## $Id$

## Set $MAKE if necessary.
@SET_MAKE@

## VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

## Python package.
steps_examplesdir = $(pythondir)/steps/examples
steps_examples_PYTHON = \
    __init__.py \
    lotka_volterra.py \
    rev_isomer.py
    
## Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak
    
## END
