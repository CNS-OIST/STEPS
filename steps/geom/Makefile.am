## ./steps/geom/Makefile.am
##
## Process this file with automake to produce Makefile.in

## $Id$

## Set $MAKE if necessary.
@SET_MAKE@

## VPATH support.
srcdir = @srcdir@
VPATH = @srcdir@

## Python modules in local packages.
python_steps_geomdir = $(pythondir)/steps/geom
python_steps_geom_PYTHON = \
    __init__.py \
    core.py 
    
## Files that should be cleaned when building a distribution.
distclean-local:
	-rm -rf .DS_Store *~ *.bak
    
## END
