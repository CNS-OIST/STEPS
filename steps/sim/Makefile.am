# ./steps/sim/Makefile.am
#
# Process this file with automake to produce Makefile.in

# $Id$

@SET_MAKE@

srcdir = @srcdir@
VPATH = @srcdir@

distclean-local:
	-rm -rf .DS_Store *~ *.bak

pyexec_steps_simdir = $(pythondir)/steps/sim
python_steps_simdir = $(pythondir)/steps/sim
header_steps_simdir = $(pkgincludedir)/steps/sim

header_steps_sim_HEADERS = \
    doxygen.hpp

AM_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) $(PYTHON_CPPFLAGS) -I$(top_srcdir)

SWIG_PYTHON_OPT = @SWIG_PYTHON_OPT@ -I$(top_srcdir) 
%_wrap.cpp: %.i
	$(SWIG) $(SWIG_PYTHON_OPT) -o $@ $<
SWIG_SOURCES = wmdirect_core.i

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# EXTENSION: WMDIRECT_CORE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

pyexec_steps_sim_LTLIBRARIES = _wmdirect_core.la
_wmdirect_core_la_LDFLAGS = -module -avoid-version $(PYTHON_LDFLAGS)

_wmdirect_core_la_SOURCES = \
    $(SWIG_SOURCES) \
    wmdirect_core_wrap.cpp

_wmdirect_core_la_LIBADD = \
    ./wmdirect/libsteps_wmdirect.la 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PYTHON PACKAGES
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

python_steps_sim_PYTHON = \
	__init__.py \
	controller.py \
	wmdirect.py \
	wmdirect_core.py
	
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

SUBDIRS = \
    shared \
    swiginf \
	wmdirect \
	.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# END
