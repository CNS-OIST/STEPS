dnl ----------------------------------------------------------------------------
dnl Synopsis:
dnl     SWIG_PYTHON([use-shadow-classes = {no, yes}])
dnl ----------------------------------------------------------------------------
dnl Description:
dnl
dnl     Checks for Python and provides the $(SWIG_PYTHON_CPPFLAGS), and 
dnl     $(SWIG_PYTHON_OPT) output variables.
dnl
dnl     $(SWIG_PYTHON_OPT) contains all necessary SWIG options to generate 
dnl     code for Python. Shadow classes are enabled unless the value of the
dnl     optional first argument is exactly 'no'. If you need multi module 
dnl     support (provided by the SWIG_MULTI_MODULE_SUPPORT macro) use 
dnl     $(SWIG_PYTHON_LIBS) to link against the appropriate library. It 
dnl     contains the SWIG Python runtime library that is needed by the type 
dnl     check system for example.
dnl ----------------------------------------------------------------------------
dnl Authors:
dnl     Sebastian Huber <sebastian-huber@web.de>
dnl     Alan W. Irwin <irwin@beluga.phys.uvic.ca>
dnl     Rafael Laboissiere <rafael@laboissiere.net>
dnl     Andrew Collier <colliera@ukzn.ac.za>
dnl ----------------------------------------------------------------------------
dnl Source:
dnl     http://autoconf-archive.cryp.to/swig_python.html
dnl ----------------------------------------------------------------------------

AC_DEFUN([SWIG_PYTHON],[
        AC_REQUIRE([AC_PROG_SWIG])
        AC_REQUIRE([AC_PYTHON_DEVEL])
        test "x$1" != "xno" || swig_shadow=" -noproxy"
        AC_SUBST([SWIG_PYTHON_OPT],[-python$swig_shadow])
        AC_SUBST([SWIG_PYTHON_CPPFLAGS],[$PYTHON_CPPFLAGS])
])
