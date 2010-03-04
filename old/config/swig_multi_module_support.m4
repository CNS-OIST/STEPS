dnl ----------------------------------------------------------------------------
dnl Synopsis:
dnl     SWIG_MULTI_MODULE_SUPPORT
dnl ----------------------------------------------------------------------------
dnl Description:
dnl     Enable support for multiple modules. This effects all invocations of 
dnl     $(SWIG). You have to link all generated modules against the appropriate 
dnl     SWIG runtime library. If you want to build Python modules for example, 
dnl     use the SWIG_PYTHON macro and link the modules against 
dnl     $(SWIG_PYTHON_LIBS).
dnl ----------------------------------------------------------------------------
dnl Authors:
dnl     Sebastian Huber <sebastian-huber@web.de>
dnl     Alan W. Irwin <irwin@beluga.phys.uvic.ca>
dnl     Rafael Laboissiere <rafael@laboissiere.net>
dnl     Andrew Collier <colliera@ukzn.ac.za>.
dnl ----------------------------------------------------------------------------
dnl Source:
dnl     http://autoconf-archive.cryp.to/swig_multi_module_support.html
dnl ----------------------------------------------------------------------------

AC_DEFUN([SWIG_MULTI_MODULE_SUPPORT],[
        AC_REQUIRE([AC_PROG_SWIG])
        SWIG="$SWIG -noruntime"
])
