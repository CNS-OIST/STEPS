dnl ----------------------------------------------------------------------------
dnl Synopsis:
dnl     SWIG_ENABLE_CXX
dnl ----------------------------------------------------------------------------
dnl Description:
dnl
dnl     Enable SWIG C++ support. This affects all invocations of $(SWIG).
dnl
dnl ----------------------------------------------------------------------------
dnl Authors:
dnl     Sebastian Huber <sebastian-huber@web.de>
dnl     Alan W. Irwin <irwin@beluga.phys.uvic.ca>
dnl     Rafael Laboissiere <rafael@laboissiere.net>
dnl     Andrew Collier <colliera@ukzn.ac.za>
dnl ----------------------------------------------------------------------------
dnl Source:
dnl     http://autoconf-archive.cryp.to/swig_enable_cxx.html
dnl ----------------------------------------------------------------------------

AC_DEFUN([SWIG_ENABLE_CXX],[
        AC_REQUIRE([AC_PROG_SWIG])
        AC_REQUIRE([AC_PROG_CXX])
        SWIG="$SWIG -c++"
])
