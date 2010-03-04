dnl ----------------------------------------------------------------------------
dnl Synopsis:
dnl     AC_CXX_NAMESPACES
dnl ----------------------------------------------------------------------------
dnl Description:
dnl     Checks for the availability of namespaces in C++.
dnl ----------------------------------------------------------------------------
dnl Authors: 
dnl     Todd Veldhuizen
dnl     Luc Maisonobe <luc@spaceroots.org>
dnl ----------------------------------------------------------------------------
dnl Source: 
dnl     http://autoconf-archive.cryp.to/ac_cxx_namespaces.html
dnl ----------------------------------------------------------------------------

AC_DEFUN([AC_CXX_NAMESPACES],[

# Actual check.
AC_CACHE_CHECK(
    [whether the compiler implements namespaces],
    [ac_cv_cxx_namespaces],
    [
        AC_LANG_PUSH(C++)
        AC_TRY_COMPILE(
            [
namespace Outer { namespace Inner { int i = 0; }}
            ],[
using namespace Outer::Inner; return i;
            ],
            [ac_cv_cxx_namespaces=yes], 
            [ac_cv_cxx_namespaces=no]
        )
        AC_LANG_POP
    ]
)

# Substitution stuff.
if test x"$ac_cv_cxx_namespaces" = xno; then
    AC_MSG_ERROR([requires namespace support])
fi
])
