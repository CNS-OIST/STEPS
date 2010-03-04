dnl ----------------------------------------------------------------------------
dnl Synopsis:
dnl     AX_CXX___ATTRIBUTE__
dnl ----------------------------------------------------------------------------
dnl Description:
dnl     Provides a test for the compiler support of __attribute__ extensions. 
dnl     defines HAVE___ATTRIBUTE__ if it is found. Adapted for C++.
dnl
dnl     Originating from the 'pork' package by Ryan McCabe <ryan@numb.org>
dnl ----------------------------------------------------------------------------
dnl Authors: 
dnl     Christian Haggstrom <chm@c00.info>
dnl ----------------------------------------------------------------------------
dnl Source: 
dnl     http://autoconf-archive.cryp.to/ax_c___attribute__.html
dnl ----------------------------------------------------------------------------

AC_DEFUN([AX_CXX___ATTRIBUTE__],[

# Actual check.
AC_MSG_CHECKING([for __attribute__ in C++])
AC_CACHE_VAL(
    [ac_cv_cxx___attribute__], 
    [
        AC_LANG_PUSH(C++)
        AC_TRY_COMPILE(
            [],[
int __attribute__ ((unused)) a;
            ],
            [ac_cv_cxx___attribute__=yes],
            [ac_cv_cxx___attribute__=no]
        )
        AC_LANG_POP
    ]
)

# Substitution stuff.
if test x"$ac_cv_cxx___attribute__" = xyes; then
    AC_DEFINE(HAVE_CXX___ATTRIBUTE__, 1, 
        [define if your C++ compiler has __attribute__])
fi
AC_MSG_RESULT($ac_cv_cxx___attribute__)

])
