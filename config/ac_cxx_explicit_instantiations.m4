dnl ----------------------------------------------------------------------------
dnl Synopsis:
dnl     AC_CXX_EXPLICIT_INSTANTIATIONS
dnl ----------------------------------------------------------------------------
dnl Description:
dnl     If the C++ compiler supports explicit instanciations syntax, 
dnl     define HAVE_EXPLICIT_INSTANTIATIONS.
dnl ----------------------------------------------------------------------------
dnl Author: 
dnl     Todd Veldhuizen
dnl     Luc Maisonobe <luc@spaceroots.org>
dnl     (Modified for STEPS)
dnl ----------------------------------------------------------------------------
dnl Source:
dnl     http://autoconf-archive.cryp.to/ac_cxx_explicit_instantiations.html
dnl ----------------------------------------------------------------------------

AC_DEFUN([AC_CXX_EXPLICIT_INSTANTIATIONS],[

# Actual check.
AC_CACHE_CHECK(
    [whether the compiler supports explicit instantiations],
    [ac_cv_cxx_explinst],
    [
        AC_LANG_PUSH(C++)
        AC_TRY_COMPILE(
            [
template <class T> class A { T t; }; template class A<int>;
            ],
            [], 
            [ac_cv_cxx_explinst=yes], 
            [ac_cv_cxx_explinst=no]
        )
        AC_LANG_POP
    ]
)

# Substitution stuff.
if test x"$ac_cv_cxx_explinst" = xyes; then
    AC_DEFINE(HAVE_EXPLICIT_INSTANTIATIONS,,
        [define if the compiler supports explicit instantiations])
else
    $1
    :
fi
])