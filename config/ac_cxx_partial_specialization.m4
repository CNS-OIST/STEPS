dnl ----------------------------------------------------------------------------
dnl Synopsis:
dnl     AC_CXX_PARTIAL_SPECIALIZATION
dnl ----------------------------------------------------------------------------
dnl Description:
dnl     If the compiler supports partial specialization, define 
dnl     HAVE_PARTIAL_SPECIALIZATION.
dnl ----------------------------------------------------------------------------
dnl Authors: 
dnl     Todd Veldhuizen
dnl     Luc Maisonobe <luc@spaceroots.org>
dnl ----------------------------------------------------------------------------
dnl Source: 
dnl     http://autoconf-archive.cryp.to/ac_cxx_partial_specialization.html
dnl ----------------------------------------------------------------------------

AC_DEFUN([AC_CXX_PARTIAL_SPECIALIZATION],[

# Actual check.
AC_CACHE_CHECK(
    [whether the compiler supports partial specialization],
    [ac_cv_cxx_partial_specialization],
    [
        AC_LANG_PUSH(C++)
        AC_TRY_COMPILE(
            [
template<class T, int N> class A            { public : enum e { z = 0 }; };
template<int N>          class A<double, N> { public : enum e { z = 1 }; };
template<class T>        class A<T, 2>      { public : enum e { z = 2 }; };
            ],[
return (A<int,3>::z == 0) && (A<double,3>::z == 1) && (A<float,2>::z == 2);
            ],
            [ac_cv_cxx_partial_specialization=yes], 
            [ac_cv_cxx_partial_specialization=no]
        )
        AC_LANG_POP
    ]
)

# Substitution stuff.
if test x"$ac_cv_cxx_partial_specialization" = xyes; then
    AC_DEFINE(HAVE_PARTIAL_SPECIALIZATION,,
        [define if the compiler supports partial specialization])
else
    $1
    :
fi
])