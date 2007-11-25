dnl ----------------------------------------------------------------------------
dnl Synopsis:
dnl     AX_NUMPY_INCLUDE
dnl ----------------------------------------------------------------------------
dnl Description:
dnl
dnl     Attempts to find the location of the NumPy headers. Requires 
dnl     AC_PYTHON_DEVEL. Defines NUMPY_CPPFLAGS.
dnl ----------------------------------------------------------------------------
dnl Author:
dnl     Stefan Wils <wils@oist.jp>
dnl ----------------------------------------------------------------------------

AC_DEFUN([AX_NUMPY_INCLUDE],
[
    AC_MSG_CHECKING(for NumPy C API include path)
    if test -z "$PYTHON_SITE_PKG"
    then
        AC_MSG_RESULT(no)
    fi
    numpy_headers_loc="$PYTHON_SITE_PKG/numpy/core/include"
    numpy_headers_arrayobject="$numpy_headers_loc/numpy/arrayobject.h"
    if test -r "$numpy_headers_arrayobject"
    then
        NUMPY_CPPFLAGS="-I$numpy_headers_loc"
        AC_SUBST(NUMPY_CPPFLAGS)
        AC_MSG_RESULT($NUMPY_CPPFLAGS)
    else
        AC_MSG_ERROR([Cannot find NumPy headers. Please specify them.])
    fi
])