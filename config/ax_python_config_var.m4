dnl ----------------------------------------------------------------------------
dnl Synopsis:
dnl     AX_PYTHON_CONFIG_VAR(PYTHON_VARIABLE, [SHELL_VARIABLE])
dnl     AX_PYTHON_CONFIG_H
dnl     AX_PYTHON_MAKEFILE
dnl ----------------------------------------------------------------------------
dnl Description:
dnl
dnl     AX_PYTHON_CONFIG_VAR:
dnl     Using the Python module distutils.sysconfig[1], return a Python 
dnl     configuration variable. PYTHON_VARIABLE is the name of the variable 
dnl     to request from Python, and SHELL_VARIABLE is the name of the shell 
dnl     variable into which the results should be deposited. If SHELL_VARIABLE 
dnl     is not specified, the macro wil prefix PY_ to the PYTHON_VARIABLE, 
dnl     e.g., LIBS -> PY_LIBS.
dnl
dnl     SHELL_VARIABLE is AC_SUBST'd. No action is taken if an error occurs. 
dnl     Note if $PYTHON is not set, AC_CHECK_PROG(PYTHON, python, python) will 
dnl     be run.
dnl
dnl     Example:
dnl         AX_PYTHON_CONFIG_VAR(LINKFORSHARED, PY_LFS)
dnl
dnl     AX_PYTHON_CONFIG_H:
dnl     Using the Python module distutils.sysconfig[1], put the full pathname 
dnl     of the config.h file used to compile Python into the shell variable 
dnl     PY_CONFIG_H. PY_CONFIG_H is AC_SUBST'd. Note if $PYTHON is not set, 
dnl     AC_CHECK_PROG(PYTHON, python, python) will be run.
dnl
dnl     AX_PYTHON_MAKEFILE:
dnl     Using the Python module distutils.sysconfig[1], put the full pathname 
dnl     of the Makefile file used to compile Python into the shell variable 
dnl     PY_MAKEFILE. PY_MAKEFILE is AC_SUBST'd. Note if $PYTHON is not set, 
dnl     AC_CHECK_PROG(PYTHON, python, python) will be run.
dnl
dnl [1] http://www.python.org/doc/current/dist/module-distutils.sysconfig.html
dnl ----------------------------------------------------------------------------
dnl Author:
dnl     Dustin Mitchell <dustin@cs.uchicago.edu>
dnl ----------------------------------------------------------------------------
dnl Sources:
dnl     http://autoconf-archive.cryp.to/ax_python_config_var.html
dnl ----------------------------------------------------------------------------

AC_DEFUN([AX_PYTHON_CONFIG_VAR],
[
 AC_MSG_CHECKING(for Python config variable $1)
 if test -z "$PYTHON"
 then
   AC_CHECK_PROG(PYTHON,python,python)
 fi
 py_error="no"
 pyval=`$PYTHON -c "from distutils import sysconfig;dnl
print sysconfig.get_config_var('$1')"` || py_error="yes"
 if test "$py_error" = "yes"
 then
   AC_MSG_RESULT(no - an error occurred)
 else
   AC_MSG_RESULT($pyval)
   m4_ifval([$2],[$2],[PY_$1])="$pyval"
   AC_SUBST(m4_ifval([$2],[$2],[PY_$1]))
 fi
])

AC_DEFUN([AX_PYTHON_CONFIG_H],
[
 AC_MSG_CHECKING(location of Python's config.h)
 if test -z "$PYTHON"
 then
   AC_CHECK_PROG(PYTHON,python,python)
 fi
 py_error="no"
 PY_CONFIG_H=`$PYTHON -c "from distutils import sysconfig;dnl
print sysconfig.get_config_h_filename()"` || py_error = "yes"
 if test "$py_error" = "yes"
 then
   AC_MSG_RESULT(no - an error occurred)
 else
   AC_MSG_RESULT($PY_CONFIG_H)
   AC_SUBST(PY_CONFIG_H)
 fi
])

AC_DEFUN([AX_PYTHON_MAKEFILE],
[
 AC_MSG_CHECKING(location of Python's Makefile)
 if test -z "$PYTHON"
 then
   AC_CHECK_PROG(PYTHON,python,python)
 fi
 py_error="no"
 PY_MAKEFILE=`$PYTHON -c "from distutils import sysconfig;dnl
print sysconfig.get_makefile_filename()"` || py_error = "yes"
 if test "$py_error" = "yes"
 then
   AC_MSG_RESULT(no - an error occurred)
 else
   AC_MSG_RESULT($PY_MAKEFILE)
   AC_SUBST(PY_MAKEFILE)
 fi
])
