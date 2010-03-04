dnl ----------------------------------------------------------------------------
dnl Synopsis:
dnl     AC_PYTHON_MODULE(modname[, fatal])
dnl ----------------------------------------------------------------------------
dnl Description:
dnl     Checks for Python module. If fatal is non-empty then absence of a 
dnl     module will trigger an error.
dnl ----------------------------------------------------------------------------
dnl Authors: 
dnl     Andrew Collier <colliera@ukzn.ac.za>
dnl ----------------------------------------------------------------------------
dnl Source: 
dnl     http://autoconf-archive.cryp.to/ac_python_module.html
dnl ----------------------------------------------------------------------------

AC_DEFUN([AC_PYTHON_MODULE],[
    if test -z $PYTHON;
    then
        PYTHON="python"
    fi
    PYTHON_NAME=`basename $PYTHON`
    AC_MSG_CHECKING($PYTHON_NAME module: $1)
        $PYTHON -c "import $1" 2>/dev/null
        if test $? -eq 0;
        then
            AC_MSG_RESULT(yes)
            eval AS_TR_CPP(HAVE_PYMOD_$1)=yes
        else
            AC_MSG_RESULT(no)
            eval AS_TR_CPP(HAVE_PYMOD_$1)=no
            #
            if test -n "$2"
            then
                AC_MSG_ERROR(failed to find required module $1)
                exit 1
            fi
        fi
])
