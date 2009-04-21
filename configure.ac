AC_PREREQ(2.50)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([test/diffusion/test_diffusion.cc])
AM_CONFIG_HEADER([config.h])
AC_SUBST([LOCAL_LIBS], '$(top_srcdir)/dumux/libdumux.la')

AC_CHECK_PROGS([TEX], [latex], [true])
AC_CHECK_PROGS([BIBTEX], [bibtex], [true])
AC_CHECK_PROGS([DVIPDF], [dvipdf], [true])
AC_CHECK_PROGS([DVIPS], [dvips], [true])
AC_CHECK_PROGS([WML], [wml], [true])
AM_CONDITIONAL([WML], [test "x$WML" != xtrue])
AC_CHECK_PROGS([TEX4HT], [tex4ht], [true])
AC_CHECK_PROGS([MK4HT], [mk4ht], [true])
AC_CHECK_PROGS([T4HT], [t4ht], [true])
AM_CONDITIONAL([TEX4HT], [test "x$TEX4HT" != xtrue])
AC_CHECK_PROGS([CONVERT], [convert], [false])
AM_CONDITIONAL([CONVERT], [test "x$CONVERT" != xfalse])
DUNE_CHECK_ALL

AC_DEFUN([SET_PARDISO], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_pardiso_ok=no

AC_ARG_WITH(pardiso,
	[AC_HELP_STRING([--with-pardiso=<lib>], [use PARDISO library <lib>])])
case $with_pardiso in
	yes | "") ;;
	no) acx_pardiso_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) PARDISO_LIBS="$with_pardiso" ;;
	*) PARDISO_LIBS="-l$with_pardiso" ;;
esac

AC_F77_FUNC(pardisoinit)

acx_pardiso_save_LIBS="$LIBS"
LIBS="$BLAS_LIBS $LIBS $FLIBS"

if test $acx_pardiso_ok = no; then
if test "x$PARDISO_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$PARDISO_LIBS $LIBS"
	AC_MSG_CHECKING([for $pardisoinit in $PARDISO_LIBS])
	AC_TRY_LINK_FUNC($pardisoinit, [acx_pardiso_ok=yes], [PARDISO_LIBS=""])
	AC_MSG_RESULT($acx_pardiso_ok)
	LIBS="$save_LIBS"
fi
fi

AC_SUBST(PARDISO_LIBS)

LIBS="$acx_pardiso_save_LIBS"

if test x"$acx_pardiso_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_PARDISO,1,[Define if you have a PARDISO library.]),[$1])
        :
else
        acx_pardiso_ok=no
        $2
fi
])dnl SET_PARDISO

SET_PARDISO

AC_SUBST(AM_CPPFLAGS, $ALL_PKG_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $ALL_PKG_LDFLAGS)
LIBS="$DUNE_LIBS $ALL_PKG_LIBS $PARDISO_LIBS"

AC_CONFIG_FILES([Makefile
    doc/Makefile 
    dumux/Makefile 
    test/Makefile
    doc/doxygen/Makefile 
    doc/handbook/Makefile
    dumux/fractionalflow/Makefile 
    dumux/io/Makefile 
    dumux/nonlinear/Makefile 
    dumux/twophase/Makefile
    test/2p2c/Makefile
    test/2p2cni/Makefile
    test/fractionalflow/Makefile
    test/transport/Makefile
    test/tutorial/Makefile
    test/twophase/Makefile
    dumux/twophase/fv/Makefile 
    dumux/twophase/problems/Makefile 
    dumux/fractionalflow/impes/Makefile 
    dune-mux.pc])


AC_OUTPUT
DUNE_SUMMARY_ALL
echo "Pardiso..........: $acx_pardiso_ok"
echo 
