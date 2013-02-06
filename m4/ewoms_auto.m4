dnl -*- autoconf -*-
AC_DEFUN([EWOMS_CHECK_AUTO],[
  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([GXX0X])
  AC_LANG_PUSH([C++])
  AC_MSG_CHECKING([whether the auto keyword is supported])
  AC_COMPILE_IFELSE([
        void f()
        {
            struct Foo
            {
                static float bar() { return 1.234; }
            };
            auto foobar = Foo::bar();
         }], [
    HAVE_AUTO=yes
    AC_MSG_RESULT(yes)], [
    HAVE_AUTO=no
    AC_MSG_RESULT(no)])
  if test "$HAVE_AUTO" = "yes"; then
    AC_DEFINE(HAVE_AUTO, 1, [Define to 1 if the auto keyword is supported])
  fi
  AC_LANG_POP
])
