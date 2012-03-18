# Check whether the variable attibute 'unused' is available and sets
# the EWOMS_UNUSED macro accordingly
AC_DEFUN([EWOMS_CHECK_UNUSED],[
  AC_MSG_CHECKING([for __attribute__((unused))])
  AC_LANG_PUSH([C++])
  AC_TRY_COMPILE([void fn(void) { int __attribute__((unused)) i; }],[],
                 [EWOMS_UNUSED="__attribute__((unused))"
                  AC_MSG_RESULT(yes)],
                 [EWOMS_UNUSED=""
                  AC_MSG_RESULT(no)])
  AC_LANG_POP([C++])
  AC_DEFINE_UNQUOTED(EWOMS_UNUSED, $DUMUX_UNUSED,
                     [USE WITH CARE: Prevents the compiler from prining a warning if it thinks that a variable might be unused])
])
