dnl Autoconf for parsem
AC_PREREQ(2.61)
AC_INIT(ParseM, 0.1, mike@mikeimelfort.com)

dnl Keep helpers in build-aux/ subdirectory to not leave too much junk.
#AC_CONFIG_AUX_DIR([build-aux])

# test for a c compiler
AC_PROG_CC
AC_PROG_CXX

# test for libm
m_flags=
AC_CHECK_LIB([m],[lgamma,log,exp],
             [m_flags="-lm"
              AC_DEFINE([HAVE_M],[1],[Defines to 1 if the m library (-lm) is found])],
             [AC_MSG_ERROR([libm not found])])


# test for libcfu
libcfu_flags=
AC_CHECK_LIB([cfu],[cfuhash_new_with_initial_size],
             [cfulib_flags="-lcfu"
              AC_DEFINE([HAVE_CFULIB],[1],[Defines to 1 if the cfu library (-lcfu) is found])],
             [AC_MSG_ERROR(l[ibcfu not found])])


AC_ARG_WITH([hts],
  [AS_HELP_STRING([--with-hts],
    [location of libhts.a, defaults to /usr/local/lib])],
  [HTS_LIB_LOC="$withval"],
  [HTS_LIB_LOC='/usr/local/lib'])
AC_SUBST([HTS_LIB_LOC])

dnl Write #defines to config.h file.
dnl You need to generate config.h.in, for example using autoheader.
AC_CONFIG_HEADERS([config.h])
dnl Generate Makefile from Makefile.in.
AC_CONFIG_FILES([Makefile])
dnl This generates fore-mentioned files.
AC_OUTPUT
