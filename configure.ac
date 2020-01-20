#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([gb], [v0.16], [ederc@mathematik.uni-kl.de])
AM_INIT_AUTOMAKE([foreign color-tests -Wall -Werror])
# AM_SILENT_RULES([yes])

# Include maintainer mode targets
AM_MAINTAINER_MODE

# set output variable INSTALL to the name of a BSD-compatible install program.
# Requires install-sh to be present as a fallback, even on systems where
# the fallback is not used.
AC_PROG_INSTALL

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT

AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS(config.h)

ac_cv_func_malloc_0_nonnull=yes
ac_cv_func_calloc_0_nonnull=yes
ac_cv_func_realloc_0_nonnull=yes

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_MKDIR_P

AM_PROG_CC_C_O

AC_PROG_CC_C99()
if test "$ac_cv_prog_cc_c99" = "no"; then
  AC_MSG_ERROR([C99 support is required but not found.])
fi

# Checks for libraries.
AC_CHECK_LIB([m], [main])

AC_CONFIG_MACRO_DIR([m4])

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

AC_CHECK_SIZEOF([int *])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday pow memset strchr strstr strtol])

# Check whether --with-gmp was given.
if test "${with_gmp+set}" = set; then :
  if test "$with_gmp" = yes ; then
	GMP_HOME_PATH="${DEFAULT_CHECKING_PATH}"
   elif test "$with_gmp" != no ; then
	GMP_HOME_PATH="$with_gmp ${DEFAULT_CHECKING_PATH}"
    fi
else
  GMP_HOME_PATH="${DEFAULT_CHECKING_PATH}"
fi

BACKUP_CFLAGS=${CFLAGS}
BACKUP_LIBS=${LIBS}

for GMP_HOME in ${GMP_HOME_PATH}
do
  if test "x$GMP_HOME" != "x/usr"; then
    GMP_CPPFLAGS="-I${GMP_HOME}/include"
    GMP_LIBS="-L${GMP_HOME}/lib -Wl,-rpath -Wl,${GMP_HOME}/lib -lgmp"
  else
    GMP_CPPFLAGS=""
    GMP_LIBS="-lgmp"
  fi
  CFLAGS="${BACKUP_CFLAGS} ${GMP_CPPFLAGS}"
  LIBS="${BACKUP_LIBS} ${GMP_LIBS}"
done
AC_SUBST(GMP_CPPFLAGS)
AC_SUBST(GMP_LIBS)
AC_CHECK_HEADERS([gmp.h], ,[AC_MSG_ERROR([GNU MP not found])])
AC_CHECK_LIB(gmp, __gmpz_init, , [AC_MSG_ERROR([GNU MP not found])])
AC_OPENMP
AX_OPENMP(
  [AC_DEFINE(HAVE_OPENMP, 1, [OpenMP enabled])],
  [AC_DEFINE(HAVE_OPENMP, 0, [OpenMP disabled])],[])
AX_CHECK_CL
AX_COUNT_CPUS
# AX_APPEND_COMPILE_FLAGS
# AX_APPEND_LINK_FLAGS
# AX_APPEND_FLAG
AX_GCC_BUILTIN

# May be used later for optimizations
#AX_EXT
#AX_LLVM
#AX_CACHE_SIZE

AC_CONFIG_FILES([Makefile test/Makefile gb.pc])
AC_OUTPUT
