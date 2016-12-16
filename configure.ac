#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([GB], [0.0.1], [ederc@mathematik.uni-kl.de])
AM_INIT_AUTOMAKE([foreign color-tests -Wall -Werror])
AM_SILENT_RULES([yes])

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

PKG_CHECK_MODULES(GBLA, gbla >= 0.0.3)

AC_OPENMP
AX_CHECK_CL
AX_COUNT_CPUS
AX_APPEND_COMPILE_FLAGS
AX_APPEND_LINK_FLAGS
AX_APPEND_FLAG
AX_GCC_BUILTIN

# May be used later for optimizations
#AX_EXT
#AX_LLVM
#AX_CACHE_SIZE

AC_CONFIG_FILES([Makefile gb.pc])
AC_OUTPUT
