AC_INIT([gb],[0.0.1])

AC_CANONICAL_HOST

AC_CONFIG_SRCDIR(src/dummy.c)

# Set up Automake
dnl foreign: do not create the GNU-specific file COPYING and do not complain
dnl   that GNU-specific files like NEWS, README, AUTHORS and ChangeLog are
dnl   missing.
dnl -Wall: set Automake to emit all warnings it can. Is NOT RELATED to setting
dnl   warnings for other tools. For example, it wil not make the compiler
dnl   get a -Wall option.
dnl subdir-objects: Put object files in a directory structure based on
dnl   the directory structure of the source files. This way, two source
dnl   files with the same name in different directories do not conflict.
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall])

# if --enable-silent-rules is passed to ./configure or if V=0 is passed
# to make, then the compilation output will be much less verbose making
# it possible to spot warnings and errors as they go by.
AM_SILENT_RULES()

# Set up the $(LN_S) macro, which creates symbolic links
AC_PROG_LN_S

LT_INIT

# set output variable INSTALL to the name of a BSD-compatible install program.
# Requires install-sh to be present as a fallback, even on systems where
# the fallback is not used.
AC_PROG_INSTALL

dnl Set the version for the library -- this concerns compatibility of the
dnl source and binary interface of the library and is not the same as the
dnl version of the project.
AC_SUBST([GB_SO_VERSION], [0:0:0])

dnl Include maintainer mode targets.
AM_MAINTAINER_MODE

dnl Needed when reconfiguring with 'autoreconf -i -s'
AC_CONFIG_MACRO_DIR([m4])

dnl Compiling with per-target flags (test_elimination.c) requires AM_PROG_CC_C_O.
AM_PROG_CC_C_O

AC_CONFIG_HEADERS(src/config.h)

AC_PROG_CC_C99()
if test "$ac_cv_prog_cc_c99" = "no"; then
  AC_MSG_ERROR([C99 support is required but not found.])
fi

# LB_CHECK_BLAS
# LB_CHECK_ATLAS
# LB_CHECK_GMP
# LB_CHECK_NTL

# SSE2 support
AC_ARG_ENABLE([sse2],
        AS_HELP_STRING([--disable-sse2], [don't use SSE2 instruction set.]),
	, [if test "$GB_wrapword" = "yes"; then enable_sse2="no"; else enable_sse2="yes"; fi])

AS_IF([test "x$enable_sse2" != "xno"], [
   if test "$GB_wrapword" = "yes"; then
      AC_MSG_ERROR([SSE2 cannot be supported when wrapping word in a C++ class.])
   fi
   case $host_cpu in i[[3456]]86*|x86_64*)
          AX_EXT()
   esac
])
if test x"$ax_cv_have_sse2_ext" = x"yes"; then
  GB_HAVE_SSE2=1
else
  GB_HAVE_SSE2=0
fi
AC_SUBST(GB_HAVE_SSE2)


AC_ARG_WITH(papi,
    AS_HELP_STRING([--with-papi@<:@=PATH@:>@], [The PAPI install prefix, if configure can't find it.]),
    [GB_config_papi=$withval])

AC_ARG_WITH(cachesize,
        AS_HELP_STRING([--with-cachesize@<:@=VALUE@:>@], [L1,L2 and L3 cache sizes in bytes, separated by a colon. Overrides cache tuning.]),[GB_config_cachesize=$withval])

AC_CHECK_HEADER([mm_malloc.h],AC_DEFINE(HAVE_MM_MALLOC,,[Support aligned allocations]),)
if test "$ac_cv_header_mm_malloc_h" = "yes"; then
  GB_HAVE_MM_MALLOC=1
else
  GB_HAVE_MM_MALLOC=0
fi
AC_SUBST(GB_HAVE_MM_MALLOC)

# Correctly working posix_memalign
AX_FUNC_POSIX_MEMALIGN
if test "$ax_cv_func_posix_memalign_works" = "yes"; then
  GB_HAVE_POSIX_MEMALIGN=1
else
  GB_HAVE_POSIX_MEMALIGN=0
fi
AC_SUBST(GB_HAVE_POSIX_MEMALIGN)

# OpenMP support
AX_OPENMP()
AC_SUBST(OPENMP_CFLAGS)
if test -n "$OPENMP_CFLAGS"; then
  GB_HAVE_OPENMP=1
else
  GB_HAVE_OPENMP=0
fi
AC_SUBST(GB_HAVE_OPENMP)

dnl Search for hwloc
AC_CHECK_LIB([hwloc], [main], [
    AC_DEFINE([HAVE_HWLOC],[1],"hwloc available.")
    LIBS="-lhwloc $LIBS"
    ], )

dnl Search for M4RI
 AC_CHECK_LIB([m4ri], [main], [
    AC_DEFINE([HAVE_M4RI],[1],"M4RI available.")
    ], )

dnl Search for Intel TBB
 AC_CHECK_LIB([tbb], [main], [
    AC_DEFINE([HAVE_TBB],[1],"TBB available.")
    ], )

dnl ----- The librt dependency
dnl On Linux TBB calls clock_gettime, which requires librt, but librt is not
dnl linked in automatically. So we need to check for that.

dnl We check if -lrt is necessary. We need librt if we are building with TBB,
dnl if linking and compiling works when linking with librt but it doesn't
dnl without linking with librt.
AC_MSG_CHECKING([if librt is needed to support TBB on this platform]);
oldLIBS=$LIBS;
LIBS="$TBB_LIBS -lrt $LIBS";
oldCFLAGS=$CFLAGS;
CFLAGS="$CFLAGS $TBB_CFLAGS";
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM(
    [[#]include [<]tbb/tbb.h[>]], [[tbb::tick_count::now();]]
  )],
  [LIBS=$oldLibs; AC_LINK_IFELSE(
    [AC_LANG_PROGRAM(
      [[#]include [<]tbb/tbb.h[>]], [[tbb::tick_count::now();]]
    )],
    [AC_MSG_RESULT([no])],
    [AC_MSG_RESULT([yes]); RT_LIBS="-lrt"]
  )],
  [AC_MSG_RESULT([no])]
)];
LIBS=$oldLIBS;
CFLAGS=$oldCFLAGS;


dnl Search for KAAPI
 AC_CHECK_LIB([kaapi], [main], [
    AC_DEFINE([HAVE_KAAPI],[1],"KAAPI available.")
    ], )

dnl Search for KAAPIC
enable_kaapic=0
AC_CHECK_LIB([kaapic], [main], [
    enable_kaapic=1
    AC_DEFINE([HAVE_KAAPIC],[1],"KAAPIC available.")
    ], )

dnl Search for KAAPI++
 AC_CHECK_LIB([kaapi++], [main], [
    AC_DEFINE([HAVE_KAAPIPP],[1],"KAAPI++ available.")
    ], )

AM_CONDITIONAL(ENABLE_KAAPIC,test "$enable_kaapic" -eq 1)

dnl Search for STARPU
enable_starpu=0
 AC_CHECK_LIB([starpu-1.2], [starpu_init], [
    AC_DEFINE([HAVE_STARPU],[1],"STARPU available.")
    enable_starpu=1
    ], )

AM_CONDITIONAL(ENABLE_STARPU,test "$enable_starpu" -eq 1)

dnl Do not silence irrelevant compiler warnings by default, since enabling this
dnl option incurs a performance penalty.
AC_ARG_ENABLE([cc-silence],
  [AS_HELP_STRING([--enable-cc-silence],
                  [Silence irrelevant compiler warnings])],
[if test "x$enable_cc_silence" = "xno" ; then
  enable_cc_silence="0"
else
  enable_cc_silence="1"
fi
],
[enable_cc_silence="0"]
)
if test "x$enable_cc_silence" = "x1" ; then
  AC_DEFINE([CC_SILENCE], [ ], "not silencing irrelevant compiler warnings by
      default, since enabling this option incurs a performance penalty")
fi

# Do not enable allocation from DSS by default.
AC_ARG_ENABLE([dss],
  [AS_HELP_STRING([--enable-dss], [Enable allocation from DSS])],
[if test "x$enable_dss" = "xno" ; then
  enable_dss="0"
else
  enable_dss="1"
fi
],
[enable_dss="0"]
)

# Check whether the BSD/SUSv1 sbrk() exists.  If not, disable DSS support.
AC_CHECK_FUNC([sbrk], [have_sbrk="1"], [have_sbrk="0"])
if test "x$have_sbrk" = "x1" ; then
  AC_DEFINE([HAVE_SBRK], [1],[sbrk available ?])
else
  enable_dss="0"
fi

if test "x$enable_dss" = "x1" ; then
  AC_DEFINE([DSS], [1], [dss available])
fi
AC_SUBST([enable_dss])

# Debugging support
AC_ARG_ENABLE([debug],
	AS_HELP_STRING([--enable-debug], [Enable assert() statements for debugging.]))

AC_ARG_ENABLE([debug-dump],
        AS_HELP_STRING([--enable-debug-dump], [Dump output at exit of every function.]))

if test "x$enable_debug_dump" = "xyes"; then
  GB_DEBUG_DUMP=1
else
  GB_DEBUG_DUMP=0
fi
AC_SUBST(GB_DEBUG_DUMP)

if test "x$enable_debug" = x"yes"; then
  DEBUG_FLAGS="-g"
  CFLAGS="-O0"
  AC_SUBST(DEBUG_FLAGS)
else
  CFLAGS="$CFLAGS -DNDEBUG"
fi

# For the testsuite. Detect if PAPI is installed. See http://icl.cs.utk.edu/papi/ .

if test -z "$GB_config_papi"; then
   AC_CHECK_LIB(papi, PAPI_start_counters,
   [
    AX_GUESS_PATH_LIB(papi)
    AX_GUESS_PATH_HEADER(papi.h)
    if test -n "$LIBPAPI_PATH"; then
       PAPI_LDFLAGS="-Wl,-rpath,$LIBPAPI_PATH"
       PAPI_LIBS="-L$LIBPAPI_PATH -lpapi"
    else
        PAPI_LIBS="-lpapi"
        if ! test -e "/usr/lib/libpapi.so"; then
           AC_MSG_WARN([Could not find libpapi.so. Use --with-papi=<install_prefix> or set LD_LIBRARY_PATH correctly before running benchmark applications.])
        fi
    fi
    if test -n "$PAPI_H_PATH"; then
       PAPI_CFLAGS="-I$PAPI_H_PATH"
       AC_DEFINE_UNQUOTED([HAVE_LIBPAPI], 1, [Define when libpapi is available.])
    else
      AC_MSG_WARN([Could not find papi.h; Use --with-papi=<install_prefix> or add -I<install_prefix>/include to either CPPFLAGS
                    or CFLAGS, or turn off papi all together by configuring with --without-papi.])
    fi
])
fi

if test x"$GB_config_papi" != x"no" && test -n "$gb_config_papi"; then
    LIBPAPI_PATH="`realpath -s $GB_config_papi/lib`"
    PAPI_H_PATH="`realpath -s $GB_config_papi/include`"
    PAPI_CFLAGS="-I$PAPI_H_PATH"
    PAPI_LDFLAGS="-Wl,-rpath,$LIBPAPI_PATH"
    PAPI_LIBS="-L$LIBPAPI_PATH -lpapi"
    AC_DEFINE_UNQUOTED([HAVE_LIBPAPI], 1, [Define when libpapi is available.])
fi

AC_SUBST(PAPI_LIBS)
AC_SUBST(PAPI_LDFLAGS)
AC_SUBST(PAPI_CFLAGS)

AC_ARG_ENABLE([cachetune],
     AS_HELP_STRING([--enable-cachetune],[calculate cache size from timing information (deprecated).]))

# Cache Sizes
if test -z $GB_config_cachesize; then
   AX_CACHE_SIZE()
   AS_IF([test "x$enable_cachetune" = "xyes"], [AC_MSG_WARN(--enable-cachetune is deprecated since it usually does not provide optimal parameters.) AX_CACHE_SIZE_TUNE()])
else
  AS_IF([test "x$enable_cachetune" = "xyes"], [AC_MSG_WARN(Ignoring cache tuning since --with-cachesize was given.)])

  ax_l1_size=`echo $GB_config_cachesize | cut -d ":" -f 1`
  ax_l2_size=`echo $GB_config_cachesize | cut -d ":" -f 2`
  ax_l3_size=`echo $GB_config_cachesize | cut -d ":" -f 3`

  GB_CPU_L1_CACHE=${ax_l1_size}
  GB_CPU_L2_CACHE=${ax_l2_size}
  GB_CPU_L3_CACHE=${ax_l3_size}
  AC_SUBST(GB_CPU_L1_CACHE)
  AC_SUBST(GB_CPU_L2_CACHE)
  AC_SUBST(GB_CPU_L3_CACHE)

fi

# PNG

have_libpng="no"
AC_ARG_ENABLE([png],
   [AC_HELP_STRING([--disable-png], [disable PNG support @<:@default=enabled@:>@])],
   [
    if test "x${enableval}" = "xyes" ; then
       want_png="yes"
    else
       want_png="no"
    fi
   ],
   [want_png="yes"])

AC_MSG_CHECKING([whether to build with PNG support])
AC_MSG_RESULT([${want_png}])
PKG_CHECK_MODULES(GBLA, gbla >= 0.0.3)
if test "x${want_png}" = "xyes" ; then
   PKG_CHECK_MODULES([PNG], [libpng],
      [have_libpng="yes"; LIBPNG_LIBADD=`pkg-config --libs libpng`],
      [have_libpng="no"])
   if ! test "x${have_libpng}" = "xyes" ; then
      AC_CHECK_LIB([png],
         [png_create_write_struct],
         [have_libpng="yes"; LIBPNG_LIBADD="-lpng"],
         [AC_CHECK_LIB([png14],
            [png_create_write_struct],
            [have_libpng="yes"; LIBPNG_LIBADD="-lpng14"],
            [AC_CHECK_LIB([png12],
               [png_create_write_struct],
               [have_libpng="yes"; LIBPNG_LIBADD="-lpng12"],
               [AC_CHECK_LIB([png10],
                  [png_create_write_struct],
                  [have_libpng="yes"; LIBPNG_LIBADD="-lpng10"],
                  [have_libpng="no"])
               ])
            ])
        ])
   fi
   if test "x${have_libpng}" = "xno" ; then
      AC_MSG_WARN([Can not find a usuable PNG library. Make sure that CPPFLAGS and LDFLAGS are correctly set.])
   fi
fi

if test "x${have_libpng}" = "xyes" ; then
   GB_HAVE_LIBPNG=1
   AC_SUBST(GB_HAVE_LIBPNG)
   AC_SUBST(LIBPNG_LIBADD)
else
   GB_HAVE_LIBPNG=0
   AC_SUBST(GB_HAVE_LIBPNG)
fi

AC_PROG_MAKE_SET

AC_CONFIG_FILES([Makefile src/gb_config.h gb.pc])
AC_OUTPUT

