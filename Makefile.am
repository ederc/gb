AUTOMAKE_OPTIONS = gnu subdir-objects color-tests
ACLOCAL_AMFLAGS	=	-I m4

AM_CFLAGS	=	-std=c11 -Wextra -Wall -Wfloat-equal -Wundef \
						-Wpointer-arith -Wcast-align -Wstrict-prototypes \
						-Wmissing-prototypes -Wstrict-overflow=1 -Wwrite-strings \
						-Waggregate-return -Wswitch-default \
						-Wswitch-enum -Wconversion -Wunreachable-code \
						-Wno-unused-result -save-temps -fverbose-asm \
						-Werror-implicit-function-declaration -pedantic \
						-Wformat=2 -Wformat-y2k -Wno-format-extra-args \
						-Wno-format-zero-length -Wformat-nonliteral \
						-Wformat-security -Wno-unknown-pragmas \
						-funroll-loops -finline-functions -fpic \
						-march=native \
						${OPENMP_CFLAGS} ${CL_CFLAGS} ${GMP_CFLAGS}

LDADD = ${GMP_LIBS} ${CL_LIBS}

lib_LTLIBRARIES = libgb.la

libgb_la_SOURCES = \
	src/gb.c

BUILT_SOURCES = config.h

pkgincludesubdir = $(includedir)/src
pkgincludesub_HEADERS = \
	src/data.h config.h

nodist_pkgincludesub_HEADERS = config.h

EXTRA_DIST=autogen.sh

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = gb.pc

libgb_la_LDFLAGS = $(LDFLAGS) -no-undefined -rdynamic -version-info 0
libgb_la_LIBADD = $(LDADD)

SUBDIRS	=	. test
