AUTOMAKE_OPTIONS = gnu subdir-objects
ACLOCAL_AMFLAGS	=	-I m4

AM_CFLAGS	=	-std=gnu99 -march=native -Wextra -Wall -Wfloat-equal -Wundef \
						-Wpointer-arith -Wcast-align -Wstrict-prototypes \
						-Wmissing-prototypes -Wstrict-overflow=5 -Wwrite-strings \
						-Waggregate-return -Wswitch-default \
						-Wswitch-enum -Wconversion -Wunreachable-code \
						-Wno-unused-result -save-temps -fverbose-asm -masm=intel \
						-Werror-implicit-function-declaration -pedantic \
						-Wformat=2 -Wformat-y2k -Wno-format-extra-args \
						-Wno-format-zero-length -Wformat-nonliteral \
						-Wformat-security -Wno-unknown-pragmas -mpopcnt \
						-funroll-loops -finline-functions -fpic \
						${OPENMP_CFLAGS} ${CL_CFLAGS}

#AM_LIBS	=	${CL_LIBS}

lib_LTLIBRARIES = libgb.la

libgb_la_SOURCES = \
	src/hash.h	src/hash.c	\
	src/poly.h	src/poly.c \
	src/spair.h src/spair.c	\
	src/symbol.h src/symbol.c	\
	src/matrix.h src/matrix.c \
	src/gbla-matrix.h src/gbla-matrix.c \
	src/gbla-elimination.h src/gbla-elimination.c \
	src/la.h src/la.c

BUILT_SOURCES = config.h


pkgincludesubdir = $(includedir)/src
pkgincludesub_HEADERS = \
	src/types.h	\
	src/hash.h	\
	src/poly.h \
	src/spair.h	\
	src/symbol.h	\
	src/matrix.h	\
	src/gbla-matrix.h	\
	src/gbla-elimination.h	\
	src/la.h

nodist_pkgincludesub_HEADERS = config.h

EXTRA_DIST=autogen.sh

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = gb.pc

libgb_la_LDFLAGS = $(LDFLAGS) -release $(RELEASE) -no-undefined -rdynamic
libgb_la_LIBADD = $(LIBPNG_LIBADD)

bin_PROGRAMS = gb

# set up the console program. Listing the headers in sources ensure that
# those files are included in distributions.
gb_SOURCES = \
	cli/gb.h	cli/gb.c	\
	cli/io.h	cli/io.c

gb_LDADD = $(top_builddir)/libgb.la $(DEPS_LIBS) ${CL_LIBS} $(GBLA_LIBS) -lm

SUBDIRS =

TESTS = cli-tests poly-tests
