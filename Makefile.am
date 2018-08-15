AUTOMAKE_OPTIONS = gnu subdir-objects color-tests
ACLOCAL_AMFLAGS	=	-I m4

AM_CFLAGS	=	-std=gnu99 -march=native -Wextra -Wall -Wfloat-equal -Wundef \
						-Wpointer-arith -Wcast-align -Wstrict-prototypes \
						-Wmissing-prototypes -Wstrict-overflow=1 -Wwrite-strings \
						-Waggregate-return -Wswitch-default \
						-Wswitch-enum -Wconversion -Wunreachable-code \
						-Wno-unused-result -save-temps -fverbose-asm -masm=intel \
						-Werror-implicit-function-declaration -pedantic \
						-Wformat=2 -Wformat-y2k -Wno-format-extra-args \
						-Wno-format-zero-length -Wformat-nonliteral \
						-Wformat-security -Wno-unknown-pragmas -mpopcnt \
						-funroll-loops -finline-functions -fpic \
						${OPENMP_CFLAGS} ${CL_CFLAGS} ${GMP_CFLAGS}

LDADD = ${GMP_LIBS} ${CL_LIBS}

lib_LTLIBRARIES = libgb.la

libgb_la_SOURCES = \
	src/gb.c

BUILT_SOURCES = config.h

pkgincludesubdir = $(includedir)/src
pkgincludesub_HEADERS = \
	src/data.h		\
	src/time.h		\
	src/stat.h		\
	src/f4.h			\
	src/hash.h		\
	src/convert.h	\
	src/order.h		\
	src/symbol.h	\
	src/update.h

nodist_pkgincludesub_HEADERS = config.h

EXTRA_DIST=autogen.sh

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = gb.pc

libgb_la_LDFLAGS = $(LDFLAGS) -release $(RELEASE) -no-undefined -rdynamic
libgb_la_LIBADD = $(LDADD)

SUBDIRS	=	. test
