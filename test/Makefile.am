# TOPSRCDIR allows to compile this testsuite for a clone (and different revision) of this
# repository. If the environment variable topsrcdir is the top source directory of the clone,
# compile this testsuite as:
#
# make TOPSRCDIR="$topsrcdir"
#
# Finally, if you also leave out TOPSRCDIR and just run 'make' then the testsuite
# is compiled against the current source tree.

AUTOMAKE_OPTIONS = foreign
TOPSRCDIR = $(srcdir)/..
TOPBUILDDIR = $(builddir)/..
DEFINES =
# include TOPBUILDIR for m4ri_config.h
AM_CFLAGS = -I$(TOPSRCDIR) -I$(TOPBUILDDIR) \
						${OPENMP_CFLAGS} ${CL_CFLAGS} ${GMP_CFLAGS}

STAGEDIR := $(realpath -s $(TOPBUILDDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) -lgb ${GMP_LIBS} ${CL_LIBS}
LDFLAGS = -no-install

TESTS =									\
	cyclic_4_drl_15_bit		\
	cyclic_5_drl_15_bit		\
	cyclic_4_drl_32_bit					

check_PROGRAMS = $(TESTS) $(XFAIL_TESTS)
