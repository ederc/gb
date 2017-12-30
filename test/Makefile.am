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
AM_CFLAGS = -I$(TOPSRCDIR) -I$(TOPBUILDDIR)

STAGEDIR := $(realpath -s $(TOPBUILDDIR)/.libs)
AM_LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) -lgb -lm
LDFLAGS = -no-install

TESTS = 											\
	test_initialize_basis 			\
	test_free_basis							\
	test_initialize_pairset			\
	test_free_pairset						\
	test_initialize_hash_table	\
	test_free_hash_table 				\
	test_divmask								\
	test_normalize_row					\
	test_julia_io								\
	test_pairset_update

check_PROGRAMS = $(TESTS)
