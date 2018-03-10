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
						${OPENMP_CFLAGS} ${CL_CFLAGS}

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
	test_input_meta_data				\
	xtest_input_meta_data				\
	test_julia_io								\
	test_pairset_update					\
	test_select_spairs					\
	test_f4											\
	test_cyclic_4_drl						\
	test_cyclic_4_lex						\
	test_cyclic_5_lex						\
	test_cyclic_7								\
	test_cyclic_9								\
	test_pearson_9							\
	test_katsura_4_lex					\
	test_katsura_5_lex					\
	test_katsura_16							\
	test_mayr_42

XFAIL_TESTS = 								\
	xtest_input_meta_data

check_PROGRAMS = $(TESTS) $(XFAIL_TESTS)
