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

TESTS =															\
	test_initialize_basis							\
	free_basis										\
	test_initialize_pairset						\
	test_free_pairset									\
	test_initialize_hash_table				\
	test_free_hash_table							\
	test_divmask											\
	test_normalize_row								\
	test_input_meta_data							\
	xtest_input_meta_data							\
	test_julia_io											\
	test_pairset_update								\
	test_select_spairs								\
	test_cyclic_4_drl_17_bit					\
	test_cyclic_4_drl_31_bit					\
	test_cyclic_4_lex_17_bit					\
	test_cyclic_5_drl_17_bit					\
	test_cyclic_5_lex_17_bit					\
	test_cyclic_6_drl_17_bit_la1_t1		\
	test_cyclic_6_drl_17_bit_la1_t2		\
	test_cyclic_6_drl_17_bit_la2_t1		\
	test_cyclic_6_drl_17_bit_la2_t2		\
	test_cyclic_6_drl_17_bit_la42_t1	\
	test_cyclic_6_drl_17_bit_la42_t2	\
	test_cyclic_6_drl_17_bit_la43_t1	\
	test_cyclic_6_drl_17_bit_la43_t2	\
	test_cyclic_6_drl_17_bit_la44_t1	\
	test_cyclic_6_drl_17_bit_la44_t2	\
	test_ideal_q											\
	test_cyclic_7											\
	test_cyclic_9											\
	test_pearson_9										\
	test_katsura_4_lex								\
	test_katsura_5_lex								\
	test_katsura_12_drl_17_bit				\
	test_katsura_13_drl_17_bit				\
	test_katsura_16										\
	test_mayr_42											\
	test_bayes_148										\
	test_cyclic_6_lex_17_bit		

XFAIL_TESTS = 											\
	xtest_input_meta_data

check_PROGRAMS = $(TESTS) $(XFAIL_TESTS)
